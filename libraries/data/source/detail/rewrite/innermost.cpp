// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/data/detail/rewrite/innermost.h"

#include "mcrl2/utilities/logger.h"
#include "mcrl2/utilities/stack_array.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/detail/rewrite/jitty_jittyc.h"
#include "mcrl2/data/detail/rewrite/substitute.h"

#include <assert.h>

// These are debug related options.

constexpr bool PrintRewriteSteps = false;
constexpr bool PrintMatchSteps   = false;

// The following options toggle tracking metrics.

/// \brief Keep track of the number of times each rule is applied (and print it after each rewrite).
constexpr bool CountRewriteSteps = false;

/// \brief Counts the number of times that the cache was used succesfull, its size and the number of insertions.
constexpr bool CountRewriteCacheMetric = false;

// The following options toggle performance features.

/// \brief Keep track of terms that are in normal form during rewriting.
constexpr bool EnableNormalForms = true;

/// \brief Enables construction stacks to reconstruct the right-hand sides bottom up.
constexpr bool EnableConstructionStack = true;

/// \brief Enable caching of rewrite results.
constexpr bool EnableCaching = true;

using namespace mcrl2;
using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

/// \returns A unique index for the head symbol that the given term starts with.
static inline std::size_t get_head_index(const data_expression& term)
{
  return core::index_traits<data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(get_nested_head(term)));
}

InnermostRewriter::InnermostRewriter(const data_specification& data_spec, const used_data_equation_selector& selector)
  : Rewriter(data_spec, selector),
    m_rewrite_cache(1024)
{
  // This probably belongs inside the Rewriter class, checks for all the selected equations whether they are valid rewrite rules.
  for (const data_equation& equation : data_spec.equations())
  {
    if (selector(equation))
    {
      try
      {
        CheckRewriteRule(equation);
      }
      catch (std::runtime_error& e)
      {
        mCRL2log(warning) << e.what() << std::endl;
        continue;
      }

      // Make sure that it is possible to insert the match data for head_index left-hand side.
      std::size_t head_index = get_head_index(equation.lhs());
      if (head_index >= m_rewrite_system.size()) { m_rewrite_system.resize(head_index + 1); }

      // Insert the left-hand side into the rewrite rule mapping and a construction stack for its right-hand side.
      m_rewrite_system[head_index].emplace_back(equation,
        ConstructionStack(equation.condition()),
        ConstructionStack(equation.rhs()));
    }
  }
}

data_expression InnermostRewriter::rewrite(const data_expression& term, substitution_type& sigma)
{
  auto result = rewrite_impl(term, sigma);
  print_rewrite_metrics();
  m_normal_forms.clear();
  return result;
}

// Private functions

data_expression InnermostRewriter::rewrite_impl(const data_expression& term, const substitution_type& sigma)
{
  // If t in variables
  if (is_variable(term))
  {
    const auto& var = static_cast<const data::variable&>(term);
    return sigma(var);
  }
  // Else if t in function_symbols (This is an extra case handled by h() in the pseudocode).
  else if (is_function_symbol(term))
  {
    const auto& function_symbol = static_cast<const data::function_symbol&>(term);
    return rewrite_function_symbol(function_symbol);
  }
  // Else if t is of the form lambda x . u
  else if (is_abstraction(term))
  {
    const auto& abstraction = static_cast<const data::abstraction&>(term);
    return rewrite_abstraction(abstraction, sigma);
  }
  // Else (t is of the form h(u_1, ..., u_n)).
  else
  {
    assert(is_application(term));
    const auto& appl = static_cast<const data::application&>(term);

    if (is_normal_form(appl))
    {
      return appl;
    }

    if (EnableCaching)
    {
      // If the cache already contains the normal form for the given term return the result immediately.
      auto it = m_rewrite_cache.find(appl);
      if (it != m_rewrite_cache.end())
      {
        if (CountRewriteCacheMetric) { m_rewrite_cache_metric.hit(); }

        if (PrintRewriteSteps)
        {
          mCRL2log(info) << "Found " << term << " to " << (*it).second << " in the rewrite cache.\n";
        }

        return (*it).second;
      }
      else if (CountRewriteCacheMetric)
      {
        m_rewrite_cache_metric.miss();
      }
    }

    return rewrite_application(appl, sigma);
  }
}

data_expression InnermostRewriter::rewrite_function_symbol(const function_symbol& symbol)
{
  data_expression rhs;
  // If R not empty, this match function already applies the substitution.
  if (match(symbol, rhs))
  {
    // Return rewrite(r^sigma', id)
    return rewrite_impl(rhs, m_identity);
  }
  else
  {
    return symbol;
  }
}

data_expression InnermostRewriter::rewrite_abstraction(const abstraction& abstraction, const substitution_type& sigma)
{
  // u' := rewrite(u, sigma[x := y]) where y are fresh variables.
  m_local_sigma.clear();
  data::variable_list new_variables = rename_bound_variables(abstraction, m_local_sigma, m_generator);

  // rewrite(u, sigma[x := y]) is equivalent to rewrite(u^[x := y], sigma);
  data_expression body_rewritten = rewrite_impl(capture_avoiding_substitution(abstraction.body(), m_local_sigma, m_generator), sigma);

  // Return lambda y . u'
  auto result = data::abstraction(abstraction.binding_operator(), new_variables, body_rewritten);

  if (PrintRewriteSteps)
  {
    mCRL2log(info) << "Applied alpha-conversion to " << abstraction << " resulting in " << result << "\n";
  }

  return static_cast<data_expression>(result);
}

data_expression InnermostRewriter::rewrite_application(const application& appl, const substitution_type& sigma)
{
  // h' := rewrite(h, sigma)
  auto head_rewritten = rewrite_impl(appl.head(), sigma);
  mark_normal_form(head_rewritten);

  // For i in {1, ..., n} do u' := rewrite(u, sigma)
  MCRL2_DECLARE_STACK_ARRAY(arguments, data_expression, appl.size());
  for (std::size_t index = 0; index < appl.size(); ++index)
  {
    arguments[index] = rewrite_impl(appl[index], sigma);
    mark_normal_form(arguments[index]);
  }

  // If h' is of the form lambda x . w
  if (is_abstraction(head_rewritten))
  {
    const auto& abstraction = static_cast<const data::abstraction&>(head_rewritten);

    // rewrite(w, sigma[x gets u']) is equivalent to rewrite(w^[x := u'], id).
    std::size_t index = 0;
    m_local_sigma.clear();

    for (auto& variable : abstraction.variables())
    {
      m_local_sigma[variable] = arguments[index];
      ++index;
    }

    // Return rewrite(w, sigma[x := u'])
    data_expression result = rewrite_impl(capture_avoiding_substitution(abstraction.body(), m_local_sigma, m_generator), m_identity);

    if (PrintRewriteSteps)
    {
      mCRL2log(info) << "Applied beta-reduction to" << appl << " resulting in " << result << "\n";
    }

    return result;
  }
  else
  {
    application new_appl(head_rewritten, arguments.begin(), arguments.end());

    // (R, sigma') := match(h'(u_1', ..., u_n')),
    data_expression rhs;

    // If R not empty, this match function already applies the substitution and rewrite steps.
    if (match(new_appl, rhs))
    {
      // Return rewrite(r^sigma', id)
      auto result = rewrite_impl(rhs, m_identity);

      if (EnableCaching)
      {
        m_rewrite_cache.emplace(new_appl, result);
      }

      return result;
    }
    else
    {
      // Return h'(u_1', ..., u_n')
      return static_cast<data_expression>(new_appl);
    }
  }
}

void InnermostRewriter::print_rewrite_metrics()
{
  if (CountRewriteSteps)
  {
    // We are going to sort them by the number of times applied.
    std::vector<std::pair<data_equation, std::size_t>> counts;
    for (auto& result : m_application_count)
    {
      counts.push_back(result);
    }

    std::sort(counts.begin(),
      counts.end(),
      [](const std::pair<data_equation, std::size_t>& left, std::pair<data_equation, std::size_t>& right) -> bool
      {
        return left.second > right.second;
      }
    );

    // Count the total number of rewrite steps applied
    std::size_t total_count = 0;
    for (auto& result : counts)
    {
      total_count += result.second;
    }

    mCRL2log(info) << "Applied " << total_count << " single rewrite steps.\n";

    for (auto& result : counts)
    {
      mCRL2log(info) << "Applied rule " << result.first << " " << result.second << " times.\n";
    }

    m_application_count.clear();
  }

  if (CountRewriteCacheMetric)
  {
    mCRL2log(info) << "Rewrite cache " << m_rewrite_cache_metric.message() << ".\n";
  }
}

bool InnermostRewriter::is_normal_form(const data_expression& term) const
{
  return EnableNormalForms && (m_normal_forms.count(term) != 0);
}

void InnermostRewriter::mark_normal_form(const data_expression& term)
{
  if (EnableNormalForms)
  {
    m_normal_forms.emplace(term);
  }
}

template<typename Substitution>
data_expression InnermostRewriter::apply_substitution(const data_expression& term, Substitution& sigma, const ConstructionStack& stack)
{
  if (EnableConstructionStack)
  {
    return stack.construct_term(sigma, m_argument_stack);
  }
  else
  {
    return capture_avoiding_substitution(term, sigma, m_generator);
  }
}

bool InnermostRewriter::match(const data_expression& term, data_expression& rhs)
{
  if (is_normal_form(term))
  {
    // By definition a normal form does not match any rewrite rule.
    return false;
  }

  std::size_t head_index = get_head_index(term);
  if (head_index >= m_rewrite_system.size())
  {
    // No left-hand side starts with this head symbol, so it cannot match.
    return false;
  }

  // Searches for a left-hand side and a substitution such that when the substitution is applied to this left-hand side it is (syntactically) equivalent
  // to the given term. However, only tries rewrite rules that start with the correct head symbol.
  for (const auto& tuple : m_rewrite_system[head_index])
  {
    const auto& equation = std::get<0>(tuple);
    const auto& condition_stack = std::get<1>(tuple);
    const auto& rhs_stack = std::get<2>(tuple);

    // Compute a matching substitution for each rule and check that the condition associated with that rule is true, either trivially or by rewrite(c^sigma, identity).
    m_local_sigma.clear();
    if (match_lhs(term, equation.lhs(), m_local_sigma))
    {
      if(PrintRewriteSteps)
      {
        mCRL2log(info) << "Matched rule " << equation << " to term " << term << "\n";
      }

      rhs = apply_substitution(equation.rhs(), m_local_sigma, rhs_stack);

      // Delaying rewriting the condition ensures that the matching substitution does not have to be saved.
      if (equation.condition() != sort_bool::true_() &&
        rewrite_impl(apply_substitution(equation.condition(), m_local_sigma, condition_stack), m_identity) != sort_bool::true_())
      {
        continue;
      }

      if (CountRewriteSteps)
      {
        ++m_application_count[equation];
      }

      return true;
    }
    else if (PrintMatchSteps)
    {
      mCRL2log(info) << "Tried rule " << equation << " to term " << term << "\n";
    }
  }

  if (PrintMatchSteps)
  {
    mCRL2log(info) << "Term " << term << " is in normal-form.\n";
  }

  return false;
}
