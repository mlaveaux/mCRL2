// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/data/detail/rewrite/innermost.h"

#include "mcrl2/atermpp/algorithm.h"
#include "mcrl2/utilities/logger.h"
#include "mcrl2/utilities/stack_array.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/replace.h"
#include "mcrl2/data/detail/rewrite/jitty_jittyc.h"
#include "mcrl2/data/detail/rewrite/substitute.h"

#include <assert.h>

// These are debug related options.

constexpr bool PrintRewriteSteps = false;

// The following options toggle tracking metrics.

/// \brief Keep track of the number of times each rule is applied (and print it after each rewrite).
constexpr bool CountRewriteSteps = false;

/// \brief Counts the number of times that the cache was used succesfull, its size and the number of insertions.
constexpr bool CountRewriteCacheMetric = false;

// The following options toggle performance features.

/// \brief Keep track of terms that are in normal form during rewriting.
constexpr bool EnableNormalForms = true;

/// \brief Enable caching of rewrite results.
constexpr bool EnableCaching = true;

/// \brief Enables construction stacks to reconstruct the right-hand sides bottom up.
constexpr bool EnableConstructionStack = false;

using namespace mcrl2;
using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

/// \returns The same data equation with all variables renamed to meta-variables that can not occur in the terms on which is matched.
template<typename Generator>
data_equation rename_meta_variables(const data_equation& equation, Generator& generator)
{
  auto variables = find_all_variables(equation);

  mutable_indexed_substitution<variable> sigma;
  for (auto& var: variables)
  {
    sigma[var] = variable(generator(), var.sort());
  }

  return replace_variables(equation, sigma);
}

/// \brief Checks every equation in the given data specification.
/// \returns A vector of equations from the data specifications that pass the given selector.
data_equation_vector filter_data_specification(const data_specification& data_spec, const used_data_equation_selector& selector)
{
  data_equation_vector equations;
  enumerator_identifier_generator generator("@");

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

      equations.emplace_back(rename_meta_variables(equation, generator));
    }
  }

  return equations;
}


InnermostRewriter::InnermostRewriter(const data_specification& data_spec, const used_data_equation_selector& selector)
  : Rewriter(data_spec, selector),
    m_rewrite_cache(1024),
    m_matcher(filter_data_specification(data_spec, selector))
{}

data_expression InnermostRewriter::rewrite(const data_expression& term, substitution_type& sigma)
{
  auto result = rewrite_impl(term, sigma);
  print_rewrite_metrics();
  m_normal_forms.clear();
  return result;
}

// Private functions

template<typename Substitution>
data_expression InnermostRewriter::apply_substitution(const data_expression& term, Substitution& sigma, const ConstructionStack& stack)
{
  if (EnableConstructionStack)
  {
    return stack.construct_term(sigma, m_generator, m_argument_stack);
  }
  else
  {
    return capture_avoiding_substitution(term, sigma, m_generator);
  }
}

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
    return rewrite_single(function_symbol);
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

data_expression InnermostRewriter::rewrite_abstraction(const abstraction& abstraction, const substitution_type& sigma)
{
  // u' := rewrite(u, sigma[x := y]) where y are fresh variables.
  m_local_sigma.clear();
  data::variable_list new_variables = rename_bound_variables(abstraction, m_local_sigma, m_generator);
  data_expression u = capture_avoiding_substitution(abstraction.body(), m_local_sigma, m_generator);

  if (PrintRewriteSteps)
  {
    mCRL2log(info) << "Applied alpha-conversion to " << abstraction << " resulting in " << data::abstraction(abstraction.binding_operator(), new_variables, u) << "\n";
  }

  // rewrite(u, sigma[x := y]) is equivalent to rewrite(u^[x := y], sigma);
  data_expression body_rewritten = rewrite_impl(u, sigma);

  // Return lambda y . u'
  auto result = data::abstraction(abstraction.binding_operator(), new_variables, body_rewritten);

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
    return rewrite_single(application(head_rewritten, arguments.begin(), arguments.end()));
  }
}

data_expression InnermostRewriter::rewrite_single(const data_expression& expression)
{
  if (is_normal_form(expression))
  {
    // By definition a normal form does not match any rewrite rule.
    return expression;
  }

  // (R, sigma') := match(h'(u_1', ..., u_n')),
  mutable_indexed_substitution<variable, data_expression> m_local_sigma;
  auto match_result = m_matcher.match(expression, m_local_sigma);

  // If R not empty
  for (const auto& match : match_result)
  {
    const auto& equation = std::get<0>(match.get());

    // Compute rhs^sigma'.
    auto rhs = apply_substitution(equation.rhs(), m_local_sigma, std::get<1>(match.get()));

    // Delaying rewriting the condition ensures that the matching substitution does not have to be saved.
    if (equation.condition() != sort_bool::true_() &&
      rewrite_impl(apply_substitution(equation.condition(), m_local_sigma, std::get<2>(match.get())), m_identity) != sort_bool::true_())
    {
      continue;
    }

    if (CountRewriteSteps)
    {
      ++m_application_count[equation];
    }

    // Return rewrite(r^sigma', id)
    auto result = rewrite_impl(rhs, m_identity);

    if (EnableCaching)
    {
      m_rewrite_cache.emplace(expression, result);
    }

    if (PrintRewriteSteps)
    {
      mCRL2log(info) << "Rewrote  " << expression << " to " << result << " using rule " << equation << "\n";
    }

    return result;
  }

  // Return h'(u_1', ..., u_n')
  return expression;
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
