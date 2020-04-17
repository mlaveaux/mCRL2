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
#include "mcrl2/data/detail/rewrite/substitute.h"

#include <assert.h>

// These are debug related options.

constexpr bool PrintRewriteSteps = false;

// Enable term rewrite engine features.

/// \brief Enable higher-order rewriting (also the head symbols).
constexpr bool EnableHigherOrder = false;

/// \brief Enable conditional rewriting, otherwise conditions cause exceptions.
constexpr bool EnableConditions = true;

/// \brief Obtain all normal forms to check that the term rewrite system is confluent.
constexpr bool EnableCheckConfluence = false;

// The following options toggle tracking metrics.

/// \brief Keep track of the number of times each rule is applied (and print it after each rewrite).
constexpr bool CountRewriteSteps = false;

/// \brief Counts the number of times that the cache was used succesfull, its size and the number of insertions.
constexpr bool CountRewriteCacheMetric = false;

// The following options toggle performance features.

/// \brief Keep track of terms that are in normal form during rewriting.
enum class NormalForm
{
  None,
  Set, /// \brief Use an unordered_set to keep track of terms in normal form.
  Tag  /// \brief Use a special term to keep track of terms in normal form.
};

constexpr NormalForm TrackNormalForms = NormalForm::Set;

/// \brief Enable caching of rewrite results.
constexpr bool EnableCaching = false;

/// \brief Enables construction stacks to reconstruct the right-hand sides bottom up.
constexpr bool EnableConstructionStack = true;

using namespace mcrl2;
using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

/// \brief Checks every equation in the given data specification.
/// \returns A vector of equations from the data specifications that pass the given selector.
data_equation_vector filter_data_specification(const data_specification& data_spec, const used_data_equation_selector& selector)
{
  data_equation_vector equations;

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

      equations.emplace_back(equation);
    }
  }

  return equations;
}


InnermostRewriter::InnermostRewriter(const data_specification& data_spec, const used_data_equation_selector& selector)
  : Rewriter(data_spec, selector),
    m_rewrite_cache(1024),
    m_matcher(filter_data_specification(data_spec, selector))
{
  mCRL2log(debug) << "InnermostRewriter (EnableHigherOrder = " << EnableHigherOrder
    << ", EnableConditions = " << EnableConditions
    << ", EnableCheckConfluence = " << EnableCheckConfluence << "):\n";

  mCRL2log(debug) << "  Performance features: TrackNormalForms = " << static_cast<std::size_t>(TrackNormalForms)
    << ", EnableCaching = " << EnableCaching
    << ", EnableConstructionStack = " << EnableConstructionStack << ".\n";
}

data_expression InnermostRewriter::rewrite(const data_expression& term, substitution_type& sigma)
{
  auto result = rewrite_impl(term, sigma);
  print_rewrite_metrics();

  // Clear the tracked normal forms to avoid unnecessary memory usage.
  if (TrackNormalForms == NormalForm::Set) { m_normal_forms.clear(); }
  return result;
}

// Private functions

template<typename Substitution, typename Transformer>
data_expression InnermostRewriter::apply_substitution(const data_expression& term, Substitution& sigma, const ConstructionStack& stack, Transformer f)
{
  if (EnableConstructionStack)
  {
    return stack.construct_term(sigma, m_generator, m_argument_stack, f);
  }
  else
  {
    return capture_avoiding_substitution(term, sigma, m_generator, f);
  }
}

data_expression InnermostRewriter::rewrite_impl(const data_expression& term, const substitution_type& sigma)
{
  if (is_normal_form(term))
  {
    // By definition a normal form does not match any rewrite rule.
    if (PrintRewriteSteps) { mCRL2log(info) << "Term " << term << " is in normal form.\n"; }

    if constexpr (TrackNormalForms == NormalForm::Tag)
    {
      return static_cast<data_expression>(term[1]);
    }
    else
    {
      return term;
    }
  }

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
  else if (is_where_clause(term))
  {
    const where_clause& w = atermpp::down_cast<where_clause>(term);
    return rewrite_where_clause(w, sigma);
  }
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
  if (data::is_lambda_binder(abstraction.binding_operator()))
  {
    // u' := rewrite(u, sigma[x := y]) where y are fresh variables.
    m_local_sigma.clear();
    data::variable_list new_variables = rename_bound_variables(abstraction, m_local_sigma, m_generator);
    data_expression u = capture_avoiding_substitution(abstraction.body(), m_local_sigma, m_generator, identity);

    // rewrite(u, sigma[x := y]) is equivalent to rewrite(u^[x := y], sigma);
    data_expression body_rewritten = rewrite_impl(u, sigma);

    // Return lambda y . u'
    auto result = data::abstraction(abstraction.binding_operator(), new_variables, body_rewritten);

    if (PrintRewriteSteps)
    {
      mCRL2log(info) << "Applied alpha-conversion to " << abstraction << " resulting in " << result << "\n";
    }

    return static_cast<data_expression>(result);
  }
  else if (data::is_exists_binder(abstraction.binding_operator()))
  {
    return existential_quantifier_enumeration(abstraction, m_identity);
  }
  else if (data::is_forall_binder(abstraction.binding_operator()))
  {
    return universal_quantifier_enumeration(abstraction, m_identity);
  }

  assert(false);
  return abstraction;
}

data_expression InnermostRewriter::rewrite_application(const application& appl, const substitution_type& sigma)
{
  // h' := rewrite(h, sigma)
  auto head_rewritten = (EnableHigherOrder ? rewrite_impl(appl.head(), sigma) : appl.head());
  if (!EnableHigherOrder)
  {
    if (!is_function_symbol(head_rewritten))
    {
      mCRL2log(error) << "Term " << appl << " is higher-order.\n";
      throw mcrl2::runtime_error("Higher-order rewriting is disabled (EnableHigherOrder = false).");
    }
  }

  // For i in {1, ..., n} do u' := rewrite(u, sigma)
  MCRL2_DECLARE_STACK_ARRAY(arguments, data_expression, appl.size());
  for (std::size_t index = 0; index < appl.size(); ++index)
  {
    arguments[index] = rewrite_impl(appl[index], sigma);
  }

  // If h' is of the form lambda x . w
  if (is_abstraction(head_rewritten))
  {
    const auto& abstraction = static_cast<const data::abstraction&>(head_rewritten);
    assert(data::is_lambda_binder(abstraction.binding_operator()));

    // rewrite(w, sigma[x gets u']) is equivalent to rewrite(w^[x := u'], id).
    m_local_sigma.clear();

    std::size_t index = 0;
    for (auto& variable : abstraction.variables())
    {
      m_local_sigma[variable] = arguments[index];
      ++index;
    }

    // Return rewrite(w, sigma[x := u'])
    data_expression result = rewrite_impl(capture_avoiding_substitution(abstraction.body(), m_local_sigma, m_generator, identity), m_identity);

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
  // (R, sigma') := match(h'(u_1', ..., u_n')),
  //std::set<data_expression> results;

  mutable_map_substitution<> matching_sigma;
  for(auto it = m_matcher.match(expression, matching_sigma); *it != nullptr; ++it)
  {
    // If R not empty
    const extended_data_equation* result = *it;
    const extended_data_equation& match = *result;

    // Compute rhs^sigma'.
    auto rhs = apply_substitution(match.equation().rhs(), matching_sigma, match.rhs_stack(),
      [this](const data_expression& expression)
      {
        return mark_normal_form(expression);
      });

    // Delaying rewriting the condition ensures that the matching substitution does not have to be saved.
    if (match.equation().condition() != sort_bool::true_())
    {
      if (EnableConditions)
      {
        if (PrintRewriteSteps) { mCRL2log(info) << "Rewriting condition " << match.equation().condition() << "\n"; }
        if (rewrite_impl(apply_substitution(match.equation().condition(), matching_sigma, match.condition_stack(), identity), m_identity) != sort_bool::true_())
        {
          continue;
        }
      }
      else
      {
        throw mcrl2::runtime_error("Conditional rewriting (EnableConditions) is disabled.");
      }
    }

    if (CountRewriteSteps) { ++m_application_count[match.equation()]; }

    if (PrintRewriteSteps) { mCRL2log(info) << "Rewrote " << expression << " to " << rhs << ".\n"; } // using rule " << match.equation() << "\n"; }

    // Return rewrite(r^sigma', id)
    if constexpr (!EnableCheckConfluence)
    {
      if (PrintRewriteSteps) { mCRL2log(info) << "Term " << expression << " is in normal form.\n"; }
      auto result = rewrite_impl(rhs, m_identity);

      if (EnableCaching) { m_rewrite_cache.emplace(expression, result); }
      return result;
    }
    else
    {
      //results.insert(rewrite_impl(rhs, m_identity));
    }
  }

  /*if constexpr (EnableCheckConfluence && !results.empty())
  {
    if constexpr (EnableCaching) { m_rewrite_cache.emplace(expression, *results.begin()); }

    if (results.size() > 1)
    {
      mCRL2log(info) << "Term rewrite system is not confluent as term " << expression << " has the following normal forms:.\n";
      for (const auto& result : results)
      {
        mCRL2log(info) << "  " << result << "\n";
      }
    }

    // Return h'(u_1', ..., u_n')
    return *results.begin();
  }*/

  return expression;
}

data_expression InnermostRewriter::rewrite_where_clause(const where_clause& clause, const substitution_type& sigma)
{
  m_local_sigma.clear();
  data::variable_list new_variables = rename_bound_variables(clause, m_local_sigma, m_generator);
  data_expression u = capture_avoiding_substitution(clause.body(), m_local_sigma, m_generator, identity);

  const data_expression result = rewrite_impl(u, sigma);
  return result;
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

// The function symbol below is used to administrate that a term is in normal form. It is put around a term.
// Terms with this auxiliary function symbol cannot be printed using the pretty printer for data expressions.
static const function_symbol& normal_form_tag()
{
  static const function_symbol this_term_is_in_normal_form(
    std::string("Rewritten@@term"),
    function_sort({ untyped_sort() },untyped_sort()));
  return this_term_is_in_normal_form;
}


data_expression add_normal_form_tag(const data_expression& term)
{
  return application(normal_form_tag(), term);
}

bool InnermostRewriter::is_normal_form(const data_expression& term) const
{
  if constexpr (TrackNormalForms == NormalForm::Set)
  {
    return (m_normal_forms.count(term) != 0);
  }
  else if constexpr (TrackNormalForms == NormalForm::Tag)
  {
    return (term[0] == normal_form_tag());
  }

  return false;
}

data_expression InnermostRewriter::mark_normal_form(const data_expression& term)
{
  if constexpr (TrackNormalForms == NormalForm::Set)
  {
    m_normal_forms.emplace(term);
    return term;
  }
  else if constexpr (TrackNormalForms == NormalForm::Tag)
  {
    return add_normal_form_tag(term);
  }

  return term;
}
