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

#include <assert.h>

constexpr bool PrintRewriteSteps = false;

constexpr bool PrintMatchSteps   = false;

/// \brief Keep track of terms that are in normal form during rewriting.
constexpr bool EnableNormalForms = true;

/// \brief Enables construction stacks to reconstruct the right-hand sides bottom up.
constexpr bool EnableConstructionStack = true;

using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

InnermostRewriter::InnermostRewriter(const data_specification& data_spec, const used_data_equation_selector& selector)
  : Rewriter(data_spec, selector)
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

      // Insert the left-hand side into the rewrite rule mapping and a construction stack for its right-hand side.
      m_rewrite_system[get_nested_head(equation.lhs())].emplace_back(equation,
        ConstructionStack(equation.condition()),
        ConstructionStack(equation.rhs()));
    }
    else
    {
      mCRL2log(info) << "Rule " << equation << " ignored by selector.\n";
    }
  }
}

data_expression InnermostRewriter::rewrite(const data_expression& term, substitution_type& sigma)
{
  auto result = rewrite_impl(term, sigma);
  m_normal_forms.clear();
  return result;
}

// Private functions

data_expression InnermostRewriter::rewrite_impl(const data_expression& term, substitution_type& sigma)
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
    return rewrite_application(appl, sigma);
  }
}

data_expression InnermostRewriter::rewrite_function_symbol(const function_symbol& symbol)
{
  data_expression rhs;
  // If R not empty, this match function already applies the substitution.
  if (match(symbol, rhs))
  {
    // Return rewrite(r^sigma', sigma)
    return rewrite_impl(rhs, m_identity);
  }

  return symbol;
}

data_expression InnermostRewriter::rewrite_abstraction(const abstraction& abstraction, substitution_type& sigma)
{
  // u' := rewrite(u, sigma[x := y]) where y are fresh variables.
  data::variable_list new_variables;
  for (const auto& var : abstraction.variables())
  {
    const variable fresh_variable(m_generator(), var.sort());
    new_variables.push_front(fresh_variable);
    sigma[var] = fresh_variable;
  }

  data_expression body_rewritten = rewrite_impl(abstraction.body(), sigma);

  // Return lambda y . u'
  auto result = data::abstraction(abstraction.binding_operator(), new_variables, body_rewritten);

  if (PrintRewriteSteps)
  {
    mCRL2log(info) << "Applied alpha-reduction to " << abstraction << " resulting in " << result << "\n";
  }

  return static_cast<data_expression>(result);
}

data_expression InnermostRewriter::rewrite_application(const application& appl, substitution_type& sigma)
{
  // h' := rewrite(h, sigma)
  auto head_rewritten = is_normal_form(appl.head()) ? appl.head() : rewrite_impl(appl.head(), sigma);
  mark_normal_form(head_rewritten);

  // For i in {1, ..., n} do u' := rewrite(u, sigma)
  MCRL2_DECLARE_STACK_ARRAY(arguments, data_expression, appl.size());
  for (std::size_t index = 0; index < appl.size(); ++index)
  {
    if (is_normal_form(appl[index]))
    {
      arguments[index] = appl[index];
    }
    else
    {
      arguments[index] = rewrite_impl(appl[index], sigma);
      mark_normal_form(arguments[index]);
    }
  }

  // If h' is of the form lambda x . u
  if (is_abstraction(head_rewritten))
  {
    const auto& abstraction = static_cast<const data::abstraction&>(head_rewritten);

    // sigma := sigma[x gets u']
    for (auto& variable : abstraction.variables())
    {
      assert(false);
      sigma[variable] = arguments[0];
    }

    data_expression result = rewrite_impl(abstraction.body(), sigma);

    if (PrintRewriteSteps)
    {
      mCRL2log(info) << "Applied beta-reduction to" << appl << " resulting in " << result << "\n";
    }

    return result;
  }
  else
  {
    application appl(head_rewritten, arguments.begin(), arguments.end());

    // (R, sigma') := match(h'(u_1', ..., u_n')),
    data_expression rhs;

    // If R not empty, this match function already applies the substitution and rewrite steps.
    if (!is_normal_form(appl) && match(appl, rhs))
    {
      // Return rewrite(r^sigma', sigma)
      return rewrite_impl(rhs, sigma);
    }
    else
    {
      // Return h'(u_1', ..., u_n')
      return static_cast<data_expression>(appl);
    }
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
data_expression InnermostRewriter::apply_substitution(const data_expression& term, const Substitution& sigma, const ConstructionStack& stack)
{
  if (EnableConstructionStack)
  {
    return stack.construct_term(sigma, m_argument_stack);
  }
  else
  {
    return capture_avoiding_substitution(term, sigma);
  }
}

bool InnermostRewriter::match(const data_expression& term, data_expression& rhs)
{
  // Searches for a left-hand side and a substitution such that when the substitution is applied to this left-hand side it is (syntactically) equivalent
  // to the given term.
  for (const auto& tuple : m_rewrite_system[get_nested_head(term)])
  {
    const auto& equation = std::get<0>(tuple);
    const auto& condition_stack = std::get<1>(tuple);
    const auto& rhs_stack = std::get<2>(tuple);

    // Compute a matching substitution for each rule and check that the condition associated with that rule is true, either trivially or by rewrite(c^sigma, identity).
    m_matching_sigma.clear();
    if (match_lhs(term, equation.lhs(), m_matching_sigma))
    {
      if(PrintMatchSteps)
      {
        mCRL2log(info) << "Matched rule " << equation << " to term " << term << "\n";
      }

      rhs = apply_substitution(equation.rhs(), m_matching_sigma, rhs_stack);

      // Delaying rewriting the conditions ensures that the matching substitution does not have to be saved.
      if (equation.condition() != sort_bool::true_() &&
        rewrite_impl(apply_substitution(equation.condition(), m_matching_sigma, condition_stack), m_identity) != sort_bool::true_())
      {
        continue;
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