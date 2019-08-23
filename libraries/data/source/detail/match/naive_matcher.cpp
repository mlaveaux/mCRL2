// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/data/detail/match/naive_matcher.h"

#include "mcrl2/data/detail/rewrite/jitty_jittyc.h"

/// \brief Print the intermediate matches that succeeded.
constexpr bool PrintMatchSteps   = false;

/// \brief When enabled, index each rewrite rule based on its head symbol for fast lookup of relevant rules.
constexpr bool EnableHeadIndexing = false;

using namespace mcrl2;
using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

/// \brief Matches a single left-hand side with the given term and creates the substitution.
template<typename Substitution>
static inline
bool match_lhs(const data_expression& term,  const data_expression& lhs, Substitution& sigma)
{
  if (term == lhs)
  {
    // If both sides are equivalent, then they match under any substitution. If one contains
    // a variable then that variable is also not bound in the other term.
    return true;
  }
  else if (is_function_symbol(lhs))
  {
    return term == lhs;
  }
  else if (is_variable(lhs))
  {
    const auto& var = static_cast<const variable&>(lhs);

    if (sigma(var) != var)
    {
      // If the variable was already assigned it must match the previously assigned value.
      return sigma(var) == term;
    }
    else
    {
      // Else substitute the current term for this variable.
      sigma[var] = term;
      return true;
    }
  }
  else
  {
    // The term and lhs are applications.
    const application& lhs_appl  = static_cast<const application&>(lhs);
    const application& term_appl = static_cast<const application&>(term);

    // Both must have the same arity, the head symbol must match and their arguments must match.
    if (lhs_appl.size() != term_appl.size())
    {
      return false;
    }

    if (!match_lhs(term_appl.head(), lhs_appl.head(), sigma))
    {
      return false;
    }

    for (std::size_t i = 0; i < term_appl.size(); i++)
    {
      if (!match_lhs(term_appl[i], lhs_appl[i], sigma))
      {
        return false;
      }
    }

    return true;
  }
}


/// \returns A unique index for the head symbol that the given term starts with.
static inline std::size_t get_head_index(const data_expression& term)
{
  return core::index_traits<data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(get_nested_head(term)));
}

NaiveMatcher::NaiveMatcher(const data_equation_vector& equations)
{
  for (auto&& equation : equations)
  {
    m_equations.emplace_back(equation,
      ConstructionStack(equation.condition()),
      ConstructionStack(equation.rhs()));
  }

  for (auto&& equation : equations)
  {
    // Make sure that it is possible to insert the match data for head_index left-hand side.
    std::size_t head_index = get_head_index(equation.lhs());
    if (head_index >= m_rewrite_system.size()) { m_rewrite_system.resize(head_index + 1); }

    // Insert the left-hand side into the rewrite rule mapping and a construction stack for its right-hand side.
    m_rewrite_system[head_index].emplace_back(equation,
      ConstructionStack(equation.condition()),
      ConstructionStack(equation.rhs()));

  }
}

std::vector<std::reference_wrapper<const data_equation_extended>> NaiveMatcher::match(const data_expression& term, mutable_indexed_substitution<>& matching_sigma)
{
  std::vector<std::reference_wrapper<const data_equation_extended>> results;

  std::size_t head_index = get_head_index(term);
  if (EnableHeadIndexing && head_index >= m_rewrite_system.size())
  {
    // No left-hand side starts with this head symbol, so it cannot match.
    return results;
  }

  // Searches for a left-hand side and a substitution such that when the substitution is applied to this left-hand side it is (syntactically) equivalent
  // to the given term. Only tries rewrite rules that start with the correct head symbol when EnableHeadIndexing is true.
  for (const auto& tuple : (EnableHeadIndexing ? m_rewrite_system[head_index] : m_equations))
  {
    const auto& equation = std::get<0>(tuple);

    // Compute a matching substitution for each rule and check that the condition associated with that rule is true, either trivially or by rewrite(c^sigma, identity).
    if (match_lhs(term, equation.lhs(), matching_sigma))
    {
      if(PrintMatchSteps)
      {
        mCRL2log(info) << "Matched rule " << equation << " to term " << term << "\n";
      }

      results.push_back(tuple);
    }
    else if (PrintMatchSteps)
    {
      mCRL2log(info) << "Tried rule " << equation << " on term " << term << "\n";
    }
  }

  if (results.empty() && PrintMatchSteps)
  {
    mCRL2log(info) << "Term " << term << " is in normal form.\n";
  }

  return results;
}
