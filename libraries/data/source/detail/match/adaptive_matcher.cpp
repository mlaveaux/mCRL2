// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/data/detail/match/adaptive_matcher.h"

#include "mcrl2/data/detail/rewrite/jitty_jittyc.h"
#include "mcrl2/utilities/stopwatch.h"
#include "mcrl2/core/index_traits.h"

#include <vector>

/// \brief Print the intermediate matches that succeeded.
constexpr bool PrintMatchSteps = true;

/// \brief Print the intermediate steps performed during construction.
constexpr bool PrintConstructionSteps = true;

using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

/// \returns A unique index for the head symbol that the given term starts with.
inline std::size_t get_head_index(const data_expression& term)
{
  return mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(get_nested_head(term)));
}

/// \brief Compute the set of positions of t such that t[p] is a variable for all p in fringe.
void fringe_impl(const atermpp::aterm_appl& appl, position current, std::set<position>& fringe)
{
  if (is_variable(appl))
  {
    // The current position has a variable.
    fringe.insert(current);
  }
  else
  {
    // Extend the position to be one deeper into the subterm.
    current.emplace_back(0);
    for (const atermpp::aterm& argument : appl)
    {
      fringe_impl(static_cast<const atermpp::aterm_appl&>(argument), current, fringe);
      ++current.back();
    }
  }

}

/// \brief Compute the set of positions of t such that t[p] is a variable for all p in fringe(t).
std::set<position> fringe(const atermpp::aterm_appl& appl)
{
  std::set<position> result;
  fringe_impl(appl, position(), result);
  return result;
}

/// \brief Decides whether the left and right terms unify, assuming that vars(left) and vars(right) are disjoint.
bool unify(const atermpp::aterm_appl& left, const atermpp::aterm_appl& right)
{
  if (left == right)
  {
    // If both sides are equivalent, then they match under any substitution. If one contains
    // a variable then that variable is also not bound in the other term.
    return true;
  }
  else if (is_variable(left) || is_variable(right))
  {
    return true;
  }
  else if (is_function_symbol(left) && is_function_symbol(right))
  {
    return left == right;
  }
  else
  {
    // The term and lhs are applications.
    const application& lhs_appl  = static_cast<const application&>(right);
    const application& term_appl = static_cast<const application&>(left);

    // Both must have the same arity, the head symbol must match and their arguments must match.
    if (lhs_appl.size() != term_appl.size())
    {
      return false;
    }

    if (!unify(term_appl.head(), lhs_appl.head()))
    {
      return false;
    }

    for (std::size_t i = 0; i < term_appl.size(); i++)
    {
      if (!unify(term_appl[i], lhs_appl[i]))
      {
        return false;
      }
    }

    return true;
  }
}

/// \brief Returns t[pos], i.e., the term at the given position using a index to keep track of the pos.
data_expression at_position_impl(const application& appl, const position& pos, std::size_t index)
{
  if ((index + 1) < pos.size())
  {
    assert(index < appl.size());
    return at_position_impl(static_cast<const application&>(appl[index]), pos, index + 1);
  }
  else
  {
    return appl[index];
  }
}

/// \brief Returns t[pos], i.e., the term at the given position.
data_expression at_position(const application& appl, const position& pos)
{
  return at_position_impl(appl, pos, 0);
}

bool compare(const position& left, const position& right)
{
  // Every element of left and right up to the length of smallest position should be less or equal.
  for (std::size_t index = 0; index <= std::min(left.size(), right.size()); ++index)
  {
    if (left[index] > right[index])
    {
      return false;
    }
  }

  return true;
}

/// \brief Replace the position variable at the given position by the expression c.
data_expression assign_at_position(const data_expression& term, const position& pos, const data_expression& c)
{
  mutable_indexed_substitution<variable, data_expression> sigma;
  sigma[position_variable(pos)] = c;
  return replace_variables(term, sigma);
}

// Public functions

template<typename Substitution>
AdaptiveMatcher<Substitution>::AdaptiveMatcher(const data_equation_vector& equations)
  : m_automaton()
{
  mcrl2::utilities::stopwatch construction;
  enumerator_identifier_generator generator("@");

  // Preprocess the term rewrite system.
  for (auto& old_equation : equations)
  {
    // Rename the variables in the equation
    std::pair<data_equation, partition> result = rename_variables_unique(make_linear(old_equation, generator));

    // Add the index of the equation
    m_linear_equations.push_back(std::make_pair(data_equation_extended(result.first), result.second));
  }

  // Determine the index of not_equal.
  m_not_equal_index = mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(not_equal()));

  // Construct the automaton.
  construct_apma(m_automaton.root(), position_variable(position()));

  mCRL2log(info) << "Adaptive Pattern Matching Automaton (states: " << m_automaton.states() << ", transitions: " << m_automaton.transitions() << ") construction took " << construction.time() << " milliseconds.\n";
}

template<typename Substitution>
void AdaptiveMatcher<Substitution>::match(const data_expression& term)
{
  m_matching_sigma.clear();
  m_match_set = nullptr;
  m_match_index = 0;

  // Start with the root state.
  std::size_t current_state = m_automaton.root();

  while (true)
  {
    // Retrieve the subterm.
    const apma_state& state = m_automaton.label(current_state);
    data_expression subterm = m_subterms[state.position];

    if (PrintMatchSteps)
    {
      mCRL2log(info) << "Matching subterm " << subterm << "\n";
    }

    // sigma := sigma \cup { (x, a) | x in L(s0) }
    if (m_automaton.label(current_state).variable.defined())
    {
      m_matching_sigma[m_automaton.label(current_state).variable] = subterm;
    }

    bool found_transition = false;
    if (is_application(subterm))
    {
      const auto& appl = static_cast<const application&>(subterm);

      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto result = m_automaton.transition(current_state, get_head_index(appl.head()));
      if (result.second)
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << current_state << " to " << result.first << " with label " << appl.head() << "\n";
        }

        found_transition = true;
        current_state = result.first;

        // Insert the subterms onto the position mapping.
        for (const atermpp::aterm& argument : appl)
        {

        }
      }
    }
    else if (is_function_symbol(subterm))
    {
      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto result = m_automaton.transition(current_state, get_head_index(subterm));
      if (result.second)
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition " << current_state << " to " << result.first << " with label " << static_cast<function_symbol>(subterm) << "\n";
        }

        found_transition = true;
        current_state = result.first;
      }
    }

    if (!found_transition)
    {
      // If delta(s0, omega) is defined then
      auto result = m_automaton.transition(current_state, m_not_equal_index);
      if (result.second)
      {
        // PMAMatch(delta(s0, omega)), t', sigma)
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << current_state << " to " << result.first << " with label omega.\n";
        }

        current_state = result.first;
      }
      // else return (emptyset, emptyset).
      else if (PrintMatchSteps)
      {
        mCRL2log(info) << "Matching failed, deadlock.\n";
        return;
      }
    }
  }

  if (PrintMatchSteps)
  {
    mCRL2log(info) << "Matching succeeded.\n";
  }

  m_match_set = &m_automaton.label(current_state).match_set;
}

template<typename Substitution>
const data_equation_extended* AdaptiveMatcher<Substitution>::next(Substitution& matching_sigma)
{
  if (m_match_set != nullptr)
  {
    while (m_match_index < m_match_set->size())
    {
      matching_sigma = m_matching_sigma;
      std::reference_wrapper<const linear_data_equation>& result = (*m_match_set)[m_match_index];
      ++m_match_index;

      if (!is_consistent(result.get().second, matching_sigma))
      {
        // This rule matched, but its variables are not consistent w.r.t. the substitution.
        continue;
      }

      return &result.get().first;
    }
  }

  return nullptr;
}

// Private functions

template<typename Substitution>
void AdaptiveMatcher<Substitution>::construct_apma(std::size_t s, data_expression pref)
{
  if (PrintConstructionSteps)
  {
    mCRL2log(info) << pref << "\n";
  }

  // L := {l in L | l unifies with pref}
  std::vector<std::reference_wrapper<const linear_data_equation>> L;
  for (const auto& equation : m_linear_equations)
  {
    if (unify(equation.first.equation().lhs(), pref))
    {
      L.emplace_back(equation);
    }
  }

  // F := fringe(pref) \ intersection_{l in L} fringe(l)
  std::set<position> F = fringe(pref);

  // if F = emptyset
  if (F.empty())
  {
    // M := M[L := L[s -> L]
    m_automaton.label(s).match_set.insert(m_automaton.label(s).match_set.begin(), L.begin(), L.end());
    return; // return M
  }
  else
  {
    // pos := select(F). For now the left-most position (left-to-right automaton).
    position pos = *std::min_element(F.begin(), F.end(), compare);

    // M := M[L := L[s -> pos]]
    m_automaton.label(s).position = m_positions.insert(pos).first;

    // for f in F s.t. exist l in L : head(l[pos]) = f
    std::set<std::pair<mcrl2::data::function_symbol, std::size_t>> symbols;
    for (const auto& equation : L)
    {
      // t := l[pos]
      data_expression t = at_position(static_cast<const application&>(equation.get().first.equation().lhs()), pos);

      if (!is_variable(t))
      {
        // head(l[pos]) in F.
        symbols.insert(std::make_pair(static_cast<const function_symbol&>(get_nested_head(t)), t.size() - 1));
      }
    }

    for (const auto& symbol : symbols)
    {
      // M := M[S := (S cup {s'})], where s' is a fresh state.
      std::size_t s_prime = m_automaton.add_state();

      // M := M[delta := (delta := delta(s, f) -> s')
      m_automaton.add_transition(s, mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(symbol.first), s_prime);

      // pref := f(omega^ar(f)), but we store in the correct position variables directly in the prefix.
      std::vector<data_expression> arguments;

      position var_position = pos;
      var_position.push_back(0);
      for (std::size_t index = 0; index < symbol.second; ++index)
      {
        var_position.back() = index;
        arguments.push_back(position_variable(var_position));
      }

      application appl(symbol.first, arguments.begin(), arguments.end());

      construct_apma(s_prime, assign_at_position(pref, pos, appl));
    }
  }
}

// Explicit instantiations.

template class mcrl2::data::detail::AdaptiveMatcher<mutable_indexed_substitution<>>;
