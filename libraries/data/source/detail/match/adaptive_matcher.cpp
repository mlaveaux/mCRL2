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
constexpr bool PrintMatchSteps = false;

/// \brief Print the intermediate steps performed during construction.
constexpr bool PrintConstructionSteps = false;

using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

/// \returns A unique index for the head symbol that the given term starts with.
inline std::size_t get_head_index(const data_expression& term)
{
  return mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(get_nested_head(term)));
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

    auto& equation = result.first;
    if (!result.second.empty())
    {
      // For now, print the now linear equations.
      mCRL2log(info) << "Renamed non-linear equation " << equation << ".\n";
    }

    // Add the index of the equation
    m_mapping.insert(std::make_pair(equation, std::make_pair(data_equation_extended(equation), result.second)));
  }

  // Determine the index of omega.
  m_not_equal_index = mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(not_equal()));

  // Construct the automaton.
  construct_apma(m_automaton.root(), position_variable(position()));

  mCRL2log(info) << "Matching automaton (states: " << m_automaton.states() << ", transitions: " << m_automaton.transitions() << ") construction took " << construction.time() << " milliseconds.\n";
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

}

// Explicit instantiations.

template class mcrl2::data::detail::AdaptiveMatcher<mutable_indexed_substitution<>>;
