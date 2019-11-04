// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_AUTOMATON_H
#define MCRL2_DATA_DETAIL_AUTOMATON_H

#include "mcrl2/utilities/unordered_map.h"

#include <deque>

namespace mcrl2
{
namespace data
{
namespace detail
{

/// \brief An automaton that uses indices for states and transition labels, optionally the states can be labelled with additional information.
template<typename StateLabel>
class IndexedAutomaton
{
public:
  IndexedAutomaton()
  {
    add_state(); // This is the zero state used to indicate no outgoing transition.
    add_state(); // This is the root state (root() == 1)
  }

  /// \brief Adds a state to the set of states and returns its index.
  std::size_t add_state()
  {
    m_states.emplace_back();
    m_transitions.emplace_back();
    return m_states.size() - 1;
  }

  /// \brief Add a transition (from, label, to) such that (to, true) = transition(from, label).
  void add_transition(std::size_t from, std::size_t label, std::size_t to)
  {
    // Ensure that the mapping for state 'from' can index the label.
    auto& mapping = m_transitions[from];
    mapping.resize(std::max(mapping.size(), label+1));

    if (mapping[label] != to)
    {
      // Only a new transition whenever its (from, label, to) tuple is different.
      ++m_noftransitions;
    }

    mapping[label] = to;
  }

  /// \returns A pair (to, true) indicating that (from, label, to) in m_transitions or (x, false) otherwise.
  std::size_t transition(std::size_t from, const std::size_t label)
  {
    const auto& mapping = m_transitions[from];
    if (label >= mapping.size())
    {
      return 0;
    }
    else
    {
      return mapping[label];
    }
  }

  /// \brief Merge this automaton with one rooted at the given root state.
  void merge(std::size_t root, const IndexedAutomaton<StateLabel>& automaton)
  {
    // Keep track of the old number of states, new states are nof_states + index in automaton.
    std::size_t nof_states = m_states.size();

    // Take the label of the other root.
    m_states[root] = automaton.label(automaton.root());

    // Add new states for every state except for the root.
    for (std::size_t state = 0; state < automaton.states() - 1; ++state)
    {
      // Add a state.
      add_state();

      // Copy the state label.
      label(nof_states + state) = automaton.label(state + 2);
    }

    // Take transitions from automaton.root() and add them for our root.
    std::size_t label = 0;
    for (std::size_t to : automaton.m_transitions[automaton.root()])
    {
      if (to != 0)
      {
        // Transition was defined, add to this automaton.
        add_transition(root, label, nof_states + to - 2);
      }
      ++label;
    }

    // For each state of this automaton except for the root.
    for (std::size_t state = 0; state < automaton.states() - 1; ++state)
    {
      std::size_t label = 0;
      for (std::size_t to : automaton.m_transitions[state + 2])
      {
        if (to != 0)
        {
          add_transition(nof_states + state, label, nof_states + to - 2);
        }
        ++label;
      }
    }
  }

  /// \returns The state label of the given state index.
  const StateLabel& label(std::size_t state) const { return m_states[state]; }
  StateLabel& label(std::size_t state) { return m_states[state]; }

  std::size_t root() const { return 1; }

  std::size_t states() const { return m_states.size() - 1; }

  std::size_t transitions() const { return m_noftransitions; }

  /// \return The size of the automaton (number of states and transitions).
  std::size_t size() const { return states() + transitions(); }

private:
  std::deque<StateLabel> m_states;

  // A mapping from states to (label, state) pairs.
  std::vector<std::vector<std::size_t>> m_transitions;

  std::size_t m_noftransitions = 0;
};

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_AUTOMATON_H
