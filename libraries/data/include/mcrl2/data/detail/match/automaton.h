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
  IndexedAutomaton() { add_state(); }

  /// \brief Adds a state to the set of states and returns its index.
  std::size_t add_state()
  {
    m_states.emplace_back();
    return m_states.size() - 1;
  }

  /// \brief Add a transition (from, label, to) such that (to, true) = transition(from, label).
  void add_transition(std::size_t from, std::size_t label, std::size_t to)
  {
    m_transitions.insert(std::make_pair(std::make_pair(from, label), to));
  }

  /// \returns A pair (to, true) indicating that (from, label, to) in m_transitions or (x, false) otherwise.
  std::pair<std::size_t, bool> transition(std::size_t from, const std::size_t label)
  {
    auto it = m_transitions.find(std::make_pair(from, label));
    if (it == m_transitions.end())
    {
      return std::make_pair(0, false);
    }
    else
    {
      return std::make_pair(it->second, true);
    }
  }

  /// \returns The state label of the given state index.
  StateLabel& label(std::size_t state) { return m_states[state]; }

  std::size_t root() const { return 0; }

  std::size_t states() const { return m_states.size(); }

  std::size_t transitions() const { return m_transitions.size(); }

private:
  std::deque<StateLabel> m_states;

  mcrl2::utilities::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t> m_transitions;
};

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_AUTOMATON_H
