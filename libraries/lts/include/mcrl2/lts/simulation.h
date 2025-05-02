// Author(s): Jan Friso Groote; Based on a simulator by Ruud Koolen and Muck van Weerdenburg.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_LPS_SIMULATION_H
#define MCRL2_LPS_SIMULATION_H

#include "mcrl2/lps/explorer.h"
#include "mcrl2/lts/trace.h"

namespace mcrl2
{

namespace lps  // TODO: This should become lts. We should not declare objects in namespace lps within the lts library. 
{

/// \brief Simulation process.
// A simulation is effectively a trace annotated with outgoing transition information
// and an operation to extend the trace with an outgoing transition from the last state.
class simulation : protected explorer<true, false, stochastic_specification>
{
  public:

    typedef transition transition_t;

    struct simulator_state_t
    {
      lps::stochastic_state source_state;
      std::size_t state_number; // This represents the number of the selected probabilistic state, or a number out of range
                                // if not source state is chosen.
      std::vector<transition_t> transitions;
      std::size_t transition_number; // This indicates the chosen transition, or a number out of range if no transition is chosen.
    };

    /// Constructor.
    simulation(const stochastic_specification& specification, data::rewrite_strategy strategy = data::rewrite_strategy());

    /// Returns the current annotated state vector.
    const std::deque<simulator_state_t>& trace() const { return m_full_trace; }

    /// Remove states from the end of the simulation, making \a state_number the last state.
    void truncate(std::size_t state_number, bool probabilistic);

    /// Choose one state among a sequence of probabilistic states. 
    void select_state(std::size_t state_number);

    /// Randomly choose one state among a sequence of probabilistic states. 
    void randomly_select_state();

    /// Choose outgoing transition indexed by transition_number and add its state to the state vector.
    void select_transition(std::size_t transition_number);

    /// Goto the state indicated by the environment, which is given by a vector of strings.
    void environment(std::vector<std::string> values);

    /// Randomly choose an outgoing transition and add its state to the state vector.
    void randomly_select_transition();

    /// Indicate whether a probabilistic state has been selected.
    bool probabilistic_state_must_be_selected(const std::size_t state_index) const
    {
      assert(state_index<trace().size());
      const simulator_state_t& current_state = trace()[state_index];
      return current_state.state_number>=current_state.source_state.size();
    }

    /// If enabled, automatically a probabilistic state is chosen using its intrinsic probability. 
    void enable_auto_select_probability(bool enable);

    /// Save the trace to a file.
    void save(const std::string& filename) const;

    /// Load a trace from a file.
    void load(const std::string& filename);

    virtual ~simulation()
    {}

  protected:
    stochastic_specification m_specification;

    bool m_auto_select_probabilistic_state=false;
    std::mt19937 m_gen; // mersenne_twister_engine seeded with rd().
    std::uniform_int_distribution<std::size_t> m_distrib; // A random generator with a uniform distribution. 

    // The complete trace.
    std::deque<simulator_state_t> m_full_trace;

    std::vector<transition_t> transitions(const lps::state& source_state);

    bool match_trace_probabilistic_state(lts::trace& trace);
    bool match_trace_transition(lts::trace& trace);
    void add_new_state(const lps::stochastic_state& s);
    void add_transitions();

};

} // namespace lps

} // namespace mcrl2

#endif
