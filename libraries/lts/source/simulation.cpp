// Author(s): Muck van Weerdenburg, Ruud Koolen; adapted by Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include "mcrl2/data//standard_numbers_utility.h"
#include "mcrl2/lps/stochastic_state.h"
#include "mcrl2/lts/simulation.h"

using namespace mcrl2;
using namespace mcrl2::lps;

std::vector<simulation::transition_t> simulation::transitions(const state& source_state)
{
  try
  {
    std::list<transition> outgoing_transitions = out_edges(source_state);
    std::vector<simulation::transition_t> output(outgoing_transitions.begin(), outgoing_transitions.end());
    return output;
  }
  catch (mcrl2::runtime_error& e)
  {
    mCRL2log(mcrl2::log::error) << "an error occurred while calculating the transitions from this state;\n" << e.what() << std::endl;
    return std::vector<simulation::transition_t>();
  }
}

simulation::simulation(const stochastic_specification& specification, data::rewrite_strategy strategy)
  : explorer(specification, explorer_options(strategy)),
    m_specification(specification),
    m_gen(),
    m_distrib(0,std::numeric_limits<std::size_t>::max())
{
  state_type initial_state;
  compute_stochastic_state(initial_state, m_initial_distribution, m_initial_state, m_global_sigma, m_global_rewr, m_global_enumerator);
  simulator_state_t state;
  state.source_state = initial_state;
  state.state_number=initial_state.size();  // This indicates that no state is selected yet.
  m_full_trace.push_back(state);

  if (initial_state.size()==1 || m_auto_select_probabilistic_state)
  {
    randomly_select_state();  // There is only one state, or we select states at random. Select the state. 
  }
}

void simulation::truncate(std::size_t state_number, bool probabilistic)
{
  assert(state_number < m_full_trace.size());
  m_full_trace.resize(state_number + 1);
  if (probabilistic && m_full_trace[state_number].source_state.size()>1)  // Indicate that this is a probabilistic state and forget about transitions.
  {
    m_full_trace[state_number].state_number=m_full_trace[state_number].source_state.size();
  }
  
}

void simulation::select_state(std::size_t state_number)
{
  simulator_state_t& state = m_full_trace.back();

  assert(state_number<state.source_state.states.size());
  state.state_number=state_number;
  state.transitions = transitions(state.source_state.states[state_number]);
}

void simulation::environment(std::vector<std::string> values)
{
  assert(values.size() == m_specification.process().process_parameters().size());
  std::vector<data::data_expression> data_values;
  for (std::size_t i = 0; i < values.size(); ++i)
  {
    data_values.push_back(data::parse_data_expression(values[i], m_specification.data()));
  }

  simulator_state_t result_state;
  compute_stochastic_state(result_state.source_state, stochastic_distribution(), data_values, m_global_sigma, m_global_rewr, m_global_enumerator);
  result_state.transitions = transitions(result_state.source_state.states[0]);
  result_state.state_number = 0;
  result_state.transition_number = 0;

  m_full_trace.clear();
  m_full_trace.emplace_back(result_state);
}

void simulation::randomly_select_state()
{
  simulator_state_t& state = m_full_trace.back();
  std::size_t state_number=0;
  
  if (state.source_state.states.size()>1)
  {
    // Generate a random size_t with random distribution. 
    double random_value=static_cast<double>(m_distrib(m_gen))/static_cast<double>(std::numeric_limits<std::size_t>::max()); 
    random_value=random_value-data::sort_real::value<double>(state.source_state.probabilities[state_number]);
    while (random_value>0)
    {
      state_number++;
      assert(state_number<state.source_state.states.size());
      random_value=random_value-data::sort_real::value<double>(state.source_state.probabilities[state_number]);
    }
  }
  state.state_number=state_number;
  state.transitions = transitions(state.source_state.states[state_number]);
}

void simulation::select_transition(std::size_t transition_number)
{
  assert(transition_number < m_full_trace.back().transitions.size());
  m_full_trace.back().transition_number = transition_number;
  simulator_state_t state;
  state.source_state = m_full_trace.back().transitions[transition_number].state;
  state.state_number=state.source_state.size();
  m_full_trace.push_back(state);

  if (state.source_state.size()==1 || m_auto_select_probabilistic_state)
  {
    randomly_select_state();
  }
}

void simulation::randomly_select_transition()
{
  std::size_t s=m_full_trace.back().transitions.size();
  std::size_t transition_number= m_distrib(m_gen) % s;
  select_transition(transition_number);
}

void simulation::enable_auto_select_probability(bool enable)
{
  m_auto_select_probabilistic_state=enable;
}

void simulation::save(const std::string& filename) const
{
  lts::trace trace;
  if (m_full_trace[0].state_number>=m_full_trace[0].source_state.size())
  {
    throw mcrl2::runtime_error("initial state is not set. Trace is not saved.");
  }
  trace.set_state(m_full_trace[0].source_state.states[m_full_trace[0].state_number]);
  for (std::size_t i = 0; i + 1 < m_full_trace.size(); i++)
  {
    assert(m_full_trace[i].state_number<m_full_trace[i].source_state.size());
    if (m_full_trace[i].transition_number<m_full_trace[i].transitions.size())
    {
      trace.add_action(m_full_trace[i].transitions[m_full_trace[i].transition_number].action);
      if (m_full_trace[i+1].state_number<m_full_trace[i+1].source_state.size())
      {
        trace.set_state(m_full_trace[i+1].source_state.states[m_full_trace[i+1].state_number]);
      }
      else
      {
        trace.set_state(m_full_trace[i+1].source_state.states[0]);
        std::cout << "the last state in the saved trace is the first concrete state from the last probabilistic state." << std::endl;
      }
    }
  }
  trace.save(filename);
}

void simulation::load(const std::string& filename)
{
  // Load the trace from file
  lts::trace trace(filename, m_specification.data(), m_specification.action_labels());
  trace.reset_position();

  // Get the first state from the generator
  m_full_trace.clear();

  state_type initial_state;
  compute_stochastic_state(initial_state, m_initial_distribution, m_initial_state, m_global_sigma, m_global_rewr, m_global_enumerator);

  add_new_state(initial_state);
  // Check that the first state (if given) matches one of the probabilistic states of the specification.
  if (trace.current_state_exists()) 
  {    
    if (std::find(initial_state.states.begin(), initial_state.states.end(), trace.current_state()) == initial_state.states.end())
    {
      throw mcrl2::runtime_error("The initial state of the trace does not match the initial state "
                                 "of this specification");
    }


    // Replay the trace using the generator.
    if (!match_trace_probabilistic_state(trace))
    {
      throw mcrl2::runtime_error("Failed to perform action " + pp(trace.current_action()) + " at position "
                                 + std::to_string(m_full_trace.size() - 1) + " from the trace loaded from " + filename + ".");
    }
  }
}

// Add a new state to m_full_trace with the indicated state. 
void simulation::add_new_state(const lps::stochastic_state& s)
{
  simulator_state_t state;
  state.source_state = s;
  state.state_number = s.size();
  state.transition_number = 0;
  m_full_trace.push_back(state);
} 

bool simulation::match_trace_probabilistic_state(lts::trace& trace)
{
  if (!trace.current_state_exists())
  {
    return true;
  } 
  assert(trace.current_state_exists());
  simulator_state_t& current = m_full_trace.back();
  for (std::size_t i = 0; i < current.source_state.size(); ++i)
  {
    if (current.source_state.states[i] == trace.current_state())
    {
      select_state(i);
      if (match_trace_transition(trace))
      {
        return true;
      }
    }
  }
  return false;
}

bool simulation::match_trace_transition(lts::trace& trace)
{
  if (!trace.current_action_exists())
  {
    return true;
  }
  simulator_state_t& current = m_full_trace.back();
  lps::multi_action action = trace.current_action();
  trace.increase_position();
  for (std::size_t i = 0; i < current.transitions.size(); ++i)
  {
    if (current.transitions[i].action == action)
    {
      add_new_state(current.transitions[i].state);
      current.transition_number = i;
      if (match_trace_probabilistic_state(trace))
      {
        return true;
      }
      m_full_trace.pop_back();
    }
  }
  return false;
}

