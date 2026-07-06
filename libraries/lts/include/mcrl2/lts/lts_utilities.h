// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

/** \file
 *
 * \brief This file contains some utility functions to manipulate lts's
 * \details The function in this file typically allow to present the
 *          transitions of a state space differently.
 * \author Jan Friso Groote
 */

#ifndef MCRL2_LTS_LTS_UTILITIES_H
#define MCRL2_LTS_LTS_UTILITIES_H

#include "mcrl2/lts/lts_lts.h"

namespace mcrl2::lts
{

/** \brief Sorts the transitions using a sort style.
 * \param[in,out] transitions A vector of transitions to be sorted. 
 * \param[in] hidden_label_set A set that tells which actions are to be interpreted as being hidden.
 *            Sorting takes place after mapping hidden actions to zero.
 * \param[in] ts The sort style to use.
 */

inline void sort_transitions(std::vector<transition>& transitions, 
                      const std::set<transition::size_type>& hidden_label_set,
                      transition_sort_style ts = src_lbl_tgt)
{
  switch (ts)
  {
    case lbl_tgt_src:
    {
      const detail::compare_transitions_lts compare(hidden_label_set);
      sort(transitions.begin(),transitions.end(),compare);
      break;
    }
    case tgt_src_lbl:
    {
      const detail::compare_transitions_tsl compare(hidden_label_set);
      sort(transitions.begin(),transitions.end(),compare);
      break;
    }
    case tgt_lbl_src:
    {
      const detail::compare_transitions_tls compare(hidden_label_set);
      sort(transitions.begin(),transitions.end(),compare);
      break;
    }
    case tgt_lbl:
    {
      const detail::compare_transitions_tl compare(hidden_label_set);
      sort(transitions.begin(),transitions.end(),compare);
      break;
    }
    case target:
    {
      const detail::compare_transitions_target compare;
      sort(transitions.begin(),transitions.end(),compare);
      break;
    }
    case src_lbl_tgt:
    default:
    {
      const detail::compare_transitions_slt compare(hidden_label_set);
      sort(transitions.begin(),transitions.end(),compare);
      break;
    }
  }
}

/** \brief Sorts the transitions on labels. Action with the tau label are put first. 
 * \details This function is linear in time, but it uses that the labels have a limited range, 
 *          i.e., it is also linear in time and space in the largest used index of a label. 
 *          This function does not employ the hidden label set. 
 * \param[in,out] transitions A vector of transitions to be sorted. 
 **/
inline std::string ptr(const transition t)
    {
      return std::to_string(t.from()) + " -" + std::to_string(t.label()) + "-> " + std::to_string(t.to());
    } 

static thread_local std::vector<std::size_t> bogus_todo_stack;

template <bool USE_STACK>
inline void group_transitions_on_label(const std::vector<transition>::iterator begin,
                                       const std::vector<transition>::iterator end,
                                       std::function<std::size_t(const transition&)> get_label, 
                                       std::vector<std::pair<std::size_t, std::size_t>>& count_sum_transitions_per_action,
                                       const std::size_t tau_label_index=0,
                                       std::vector<std::size_t>& todo_stack=bogus_todo_stack)
{
#ifndef NDEBUG
  for(std::pair<std::size_t, std::size_t> p: count_sum_transitions_per_action){ assert(p.first==0); }
  assert(todo_stack.empty());
#endif

  if (USE_STACK)
  {
    todo_stack.push_back(tau_label_index);
  }
  for(std::vector<transition>::iterator ti=begin; ti!=end; ++ti)
  {
    const transition& t=*ti;
    std::size_t label=get_label(t);
    if (USE_STACK && label!=tau_label_index && count_sum_transitions_per_action[label].first==0)
    {
      todo_stack.push_back(label);
    }
    count_sum_transitions_per_action[label].first++;
  }

  // Sum the number of transitions per label. By summing these numbers up, count_sum_transitions_per_action[a].second indicates the starting
  // position where a transition with label a must be placed. 

  std::size_t sum=0;
  if (USE_STACK) 
  {
    for(std::size_t a: todo_stack)
    {
      sum=sum+count_sum_transitions_per_action[a].first;
      count_sum_transitions_per_action[a].second=sum;
    }
  }
  else
  {
    for(std::pair<std::size_t, std::size_t>& p: count_sum_transitions_per_action)
    {
      sum=sum+p.first;
      p.second=sum;
    }
  }

  // Move from left to right through the transitions and move them to the required place. 
  std::vector<transition>::iterator current_leftmost_position=begin;
  std::vector<std::size_t>::iterator current_leftmost_label=todo_stack.begin();
  std::size_t current_label=USE_STACK?*current_leftmost_label:0;
  while (current_leftmost_position!=end)
  {
    if (USE_STACK)
    {
      assert(current_leftmost_label!=todo_stack.end());
      while (count_sum_transitions_per_action[current_label].first==0)
      {
        current_leftmost_label++;
        if (current_leftmost_label==todo_stack.end())
        {
          break;
        }
        current_label=*current_leftmost_label;
      }
      // Move to the first position with the current label that is potentially not correct.
      if (current_leftmost_label==todo_stack.end())
      {
        break;
      }
    }
    else
    {
      while (current_label<count_sum_transitions_per_action.size() &&
             count_sum_transitions_per_action[current_label].first==0)
      {
        current_label++;
      }
      // Move to the first position with the current label that is potentially not correct.
      if (current_label==count_sum_transitions_per_action.size())
      {
        break;
      }
    }
    current_leftmost_position=begin+
                              static_cast<std::ptrdiff_t>(count_sum_transitions_per_action[current_label].second-
                              count_sum_transitions_per_action[current_label].first);
    
// std::cerr << "current_leftmost_position: " << &*current_leftmost_position << "\n";
// std::cerr << "BBBBB==============\n"; for(auto t=begin; t!=end; ++t){ std::cerr << ptr(*t) << "\n";} std::cerr << "---------\n";
    // Shift the current leftmost position until a transition is found that must be moved.
    while (current_leftmost_position!=end && get_label(*current_leftmost_position)==current_label)
    {
      current_leftmost_position++;
      count_sum_transitions_per_action[current_label].first--;
    }
    
    // Check whether a transition with the current label is expected at the current_leftmost_position.
    if (count_sum_transitions_per_action[current_label].first>0)
    {
      // At the current position there is a transition that needs to be moved. Move transitions out of the way, until
      // the transitions that belongs here is found and is moved in place. 
      transition transition_on_the_move= *current_leftmost_position;
      do
      {
        std::size_t label=get_label(transition_on_the_move);
        std::size_t new_position=count_sum_transitions_per_action[label].second-count_sum_transitions_per_action[label].first;
        assert(0<std::distance(begin,end)-new_position);
        count_sum_transitions_per_action[label].first--;
        // Search for a target position with a non-matching label. 
        while (get_label(*(begin+static_cast<std::ptrdiff_t>(new_position)))==label)
        {
          new_position++;
          count_sum_transitions_per_action[label].first--;
        }
        assert(get_label(transition_on_the_move)!=get_label(*(begin+new_position)));
        std::swap(transition_on_the_move,*(begin+static_cast<std::ptrdiff_t>(new_position)));
      }
      while (get_label(transition_on_the_move)!=current_label);
      // We found the transition that we must put at the current place. 
      *current_leftmost_position=transition_on_the_move;
      current_leftmost_position++;
      count_sum_transitions_per_action[current_label].first--;
    }
  }
}

inline void group_transitions_on_label(std::vector<transition>& transitions, 
            std::function<std::size_t(const transition&)> get_label, 
            const std::size_t number_of_labels, 
            const std::size_t tau_label_index)
{
  std::vector<std::pair<std::size_t, std::size_t>> count_sum_transitions_per_action(number_of_labels, {0,0});
  std::vector<std::size_t> todo_stack;
  group_transitions_on_label<true>(transitions.begin(), transitions.end(), get_label, 
                                   count_sum_transitions_per_action, tau_label_index, todo_stack);
}
      



inline void group_transitions_on_label_tgt(std::vector<transition>& transitions, const std::size_t number_of_labels, const std::size_t tau_label_index,
                                          const std::size_t number_of_states)
{
  std::vector<std::pair<std::size_t, std::size_t>> count_sum_transitions_per_action(number_of_labels, {0,0});
  std::vector<std::size_t> todo_stack;

  group_transitions_on_label<true>(transitions.begin(), 
                                   transitions.end(),
                                   [](const transition& t){ return t.label(); }, 
                                   count_sum_transitions_per_action, 
                                   tau_label_index, 
                                   todo_stack);
  std::vector<std::size_t> todo_stack_target;
  std::vector<std::pair<std::size_t, std::size_t>> value_sum_counter(number_of_states, {0,0});
  std::vector<transition>::iterator begin=transitions.begin();
  for(std::size_t label: todo_stack)
  {
    todo_stack_target.clear();
    std::vector<transition>::iterator end=transitions.begin()+static_cast<std::ptrdiff_t>(count_sum_transitions_per_action[label].second);
    group_transitions_on_label<true>(begin, end, [](const std::vector<transition>::value_type& t){ return t.to(); }, value_sum_counter, 0, todo_stack_target); 
    begin=end;
  }
} 


inline void group_transitions_on_tgt_label(std::vector<transition>& transitions, 
                                           const std::size_t number_of_labels, 
                                           const std::size_t tau_label_index,
                                           const std::size_t number_of_states)
{
  constexpr std::size_t log2cache_size = 16; // This is an approximation of the 2log of the cache size in pointers. 
  constexpr std::size_t one=1;
  assert(number_of_states>0);
  assert(number_of_states< (one << 2*log2cache_size));
  size_t relevant_bits=(number_of_states<2?1:static_cast<std::size_t>(std::log2(number_of_states-1)+1));
  assert((one<<relevant_bits)>=number_of_states);
  size_t mask = 0;
  size_t shift = 0;
  if (relevant_bits>log2cache_size)
  {
    // Sort in two phases, with the relevant bits evenly divided. 
    mask=(one<<(relevant_bits/2))-1;
    shift=relevant_bits-(relevant_bits/2);
  } 
  else
  {
    mask=(one<<relevant_bits)-1;
    shift=0;
  }
  relevant_bits=shift;

  std::vector<std::pair<std::size_t, std::size_t>> count_sum_transitions_per_tgt1(mask+1, {0,0});
  std::vector<std::size_t> todo_stack;
  
  group_transitions_on_label<false>(transitions.begin(),
                                    transitions.end(),
                                    [shift](const transition& t){ return (t.to()>>shift); },
                                    count_sum_transitions_per_tgt1);
  
  if (relevant_bits>0)
  {
    mask=(one<<shift)-1;
    std::vector<std::pair<std::size_t, std::size_t>> count_sum_transitions_per_tgt2(mask+1, {0,0});
    std::vector<std::size_t> todo_stack_target;
    std::vector<std::pair<std::size_t, std::size_t>> value_sum_counter(number_of_labels, {0,0});

    std::vector<transition>::iterator begin=transitions.begin();
    for(std::pair<std::size_t, std::size_t>& p1: count_sum_transitions_per_tgt1)
    {
      std::vector<transition>::iterator end=transitions.begin()+static_cast<std::ptrdiff_t>(p1.second);
      group_transitions_on_label<false>(begin, 
                                        end, 
                                        [mask](const transition& t){ return (t.to() & mask); }, 
                                        count_sum_transitions_per_tgt2);
      std::vector<transition>::iterator begina=begin;
      for(std::pair<std::size_t, std::size_t>& p2: count_sum_transitions_per_tgt2)
      {

        std::vector<transition>::iterator enda=begin+static_cast<std::ptrdiff_t>(p2.second);
        group_transitions_on_label<true>(begina, 
                                         enda, 
                                         [](const std::vector<transition>::value_type& t){ return t.label(); }, 
                                         value_sum_counter, 
                                         tau_label_index, 
                                         todo_stack_target);
        for(std::size_t i:todo_stack)
        {
          if (value_sum_counter[i].first != 0)
          {
            throw mcrl2::runtime_error("BLAH " + std::to_string(value_sum_counter[i].first));
          };
          value_sum_counter[i].second=0;
        }
        todo_stack_target.clear();
        begina=enda;
        
      }
      begin=end;
    }
  }
  else 
  {
    std::vector<std::size_t> todo_stack_target;
    std::vector<std::pair<std::size_t, std::size_t>> value_sum_counter(number_of_labels, {0,0});
    std::vector<transition>::iterator begin=transitions.begin();
    for(std::pair<std::size_t, std::size_t>& p: count_sum_transitions_per_tgt1)
    {
      std::vector<transition>::iterator end=transitions.begin()+static_cast<std::ptrdiff_t>(p.second);
     group_transitions_on_label<true>(begin, 
                                      end, 
                                      [](const transition& t){ return t.label(); }, 
                                      value_sum_counter, 
                                      0, 
                                      todo_stack_target);
      for(std::size_t i:todo_stack)
      {
        value_sum_counter[i]={0,0};
      }
      todo_stack_target.clear(); 
      begin=end;
    }
  }
#ifndef NDEBUG
  std::size_t last_state=0;
  std::size_t last_label=tau_label_index;
  std::unordered_set<std::size_t> seen_labels;
  for(const transition& t: transitions)
  {
    if (last_state<t.to())
    {
      last_label=tau_label_index;
      last_state=t.to();
      seen_labels.clear();
    }
    else
    {
      assert(last_state==t.to());
      if (last_label!=t.label())
      {
        seen_labels.insert(last_label);
        assert(seen_labels.count(t.label())==0);
        last_label=t.label();
      }
    }
  }
#endif
} 

template <class LTS_TYPE>
inline void group_transitions_on_tgt_label(LTS_TYPE& aut)
{
  group_transitions_on_tgt_label(aut.get_transitions(), aut.num_action_labels(), aut.tau_label_index(), aut.num_states());
}

/** \brief Sorts the transitions using a sort style.
 * \param[in,out] transitions A vector of transitions to be sorted. 
 * \param[in] ts The sort style to use.
 */

inline void sort_transitions(std::vector<transition>& transitions, transition_sort_style ts = src_lbl_tgt)
{
  sort_transitions(transitions, std::set<transition::size_type>(), ts);
}


namespace detail
{

// An indexed sorted vector below contains the outgoing or incoming transitions per state,
// grouped per state. The input consists of a vector of transitions. The incoming/outcoming
// tau transitions are grouped by state in the m_states_with_outgoing_or_incoming_transition. 
// It is as long as the lts aut has transitions. The vector m_indices is as long as the number
// of states plus 1. For each state it contains the place in the other vector where its tau transitions
// start. So, the tau transitions reside at position indices[s] to indices[s+1]. These indices
// can be acquired using the functions lowerbound and upperbound. 
// This data structure is chosen due to its minimal memory and time footprint. 
template <class CONTENT>
class indexed_sorted_vector_for_transitions
{
  protected:
    using state_type = std::size_t;
    using label_type = std::size_t;
    using label_state_pair = std::pair<label_type, state_type>;

    std::vector < CONTENT > m_states_with_outgoing_or_incoming_transition;
    std::vector <size_t> m_indices;

  public:

    indexed_sorted_vector_for_transitions(const std::vector < transition >& transitions , state_type num_states, bool outgoing)
     : m_indices(num_states+1,0)
    {
      // First count the number of outgoing transitions per state and put it in indices.
      for(const transition& t: transitions)
      {
        m_indices[outgoing?t.from():t.to()]++;
      }

      // Calculate the m_indices where the states with outgoing/incoming tau transition must be placed.
      // Put the starting index for state i at position i-1. When placing the transitions these indices
      // are decremented properly. 
      
      size_t sum=0;
      for(state_type& i: m_indices)  // The vector is changed. This must be a reference. 
      {
        sum=sum+i;
        i=sum;
      }

      // Now declare enough space for all transitions and store them in reverse order, while
      // at the same time decrementing the indices in m_indices. 
      m_states_with_outgoing_or_incoming_transition.resize(sum);
      for(const transition& t: transitions)
      {
        if (outgoing)
        {
          assert(t.from()<m_indices.size());
          assert(m_indices[t.from()]>0);
          m_indices[t.from()]--;
          assert(m_indices[t.from()] < m_states_with_outgoing_or_incoming_transition.size());
          m_states_with_outgoing_or_incoming_transition[m_indices[t.from()]]=label_state_pair(t.label(), t.to());
        }
        else
        {
          assert(t.to()<m_indices.size());
          assert(m_indices[t.to()]>0);
          m_indices[t.to()]--;
          assert(m_indices[t.to()] < m_states_with_outgoing_or_incoming_transition.size());
          m_states_with_outgoing_or_incoming_transition[m_indices[t.to()]]=label_state_pair(t.label(), t.from());
        }
      }
      assert(m_indices.at(num_states)==m_states_with_outgoing_or_incoming_transition.size());
    }

    // Get the indexed transitions. 
    const std::vector<CONTENT>& get_transitions() const
    {
      return m_states_with_outgoing_or_incoming_transition;
    }
  
    // Get the lowest index of incoming/outging transitions stored in m_states_with_outgoing_or_incoming_transition.
    size_t lowerbound(const state_type s) const
    {
      assert(s+1<m_indices.size());
      return m_indices[s];
    }

    // Get 1 beyond the higest index of incoming/outging transitions stored in m_states_with_outgoing_or_incoming_transition.
    size_t upperbound(const state_type s) const
    {
      assert(s+1<m_indices.size());
      return m_indices[s+1];
    }

    // Drastically clear the vectors by resetting its memory usage to minimal. 
    void clear()   
    {
      std::vector <state_type>().swap(m_states_with_outgoing_or_incoming_transition);
      std::vector <size_t>().swap(m_indices);
      
    }
};

} // end namespace detail 

//
/// \brief Type for exploring transitions per state.
using outgoing_pair_t = std::pair<transition::size_type, transition::size_type>;

using outgoing_transitions_per_state_t = detail::indexed_sorted_vector_for_transitions<outgoing_pair_t>;

/// \brief Label of a pair of a label and target state. 
inline std::size_t label(const outgoing_pair_t& p)
{
  return p.first;
}

/// \brief Target state of a label state pair. 
inline std::size_t to(const outgoing_pair_t& p)
{
  return p.second;
}

/// \brief Type for exploring transitions per state and action.
// It can be considered to replace this function with an unordered_multimap.
// This may increase memory requirements, but would allow for constant versus logarithmic access times
// of elements.
using outgoing_transitions_per_state_action_t
    = std::multimap<std::pair<transition::size_type, transition::size_type>, transition::size_type>;

/// \brief From state of an iterator exploring transitions per outgoing state and action.
inline std::size_t from(const outgoing_transitions_per_state_action_t::const_iterator& i)
{
  return i->first.first;
}

/// \brief Label of an iterator exploring transitions per outgoing state and action.
inline std::size_t label(const outgoing_transitions_per_state_action_t::const_iterator& i)
{
  return i->first.second;
}

/// \brief To state of an iterator exploring transitions per outgoing state and action.
inline std::size_t to(const outgoing_transitions_per_state_action_t::const_iterator& i)
{
  return i->second;
}

/// \brief Provide the transitions as a multimap accessible per from state and label.
inline outgoing_transitions_per_state_action_t transitions_per_outgoing_state_action_pair(const std::vector<transition>& trans)
{
  outgoing_transitions_per_state_action_t result;
  for (const transition& t: trans)
  {
    result.insert(std::pair<std::pair<transition::size_type, transition::size_type>, transition::size_type>(
                    std::pair<transition::size_type, transition::size_type>(t.from(), t.label()), t.to()));
  }
  return result;
}

/// \brief Provide the transitions as a multimap accessible per from state and label.
inline outgoing_transitions_per_state_action_t transitions_per_outgoing_state_action_pair(
                 const std::vector<transition>& trans, 
                 const std::set<transition::size_type>& hide_label_set)
{
  outgoing_transitions_per_state_action_t result;
  for (const transition& t: trans)
  {
    result.insert(std::pair<std::pair<transition::size_type, transition::size_type>, transition::size_type>(
                    std::pair<transition::size_type, transition::size_type>(t.from(), detail::apply_hidden_labels(t.label(),hide_label_set)), t.to()));
  }
  return result;
} 

/// \brief Provide the transitions as a multimap accessible per from state and label, ordered backwardly.
inline outgoing_transitions_per_state_action_t transitions_per_outgoing_state_action_pair_reversed(const std::vector<transition>& trans)
{
  outgoing_transitions_per_state_action_t result;
  for (const transition& t: trans)
  {
    result.insert(std::pair<std::pair<transition::size_type, transition::size_type>, transition::size_type>(
                    std::pair<transition::size_type, transition::size_type>(t.to(), t.label()), t.from()));
  }
  return result;
}

/// \brief Provide the transitions as a multimap accessible per from state and label, ordered backwardly.
inline outgoing_transitions_per_state_action_t transitions_per_outgoing_state_action_pair_reversed(
                const std::vector<transition>& trans,
                const std::set<transition::size_type>& hide_label_set)
{
  outgoing_transitions_per_state_action_t result;
  for (const transition& t: trans)
  {
    result.insert(std::pair<std::pair<transition::size_type, transition::size_type>, transition::size_type>(
                    std::pair<transition::size_type, transition::size_type>(t.to(), detail::apply_hidden_labels(t.label(),hide_label_set)), t.from()));
  }
  return result;
} 

namespace detail
{
// Yields a label with an obscure name referring to divergence.

template < class LABEL_TYPE >
// LABEL_TYPE make_divergence_label(const std::string& s, const LABEL_TYPE& l)
LABEL_TYPE make_divergence_label(const std::string& s)
{
  return LABEL_TYPE(s);
}

template <>
inline mcrl2::lts::action_label_lts make_divergence_label<mcrl2::lts::action_label_lts>(const std::string& s)
{
  return action_label_lts(lps::multi_action(process::action(process::action_label(core::identifier_string(s),
                                                                                  data::sort_expression_list()),
                                                            data::data_expression_list())));
}

// Make a new divergent_transition_label and replace each self loop with it.
// Return the number of the divergent transition label.
template < class LTS_TYPE >
std::size_t mark_explicit_divergence_transitions(LTS_TYPE& l)
{
  // Below we create an odd action label, representing divergence.
  const typename LTS_TYPE::action_label_t lab=make_divergence_label<typename LTS_TYPE::action_label_t>("!@*&divergence&*@!"); 
  std::size_t divergent_transition_label=l.add_action(lab);
  assert(divergent_transition_label+1==l.num_action_labels());
  for(transition& t: l.get_transitions())
  {
    if (l.is_tau(l.apply_hidden_label_map(t.label())) && t.to()==t.from())
    {
      t = transition(t.to(),divergent_transition_label,t.to());
    }
  }
  return divergent_transition_label;
}

// Replace each transition in a state that is an outgoing divergent_transition with a tau_loop in that state.
template < class LTS_TYPE >
void unmark_explicit_divergence_transitions(LTS_TYPE& l, const std::size_t divergent_transition_label)
{
  for(transition& t: l.get_transitions())
  {
    if (t.label()==divergent_transition_label)
    { 
      t = transition(t.from(),l.tau_label_index(),t.from());
    }
  }
}

} // namespace detail

}

#endif // MCRL2_LTS_LTS_UTILITIES_H
