// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts/detail/liblts_failures_refinement.h

// This file contains an implementation of 
// T. Wang, S. Song, J. Sun, Y. Liu, J.S. Dong, X. Wang and S. Li.
// More Anti-Chain Based Refinement Checking. In proceedings ICFEM 2012, editors T. Aoki and K. Tagushi,
// Lecture Notes in Computer Science no 7635, pages 364-380, 2012.
//
// There are six algorithms. One for trace inclusion, one for failures inclusion and one for failures-divergence inclusion.
// All algorithms come in a variant with and without internal steps. 
// It is possible to generate a counter transition system in case the inclusion is answered by no.

#ifndef _LIBLTS_FAILURES_REFINEMENT_H
#define _LIBLTS_FAILURES_REFINEMENT_H

#include "unordered_set"
#include "mcrl2/lts/detail/counter_example.h"
#include "mcrl2/lts/detail/exploration_strategy.h"
#include "mcrl2/lts/detail/liblts_bisim_gjkw.h"

namespace mcrl2
{
namespace lts
{
namespace detail
{
  typedef std::size_t state_type;
  typedef std::size_t label_type;
  typedef std::set<state_type> set_of_states;
  typedef std::multimap <detail::state_type,detail::set_of_states> anti_chain_type;
  typedef std::set < label_type > action_label_set;

  template < class COUNTER_EXAMPLE_CONSTRUCTOR >
  class state_states_counter_example_index_triple
  {
    protected:
      detail::state_type m_state;
      detail::set_of_states m_states;
      typename COUNTER_EXAMPLE_CONSTRUCTOR::index_type m_counter_example_index;

    public:
      state_states_counter_example_index_triple()
      {}

      /// \brief Constructor.
      state_states_counter_example_index_triple(
              const state_type state, 
              const set_of_states& states, 
              const typename COUNTER_EXAMPLE_CONSTRUCTOR::index_type& counter_example_index)
       : m_state(state),
         m_states(states),
         m_counter_example_index(counter_example_index)
      {}

      /// \brief Get the state.
      state_type state() const
      {
        return m_state;
      }

      void swap(state_states_counter_example_index_triple& other)
      {
        std::swap(m_state,other.m_state);
        std::swap(m_states,other.m_states);
        std::swap(m_counter_example_index,other.m_counter_example_index);
      }

      /// \brief Get the set of states.
      const set_of_states& states() const
      {
        return m_states;
      }

      /// \brief Get the counter example index.
      const typename COUNTER_EXAMPLE_CONSTRUCTOR::index_type& counter_example_index() const
      {
        return m_counter_example_index;
      }
  };
 
  template <  class COUNTER_EXAMPLE_CONSTRUCTOR >
  inline bool antichain_insert(
                  anti_chain_type& anti_chain, 
                  const state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>& impl_spec);

  // The class below recalls what the stable states and the states with a divergent
  // self loop of a transition system are, such that it does not have to be recalculated each time again.
  template < class LTS_TYPE >
  class lts_cache
  {
    protected:
      LTS_TYPE& m_l;
      std::vector<std::vector<state_type> > m_tau_reachable_states;
      std::vector<std::vector<transition> > m_sorted_transitions;
      std::vector<bool> m_divergent;
      std::vector<action_label_set> m_enabled_actions;
      
    private:
      void calculate_weak_property_cache(const bool weak_reduction)
      {
        scc_partitioner<LTS_TYPE> strongly_connected_component_partitioner(m_l);
        for(const transition& t: m_l.get_transitions()) 
        {
          assert(t.from()<m_l.num_states());
          if (m_l.is_tau(m_l.apply_hidden_label_map(t.label())) && weak_reduction)
          {
            m_tau_reachable_states[t.from()].push_back(t.to());  // There is an outgoing tau. 
          }
          m_sorted_transitions[t.from()].push_back(t);
         
          if (weak_reduction && m_l.is_tau(m_l.apply_hidden_label_map(t.label())) && 
              strongly_connected_component_partitioner.in_same_class(t.from(),t.to()))
          {
            m_divergent[t.from()]=true;  // There is a self loop.
          }
          m_enabled_actions[t.from()].insert(m_l.apply_hidden_label_map(t.label()));
        }
      }

    public:
      
      lts_cache(/* const */ LTS_TYPE& l, const bool weak_reduction)    /* l is not changed, but the use of the scc partitioner requires l to be non const */
        : m_l(l),
          m_tau_reachable_states(l.num_states()),
          m_sorted_transitions(l.num_states()),
          m_divergent(l.num_states(),false),
          m_enabled_actions(l.num_states())
      { 
        calculate_weak_property_cache(weak_reduction);
      }

      bool stable(const state_type s) const
      {
        return m_tau_reachable_states[s].size()==0;
      }

      const std::vector<state_type>& tau_reachable_states(const state_type s) const
      {
        return m_tau_reachable_states[s];
      }

      const std::vector<transition>& transitions(const state_type s) const
      {
        assert(s<m_sorted_transitions.size());
        return m_sorted_transitions[s];
      }

      bool diverges(const state_type s) const
      {
        return m_divergent[s];
      }

      const action_label_set& action_labels(const state_type s) const
      {
        return m_enabled_actions[s];
      }
  };

  template < class LTS_TYPE >
  set_of_states collect_reachable_states_via_taus(
                 const state_type s, 
                 const lts_cache<LTS_TYPE>& weak_property_cache,
                 const bool weak_reduction);

  template < class LTS_TYPE >
  set_of_states collect_reachable_states_via_an_action(
                 const state_type s, 
                 const label_type e,
                 const lts_cache<LTS_TYPE>& weak_property_cache,
                 const bool weak_reduction,
                 const LTS_TYPE& l);

  template < class LTS_TYPE >
  bool refusals_contained_in(
              const state_type impl, 
              const set_of_states& spec, 
              const lts_cache<LTS_TYPE>& weak_property_cache,
              label_type& culprit,
              const LTS_TYPE& l,
              const bool provide_a_counter_example);

  /* Construct a path to state s using the backward map, and return it in result */
  inline
  void reconstruct_path(const state_type s, 
                        const std::map<state_type, std::pair<label_type,state_type> >& backward_map, 
                        std::vector<label_type>& result)
  {
    if (backward_map.count(s)==0)
    {
      return;
    }
    const std::pair<label_type,state_type> p=backward_map.at(s);
    {
      reconstruct_path(p.second, backward_map, result);
      result.push_back(p.first);
    }
  }

  template < class LTS_TYPE >
  std::vector<label_type> find_path_to_stable_state_without_action_in_impl(const label_type offending_action, 
                                                                           const state_type s, 
                                                                           const lts_cache<LTS_TYPE>& lts_cache,
                                                                           const LTS_TYPE& l,
                                                                           const bool find_trace_with_taus)
  {
    set_of_states visited;
    visited.insert(s);
    std::deque<state_type> todo_stack={s};
    std::map<state_type, std::pair<label_type,state_type> > backward_map;
    while (todo_stack.size()>0)  
    {
      state_type current_state=todo_stack.front();
      todo_stack.pop_front();

      if (lts_cache.stable(current_state))  
      {
        // Check that the current state does not occur in all outgoing transitions. 
        bool found=false;
        for(const transition& t: lts_cache.transitions(current_state))
        {
          if (t.label()==offending_action)
          {
            found=true;
            break;
          }
        }
        if (!found)
        { 
          std::vector<label_type> resulting_path;
          reconstruct_path(current_state, backward_map, resulting_path);
          assert(find_trace_with_taus || resulting_path.empty()); // There are no tau's in the path if taus were not requested. 
          return resulting_path;
        }
      }
      else // The current state is not stable. 
      {
        for(const transition& t: lts_cache.transitions(current_state))
        {
          if (find_trace_with_taus && l.is_tau(l.apply_hidden_label_map(t.label())))
          {
            if (visited.insert(t.to()).second)  // The state to() was not yet explored.
            {
              todo_stack.push_back(t.to());
              backward_map[t.to()]=std::pair<label_type,state_type>(t.label(), t.from()); // Store how the to state could be reached.
            }
          }
        }
      }
    }
    // The action was not found. Do not return a path. 
    return std::vector<label_type>();
  }

} // namespace detail

enum refinement_type { trace, failures, failures_divergence };

/// \brief This function checks using algorithms in the paper mentioned above that
/// whether transition system l1 is included in transition system l2, in the
/// sense of trace inclusions, failures inclusion and divergence failures
/// inclusion.
/// \param weak_reduction If set, it will do so where tau's are included.
/// \param generate_counter_example If set, a labelled transition system is generated
///        that can act as a counterexample. It consists of a trace, followed by
///        outgoing transitions representing a refusal set.
/// \param strategy Choose between breadth and depth first.
template < class LTS_TYPE, class COUNTER_EXAMPLE_CONSTRUCTOR = detail::dummy_counter_example_constructor >
bool destructive_refinement_checker(
                        LTS_TYPE& l1, 
                        LTS_TYPE& l2, 
                        const refinement_type refinement, 
                        const bool weak_reduction,
                        const exploration_strategy strategy,
                        COUNTER_EXAMPLE_CONSTRUCTOR generate_counter_example = detail::dummy_counter_example_constructor())
{
  if (!generate_counter_example.is_dummy())  // Counter example is requested, apply bisimulation to l2.
  {
    const bool preserve_divergence=weak_reduction && (refinement!=trace);
    l2.clear_state_labels();
    if (weak_reduction)
    {
      detail::scc_partitioner<LTS_TYPE> scc_part(l2);
      scc_part.replace_transition_system(preserve_divergence);
    }

    detail::bisim_partitioner_gjkw<LTS_TYPE> bisim_part(l2,weak_reduction,preserve_divergence);
    // Assign the reduced LTS, and set init_l2.
    bisim_part.replace_transition_system(weak_reduction,preserve_divergence);
  }

  std::size_t init_l2 = l2.initial_state() + l1.num_states();
  mcrl2::lts::detail::merge(l1,l2);
  l2.clear(); // No use for l2 anymore.
  // For weak-failures and failures-divergence, the existence of tau loops make a difference.
  // Therefore, we apply bisimulation reduction preserving divergences.
  // A typical example is a.(b+c) which is not weak-failures included n a.tau*.(b+c). The lhs has failure pairs
  // <a,{a}>, <a,{}> while the rhs has only failure pairs <a,{}>, as the state after the a is not stable.
  
  if (generate_counter_example.is_dummy())  // No counter example is requested. We can use bisimulation preprocessing.
  {
    const bool preserve_divergence=weak_reduction && (refinement!=trace);
    l1.clear_state_labels(); 
    if (weak_reduction)
    {
      detail::scc_partitioner<LTS_TYPE> scc_part(l1);
      init_l2=scc_part.get_eq_class(init_l2);
      scc_part.replace_transition_system(preserve_divergence);
    }

    detail::bisim_partitioner_gjkw<LTS_TYPE> bisim_part(l1,weak_reduction,preserve_divergence);
    // Assign the reduced LTS, and set init_l2.
    init_l2=bisim_part.get_eq_class(init_l2);
    bisim_part.replace_transition_system(weak_reduction,preserve_divergence);

    if (bisim_part.in_same_class(init_l2, l1.initial_state()))
    {
      mCRL2log(log::debug) << "Both LTSs are (divergence-preserving) branching bisimular, so no need to check refinement relation.\n";
      return true;
    }
  }

  const detail::lts_cache<LTS_TYPE> weak_property_cache(l1,weak_reduction);
  std::deque< detail::state_states_counter_example_index_triple < COUNTER_EXAMPLE_CONSTRUCTOR > > 
              working(  // let working be a stack containg the triple (init1,{s|init2-->s},root_index);
                    { detail::state_states_counter_example_index_triple< COUNTER_EXAMPLE_CONSTRUCTOR >(
                                  l1.initial_state(), 
                                  detail::collect_reachable_states_via_taus(init_l2,weak_property_cache,weak_reduction),
                                  generate_counter_example.root_index() ) });
                                                      // let antichain := emptyset;
  detail::anti_chain_type anti_chain;
  detail::antichain_insert(anti_chain, working.front());   // antichain := antichain united with (impl,spec); 
                                                           // This line occurs at another place in the code than in 
                                                           // the original algorithm, where insertion in the anti-chain
                                                           // was too late, causing too many impl-spec pairs to be investigated.
  std::size_t antichain_misses = 0; // Number of times a pair was inserted into the antichain.
  std::size_t antichain_inserts = 0; // Number of times antichain_insert was called.

  while (working.size()>0)                            // while working!=empty
  {
    detail::state_states_counter_example_index_triple < COUNTER_EXAMPLE_CONSTRUCTOR > impl_spec;   // pop (impl,spec) from working;
    impl_spec.swap(working.front());  
    working.pop_front();     // At this point it could be checked whether impl_spec still exists in anti_chain. 
                             // Small scale experiments show that this is a little bit more expensive than doing the explicit check below.
    
    if (refinement==failures_divergence && weak_property_cache.diverges(impl_spec.state()))
                                                      // if impl diverges
    {
      bool spec_diverges=false;
      for(const detail::state_type s: impl_spec.states())       // if spec does not diverge
      {
        if (weak_property_cache.diverges(s))
        {
          spec_diverges=true;
          break;
        }
      }
      if (!spec_diverges)
      {
        generate_counter_example.save_counter_example(impl_spec.counter_example_index(),l1);
        return false;                                 // return false; 
      }
    }
    else 
    {
      if (refinement==failures || refinement==failures_divergence)
      { 
        detail::label_type offending_action=std::size_t(-1);
        // if refusals(impl) not contained in refusals(spec) then 
        if (!detail::refusals_contained_in(impl_spec.state(),
                                           impl_spec.states(),
                                           weak_property_cache,
                                           offending_action,
                                           l1,
                                           !generate_counter_example.is_dummy()))   
        {
          std::vector<detail::label_type> counter_example_extension;
          if (offending_action!=std::size_t(-1))
          { 
            counter_example_extension = 
                       detail::find_path_to_stable_state_without_action_in_impl(offending_action, 
                                                                                impl_spec.state(),
                                                                                weak_property_cache,
                                                                                l1, 
                                                                                failures_divergence || weak_reduction);
          }
          generate_counter_example.save_counter_example(impl_spec.counter_example_index(),l1, counter_example_extension);
          return false;                               // return false; 
        }
      }
      
      for(const transition& t: weak_property_cache.transitions(impl_spec.state()))
      {
        const typename COUNTER_EXAMPLE_CONSTRUCTOR::index_type new_counterexample_index=
               generate_counter_example.add_transition(t.label(),impl_spec.counter_example_index());
        detail::set_of_states spec_prime;
        if (l1.is_tau(l1.apply_hidden_label_map(t.label())) && weak_reduction)                   // if e=tau then
        {
          spec_prime=impl_spec.states();        // spec' := spec;
        }
        else
        {                                           // spec' := {s' | exists s in spec. s-e->s'};
          for(const detail::state_type s: impl_spec.states())  
          {
            detail::set_of_states reachable_states_from_s_via_e=
                    detail::collect_reachable_states_via_an_action(s,l1.apply_hidden_label_map(t.label()),weak_property_cache,weak_reduction,l1);
            spec_prime.insert(reachable_states_from_s_via_e.begin(),reachable_states_from_s_via_e.end());
          }
        }
        if (spec_prime.empty())                     // if spec'={} then
        {
          generate_counter_example.save_counter_example(new_counterexample_index,l1);
          return false;                             //    return false;  
        }
                                                    // if (impl',spec') in antichain is not true then
        ++antichain_inserts;
        const detail::state_states_counter_example_index_triple < COUNTER_EXAMPLE_CONSTRUCTOR > 
                          impl_spec_counterex(t.to(),spec_prime,new_counterexample_index);
        if (detail::antichain_insert(anti_chain, impl_spec_counterex))   
        {
          ++antichain_misses;
          if (strategy == exploration_strategy::es_breadth)
          {
            working.push_back(impl_spec_counterex);   // push(impl,spec') into working;
          }
          else if (strategy == exploration_strategy::es_depth)
          {
            working.push_front(impl_spec_counterex);   // push(impl,spec') into working;
          }
        }
      }
    }
    
  }

  mCRL2log(log::debug) << antichain_inserts - antichain_misses << " pairs already included in the antichain, "
                       << antichain_misses << " pairs inserted in the antichain, "
                       << anti_chain.size() << " pairs in the antichain.\n";
  return true;                                      // return true;
}


namespace detail
{
  /* This function generates the set of states reachable from s within labelled
     transition system l1 by internal transitions, provided weak_reduction is true.
     Otherwise it generates a set with only state s in it.
  */
  template < class LTS_TYPE >
  set_of_states collect_reachable_states_via_taus(
              const set_of_states& s, 
              const lts_cache<LTS_TYPE>& weak_property_cache, 
              const bool weak_reduction)
  {
    set_of_states result(s);
    if (!weak_reduction)
    {
      return result;
    }
    set_of_states visited(s);
    std::deque<state_type> todo_stack(s.begin(),s.end());
    while (todo_stack.size()>0)
    {
      state_type current_state=todo_stack.front();
      todo_stack.pop_front();
      for(const state_type s: weak_property_cache.tau_reachable_states(current_state)) 
      {
        if (visited.insert(s).second)  // The element has been inserted.
        {
          todo_stack.push_back(s);
          result.insert(s);
        }
      }
    }
    
    return result;
  }

  template < class LTS_TYPE >
  set_of_states collect_reachable_states_via_taus(
                  const state_type s, 
                  const lts_cache<LTS_TYPE>& weak_property_cache, 
                  const bool weak_reduction)
  {
    set_of_states set_with_s({s});
    return collect_reachable_states_via_taus(set_with_s, weak_property_cache, weak_reduction);
  }

  template < class LTS_TYPE >
  set_of_states collect_reachable_states_via_an_action(
                 const state_type s,
                 const label_type e,  // This is already the hidden action.
                 const lts_cache<LTS_TYPE>& weak_property_cache,
                 const bool weak_reduction,
                 const LTS_TYPE& l)
  {
    const set_of_states set_before_action_e=collect_reachable_states_via_taus(s,weak_property_cache,weak_reduction);
    set_of_states states_reachable_via_e;
    for(const state_type s: set_before_action_e)
    {
      for(const transition& t: weak_property_cache.transitions(s))
      {
        {
          if (l.apply_hidden_label_map(t.label())==e)
          { 
            assert(set_before_action_e.count(t.from())>0);
            states_reachable_via_e.insert(t.to());
          }
        }
      }
    }
    return collect_reachable_states_via_taus(states_reachable_via_e, weak_property_cache, weak_reduction);
  }

  /* This function implements the insertion of <p,state(), p.states()> in the anti_chain.
     Concretely, this means that p.states() is inserted among the sets s1,...,sn associated to p.state().
     It is important that an anti_chain contains for each state a set of states of which 
     no set is a subset of another. The idea is that if such two sets occur, it is enough
     to keep the smallest.
     If p.states() is smaller than a set si associated to p.state(), this set is removed.
     If p.states() is larger than a set si, there is no need to add p.states(), as a better candidate
     is already there. 
     This function returns true if insertion was succesful, and false otherwise.
   */
  template <  class COUNTER_EXAMPLE_CONSTRUCTOR >
  inline bool antichain_insert(
                  anti_chain_type& anti_chain, 
                  const state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>& impl_spec)
  {
    // First check whether there is a set in the antichain for impl_spec.state() which is smaller than impl_spec.states().
    // If so, impl_spec.states() does not have to be inserted in the anti_chain.
    for(anti_chain_type::const_iterator i=anti_chain.lower_bound(impl_spec.state()); i!=anti_chain.upper_bound(impl_spec.state()); ++i)
    {
      const set_of_states s=i->second;
      // If s is included in impl_spec.states() 
      if (std::includes(impl_spec.states().begin(), impl_spec.states().end(), s.begin(),s.end()))
      {
        return false;
      }
    }

    // Here impl_spec.states() must be inserted in the antichain. Moreover, all sets in the antichain that 
    // are a superset of impl_spec.states() must be removed.
    for(anti_chain_type::iterator i=anti_chain.lower_bound(impl_spec.state()); i!=anti_chain.upper_bound(impl_spec.state()); )
    {
      const set_of_states s=i->second;
      // if s is a superset of impl_spec.states() 
      // if (std::includes(impl_spec.states().begin(),impl_spec.states().end(),s.begin(),s.end()))  
      if (std::includes(s.begin(), s.end(), impl_spec.states().begin(),impl_spec.states().end()))  
      {
        // set s must be removed. 
        i=anti_chain.erase(i);
      }
      else
      {
        ++i;
      }
    }
    anti_chain.insert(std::pair<detail::state_type,detail::set_of_states>(impl_spec.state(),impl_spec.states()));
    return true;
  }
  
  /* Calculate the states that are stable and reachable through tau-steps */
  template < class LTS_TYPE >
  const set_of_states& calculate_tau_reachable_states(
        const set_of_states& states, 
        const lts_cache<LTS_TYPE>& weak_property_cache)
  {
    static std::map<set_of_states, set_of_states> cache;
    const std::map<set_of_states, set_of_states>::const_iterator i=cache.find(states);
    if (i!=cache.end())
    {
      return i->second;
    }

    static set_of_states visited;
    assert(visited.empty());
    static std::stack < state_type > todo_stack;
    assert(todo_stack.empty());
    set_of_states result;

    for(const state_type s: states)
    {
      if (weak_property_cache.stable(s))
      { 
        // Put the outgoing action labels in a set and put these in the result.
        result.insert(s);
      }
      else
      {
        visited.insert(s);
        todo_stack.push(s);
      }
    }
    
    while (todo_stack.size()>0)
    {
      const state_type current_state=todo_stack.top();
      todo_stack.pop();
      // Consider the states reachable in one tau step. 
      for(const state_type s: weak_property_cache.tau_reachable_states(current_state)) 
      {
        if (weak_property_cache.stable(s))
        { 
          // Put the outgoing action labels in a set and put these in the result.
          result.insert(s);
        }
        else
        {
          if (visited.insert(s).second) // t.to() is a new state.
          { 
            todo_stack.push(s);
          }
        }
      }
    }
    cache[states]=result;
    visited.clear();
    return cache[states];
  }

  /* This function checks that the refusals(impl) are contained in the refusals of spec, where
     the refusals of spec are defined by { r | exists s in spec. r in refusals(s) }.
     This is equivalent to saying that for all stable states s' reachable via tau's from impl there is a stable t'
     reachable via tau's from some t in spec, enable(t') is contained in enable(s'). The last expression is calculated below. 
     If enable(t') is not included in enable(s'), their is a problematic action a. This action is returned as "culprit".
     It can be used to construct an extended counterexample. 
  */
  template < class LTS_TYPE >
  bool refusals_contained_in(
              const state_type impl, 
              const set_of_states& spec, 
              const lts_cache<LTS_TYPE>& weak_property_cache,
              label_type& culprit,
              const LTS_TYPE& l,
              const bool provide_a_counter_example)
  {
    if (!weak_property_cache.stable(impl)) return true; // Checking in case of instability is not necessary, but rather time consuming. 

    // This function calculates whether refusals(impl) are not included in the refusals(spec).
    // This is equivalent to:
    // There is a tau-reachable stable state s' from impl, such that for each tau-reachable stable state s''
    // from any of the states in spec: enable(s'')\enable(s') is not empty.

    // First calculate the refusal sets reachable from spec.
    const set_of_states& tau_reachable_states_of_the_specification=calculate_tau_reachable_states(spec,weak_property_cache);

    // Now walk through the tau-reachable stable states s' of impl.
    static std::unordered_set<state_type> visited;
    assert(visited.empty());
    visited.insert(impl);
    static std::stack < state_type > todo_stack;
    assert(todo_stack.empty());
    todo_stack.push(impl);

    while (todo_stack.size()>0)
    {
      const state_type current_state=todo_stack.top();
      todo_stack.pop();
      if (weak_property_cache.stable(current_state))
      {
        // Put the outgoing action labels in a set and put these in the result.
       
        const action_label_set& impl_enabled_action_set=weak_property_cache.action_labels(current_state);

        bool success=false;
        // Compare the obtained enable set of s' with all those of the specification.
        // for(const action_label_set& spec_action_labels: enabled_stable_sets_of_specification)
        for(const state_type s: tau_reachable_states_of_the_specification)
        {
          // Check whether the enabled actions of spec are included in the enabled actions of impl.
          // This is equivalent to checking that all spec_action_labels are in the impl_enabled_action_set.
          const action_label_set& spec_action_labels=weak_property_cache.action_labels(s);
          // Warning: std::includes checks whether the second range is included in the first. 
          bool inclusion_success=std::includes(impl_enabled_action_set.begin(), impl_enabled_action_set.end(),
                                               spec_action_labels.begin(), spec_action_labels.end());
          if (inclusion_success)
          {
            success=true;
            break;
          }
          else 
          { 
            // Find the offending action. 
            for(const label_type a: spec_action_labels)
            {
              if (impl_enabled_action_set.count(a)==0) // We want to know which action caused the problem. 
              {
                culprit=a;
              }
            }
          }
        }
        if (!success)
        {
          if (provide_a_counter_example)
          {
            // Print the acceptance set of the implementation 
            if (impl_enabled_action_set.empty())
            {
              mCRL2log(log::verbose) << "The acceptance of the left process is empty.\n";
            }
            else
            {
              mCRL2log(log::verbose) << "A stable acceptance set of the left process is:\n";
              for(const label_type a: impl_enabled_action_set)
              {
                mCRL2log(log::verbose) << l.action_label(a) << "\n";
              }
            }
            // Print the acceptance sets of the specification. 
            // if (enabled_stable_sets_of_specification.empty())
            if (tau_reachable_states_of_the_specification.empty())
            {
              mCRL2log(log::verbose) << "The process at the right has no acceptance sets.\n";
            }
            else 
            {
              mCRL2log(log::verbose) << "Below all corresponding stable acceptance sets of the right process are provided:\n";
              for(const state_type s: tau_reachable_states_of_the_specification)
              {
                const action_label_set& spec_action_labels=weak_property_cache.action_labels(s);
                mCRL2log(log::verbose) << "An acceptance set of the right process is:\n";
                for(const label_type a: spec_action_labels)
                {
                  mCRL2log(log::verbose) << l.action_label(a) << "\n";
                }
              }
            }
            mCRL2log(log::verbose) << "Finished printing acceptance sets.\n";          
            // Ready printing acceptance sets. 
          }
          
          visited.clear();
          todo_stack=std::stack < state_type >();
          return false;
        }
      }
      else
      {
        // Put the states reachable in one tau step onto the todo stack, if they have not 
        // been visited yet. 
        for(const state_type s: weak_property_cache.tau_reachable_states(current_state)) 
        {
          if (visited.insert(s).second) // s is a new state.
          {
            todo_stack.push(s);
          }
        }

      }
    }
    visited.clear();
    return true;   
  }
  
  
} // namespace detail
} // namespace lts
} // namespace mcrl2

#endif //  _LIBLTS_FAILURES_REFINEMENT_H


