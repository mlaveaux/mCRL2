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
// M. Laveaux, J.F. Groote and T.A.C. Willemse
// Correct and Efficient Antichain Algorithms for Refinement Checking. Logical Methods in Computer Science 17(1) 2021
//
// There are six algorithms. One for trace inclusion, one for failures inclusion and one for failures-divergence
// inclusion. All algorithms come in a variant with and without internal steps. It is possible to generate a counter
// transition system in case the inclusion is answered by no.

#ifndef LIBLTS_FAILURES_REFINEMENT_H
#define LIBLTS_FAILURES_REFINEMENT_H

#include "mcrl2/lps/exploration_strategy.h"
#include "mcrl2/lts/detail/counter_example.h"
#include "mcrl2/lts/detail/liblts_bisim_dnj.h"

namespace mcrl2::lts 
{
  
namespace detail
{

typedef std::size_t state_type;
typedef std::size_t label_type;
typedef std::set<state_type> set_of_states;
typedef std::multimap<detail::state_type, detail::set_of_states> anti_chain_type;
typedef std::set<label_type> action_label_set;

template <class COUNTER_EXAMPLE_CONSTRUCTOR>
class state_states_counter_example_index_triple
{
protected:
  detail::state_type m_state;
  detail::set_of_states m_states;
  typename COUNTER_EXAMPLE_CONSTRUCTOR::index_type m_counter_example_index;

public:
  state_states_counter_example_index_triple() {}

  /// \brief Constructor.
  state_states_counter_example_index_triple(const state_type state,
      const set_of_states& states,
      const typename COUNTER_EXAMPLE_CONSTRUCTOR::index_type& counter_example_index)
      : m_state(state),
        m_states(states),
        m_counter_example_index(counter_example_index)
  {}

  /// \brief Get the state.
  state_type state() const { return m_state; }

  void swap(state_states_counter_example_index_triple& other)
  {
    std::swap(m_state, other.m_state);
    std::swap(m_states, other.m_states);
    std::swap(m_counter_example_index, other.m_counter_example_index);
  }

  /// \brief Get the set of states.
  const set_of_states& states() const { return m_states; }

  /// \brief Get the counter example index.
  const typename COUNTER_EXAMPLE_CONSTRUCTOR::index_type& counter_example_index() const
  {
    return m_counter_example_index;
  }
};

inline bool antichain_include(anti_chain_type& anti_chain,
  const detail::state_type& impl,
  const detail::set_of_states& spec);

inline bool antichain_insert(anti_chain_type& anti_chain,
  const detail::state_type& impl,
  const detail::set_of_states& spec);

// The class below recalls what the stable states and the states with a divergent
// self loop of a transition system are, such that it does not have to be recalculated each time again.
template <class LTS_TYPE>
class lts_cache
{
protected:
  LTS_TYPE& m_l;
  std::vector<std::vector<state_type>> m_tau_reachable_states;
  std::vector<std::vector<transition>> m_sorted_transitions;
  std::vector<bool> m_divergent;
  std::vector<action_label_set> m_enabled_actions;

private:
  void calculate_weak_property_cache(const bool weak_reduction)
  {
    scc_partitioner<LTS_TYPE> strongly_connected_component_partitioner(m_l);
    for (const transition& t : m_l.get_transitions())
    {
      assert(t.from() < m_l.num_states());
      if (m_l.is_tau(m_l.apply_hidden_label_map(t.label())) && weak_reduction)
      {
        m_tau_reachable_states[t.from()].push_back(t.to()); // There is an outgoing tau.
      }
      m_sorted_transitions[t.from()].push_back(t);

      if (weak_reduction && m_l.is_tau(m_l.apply_hidden_label_map(t.label()))
          && strongly_connected_component_partitioner.in_same_class(t.from(), t.to()))
      {
        m_divergent[t.from()] = true; // There is a self loop.
      }
      m_enabled_actions[t.from()].insert(m_l.apply_hidden_label_map(t.label()));
    }
  }

public:

  lts_cache(/* const */ LTS_TYPE& l,
      const bool weak_reduction) /* l is not changed, but the use of the scc partitioner requires l to be non const */
      : m_l(l),
        m_tau_reachable_states(l.num_states()),
        m_sorted_transitions(l.num_states()),
        m_divergent(l.num_states(), false),
        m_enabled_actions(l.num_states())
  {
    calculate_weak_property_cache(weak_reduction);
  }

  bool stable(const state_type s) const { return m_tau_reachable_states[s].size() == 0; }

  const std::vector<state_type>& tau_reachable_states(const state_type s) const { return m_tau_reachable_states[s]; }

  const std::vector<transition>& transitions(const state_type s) const
  {
    assert(s < m_sorted_transitions.size());
    return m_sorted_transitions[s];
  }

  bool diverges(const state_type s) const { return m_divergent[s]; }

  const action_label_set& action_labels(const state_type s) const { return m_enabled_actions[s]; }
};

template <class LTS_TYPE>
set_of_states collect_reachable_states_via_taus(const state_type s,
    const lts_cache<LTS_TYPE>& weak_property_cache,
    const bool weak_reduction);

template <class LTS_TYPE>
set_of_states collect_reachable_states_via_an_action(const state_type s,
    const label_type e,
    const lts_cache<LTS_TYPE>& weak_property_cache,
    const bool weak_reduction,
    const LTS_TYPE& l);

template <class LTS_TYPE>
bool refusals_contained_in(const state_type impl,
    const set_of_states& spec,
    const lts_cache<LTS_TYPE>& weak_property_cache,
    label_type& culprit,
    const LTS_TYPE& l,
    const bool provide_a_counter_example,
    const bool structured_output);

} // namespace detail

enum class refinement_type
{
  trace,
  failures,
  failures_divergence
};

template <typename T>
struct refinement_statistics
{
  refinement_statistics(detail::anti_chain_type& antichain, std::deque<T>& working)
      : antichain(antichain),
        working(working)
  {}

  detail::anti_chain_type& antichain;
  std::deque<T>& working;
  std::size_t max_working = 0;       // The largest size of working.
  std::size_t max_antichain = 0;     // The largest size of the antichain.
  std::size_t antichain_misses = 0;  // Number of times a pair was inserted into the antichain.
  std::size_t antichain_inserts = 0; // Number of times antichain_insert was called.
};

/// \brief Print a message to debugging containing information about the given statistics.
template <typename T>
void report_statistics(refinement_statistics<T>& stats)
{
  mCRL2log(log::debug) << "working (current: " << stats.working.size() << ", max: " << stats.max_working << ").\n";
  mCRL2log(log::debug) << "antichain (hits: " << stats.antichain_inserts - stats.antichain_misses
                       << ", misses: " << stats.antichain_misses << ", size: " << stats.antichain.size()
                       << ", max: " << stats.max_antichain << ")\n";
}

/// \brief Preprocess the LTS for destructive refinement checking.
/// \param lts The lts to preprocess.
/// \param init The initial state of the right LTS that was merged.
/// \return A pair where the first element is the state number of init in the reduced
///         lts and the second value indicate whether this state in equal to lts.initial_state.
template <typename LTS_TYPE>
std::pair<std::size_t, bool>
reduce(LTS_TYPE& lts, const bool weak_reduction, const bool preserve_divergence, std::size_t l2_init)
{
  lts.clear_state_labels();
  if (weak_reduction)
  {
    // Remove inert tau loops when requested, but preserve divergences for failures and failures-divergence.
    detail::scc_partitioner<LTS_TYPE> scc_part(lts);
    l2_init = scc_part.get_eq_class(l2_init);
    scc_part.replace_transition_system(preserve_divergence);
  }

  detail::bisim_partitioner_dnj<LTS_TYPE> bisim_part(lts, weak_reduction, preserve_divergence);
  // Assign the reduced LTS, and set init_l2.
  l2_init = bisim_part.get_eq_class(l2_init);
  bisim_part.finalize_minimized_LTS();

  return std::make_pair(l2_init, l2_init == lts.initial_state());
}

/// \brief This function checks using algorithms in the paper mentioned above
/// whether transition system l1 is included in transition system l2, in the
/// sense of trace inclusions, failures inclusion and divergence failures
/// inclusion.
/// \param weak_reduction Remove inert tau loops.
/// \param strategy Choose between breadth and depth first.
/// \param preprocess Uses (divergence preserving) branching bisimulation and tau scc reduction to reduce the input
/// LTSs. \param generate_counter_example If set, a labelled transition system is generated
///        that can act as a counterexample. It consists of a trace, followed by
///        outgoing transitions representing a refusal set.
template <class LTS_TYPE, class COUNTER_EXAMPLE_CONSTRUCTOR = detail::dummy_counter_example_constructor>
bool destructive_refinement_checker(LTS_TYPE& l1,
    LTS_TYPE& l2,
    const refinement_type refinement,
    const bool weak_reduction,
    const lps::exploration_strategy strategy,
    const bool preprocess = true,
    COUNTER_EXAMPLE_CONSTRUCTOR generate_counter_example = detail::dummy_counter_example_constructor())
{
  assert(strategy == lps::exploration_strategy::es_breadth
         || strategy == lps::exploration_strategy::es_depth); // Need a valid strategy.

  // For weak-failures and failures-divergence, the existence of tau loops make a difference.
  // Therefore, we apply bisimulation reduction preserving divergences.
  // A typical example is a.(b+c) which is not weak-failures included n a.tau*.(b+c). The lhs has failure pairs
  // <a,{a}>, <a,{}> while the rhs has only failure pairs <a,{}>, as the state after the a is not stable.
  const bool preserve_divergence = weak_reduction && (refinement != refinement_type::trace);

  if (!generate_counter_example.is_dummy() && preprocess)
  {
    // Counter example is requested, apply bisimulation to l2.
    reduce(l2, weak_reduction, preserve_divergence, l2.initial_state());
  }

  std::size_t init_l2 = l2.initial_state() + l1.num_states();
  mcrl2::lts::detail::merge(l1, l2);
  l2.clear(); // No use for l2 anymore.

  if (generate_counter_example.is_dummy() && preprocess)
  {
    // No counter example is requested. We can use bisimulation preprocessing.
    bool initial_equal = false;
    std::tie(init_l2, initial_equal) = reduce(l1, weak_reduction, preserve_divergence, init_l2);

    if (initial_equal && weak_reduction)
    {
      mCRL2log(log::verbose) << "The two LTSs are";
      if (preserve_divergence)
      {
        mCRL2log(log::verbose) << " divergence-preserving";
      }
      mCRL2log(log::verbose) << " branching bisimilar, so there is no need to check the refinement relation.\n";
      return true;
    }
  }

  std::deque<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> working;
  detail::anti_chain_type anti_chain;
  detail::anti_chain_type positive_anti_chain;
  refinement_statistics<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> stats(
      anti_chain,
      working);
  const detail::lts_cache<LTS_TYPE> weak_property_cache(l1, weak_reduction);

  bool result = check_refinement(l1,
      weak_property_cache,
      working,
      stats,
      anti_chain,
      positive_anti_chain,
      generate_counter_example,
      l1.initial_state(),
      init_l2,
      refinement,
      weak_reduction,
      strategy);
  report_statistics(stats);
  return result;
}

template <typename LTS_TYPE, typename COUNTER_EXAMPLE_CONSTRUCTOR>
bool check_refinement(LTS_TYPE& l1,
    const detail::lts_cache<LTS_TYPE>& weak_property_cache,
    std::deque<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>>& working,
    refinement_statistics<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>>& stats,
    detail::anti_chain_type& anti_chain,
    detail::anti_chain_type& anti_chain_positive,
    COUNTER_EXAMPLE_CONSTRUCTOR& generate_counter_example,
    detail::state_type init_l1,
    detail::state_type init_l2,
    const refinement_type refinement,
    bool weak_reduction,
    const lps::exploration_strategy strategy)
{
  // let working be a stack containg the triple (init1,{s|init2-->s},root_index);
  working.clear();
  working.push_back({detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>(init_l1,
      detail::collect_reachable_states_via_taus(init_l2, weak_property_cache, weak_reduction),
      generate_counter_example.root_index())});

  // let antichain := emptyset;
  anti_chain.clear();
  detail::antichain_insert(anti_chain,
      working.front().state(),
      working.front().states()); // antichain := antichain united with (impl,spec);
                        // This line occurs at another place in the code than in
                        // the original algorithm, where insertion in the anti-chain
                        // was too late, causing too many impl-spec pairs to be investigated.

  while (!working.empty()) // while working!=empty
  {
    detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>
        impl_spec; // pop (impl,spec) from working;
    impl_spec.swap(working.front());
    stats.max_working = std::max(working.size(), stats.max_working);
    stats.max_antichain = std::max(anti_chain.size(), stats.max_antichain);
    working.pop_front(); // At this point it could be checked whether impl_spec still exists in anti_chain.
                         // Small scale experiments show that this is a little bit more expensive than doing the
                         // explicit check below.

    bool spec_diverges = false;
    if (refinement == refinement_type::failures_divergence)
    {
      // Only compute when the result is required.
      for (detail::state_type s : impl_spec.states())
      {
        if (weak_property_cache.diverges(s))
        {
          spec_diverges = true;
          break;
        }
      }
    }

    // if not diverges(spec) or not CheckDiv (refinement == failures_divergence_preorder)
    if (!spec_diverges || refinement != refinement_type::failures_divergence)
    {
      if (weak_property_cache.diverges(impl_spec.state())
          && refinement == refinement_type::failures_divergence) // if impl diverges and CheckDiv
      {
        generate_counter_example.save_counter_example(impl_spec.counter_example_index(), l1);
        report_statistics(stats);
        return false; // return false;
      }

      if (refinement == refinement_type::failures || refinement == refinement_type::failures_divergence)
      {
        detail::label_type offending_action = std::size_t(-1);
        // if refusals(impl) not contained in refusals(spec) then
        if (!detail::refusals_contained_in(impl_spec.state(),
                impl_spec.states(),
                weak_property_cache,
                offending_action,
                l1,
                !generate_counter_example.is_dummy(),
                generate_counter_example.is_structured()))
        {
          generate_counter_example.save_counter_example(impl_spec.counter_example_index(), l1);
          report_statistics(stats);
          return false; // return false;
        }
      }

      for (const transition& t : weak_property_cache.transitions(impl_spec.state()))
      {
        const typename COUNTER_EXAMPLE_CONSTRUCTOR::index_type new_counterexample_index
            = generate_counter_example.add_transition(t.label(), impl_spec.counter_example_index());

        detail::set_of_states spec_prime;
        if (l1.is_tau(l1.apply_hidden_label_map(t.label())) && weak_reduction) // if e=tau then
        {
          spec_prime = impl_spec.states(); // spec' := spec;
        }
        else
        { // spec' := {s' | exists s in spec. s-e->s'};
          for (const detail::state_type s : impl_spec.states())
          {
            detail::set_of_states reachable_states_from_s_via_e = detail::collect_reachable_states_via_an_action(s,
                l1.apply_hidden_label_map(t.label()),
                weak_property_cache,
                weak_reduction,
                l1);
            spec_prime.insert(reachable_states_from_s_via_e.begin(), reachable_states_from_s_via_e.end());
          }
        }

        if (spec_prime.empty()) // if spec'={} then
        {
          generate_counter_example.save_counter_example(new_counterexample_index, l1);
          report_statistics(stats);
          return false; //    return false;
        }
        
        // if (impl',spec') in antichain is not true then
        ++stats.antichain_inserts;
        const detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR> impl_spec_counterex(t.to(),
            spec_prime,
            new_counterexample_index);
        if (!detail::antichain_include(anti_chain_positive, t.to(), spec_prime) && detail::antichain_insert(anti_chain, t.to(), spec_prime))
        {
          ++stats.antichain_misses;
          if (strategy == lps::exploration_strategy::es_breadth)
          {
            working.push_back(impl_spec_counterex); // add(impl,spec') at the bottom of the working;
          }
          else if (strategy == lps::exploration_strategy::es_depth)
          {
            working.push_front(impl_spec_counterex); // push(impl,spec') into working;
          }
        }
      }
    }
  }

  for (const auto& [impl, spec] : anti_chain)
  {
    detail::antichain_insert(anti_chain_positive, impl, spec);
  }

  return true; // return true;
}

namespace detail
{
/* This function generates the set of states reachable from s within labelled
   transition system l1 by internal transitions, provided weak_reduction is true.
   Otherwise it generates a set with only state s in it.
*/
template <class LTS_TYPE>
set_of_states collect_reachable_states_via_taus(const set_of_states& s,
    const lts_cache<LTS_TYPE>& weak_property_cache,
    const bool weak_reduction)
{
  set_of_states result(s);
  if (!weak_reduction)
  {
    return result;
  }
  set_of_states visited(s);
  std::deque<state_type> todo_stack(s.begin(), s.end());
  while (todo_stack.size() > 0)
  {
    state_type current_state = todo_stack.front();
    todo_stack.pop_front();
    for (const state_type s : weak_property_cache.tau_reachable_states(current_state))
    {
      if (visited.insert(s).second) // The element has been inserted.
      {
        todo_stack.push_back(s);
        result.insert(s);
      }
    }
  }

  return result;
}

template <class LTS_TYPE>
set_of_states collect_reachable_states_via_taus(const state_type s,
    const lts_cache<LTS_TYPE>& weak_property_cache,
    const bool weak_reduction)
{
  set_of_states set_with_s({s});
  return collect_reachable_states_via_taus(set_with_s, weak_property_cache, weak_reduction);
}

template <class LTS_TYPE>
set_of_states collect_reachable_states_via_an_action(const state_type s,
    const label_type e, // This is already the hidden action.
    const lts_cache<LTS_TYPE>& weak_property_cache,
    const bool weak_reduction,
    const LTS_TYPE& l)
{
  const set_of_states set_before_action_e = collect_reachable_states_via_taus(s, weak_property_cache, weak_reduction);
  set_of_states states_reachable_via_e;
  for (const state_type s : set_before_action_e)
  {
    for (const transition& t : weak_property_cache.transitions(s))
    {
      {
        if (l.apply_hidden_label_map(t.label()) == e)
        {
          assert(set_before_action_e.count(t.from()) > 0);
          states_reachable_via_e.insert(t.to());
        }
      }
    }
  }
  return collect_reachable_states_via_taus(states_reachable_via_e, weak_property_cache, weak_reduction);
}

inline bool antichain_include(anti_chain_type& anti_chain,
  const detail::state_type& impl,
  const detail::set_of_states& spec)
{
  // First check whether there is a set in the antichain for impl_spec.state() which is smaller than impl_spec.states().
  // If so, impl_spec.states() is included in the antichain.
  for (anti_chain_type::const_iterator i = anti_chain.lower_bound(impl);
       i != anti_chain.upper_bound(impl);
       ++i)
  {
    const set_of_states s = i->second;
    // If s is included in impl_spec.states()
    if (std::includes(spec.begin(), spec.end(), s.begin(), s.end()))
    {
      return true;
    }
  }

  return false;
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
inline bool antichain_insert(anti_chain_type& anti_chain,
  const detail::state_type& impl,
  const detail::set_of_states& spec)
{
  if (antichain_include(anti_chain, impl, spec))
  {
    return false;
  }

  // Here impl_spec.states() must be inserted in the antichain. Moreover, all sets in the antichain that
  // are a superset of impl_spec.states() must be removed.
  for (anti_chain_type::iterator i = anti_chain.lower_bound(impl);
       i != anti_chain.upper_bound(impl);)
  {
    const set_of_states s = i->second;
    // if s is a superset of impl_spec.states()
    // if (std::includes(impl_spec.states().begin(),impl_spec.states().end(),s.begin(),s.end()))
    if (std::includes(s.begin(), s.end(), spec.begin(), spec.end()))
    {
      // set s must be removed.
      i = anti_chain.erase(i);
    }
    else
    {
      ++i;
    }
  }
  anti_chain.emplace(impl, spec);
  return true;
}

/* Calculate the states that are stable and reachable through tau-steps */
template <class LTS_TYPE>
const set_of_states& calculate_tau_reachable_states(const set_of_states& states,
    const lts_cache<LTS_TYPE>& weak_property_cache)
{
  static std::map<set_of_states, set_of_states> cache;
  const std::map<set_of_states, set_of_states>::const_iterator i = cache.find(states);
  if (i != cache.end())
  {
    return i->second;
  }

  static set_of_states visited;
  assert(visited.empty());
  static std::stack<state_type> todo_stack;
  assert(todo_stack.empty());
  set_of_states result;

  for (const state_type s : states)
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

  while (todo_stack.size() > 0)
  {
    const state_type current_state = todo_stack.top();
    todo_stack.pop();
    // Consider the states reachable in one tau step.
    for (const state_type s : weak_property_cache.tau_reachable_states(current_state))
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
  cache[states] = result;
  visited.clear();
  return cache[states];
}

/// \brief This function checks that the refusals(impl) are contained in the refusals of spec, where
///        the refusals of spec are defined by { r | exists s in spec. r in refusals(s) and stable(r) }.
/// \details This is equivalent to saying that for all enabled actions of impl it must be contained in the enabled
/// actions
///          of every stable state in spec.
///          If enable(t') is not included in enable(s'), their is a problematic action a. This action is returned as
///          "culprit". It can be used to construct an extended counterexample.
template <class LTS_TYPE>
bool refusals_contained_in(const state_type impl,
    const set_of_states& spec,
    const lts_cache<LTS_TYPE>& weak_property_cache,
    label_type& culprit,
    const LTS_TYPE& l,
    const bool provide_a_counter_example,
    const bool structured_output)
{
  if (!weak_property_cache.stable(impl))
    return true; // Checking in case of instability is not necessary, but rather time consuming.

  const action_label_set& impl_action_labels = weak_property_cache.action_labels(impl);
  bool success = false;

  // Compare the obtained enable set of s' with all those of the specification.
  for (const state_type s : spec)
  {
    // Only stable states in this set should be checked.
    if (!weak_property_cache.stable(s))
    {
      continue;
    }

    // Check whether the enabled actions of spec are included in the enabled actions of impl.
    // This is equivalent to checking that all spec_action_labels are in the impl_enabled_action_set.
    const action_label_set& spec_action_labels = weak_property_cache.action_labels(s);

    // Warning: std::includes checks whether the second range is included in the first.
    bool inclusion_success = std::includes(impl_action_labels.begin(),
        impl_action_labels.end(),
        spec_action_labels.begin(),
        spec_action_labels.end());
    if (inclusion_success)
    {
      success = true;
      break;
    }
    else
    {
      // Find the offending action.
      for (const label_type a : spec_action_labels)
      {
        if (impl_action_labels.count(a) == 0) // We want to know which action caused the problem.
        {
          culprit = a;
          break;
        }
      }
    }
  }

  if (!success && provide_a_counter_example)
  {
    // Print the acceptance set of the implementation
    if (impl_action_labels.empty())
    {
      if (structured_output)
      {
        std::cout << "left-acceptance:\n";
      }
      else
      {
        mCRL2log(log::verbose) << "The acceptance of the left process is empty.\n";
      }
    }
    else
    {
      if (structured_output)
      {
        std::cout << "left-acceptance: ";
      }
      else
      {
        mCRL2log(log::verbose) << "A stable acceptance set of the left process is:\n";
      }
      std::string sep = "";
      for (const label_type a : impl_action_labels)
      {
        if (structured_output)
        {
          std::cout << sep << l.action_label(a);
          sep = ";";
        }
        else
        {
          mCRL2log(log::verbose) << l.action_label(a) << "\n";
        }
      }
      if (structured_output)
      {
        std::cout << "\n";
      }
    }

    // Print the acceptance sets of the specification.
    if (spec.empty())
    {
      if (structured_output)
      {
        std::cout << "right-acceptance-sets: 0\n";
      }
      else
      {
        mCRL2log(log::verbose) << "The process at the right has no acceptance sets.\n";
      }
    }
    else
    {
      set_of_states stable;
      // Only the stable specification states contributed to this counter example.
      std::copy_if(spec.begin(),
          spec.end(),
          std::inserter(stable, stable.end()),
          [=](const state_type s) { return weak_property_cache.stable(s); });

      if (structured_output)
      {
        std::cout << "right-acceptance-sets: " << stable.size() << "\n";
      }
      else
      {
        mCRL2log(log::verbose) << "Below all corresponding stable acceptance sets of the right process are provided:\n";
      }
      for (const state_type s : stable)
      {
        const action_label_set& spec_action_labels = weak_property_cache.action_labels(s);
        if (structured_output)
        {
          std::cout << "right-acceptance: ";
        }
        else
        {
          mCRL2log(log::verbose) << "An acceptance set of the right process is:\n";
        }
        std::string sep = "";
        for (const label_type a : spec_action_labels)
        {
          if (structured_output)
          {
            std::cout << sep << l.action_label(a);
            sep = ";";
          }
          else
          {
            mCRL2log(log::verbose) << l.action_label(a) << "\n";
          }
        }
        if (structured_output)
        {
          std::cout << "\n";
        }
      }
    }
    if (!structured_output)
    {
      mCRL2log(log::verbose) << "Finished printing acceptance sets.\n";
    }
    // Done printing acceptance sets.
  }

  return success;
}

} // namespace detail

} // namespace mcrl2::lts

#endif //  LIBLTS_FAILURES_REFINEMENT_H
