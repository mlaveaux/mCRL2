// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts/detail/liblts_failures_refinement.h

#ifndef LIBLTS_IMPOSSIBLE_FUTURES_H
#define LIBLTS_IMPOSSIBLE_FUTURES_H

#include <boost/container/vector.hpp>
#include <deque>

#include "mcrl2/lts/detail/counter_example.h"
#include "mcrl2/lts/detail/liblts_failures_refinement.h"

namespace mcrl2::lts
{
  
template <typename LTS_TYPE, typename COUNTER_EXAMPLE_CONSTRUCTOR>
bool check_trace_inclusion(LTS_TYPE& l1,
    const detail::lts_cache<LTS_TYPE>& weak_property_cache,
    std::deque<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>>& working,
    refinement_statistics<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>>& stats,
    detail::anti_chain_type& anti_chain,
    detail::anti_chain_type& anti_chain_positive,
    COUNTER_EXAMPLE_CONSTRUCTOR& generate_counter_example,
    detail::state_type init_l1,
    detail::state_type init_l2,
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

  for (const auto& [impl, spec] : anti_chain)
  {
    detail::antichain_insert(anti_chain_positive, impl, spec);
  }

  return true; // return true;
}

namespace detail {

/// \brief This function checks using algorithms in the paper mentioned above
/// whether transition system l1 is included in transition system l2, in the
/// sense of trace inclusions, failures inclusion and divergence failures
/// inclusion.
template <typename LTS,
  typename COUNTER_EXAMPLE_CONSTRUCTOR = detail::dummy_counter_example_constructor>
bool destructive_impossible_futures(LTS& l1, LTS& l2, const lps::exploration_strategy strategy)
{
  // Remove tau-loops from l1 to allow the (impl, spec) => (impl', spec) optimisation.
  scc_partitioner<LTS> scc_partitioner(l1);
  scc_partitioner.replace_transition_system(false);

  std::size_t init_l2 = l2.initial_state() + l1.num_states();
  mcrl2::lts::detail::merge(l1, l2);

  const detail::lts_cache<LTS> weak_property_cache(l1, true);

  std::deque<state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> working = std::deque(
      {state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>(l1.initial_state(),
          detail::collect_reachable_states_via_taus(init_l2, weak_property_cache, true),
          COUNTER_EXAMPLE_CONSTRUCTOR())});
  detail::anti_chain_type anti_chain;
  detail::antichain_insert(anti_chain, working.front().state(), working.front().states()); // antichain := antichain united with (impl,spec);
  refinement_statistics<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> stats(anti_chain, working);

  // Used for the weak trace refinement checks
  detail::anti_chain_type positive_anti_chain;
  detail::anti_chain_type inner_anti_chain;
  std::deque<state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> inner_working;
  refinement_statistics<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> inner_stats(inner_anti_chain, inner_working);
  
  COUNTER_EXAMPLE_CONSTRUCTOR inner_generate_counterexample = COUNTER_EXAMPLE_CONSTRUCTOR();

  while (!working.empty())
  {
    // pop(impl,spec) from working;
    const auto front = working.front();
    working.pop_front();

    const detail::state_type impl = front.state();
    const detail::set_of_states& spec = front.states();

    if (weak_property_cache.stable(impl) && !std::any_of(spec.begin(),
            spec.end(),
            [&](const auto& t)
            {
              // Print the current (impl,spec) pair being inspected
              std::cout << "Checking (" << impl << ", " << t << ")" << std::endl;

              return check_trace_inclusion(l1,
                  weak_property_cache,
                  inner_working,
                  inner_stats,
                  inner_anti_chain,
                  positive_anti_chain,
                  inner_generate_counterexample,
                  t,
                  impl,
                  true,
                  strategy);
            }))
    {
      return false;
    }

    for (const transition& t : weak_property_cache.transitions(impl))
    {
      detail::set_of_states spec_prime;
      if (l1.is_tau(l1.apply_hidden_label_map(t.label())))
      {
        spec_prime = spec;
      }
      else
      {
        for (const detail::state_type s : spec)
        {
          detail::set_of_states reachable_states_from_s_via_e = detail::collect_reachable_states_via_an_action(s,
              l1.apply_hidden_label_map(t.label()),
              weak_property_cache,
              true,
              l1);
          spec_prime.insert(reachable_states_from_s_via_e.begin(), reachable_states_from_s_via_e.end());
        }
      }

      if (spec_prime.empty())
      {
        return false;
      }

      auto impl_spec_counterex
          = state_states_counter_example_index_triple<detail::dummy_counter_example_constructor>(t.to(),
              spec_prime,
              detail::dummy_counter_example_constructor());

      ++stats.antichain_inserts;
      if (detail::antichain_insert(anti_chain, t.to(), spec_prime))
      {
        ++stats.antichain_misses;
        if (strategy == lps::exploration_strategy::es_breadth)
        {
          working.push_back(impl_spec_counterex);
        }
        else if (strategy == lps::exploration_strategy::es_depth)
        {
          working.push_front(impl_spec_counterex);
        }
      }
    }
  }

  report_statistics(inner_stats);
  report_statistics(stats);
  return true;
}

}

} // namespace mcrl2::lts::detail

#endif // LIBLTS_IMPOSSIBLE_FUTURES_H