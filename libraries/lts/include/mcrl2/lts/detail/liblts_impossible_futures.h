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
bool check_trace_inclusion_naive(LTS_TYPE& l1,
    const detail::lts_cache<LTS_TYPE>& weak_property_cache,
    std::deque<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>>& working,
    std::multiset<std::pair<detail::state_type, detail::set_of_states>>& discovered,
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
  discovered.clear();
  discovered.insert({working.front().state(), working.front().states()});

  while (!working.empty()) // while working!=empty
  {
    // pop (impl,spec) from working;
    detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>
        impl_spec = working.front();
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
        return false; //    return false;
      }
      
      // if (impl',spec') in antichain is not true then
      const detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR> impl_spec_counterex(t.to(),
          spec_prime,
          new_counterexample_index);
      if (discovered.find({t.to(), spec_prime}) == discovered.end())
      {
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
  std::size_t init_l2 = l2.initial_state() + l1.num_states();
  mcrl2::lts::detail::merge(l1, l2);

  const detail::lts_cache<LTS> weak_property_cache(l1, true);

  std::deque<state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> working = std::deque(
      {state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>(l1.initial_state(),
          detail::collect_reachable_states_via_taus(init_l2, weak_property_cache, true),
          COUNTER_EXAMPLE_CONSTRUCTOR())});
          
  std::multiset<std::pair<detail::state_type, detail::set_of_states>> discovered;
  discovered.insert({working.front().state(), working.front().states()}); // antichain := antichain united with (impl,spec);

  // Used for the weak trace refinement checks
  std::deque<state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> inner_working;
  std::multiset<std::pair<detail::state_type, detail::set_of_states>> inner_discovered;
  
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
              return check_trace_inclusion_naive(l1,
                  weak_property_cache,
                  inner_working,
                  inner_discovered,
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

      if (discovered.find({t.to(), spec_prime}) == discovered.end())
      {
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

  return true;
}

}

} // namespace mcrl2::lts::detail

#endif // LIBLTS_IMPOSSIBLE_FUTURES_H