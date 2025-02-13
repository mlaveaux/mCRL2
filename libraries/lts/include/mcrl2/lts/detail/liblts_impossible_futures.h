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

namespace mcrl2::lts::detail
{

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
  detail::anti_chain_type anti_chain;
  detail::antichain_insert(anti_chain, working.front()); // antichain := antichain united with (impl,spec);

  // Used for the weak trace refinement checks
  detail::anti_chain_type inner_anti_chain;
  detail::anti_chain_type inner_anti_chain_updated;
  std::deque<state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> inner_working;
  refinement_statistics<detail::state_states_counter_example_index_triple<COUNTER_EXAMPLE_CONSTRUCTOR>> stats(inner_anti_chain, inner_working);
  
  COUNTER_EXAMPLE_CONSTRUCTOR inner_generate_counterexample = COUNTER_EXAMPLE_CONSTRUCTOR();

  while (!working.empty())
  {
    // pop(impl,spec) from working;
    const auto front = working.front();
    working.pop_front();

    const detail::state_type impl = front.state();
    const detail::set_of_states& spec = front.states();

    if (!std::any_of(spec.begin(),
            spec.end(),
            [&](const auto& t)
            {
              // Print the current (impl,spec) pair being inspected
              std::cout << "Checking (" << impl << ", " << t << ")" << std::endl;

              inner_anti_chain = inner_anti_chain_updated;
              if (!check_refinement(l1,
                  weak_property_cache,
                  inner_working,
                  stats,
                  inner_anti_chain_updated,
                  inner_generate_counterexample,
                  t,
                  impl,
                  refinement_type::trace,
                  true,
                  strategy))
                {
                  // Reset for failing inclusion checks.
                  std::swap(inner_anti_chain_updated, inner_anti_chain);
                  return false;
                }

                return true;
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

      if (detail::antichain_insert(anti_chain, impl_spec_counterex))
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

} // namespace mcrl2::lts::detail

#endif // LIBLTS_IMPOSSIBLE_FUTURES_H