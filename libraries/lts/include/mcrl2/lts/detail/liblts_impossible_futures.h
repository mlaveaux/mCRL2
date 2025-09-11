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

#include "mcrl2/lts/detail/check_complexity.h"
#include "mcrl2/lts/detail/counter_example.h"
#include "mcrl2/lts/detail/liblts_failures_refinement.h"

namespace mcrl2::lts
{

using state_type_if = detail::state_states_counter_example_index_triple<detail::counter_example_constructor>;

template <typename LTS_TYPE>
std::pair<bool, trace> check_trace_inclusion(LTS_TYPE& l1,
    const detail::lts_cache<LTS_TYPE>& weak_property_cache,
    std::deque<state_type_if>& working,
    refinement_statistics<state_type_if>& stats,
    detail::anti_chain_type& anti_chain,
    detail::anti_chain_type& anti_chain_positive,
    detail::anti_chain_type& anti_chain_negative,
    detail::state_type init_l1,
    detail::state_type init_l2,
    bool weak_reduction,
    const lps::exploration_strategy strategy)
{
  detail::counter_example_constructor generate_counter_example("trace", "unused", false);

  // let working be a stack containg the triple (init1,{s|init2-->s},root_index);
  working.clear();
  working.push_back({state_type_if(init_l1,
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
    // pop (impl,spec) from working;
    state_type_if impl_spec = working.front();
    stats.max_working = std::max(working.size(), stats.max_working);
    stats.max_antichain = std::max(anti_chain.size(), stats.max_antichain);
    working.pop_front(); // At this point it could be checked whether impl_spec still exists in anti_chain.
                         // Small scale experiments show that this is a little bit more expensive than doing the
                         // explicit check below.

    if (detail::antichain_include_inverse(anti_chain_negative, impl_spec.state(), impl_spec.states()))
    {
      return std::make_pair(false,
          generate_counter_example.get_trace(l1, impl_spec.counter_example_index())); //    return false;
    }

    for (const transition& t : weak_property_cache.transitions(impl_spec.state()))
    {
      const typename detail::counter_example_constructor::index_type new_counterexample_index
          = generate_counter_example.add_transition(t.label(), impl_spec.counter_example_index());

      detail::set_of_states spec_prime;
      if (l1.is_tau(l1.apply_hidden_label_map(t.label())) && weak_reduction) // if e=tau then
      {
        spec_prime = impl_spec.states(); // spec' := spec;
      }
      else
      { // spec' := {s' | exists s in spec. s-e->s'};
        for (const auto& s : impl_spec.states())
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
        detail::antichain_insert(anti_chain_negative,
            init_l1,
            detail::collect_reachable_states_via_taus(init_l2, weak_property_cache, weak_reduction));
        return std::make_pair(false,
            generate_counter_example.get_trace(l1, new_counterexample_index)); //    return false;
      }

      // if (impl',spec') in antichain is not true then
      ++stats.antichain_inserts;
      const state_type_if impl_spec_counterex(t.to(), spec_prime, new_counterexample_index);
      if (!detail::antichain_include(anti_chain_positive, t.to(), spec_prime)
          && detail::antichain_insert(anti_chain, t.to(), spec_prime))
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

  return std::make_pair(true, trace()); // return true;
}

/// Generates a counter example of the shape AG EF, which is a modal formula [\sigma]\bigvee{i in I} <rho_i>true.
///
/// This counter examples that there is an implementation state (reached from the initial state by \sigma) whose
/// language is not included in any of the specification states reachable by \sigma.
inline void
write_counter_example(std::ostream& stream, const trace& initial_trace, const std::vector<trace>& inner_traces)
{
  if (!initial_trace.actions().empty())
  {
    stream << "<";

    bool first_action = true;
    for (const auto& action : initial_trace.actions())
    {
      if (!first_action)
      {
        stream << " . ";
      }

      stream << action;

      first_action = false;
    }

    stream << ">(";
  }

  bool first_trace = true;
  for (const auto& inner_trace : inner_traces)
  {
    if (!first_trace)
    {
      stream << " && ";
    }

    stream << "[";

    bool first_action = true;
    for (const auto& action : inner_trace.actions())
    {
      if (!first_action)
      {
        stream << " . tau* . ";
      }
      else
      {
        stream << "tau* . ";
        first_action = false;
      }

      stream << action; 
    }

    // Deal with the empty case.
    if (first_action)
    {
      stream << "tau*";
    }

    stream << "]false";
    first_trace = false;
  }

  if (!initial_trace.actions().empty())
  {
    stream << ")";
  }
}

namespace detail
{

/// \brief Checks impossible futures refinement for the given LTSs.
///
/// Impossible futures are defined in the article:
///
/// Marc Voorhoeve, Sjouke Mauw. Impossible futures and determinism, Inf. Process. Lett. 80, 2001.
template <typename LTS>
bool destructive_impossible_futures(LTS& l1,
    LTS& l2,
    const lps::exploration_strategy strategy,
    std::string counter_example_file,
    bool generate_counter_example,
    bool structured_output)
{
  if (counter_example_file.empty())
  {
    counter_example_file = "counter_example_impossible_futures.mcf";
  }

  // Merge this LTS and l and store the result in this LTS.
  // In the resulting LTS, the initial state i of l will have the
  // state number i + N where N is the number of states in this LTS (before the merge).
  // The initial state of l1 is 0, and the initial state of l2 is l1.num_states().

  // Remove tau-loops from l1 to allow the (impl, spec) => (impl', spec) optimisation.
  scc_partitioner<LTS> scc_partitioner(l1);
  scc_partitioner.replace_transition_system(false);

  std::size_t init_l2 = l2.initial_state() + l1.num_states();
  mcrl2::lts::detail::merge(l1, l2);

  const detail::lts_cache<LTS> weak_property_cache(l1, true);

  // The name and output are not used anyway.
  detail::counter_example_constructor ce_constructor("impossible_futures", counter_example_file, structured_output);

  std::deque<state_type_if> working = std::deque({state_type_if(l1.initial_state(),
      detail::collect_reachable_states_via_taus(init_l2, weak_property_cache, true),
      ce_constructor.root_index())});
  detail::anti_chain_type anti_chain;
  detail::antichain_insert(anti_chain,
      working.front().state(),
      working.front().states()); // antichain := antichain united with (impl,spec);
  refinement_statistics<state_type_if> stats(anti_chain, working);

  // Used for the weak trace refinement checks
  detail::anti_chain_type positive_anti_chain;
  detail::anti_chain_type negative_anti_chain;
  detail::anti_chain_type inner_anti_chain;
  std::deque<state_type_if> inner_working;
  refinement_statistics<state_type_if> inner_stats(inner_anti_chain, inner_working);

  std::vector<trace> inner_counter_examples;

  while (!working.empty())
  {
    // pop(impl,spec) from working;
    const auto front = working.front();
    working.pop_front();

    const detail::state_type impl = front.state();
    const detail::set_of_states& spec = front.states();

    inner_counter_examples.clear();
    mCRL2log(log::debug) << "Checking impl: " << impl << ", spec: ";
    for (const auto& state : spec) 
    {
      mCRL2log(log::debug) << "  " << state;
    }
    mCRL2log(log::debug)  << "\n";

    if (weak_property_cache.stable(impl)
        && !std::any_of(spec.begin(),
            spec.end(),
            [&](const auto& t)
            {
              auto [result, trace] = check_trace_inclusion(l1,
                  weak_property_cache,
                  inner_working,
                  inner_stats,
                  inner_anti_chain,
                  positive_anti_chain,
                  negative_anti_chain,
                  t,
                  impl,
                  true,
                  strategy);

              if (!result)
              {
                inner_counter_examples.emplace_back(trace);
              }

              return result;
            }))
    {
      if (generate_counter_example)
      {
        // Construct the counter example.
        if (structured_output)
        {
          write_counter_example(std::cout,
              ce_constructor.get_trace(l1, front.counter_example_index()),
              inner_counter_examples);
        }
        else
        {
          std::ofstream file(counter_example_file);
          if (!file)
          {
            throw mcrl2::runtime_error("Could not open file " + counter_example_file);
          }

          write_counter_example(file,
              ce_constructor.get_trace(l1, front.counter_example_index()),
              inner_counter_examples);
        }
      }

      return false;
    }

    for (const transition& t : weak_property_cache.transitions(impl))
    {
      const typename counter_example_constructor::index_type new_counterexample_index
          = ce_constructor.add_transition(t.label(), front.counter_example_index());

      mCRL2log(log::debug) << "Taking transition: " << t.label() << " from " << impl << " to " << t.to() << std::endl;

      detail::set_of_states spec_prime;
      if (l1.is_tau(l1.apply_hidden_label_map(t.label())))
      {
        spec_prime = spec;
      }
      else
      {
        for (const detail::state_type& s : spec)
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
        if (generate_counter_example)
        {
          // Construct the counter example.
          if (structured_output)
          {
            write_counter_example(std::cout,
                ce_constructor.get_trace(l1, new_counterexample_index),
                std::vector<trace>());
          }
          else
          {
            std::ofstream file(counter_example_file);
            if (!file)
            {
              throw mcrl2::runtime_error("Could not open file " + counter_example_file);
            }

            write_counter_example(file,
                ce_constructor.get_trace(l1, new_counterexample_index),
                std::vector<trace>());
          }
        }

        return false;
      }

      state_states_counter_example_index_triple<counter_example_constructor> impl_spec_counterex
          = detail::state_states_counter_example_index_triple<detail::counter_example_constructor>(t.to(),
              spec_prime,
              new_counterexample_index);

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

} // namespace detail

} // namespace mcrl2::lts

#endif // LIBLTS_IMPOSSIBLE_FUTURES_H