// Author(s): Wieger Wesselink 
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org./mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/explorer.h
/// \brief Generic state space explorer for (stochastic) linear process specifications.

#ifndef MCRL2_LPS_EXPLORER_H
#define MCRL2_LPS_EXPLORER_H

#include <random>
#include <thread>
#include <type_traits>
#include <condition_variable>
#include "mcrl2/utilities/detail/io.h"
#include "mcrl2/utilities/skip.h"
#include "mcrl2/atermpp/standard_containers/vector.h"
#include "mcrl2/atermpp/standard_containers/indexed_set.h"
#include "mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h"
#include "mcrl2/data/consistency.h"
#include "mcrl2/data/enumerator.h"
#include "mcrl2/data/find_quantifier_variables.h"
#include "mcrl2/data/detail/enumerator_iteration_limit.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/lps/detail/instantiate_global_variables.h"
#include "mcrl2/lps/explorer_options.h"
#include "mcrl2/lps/explorer_todo_set.h"
#include "mcrl2/lps/explorer_utilities.h"
#include "mcrl2/lps/find_representative.h"
#include "mcrl2/lps/one_point_rule_rewrite.h"
#include "mcrl2/lps/order_summand_variables.h"
#include "mcrl2/lps/replace_constants_by_variables.h"
#include "mcrl2/lps/resolve_name_clashes.h"
#include "mcrl2/lps/stochastic_state.h"
#include "mcrl2/lps/explorer_projections.h"

#ifdef MCRL2_USE_CONTROL_FLOW
#include <boost/container/small_vector.hpp>
#include "mcrl2/lps/explorer_control_flow.h"
#endif

namespace mcrl2::lps {

template <bool Stochastic, bool Timed, typename Specification>
class explorer: public abortable
{
  public:
    using state_type = std::conditional_t<Stochastic, stochastic_state, state>;
    using state_index_type = std::conditional_t<Stochastic, std::list<std::size_t>, std::size_t>;
    static constexpr bool is_stochastic = Stochastic;
    static constexpr bool is_timed = Timed;

    using indexed_set_for_states_type = atermpp::indexed_set<state, mcrl2::utilities::detail::GlobalThreadSafe>;

    struct transition
    {
      lps::multi_action action;
      state_type state;

      transition(lps::multi_action  action_, const state_type& state_)
       : action(std::move(action_)), state(state_)
      {}
    };

  protected:
    using enumerator_element = data::enumerator_list_element_with_substitution<>;
    const explorer_options m_options;  // must not be a reference.

    // The four data structures that must be separate per thread.
    mutable data::mutable_indexed_substitution<> m_global_sigma;
    data::rewriter m_global_rewr;
    data::enumerator_algorithm<> m_global_enumerator;
    data::enumerator_identifier_generator m_global_id_generator;

    Specification m_global_lpsspec;
    // Mutexes
    std::mutex m_exclusive_state_access;
    // Mutex locked if a process finds the global buffer empty, and unlocked if it is (attempted to be) refilled
    std::mutex m_global_todo_buffer_mutex;
    std::condition_variable m_signal_global_todo_buffer_filled;

    std::vector<data::variable> m_process_parameters;
    std::size_t m_n; // m_n = m_process_parameters.size()
    data::data_expression_list m_initial_state;
    lps::stochastic_distribution m_initial_distribution;
    bool m_recursive = false;
    std::vector<explorer_summand> m_regular_summands;
    std::vector<explorer_summand> m_confluent_summands;

    volatile std::atomic<bool> m_must_abort = false;

    // N.B. The keys are stored in term_appl instead of data_expression_list for performance reasons.
    summand_cache_map global_cache;

    indexed_set_for_states_type m_discovered;

    // used by make_timed_state, to avoid needless creation of vectors
    mutable std::vector<data::data_expression> timed_state;

#ifdef MCRL2_USE_CONTROL_FLOW
    std::vector<cf_graph> m_control_flow_graphs;
#endif

    Specification preprocess(const Specification& lpsspec)
    {
      Specification result = lpsspec;
      detail::instantiate_global_variables(result);
      lps::order_summand_variables(result);
      resolve_summand_variable_name_clashes(result); // N.B. This is a required preprocessing step.
      if (m_options.one_point_rule_rewrite)
      {
        one_point_rule_rewrite(result);
      }
      if (m_options.replace_constants_by_variables)
      {
        replace_constants_by_variables(result, m_global_rewr, m_global_sigma);
      }
      return result;
    }

    // Evaluates whether t0 <= t1
    bool less_equal(const data::data_expression& t0, 
                    const data::data_expression& t1, 
                    data::mutable_indexed_substitution<>& sigma,
                    data::rewriter& rewr) const
    {
      return rewr(data::less_equal(t0, t1),sigma) == data::sort_bool::true_();
    }

    // Find a unique representative in the confluent tau-graph reachable from u0.
    template <typename SummandSequence>
    state find_representative(state& u0, 
                              const SummandSequence& summands, 
                              data::mutable_indexed_substitution<>& sigma,
                              data::rewriter& rewr, 
                              data::enumerator_algorithm<>& enumerator,
                              data::enumerator_identifier_generator& id_generator)
    {
      bool recursive_undo = m_recursive;
      m_recursive = true;
      data::data_expression_list process_parameter_undo = process_parameter_values(sigma);
      state result = lps::find_representative(u0, 
                                              [&](const state& u) 
                                              { return generate_successors(u, summands, sigma, rewr, enumerator, id_generator); });
      set_process_parameter_values(process_parameter_undo, sigma);
      m_recursive = recursive_undo;
      return result;
    }

    template <typename DataExpressionSequence>
    void compute_state(state& result,
                       const DataExpressionSequence& v,
                       data::mutable_indexed_substitution<>& sigma,
                       const data::rewriter& rewr) const
    {
      lps::make_state(result, 
                      v.begin(), 
                      m_n, 
                      [&](data::data_expression& result, const data::data_expression& x) { rewr(result, x, sigma); return; });
    }

    template <typename DataExpressionSequence>
    void compute_stochastic_state(stochastic_state& result,
                                  const stochastic_distribution& distribution, 
                                  const DataExpressionSequence& next_state,
                                  data::mutable_indexed_substitution<>& sigma,
                                  const data::rewriter& rewr,
                                  const data::enumerator_algorithm<>& enumerator) const
    {                                              
      result.clear();
      if (distribution.is_defined())
      {
        enumerator.enumerate<enumerator_element>(
                    distribution.variables(), 
                    distribution.distribution(),
                    sigma,
                    [&](const enumerator_element& p) 
                    {
                      p.add_assignments(distribution.variables(), sigma, rewr);
                      result.probabilities.push_back(p.expression());
                      result.states.emplace_back();
                      compute_state(result.states.back(), next_state, sigma, rewr);
                      return false;
                    },
                    [](const data::data_expression& x) { return x == data::sort_real::real_zero(); }
        );
        data::remove_assignments(sigma, distribution.variables());
        if (m_options.check_probabilities)
        {
          check_stochastic_state(result, rewr);
        }
      }
      else
      {
        result.probabilities.push_back(data::sort_real::real_one());
        result.states.emplace_back();
        compute_state(result.states.back(),next_state,sigma,rewr);
      }
    }
    
    template <typename DataExpressionSequence, typename IndexSequence>
    void compute_projected_state(
        lps::state& result,
        const DataExpressionSequence& next_state,
        const IndexSequence& I_w,
        data::mutable_indexed_substitution<>& sigma,
        const data::rewriter& rewr) const
    {
      lps::make_state(
        result,
        I_w.begin(),
        I_w.size(),
        [&](data::data_expression& target, std::size_t i)
        {
          rewr(target, next_state[i], sigma);
        }
      );
    }

    template <typename DataExpressionSequence, typename IndexSequence>
    void compute_stochastic_projected_state(
        stochastic_state& result,
        const stochastic_distribution& distribution,
        const DataExpressionSequence& next_state,
        const IndexSequence& I_w,
        data::mutable_indexed_substitution<>& sigma,
        const data::rewriter& rewr,
        data::enumerator_algorithm<>& enumerator) const
    {
      result.clear();
      if (distribution.is_defined())
      {
        enumerator.enumerate<enumerator_element>(
          distribution.variables(),
          distribution.distribution(),
          sigma,
          [&](const enumerator_element& p)
          {
            p.add_assignments(distribution.variables(), sigma, rewr);
            result.probabilities.push_back(p.expression());
            result.states.emplace_back();
            compute_projected_state(result.states.back(), next_state, I_w, sigma, rewr);
            return false;
          },
          [](const data::data_expression& x) { return x == data::sort_real::real_zero(); }
        );
        data::remove_assignments(sigma, distribution.variables());

        if (m_options.check_probabilities)
        {
          check_stochastic_state(result, rewr);
        }
      }
      else
      {
        result.probabilities.push_back(data::sort_real::real_one());
        result.states.emplace_back();
        compute_projected_state(result.states.back(), next_state, I_w, sigma, rewr);
      }
    }

    /// Rewrite action a, and put it back in place.
    lps::multi_action rewrite_action(
                             const lps::multi_action& a,
                             data::mutable_indexed_substitution<>& sigma,
                             data::rewriter& rewr) const
    {
      const process::action_list& actions = a.actions();
      const data::data_expression& time = a.time();
      return lps::multi_action(
          process::action_list(
            actions.begin(),
            actions.end(),
            [&](process::action& result, const process::action& a)
            {
              const data::data_expression_list& args = a.arguments();
              process::make_action(result, 
                                          a.label(), 
                                          data::data_expression_list(args.begin(), 
                                                                     args.end(), 
                                                                     [&](data::data_expression& result, 
                                                                         const data::data_expression& x) -> void
                                                                                 { rewr(result, x, sigma); }));
            }
          ),
          a.has_time() ? rewr(time, sigma) : time
        );
    }

    void check_enumerator_solution(const data::data_expression& p_expression, // WAS: const enumerator_element& p, 
                                   const explorer_summand& summand,
                                   data::mutable_indexed_substitution<>& sigma,
                                   data::rewriter& rewr) const
    {
      if (p_expression != data::sort_bool::true_())
      {
        std::string printed_condition = data::pp(p_expression);
        data::remove_assignments(sigma, m_process_parameters);
        data::remove_assignments(sigma, summand.variables);
        data::data_expression reduced_condition = rewr(summand.condition, sigma);

        if (!data::find_quantifier_variables(p_expression).empty()) 
        {
          mCRL2log(log::info) << "The condition contains quantifiers, which means that rewriting to normal form could fail when the quantifier enumeration limit is too low. Use ---qlimit/-Q to increase the limit." << std::endl;
        }

        throw data::enumerator_error("Condition " + data::pp(reduced_condition) +
                                     " does not rewrite to true or false. \nCulprit: "
                                     + printed_condition.substr(0,2000)
                                     + (printed_condition.size() > 2000 ? "..." : ""));
      }
    }

    // Generic enumeration primitive.
    // The callback determines whether solutions are:
    //  - consumed immediately,
    //  - stored in an enumeration cache, or
    //  - transformed into projected transitions.
    template <typename EmitSolution>
    void enumerate_solutions(
      const explorer_summand& summand,
      data::mutable_indexed_substitution<>& sigma,
      data::rewriter& rewr,
      data::data_expression& condition,
      data::enumerator_algorithm<>& enumerator,
      EmitSolution emit_solution
    )
    {
      rewr(condition, summand.condition, sigma);
      if (data::is_false(condition))
      {
        return;
      }

      if (summand.variables.empty())
      {
        // No enumeration needed
        check_enumerator_solution(condition, summand, sigma, rewr);
        emit_solution(data::data_expression_list{});
      }
      else
      {
        enumerator.enumerate<enumerator_element>(
          summand.variables,
          condition,
          sigma,
          [&](const enumerator_element& p)
          {
            check_enumerator_solution(p.expression(), summand, sigma, rewr);
            emit_solution(p.assign_expressions(summand.variables, rewr));
            return false;
          },
          data::is_false
        );
      }
    }

    // Projects the current state `x` on the read parameters `I_r`.
    // * The current state x is encoded in the substitution sigma.
    // * As a side effect the projected state is assigned to `result`.
    void project(lps::state& result, const std::vector<std::size_t>& I_r, const data::mutable_indexed_substitution<>& sigma)
    {
      lps::make_state(
        result,
        I_r.begin(),
        I_r.size(),
        [&](data::data_expression& target, std::size_t j)
        {
          sigma.apply(m_process_parameters[j], target);
        }
      );
    }

    // Generates outgoing transitions for a summand, and reports them via the callback function report_transition.
    // It is assumed that the substitution sigma contains the assignments corresponding to the current state.
    template <bool ReportActions = true, typename SummandSequence, typename ReportTransition = utilities::skip>
    void generate_transitions(
      const explorer_summand& summand,
      const SummandSequence& confluent_summands,
      data::mutable_indexed_substitution<>& sigma,
      data::rewriter& rewr,
      data::data_expression& condition,                // These three variables are passed on such
      state_type& s1,                                  // that they don't have to be declared often.
      atermpp::aterm key,
      data::enumerator_algorithm<>& enumerator,
      data::enumerator_identifier_generator& id_generator,
#ifdef MCRL2_USE_CONTROL_FLOW
      const boost::container::small_vector<std::size_t, 64>& active_cfg_vertices,
#endif
      ReportTransition report_transition = ReportTransition()
    )
    {
      bool variables_are_assigned_to_sigma=false;

#ifdef MCRL2_USE_CONTROL_FLOW
      if (m_options.use_control_flow && !m_control_flow_graphs.empty() && !cfg_allows_summand(summand.index, m_control_flow_graphs, active_cfg_vertices))
      {
        return;
      }
#endif

      auto process_transition =
        [&](const data::data_expression_list* assignments)
        {
          if (assignments)
          {
            data::add_assignments(sigma, summand.variables, *assignments);
            variables_are_assigned_to_sigma = true;
          }

          if constexpr (Stochastic)
          {
            compute_stochastic_state(s1, summand.distribution, summand.next_state, sigma, rewr, enumerator);
          }
          else
          {
            compute_state(s1, summand.next_state, sigma, rewr);
            if (!confluent_summands.empty())
            {
              s1 = find_representative(s1, confluent_summands, sigma, rewr, enumerator, id_generator);
            }
          }

          if constexpr (ReportActions)
          {
            if (m_options.rewrite_actions)
            {
              lps::multi_action a = rewrite_action(summand.multi_action, sigma, rewr);
              report_transition(a, s1);
            }
            else
            {
              report_transition(summand.multi_action, s1);
            }
          }
          else
          {
            report_transition(s1);
          }

          if (m_recursive && variables_are_assigned_to_sigma)
          {
            data::remove_assignments(sigma, summand.variables);
            variables_are_assigned_to_sigma = false;
          }
        };

      auto enumerate_projected_transitions =
        [&]()
        {
          std::vector<projected_transition> result;
          // Enumerate all satisfying valuations and compute projected transitions (action, projected state).
          // Results are stored in the projection cache.
          enumerate_solutions(
            summand, sigma, rewr, condition, enumerator,
            [&](const data::data_expression_list& e)
            {
              data::add_assignments(sigma, summand.variables, e);
              lps::multi_action a = m_options.rewrite_actions ? rewrite_action(summand.multi_action, sigma, rewr) : summand.multi_action;

              if constexpr (Stochastic)
              {
                // lps::stochastic_state q;
                // compute_stochastic_projected_state(q, summand.distribution, summand.next_state, summand.I_w, sigma, rewr, enumerator);
                // out.emplace_back(a, q);
                throw mcrl2::runtime_error("A stochastic projection cache has not been implemented yet");
              }
              else
              {
                lps::state q;
                compute_projected_state(q, summand.next_state, summand.I_w, sigma, rewr);
                result.emplace_back(a, q);
              }

              data::remove_assignments(sigma, summand.variables);
            }
          );
          return result;
        };

      auto consume_projected_transitions =
        [&](const projected_transition& t)
        {
          const lps::multi_action& a = t.first;
          const lps::state& q = t.second;

          // reconstruct full state from sigma and the projected write values q
          state s_y;
          std::size_t j = 0;
          lps::make_state(
            s_y,
            counting_iterator(0),
            m_n, // total number of process parameters
            [&](data::data_expression& r, std::size_t i)
            {
              if (j < summand.I_w.size() && summand.I_w[j] == i)
              {
                r = q[j++];
              }
              else
              {
                // reconstruct x[i] from the substitution sigma
                rewr(r, m_process_parameters[i], sigma);
              }
            }
          );

          if constexpr (!Stochastic)
          {
            state_type y = s_y;
            if (!confluent_summands.empty())
            {
              y = find_representative(y, confluent_summands, sigma, rewr, enumerator, id_generator);
            }
            if constexpr (utilities::is_applicable<ReportTransition, state_type, void>::value)
            {
              report_transition(y);
            }
            else
            {
              report_transition(a, y);
            }
          }
          else
          {
            // Convert reconstructed state to a Dirac stochastic state.
            state_type y;
            y.clear();
            y.probabilities.push_back(data::sort_real::real_one());
            y.states.push_back(s_y);
            if constexpr (utilities::is_applicable<ReportTransition, state_type, void>::value)
            {
              report_transition(y);
            }
            else
            {
              report_transition(a, y);
            }
          }
        };

      if (!m_recursive)
      {
        id_generator.clear();
      }

      if (m_options.use_projections && !summand.I_r.empty())
      {
        lps::state key;
        project(key, summand.I_r, sigma);
        utilities::shared_guard g = atermpp::detail::g_thread_term_pool().lock_shared();

        auto it = summand.projection_cache.find(key);
        if (it == summand.projection_cache.end())
        {
          g.unlock();

          std::vector<projected_transition> projected = enumerate_projected_transitions();
          it = summand.projection_cache.emplace(key, std::move(projected)).first;
        }
        else
        {
          g.unlock();
        }

        for (const projected_transition& t : it->second)
        {
          consume_projected_transitions(t);
        }
      }
      else if (summand.cache_strategy == caching::none)
      {
        // Enumerate and immediately consume all satisfying valuations (no caching, no projection).
        enumerate_solutions(
          summand, sigma, rewr, condition, enumerator,
          [&](const data::data_expression_list& e)
          {
            process_transition(&e);
          }
        );
      }
      else
      {
        summand_cache_map& cache = summand.cache_strategy == caching::global ? global_cache : summand.local_cache;
        utilities::shared_guard g = atermpp::detail::g_thread_term_pool().lock_shared();

        auto q = cache.find(detail::cheap_cache_key(sigma, summand.gamma));
        if (q == cache.end())
        {
          g.unlock();

          atermpp::term_list<data::data_expression_list> solutions;
          // Enumerate all satisfying valuations for this summand and store them in the cache.
          enumerate_solutions(
            summand, sigma, rewr, condition, enumerator,
            [&](const data::data_expression_list& e)
            {
              solutions.push_front(e);
            }
          );
          summand.compute_key(key, sigma);
          q = cache.insert({key, solutions}).first;
        }
        else
        {
          g.unlock();
        }

        for (const auto& e : static_cast<atermpp::term_list<data::data_expression_list>&>(q->second))
        {
          process_transition(e.empty() ? nullptr : &e);
        }
      }

      if (!m_recursive && variables_are_assigned_to_sigma)
      {
        data::remove_assignments(sigma, summand.variables);
      }
    }

    template <typename SummandSequence>
    std::list<transition> out_edges(const state& s, 
                                    const SummandSequence& regular_summands, 
                                    const SummandSequence& confluent_summands,
                                    data::mutable_indexed_substitution<>& sigma,
                                    data::rewriter& rewr,
                                    data::enumerator_algorithm<>& enumerator,
                                    data::enumerator_identifier_generator& id_generator
                                   )
    {
      data::data_expression condition;   // This variable is used often, and it is time consuming to declare it too often.
      state_type state_;                  // The same holds for this variable. 
      atermpp::aterm key;
      std::list<transition> transitions;
      data::add_assignments(sigma, m_process_parameters, s);
#ifdef MCRL2_USE_CONTROL_FLOW
      auto active_cfg_vertices = compute_active_cfg_vertices(sigma, m_process_parameters, m_control_flow_graphs);
#endif
      for (const explorer_summand& summand: regular_summands)
      {
        generate_transitions(
          summand,
          confluent_summands,
          sigma,
          rewr,
          condition,
          state_,
          key,
          enumerator,
          id_generator,
#ifdef MCRL2_USE_CONTROL_FLOW
          active_cfg_vertices,
#endif
          [&](const lps::multi_action& a, const state_type& s1)
          {
            if constexpr (Timed)
            {
              const data::data_expression& t = s[m_n];
              if (a.has_time() && less_equal(a.time(), t, sigma, rewr))
              {
                return;
              }
              const data::data_expression& t1 = a.has_time() ? a.time() : t;
              make_timed_state(state_, s1, t1);
              transitions.emplace_back(a, state_);
            }
            else
            {
              transitions.emplace_back(a, s1);
            }
          }
        );
      }
      return transitions;
    }

    // pre: s0 is in normal form
    template <typename SummandSequence>
    std::vector<state> generate_successors(
      const state& s0,
      const SummandSequence& summands,
      data::mutable_indexed_substitution<>& sigma,
      data::rewriter& rewr,
      data::enumerator_algorithm<>& enumerator,
      data::enumerator_identifier_generator& id_generator,
      const SummandSequence& confluent_summands = SummandSequence()
    )
    {
      data::data_expression condition; 
      state_type state_;
      atermpp::aterm key;
      std::vector<state> result;
      data::add_assignments(sigma, m_process_parameters, s0);
#ifdef MCRL2_USE_CONTROL_FLOW
      auto active_cfg_vertices = compute_active_cfg_vertices(sigma, m_process_parameters, m_control_flow_graphs);
#endif
      for (const explorer_summand& summand: summands)
      {
        generate_transitions<false>(
          summand,
          confluent_summands,
          sigma,
          rewr,
          condition,
          state_,
          key,
          enumerator,
          id_generator,
#ifdef MCRL2_USE_CONTROL_FLOW
          active_cfg_vertices,
#endif
          // [&](const lps::multi_action& /* a */, const state& s1) OLD. Calculates transitions, that are not used. 
          [&](const state& s1)
          {
            result.push_back(s1);
          }
        );
        data::remove_assignments(sigma, summand.variables);
      }
      return result;
    }

    // Add operations on reals that are needed for the exploration.
    std::set<data::function_symbol> add_real_operators(std::set<data::function_symbol> s) const
    {
      std::set<data::function_symbol> result = std::move(s);
      result.insert(data::less_equal(data::sort_real::real_()));
      result.insert(data::greater_equal(data::sort_real::real_()));
      result.insert(data::sort_real::plus(data::sort_real::real_(), data::sort_real::real_()));
      return result;
    }

    std::unique_ptr<todo_set> make_todo_set(const state& init)
    {
      switch (m_options.search_strategy)
      {
        case lps::es_breadth: return std::make_unique<breadth_first_todo_set>(init);
        case lps::es_depth: return std::make_unique<depth_first_todo_set>(init);
        case lps::es_highway: return std::make_unique<highway_todo_set>(init, m_options.highway_todo_max);
        default: throw mcrl2::runtime_error("unsupported search strategy");
      }
    }

    template <typename ForwardIterator>
    std::unique_ptr<todo_set> make_todo_set(ForwardIterator first, ForwardIterator last)
    {
      switch (m_options.search_strategy)
      {
        case lps::es_breadth: return std::make_unique<breadth_first_todo_set>(first, last);
        case lps::es_depth: return std::make_unique<depth_first_todo_set>(first, last);
        case lps::es_highway: return std::make_unique<highway_todo_set>(first, last, m_options.highway_todo_max);
        default: throw mcrl2::runtime_error("unsupported search strategy");
      }
    }

  private:
    bool is_confluent_tau(const multi_action& a)
    {
      if (a.actions().empty())
      {
        return m_options.confluence_action == "tau";
      }
      else if (a.actions().size() == 1)
      {
        return std::string(a.actions().front().label().name()) == m_options.confluence_action;
      }
      return false;
    }

  public:
    explorer(const Specification& lpsspec, const explorer_options& options_, const data::rewriter& rewr
#ifdef MCRL2_USE_CONTROL_FLOW
            , const pbes_system::pbesstategraph_options& stategraph_options = {}
            , bool compute_marking = false
#endif
            )
      : m_options(options_),
        m_global_rewr(rewr),
        m_global_enumerator(m_global_rewr, lpsspec.data(), m_global_rewr, m_global_id_generator, false),
        m_global_lpsspec(preprocess(lpsspec)),
        m_discovered(m_options.number_of_threads)
    {
#ifdef MCRL2_USE_CONTROL_FLOW
      if constexpr (!Stochastic)
      {
        lps::specification control_flow_lpsspec = lpsspec;
        m_control_flow_graphs = compute_control_flow_graphs(control_flow_lpsspec, stategraph_options, compute_marking);
        if (mCRL2logEnabled(log::debug))
        {
          for (const cf_graph& G: m_control_flow_graphs)
          {
            print_graph(G);
          }
        }
      }
#endif

      const data::variable_list& params = m_global_lpsspec.process().process_parameters();
      m_process_parameters = std::vector<data::variable>(params.begin(), params.end());
      m_n = m_process_parameters.size();
      timed_state.resize(m_n + 1);
      m_initial_state = m_global_lpsspec.initial_process().expressions();
      m_initial_distribution = initial_distribution(m_global_lpsspec);

      // Split the summands in regular and confluent summands
      const auto& lpsspec_summands = m_global_lpsspec.process().action_summands();
      for (std::size_t i = 0; i < lpsspec_summands.size(); i++)
      {
        const auto& summand = lpsspec_summands[i];
        caching cache_strategy = m_options.cached ? (m_options.global_cache ? lps::caching::global : lps::caching::local) : lps::caching::none;
        if (is_confluent_tau(summand.multi_action()))
        {
          m_confluent_summands.emplace_back(summand, i, m_global_lpsspec.process().process_parameters(), cache_strategy);
        }
        else
        {
          m_regular_summands.emplace_back(summand, i, m_global_lpsspec.process().process_parameters(), cache_strategy);
        }
      }

      // Initialize the read/write projections
      if (m_options.use_projections)
      {
        auto [R, W] = compute_read_write_patterns_separated(lpsspec);
        for (std::size_t i = 0; i < m_regular_summands.size(); ++i)
        {
          m_regular_summands[i].set_projection_attributes(R[i], W[i]);
        }
      }
    }

    ~explorer() override = default;

    // Get the initial state of the specification. 
    const data::data_expression_list& initial_state() const
    {
      return m_initial_state;
    }

    // Make the rewriter available to be used in a class that uses this explorer class.
    const data::rewriter& get_rewriter() const
    {
      return m_global_rewr;
    }

    // Compute the initial stochastic state
    void compute_initial_stochastic_state(stochastic_state& result) const
    {
      compute_stochastic_state(result, m_initial_distribution, m_initial_state, m_global_sigma, m_global_rewr, m_global_enumerator);
    }

    // Convenience overload: use internal sigma/rewriter/enumerator
    template <typename DataExpressionSequence>
    void compute_stochastic_state(stochastic_state& result,
                                  const stochastic_distribution& distribution,
                                  const DataExpressionSequence& next_state) const
    {
      compute_stochastic_state(result, distribution, next_state, m_global_sigma, m_global_rewr, m_global_enumerator);
    }

    // Utility function to obtain the outgoing transitions of the current state.
    // Should not be used concurrently. 
    std::list<transition> out_edges(const state& s)
    {
      return out_edges(s, m_regular_summands, m_confluent_summands, m_global_sigma, m_global_rewr, m_global_enumerator, m_global_id_generator);
    }

    // Utility function to obtain the outgoing transitions of the current state.
    // Only transitions for the summand with the indicated index are generated.
    // Should not be used concurrently. 
    std::list<transition> out_edges(const state& s, const std::size_t summand_index)
    {
      assert(summand_index<m_regular_summands.size());
      return out_edges(s, 
                       std::vector(1, m_regular_summands[summand_index]), 
                       m_confluent_summands, 
                       m_global_sigma, m_global_rewr, m_global_enumerator, m_global_id_generator);
    }

    // Returns the concatenation of s and [t]
    void make_timed_state(state& result, const state& s, const data::data_expression& t) const
    {
      std::copy(s.begin(), s.end(), timed_state.begin());
      timed_state.back() = t;
      lps::make_state(result, timed_state.begin(), m_n + 1);
    }

    state_type make_state(const state& s) const
    {
      if constexpr (Stochastic)
      {
        return stochastic_state(s);
      }
      else
      {
        return s;
      }
    }

    const state_type& make_state(const stochastic_state& s) const
    {
      return s;
    }

    template<typename StateType,
      typename SummandSequence,
      typename DiscoverState = utilities::skip,
      typename ExamineTransition = utilities::skip,
      typename StartState = utilities::skip,
      typename FinishState = utilities::skip,
      typename DiscoverInitialState = utilities::skip>
    void generate_state_space_thread(std::unique_ptr<todo_set>& todo,
      std::size_t thread_index,
      std::atomic<std::size_t>& number_of_active_processes,
      const SummandSequence& regular_summands,
      const SummandSequence& confluent_summands,
      indexed_set_for_states_type& discovered,
      DiscoverState discover_state,
      ExamineTransition examine_transition,
      StartState start_state,
      FinishState finish_state,
      data::rewriter thread_rewr,
      data::mutable_indexed_substitution<> thread_sigma // This is intentionally a copy.
    );

    // pre: s0 is in normal form
    template <
      typename StateType,
      typename SummandSequence,
      typename DiscoverState = utilities::skip,
      typename ExamineTransition = utilities::skip,
      typename StartState = utilities::skip,
      typename FinishState = utilities::skip,
      typename DiscoverInitialState = utilities::skip
    >
    void generate_state_space(
      bool recursive,
      const StateType& s0,
      const SummandSequence& regular_summands,
      const SummandSequence& confluent_summands,
      indexed_set_for_states_type& discovered,
      DiscoverState discover_state = DiscoverState(),
      ExamineTransition examine_transition = ExamineTransition(),
      StartState start_state = StartState(),
      FinishState finish_state = FinishState(),
      [[maybe_unused]]
      DiscoverInitialState discover_initial_state = DiscoverInitialState()
    );

    /// \brief Generates the state space, and reports all discovered states and transitions by means of callback
    /// functions.
    /// \param discover_state Is invoked when a state is encountered for the first time.
    /// \param examine_transition Is invoked on every transition.
    /// \param start_state Is invoked on a state right before its outgoing transitions are being explored.
    /// \param finish_state Is invoked on a state after all of its outgoing transitions have been explored.
    template <
      typename DiscoverState = utilities::skip,
      typename ExamineTransition = utilities::skip,
      typename StartState = utilities::skip,
      typename FinishState = utilities::skip,
      typename DiscoverInitialState = utilities::skip
    >
    void generate_state_space(
      bool recursive,
      DiscoverState discover_state = DiscoverState(),
      ExamineTransition examine_transition = ExamineTransition(),
      StartState start_state = StartState(),
      FinishState finish_state = FinishState(),
      DiscoverInitialState discover_initial_state = DiscoverInitialState()
    );

    /// \brief Generates outgoing transitions for a given state.
    std::vector<std::pair<lps::multi_action, state_type>> generate_transitions(
                   const state& d0,
                   data::mutable_indexed_substitution<>& sigma,
                   data::rewriter& rewr,
                   data::enumerator_algorithm<>& enumerator,
                   data::enumerator_identifier_generator& id_generator)
    {
      data::data_expression_list process_parameter_undo = process_parameter_values(sigma);
      std::vector<std::pair<lps::multi_action, state_type>> result;
      data::add_assignments(sigma, m_process_parameters, d0);
      data::data_expression condition;
      atermpp::aterm key;
      state_type state;
#ifdef MCRL2_USE_CONTROL_FLOW
      auto active_cfg_vertices = compute_active_cfg_vertices(sigma, m_process_parameters, m_control_flow_graphs);
#endif
      for (const explorer_summand& summand: m_regular_summands)
      {
        generate_transitions(
          summand,
          m_confluent_summands,
          sigma,
          rewr,
          condition,
          state,
          key,
          enumerator,
          id_generator,
#ifdef MCRL2_USE_CONTROL_FLOW
          active_cfg_vertices,
#endif
          [&](const lps::multi_action& a, const state_type& d1)
          {
            result.emplace_back(lps::multi_action(a.actions(), a.time()), d1);
          }
        );
      }
      set_process_parameter_values(process_parameter_undo, sigma);
      return result;
    }

    /// \brief Generates outgoing transitions for a given state, using the global substitution, rewriter, enumerator and id_generator.
    /// \details This function is not suitable to be used in parallel threads, but can only be used for pre or post processing. 
    std::vector<std::pair<lps::multi_action, state_type>> generate_transitions(
                   const state& d0)
    {
      assert(m_options.number_of_threads==1); // A global rewriter is invoked, and this can only happen in a single threaded setting. 
      return generate_transitions(d0, m_global_sigma, m_global_rewr, m_global_enumerator, m_global_id_generator);
    }

    /// \brief Generates outgoing transitions for a given state.
    std::vector<std::pair<lps::multi_action, state>> generate_transitions(
              const data::data_expression_list& init,
              data::mutable_indexed_substitution<>& sigma,
              data::rewriter& rewr)
    {
      state d0;
      compute_state(d0,init,sigma,rewr);
      return generate_transitions(d0);
    }

    /// \brief Generates outgoing transitions for a given state, reachable via the summand with index i.
    std::vector<std::pair<lps::multi_action, state_type>> generate_transitions(
                const data::data_expression_list& init, 
                std::size_t i,
                data::mutable_indexed_substitution<>& sigma,
                data::rewriter& rewr,
                data::data_expression& condition,
                state_type& d0,
                data::enumerator_algorithm<>& enumerator,
                data::enumerator_identifier_generator& id_generator)
    {
      data::data_expression_list process_parameter_undo = process_parameter_values(sigma);
      compute_state(d0,init,sigma,rewr);
      std::vector<std::pair<lps::multi_action, state_type>> result;
      data::add_assignments(sigma, m_process_parameters, d0);
#ifdef MCRL2_USE_CONTROL_FLOW
      auto active_cfg_vertices = compute_active_cfg_vertices(sigma, m_process_parameters, m_control_flow_graphs);
#endif
      generate_transitions(
        m_regular_summands[i],
        m_confluent_summands,
        sigma,
        rewr,
        condition,
        d0,
        enumerator,
        id_generator,
#ifdef MCRL2_USE_CONTROL_FLOW
        active_cfg_vertices,
#endif
        [&](const lps::multi_action& a, const state_type& d1)
        {
          result.emplace_back(lps::multi_action(a), d1);
        }
      );
      data::remove_assignments(sigma, m_regular_summands[i].variables);
      set_process_parameter_values(process_parameter_undo, sigma);
      return result;
    }

    // --- Recursive DFS (core) ---------------------------------------------
    // pre: s0 is in normal form
    // N.B. Does not support stochastic specifications!
    template <
      typename SummandSequence,
      typename DiscoverState = utilities::skip,
      typename ExamineTransition = utilities::skip,
      typename TreeEdge = utilities::skip,
      typename BackEdge = utilities::skip,
      typename ForwardOrCrossEdge = utilities::skip,
      typename FinishState = utilities::skip
    >
    void generate_state_space_dfs_recursive(
      const state& s0,
      std::unordered_set<state> gray,
      std::unordered_set<state>& discovered,
      const SummandSequence& regular_summands,
      const SummandSequence& confluent_summands,
      DiscoverState discover_state = DiscoverState(),
      ExamineTransition examine_transition = ExamineTransition(),
      TreeEdge tree_edge = TreeEdge(),
      BackEdge back_edge = BackEdge(),
      ForwardOrCrossEdge forward_or_cross_edge = ForwardOrCrossEdge(),
      FinishState finish_state = FinishState()
    );

    // --- Recursive DFS (entry point) ---------------------------------------
    template <
      typename DiscoverState = utilities::skip,
      typename ExamineTransition = utilities::skip,
      typename TreeEdge = utilities::skip,
      typename BackEdge = utilities::skip,
      typename ForwardOrCrossEdge = utilities::skip,
      typename FinishState = utilities::skip
    >
    void generate_state_space_dfs_recursive(
      bool recursive,
      DiscoverState discover_state = DiscoverState(),
      ExamineTransition examine_transition = ExamineTransition(),
      TreeEdge tree_edge = TreeEdge(),
      BackEdge back_edge = BackEdge(),
      ForwardOrCrossEdge forward_or_cross_edge = ForwardOrCrossEdge(),
      FinishState finish_state = FinishState()
    );

    // --- Iterative DFS (core) ----------------------------------------------
    // pre: s0 is in normal form
    // N.B. Does not support stochastic specifications!
    template <
      typename SummandSequence,
      typename DiscoverState = utilities::skip,
      typename ExamineTransition = utilities::skip,
      typename TreeEdge = utilities::skip,
      typename BackEdge = utilities::skip,
      typename ForwardOrCrossEdge = utilities::skip,
      typename FinishState = utilities::skip
    >
    void generate_state_space_dfs_iterative(
      const state& s0,
      std::unordered_set<state>& discovered,
      const SummandSequence& regular_summands,
      const SummandSequence& confluent_summands,
      DiscoverState discover_state = DiscoverState(),
      ExamineTransition examine_transition = ExamineTransition(),
      TreeEdge tree_edge = TreeEdge(),
      BackEdge back_edge = BackEdge(),
      ForwardOrCrossEdge forward_or_cross_edge = ForwardOrCrossEdge(),
      FinishState finish_state = FinishState()
    );

    // --- Iterative DFS (entry point) ---------------------------------------
    template <
      typename DiscoverState = utilities::skip,
      typename ExamineTransition = utilities::skip,
      typename TreeEdge = utilities::skip,
      typename BackEdge = utilities::skip,
      typename ForwardOrCrossEdge = utilities::skip,
      typename FinishState = utilities::skip
    >
    void generate_state_space_dfs_iterative(
      bool recursive,
      DiscoverState discover_state = DiscoverState(),
      ExamineTransition examine_transition = ExamineTransition(),
      TreeEdge tree_edge = TreeEdge(),
      BackEdge back_edge = BackEdge(),
      ForwardOrCrossEdge forward_or_cross_edge = ForwardOrCrossEdge(),
      FinishState finish_state = FinishState()
    );

    /// \brief Abort the state space generation
    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    void abort() override
    {
      m_must_abort = true;
    }

    /// \brief Returns a mapping containing all discovered states.
    const indexed_set_for_states_type& state_map() const
    {
      return m_discovered;
    }

    const std::vector<explorer_summand>& regular_summands() const
    {
      return m_regular_summands;
    }

    const std::vector<explorer_summand>& confluent_summands() const
    {
      return m_confluent_summands;
    }

    const std::vector<data::variable>& process_parameters() const
    {
      return m_process_parameters;
    }

    data::data_expression_list process_parameter_values(data::mutable_indexed_substitution<>& sigma) const
    {
      return data::data_expression_list{m_process_parameters.begin(), m_process_parameters.end(), [&](const data::variable& x) { return sigma(x); }};
    }

    /// \brief Process parameter values for use in a single thread. 
    data::data_expression_list process_parameter_values() const
    {
      assert(m_options.number_of_threads==1); // Using a global sigma is not thread safe. 
      return process_parameter_values(m_global_sigma);
    }
    void set_process_parameter_values(const data::data_expression_list& values, data::mutable_indexed_substitution<>& sigma)
    {
      data::add_assignments(sigma, m_process_parameters, values);
    }

    void set_process_parameter_values(const data::data_expression_list& values)
    {
       assert(m_options.number_of_threads==1); // Using a global sigma is not thread safe. 
       set_process_parameter_values(values, m_global_sigma);
    }
};

} // namespace mcrl2::lps

#endif // MCRL2_LPS_EXPLORER_H

// implementation
#include "mcrl2/lps/explorer_bfs.h"
#include "mcrl2/lps/explorer_dfs.h"
