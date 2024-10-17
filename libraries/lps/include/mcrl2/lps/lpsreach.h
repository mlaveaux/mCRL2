// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/lpsreach.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_LPSREACH_H
#define MCRL2_LPS_LPSREACH_H

#ifdef MCRL2_ENABLE_SYLVAN

#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/detail/instantiate_global_variables.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/one_point_rule_rewrite.h"
#include "mcrl2/lps/order_summand_variables.h"
#include "mcrl2/lps/replace_constants_by_variables.h"
#include "mcrl2/lps/resolve_name_clashes.h"
#include "mcrl2/lps/symbolic_lts.h"
#include "mcrl2/symbolic/ordering.h"
#include "mcrl2/symbolic/print.h"
#include "mcrl2/symbolic/symbolic_reachability.h"
#include "mcrl2/utilities/parse_numbers.h"
#include "mcrl2/utilities/stack_array.h"
#include "mcrl2/utilities/stopwatch.h"

#include <sylvan_ldd.hpp>

#include <chrono>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>

namespace mcrl2::lps {

class lpsreach_algorithm
{
    using ldd = sylvan::ldds::ldd;
    using enumerator_element = data::enumerator_list_element_with_substitution<>;

    template <typename Context, bool ActionLabel>
    friend void symbolic::learn_successors_callback(WorkerP*, Task*, std::uint32_t* v, std::size_t n, void* context);

  protected:
    const symbolic::symbolic_reachability_options& m_options;
    data::rewriter m_rewr;
    data::mutable_indexed_substitution<> m_sigma;
    data::enumerator_identifier_generator m_id_generator;
    data::enumerator_algorithm<> m_enumerator;
    std::vector<boost::dynamic_bitset<>> m_summand_patterns;
    std::vector<boost::dynamic_bitset<>> m_group_patterns;
    std::vector<std::size_t> m_variable_order;
    symbolic_lts m_lts;
    
    /// \brief Rewrites all arguments of the given action.
    template<typename Rewriter, typename Substitution>
    lps::multi_action rewrite_action(const lps::multi_action& a, const Rewriter& rewr, const Substitution& sigma)
    {
      const process::action_list& actions = a.actions();
      const data::data_expression& time = a.time();
      return
        lps::multi_action(
          process::action_list(
            actions.begin(),
            actions.end(),
            [&](const process::action& a)
            {
              const auto& args = a.arguments();
              return process::action(a.label(), data::data_expression_list(args.begin(), args.end(), [&](const data::data_expression& x) { return rewr(x, sigma); }));
            }
          ),
          a.has_time() ? rewr(time, sigma) : time
        );
    }

    // R.L := R.L U {(x,y) in R | x in X}
    void learn_successors(std::size_t i, symbolic::summand_group& R, const ldd& X)
    {
      mCRL2log(log::trace) << "learn successors of summand group " << i << " for X = " << print_states(m_lts.data_index, X, R.read) << std::endl;

      using namespace sylvan::ldds;
      std::pair<lpsreach_algorithm&, symbolic::summand_group&> context{*this, R};
      sat_all_nopar(X, symbolic::learn_successors_callback<std::pair<lpsreach_algorithm&, lps_summand_group&>, true>, &context);
    }

    template <typename Specification>
    Specification preprocess(const Specification& lpsspec)
    {
      Specification result = lpsspec;
      lps::detail::instantiate_global_variables(result);
      lps::order_summand_variables(result);
      resolve_summand_variable_name_clashes(result); // N.B. This is a required preprocessing step.
      if (m_options.one_point_rule_rewrite)
      {
        one_point_rule_rewrite(result);
      }
      if (m_options.replace_constants_by_variables)
      {
        replace_constants_by_variables(result, m_rewr, m_sigma);
      }

      return result;
    }

    std::string print_size(const sylvan::ldds::ldd& L)
    {
      return symbolic::print_size(L, m_options.print_exact, m_options.print_nodesize);
    }

  public:
    lpsreach_algorithm(const lps::specification& lpsspec, const symbolic::symbolic_reachability_options& options_)
      : m_options(options_),
        m_rewr(symbolic::construct_rewriter(lpsspec.data(), m_options.rewrite_strategy, lps::find_function_symbols(lpsspec), m_options.remove_unused_rewrite_rules)),
        m_enumerator(m_rewr, lpsspec.data(), m_rewr, m_id_generator, false)
    {
      using utilities::detail::as_vector;

      lps::specification lpsspec_ = preprocess(lpsspec);
      m_lts.process_parameters = lpsspec_.process().process_parameters();

      // Rewrite the initial expressions to normal form,
      std::vector<data::data_expression> initial_values;
      for (const data::data_expression& expression : lpsspec_.initial_process().expressions())
      {
        initial_values.push_back(m_rewr(expression));
      }

      data::data_expression_list initial_state(initial_values.begin(), initial_values.end());

      m_summand_patterns = compute_read_write_patterns(lpsspec_);
      symbolic::adjust_read_write_patterns(m_summand_patterns, m_options);

      m_variable_order = symbolic::compute_variable_order(m_options.variable_order, m_lts.process_parameters.size(), m_summand_patterns);
      mCRL2log(log::debug) << "variable order = " << core::detail::print_list(m_variable_order) << std::endl;
      m_summand_patterns = symbolic::reorder_read_write_patterns(m_summand_patterns, m_variable_order);
      mCRL2log(log::debug) << symbolic::print_read_write_patterns(m_summand_patterns);

      m_lts.process_parameters = symbolic::permute_copy(m_lts.process_parameters, m_variable_order);
      for (const data::variable& param : m_lts.process_parameters)
      {
        m_lts.data_index.push_back(symbolic::data_expression_index(param.sort()));
      }

      m_lts.initial_state = symbolic::state2ldd(symbolic::permute_copy(initial_state, m_variable_order), m_lts.data_index);
      mCRL2log(log::debug) << "process parameters = " << core::detail::print_list(m_lts.process_parameters) << std::endl;

      std::vector<std::set<std::size_t>> groups = symbolic::compute_summand_groups(m_options.summand_groups, m_summand_patterns);
      for (const auto& group: groups)
      {
        mCRL2log(log::debug) << "group " << core::detail::print_set(group) << std::endl;
      }
      m_group_patterns = symbolic::compute_summand_group_patterns(m_summand_patterns, groups);
      for (std::size_t j = 0; j < m_group_patterns.size(); j++)
      {
        m_lts.summand_groups.emplace_back(lpsspec_, m_lts.process_parameters, groups[j], m_group_patterns[j], m_summand_patterns, m_variable_order);
      }

      for (std::size_t i = 0; i < m_lts.summand_groups.size(); i++)
      {
        mCRL2log(log::debug) << "=== summand group " << i << " ===\n" << m_lts.summand_groups[i] << std::endl;
      }
    }

    /// \brief Computes relprod(U, group).
    ldd relprod_impl(const ldd& U, const lps_summand_group& group, std::size_t i)
    {
      if (m_options.no_relprod)
      {
        ldd z = symbolic::alternative_relprod(U, group);
        mCRL2log(log::trace) << "relprod(" << i << ", todo) = " << print_states(m_lts.data_index, z) << std::endl;
        return z;
      }
      else
      {
        ldd z = relprod(U, group.L, group.Ir);
        mCRL2log(log::trace) << "relprod(" << i << ", todo) = " << print_states(m_lts.data_index, z) << std::endl;
        return z;
      }
    }

    /// \brief Perform a single breadth first step.
    /// \returns The tuple <visited, todo, deadlocks>
    std::tuple<ldd, ldd, ldd> step(const ldd& visited, const ldd& todo, bool learn_transitions = true, bool detect_deadlocks = false)
    {
      using namespace sylvan::ldds;
      auto& R = m_lts.summand_groups;

      ldd todo1 = empty_set();
      ldd potential_deadlocks = detect_deadlocks ? todo : empty_set();

      if (!m_options.saturation)
      {
        // regular and chaining.
        todo1 = m_options.chaining ? todo : empty_set();

        for (std::size_t i = 0; i < R.size(); i++)
        {
          if (learn_transitions)
          {
            ldd proj = project(m_options.chaining ? todo1 : todo, R[i].Ip);
            learn_successors(i, R[i], m_options.cached ? minus(proj, R[i].Ldomain) : proj);

            mCRL2log(log::trace) << "L =\n" << print_relation(m_lts.data_index, R[i].L, R[i].read, R[i].write) << std::endl;
          }

          todo1 = union_(todo1, relprod_impl(m_options.chaining ? todo1 : todo, R[i], i));

          if (detect_deadlocks)
          {
            potential_deadlocks = minus(potential_deadlocks, relprev(todo1, R[i].L, R[i].Ir, potential_deadlocks));
          }
        }
      }
      else
      {
        // saturation and chaining
        todo1 = todo;
        ldd todo1_old; // the old todo set.

        for (std::size_t i = 0; i < R.size(); i++)
        {
          if (learn_transitions)
          {
            ldd proj = project(todo1, R[i].Ip);
            learn_successors(i, R[i], m_options.cached ? minus(proj, R[i].Ldomain) : proj);

            mCRL2log(log::trace) << "L =\n" << print_relation(m_lts.data_index, R[i].L, R[i].read, R[i].write) << std::endl;
          }

          // Apply one transition relation repeatedly.
          do
          {
            todo1_old = todo1;
            todo1 = union_(todo1, relprod_impl(todo1, R[i], i));
          }
          while (todo1 != todo1_old);

          if (detect_deadlocks)
          {
            potential_deadlocks = minus(potential_deadlocks, relprev(todo1, R[i].L, R[i].Ir, potential_deadlocks));
          }

          // Apply all previously learned transition relations repeatedly.
          if (m_options.chaining)
          {
            do
            {
              todo1_old = todo1;
              for (std::size_t j = 0; j <= i; j++)
              {
                todo1 = union_(todo1, relprod_impl(todo1, R[j], j));
              }
            }
            while (todo1 != todo1_old);
          }
        }
      }

      // after all transition groups are applied the remaining potential deadlocks are actual deadlocks.
      return std::make_tuple(union_(visited, todo), minus(todo1, visited), potential_deadlocks);
    }

    ldd run()
    {
      using namespace sylvan::ldds;
      auto& R = m_lts.summand_groups;
      std::size_t iteration_count = 0;

      mCRL2log(log::trace) << "initial state = " << core::detail::print_list(m_lts.initial_state) << std::endl;

      auto start = std::chrono::steady_clock::now();
      ldd x = m_lts.initial_state;
      std::chrono::duration<double> elapsed_seconds = std::chrono::steady_clock::now() - start;
      ldd visited = empty_set();
      ldd todo = x;
      ldd deadlocks = empty_set();

      while (todo != empty_set() && (m_options.max_iterations == 0 || iteration_count < m_options.max_iterations))
      {
        stopwatch loop_start;
        iteration_count++;
        mCRL2log(log::trace) << "--- iteration " << iteration_count << " ---" << std::endl;
        mCRL2log(log::trace) << "todo = " << print_states(m_lts.data_index, todo) << std::endl;

        std::tie(visited, todo, deadlocks) = step(visited, todo, true, m_options.detect_deadlocks);

        mCRL2log(log::verbose) << "explored " << std::setw(12) << print_size(union_(visited, todo)) << " states after "
                               << std::setw(3) << iteration_count << " iterations (time = " << std::setprecision(2)
                               << std::fixed << loop_start.seconds() << "s)" << std::endl;
        if (m_options.detect_deadlocks)
        {
          mCRL2log(log::verbose) << "found " << std::setw(12) << print_size(deadlocks) << " deadlocks" << std::endl;
        }

        sylvan::sylvan_stats_report(stderr);
      }

      elapsed_seconds = std::chrono::steady_clock::now() - start;
      std::cout << "number of states = " << print_size(visited) << " (time = " << std::setprecision(2) << std::fixed << elapsed_seconds.count() << "s)" << std::endl;
      mCRL2log(log::verbose) << "used variable order = " << core::detail::print_list(m_variable_order) << std::endl;

      double total_time = 0.0;
      for (std::size_t i = 0; i < R.size(); i++)
      {
        mCRL2log(log::verbose) << "group " << std::setw(4) << i << " contains " << std::setw(10) << print_size(R[i].L) << " transitions (learn time = "
                               << std::setw(5) << std::setprecision(2) << std::fixed << R[i].learn_time << "s with " << std::setw(9) << R[i].learn_calls 
                               << " calls, cached " << print_size(R[i].Ldomain) << " values" << ")" << std::endl;

        total_time += R[i].learn_time;
      }
      mCRL2log(log::verbose) << "learning transitions took " << total_time << "s" << std::endl;

      std::size_t i = 0;
      for (const auto& param : m_lts.process_parameters)
      {
        auto& table = m_lts.data_index[i];

        mCRL2log(log::verbose) << "Parameter " << i << " (" << param << ")" << " has " << table.size() << " values."<< std::endl;
        for (const auto& data : table)
        {
          mCRL2log(log::debug) << table.index(data) << ": " << data << std::endl;
        }

        ++i;
      }

      mCRL2log(log::verbose) << "There are " << m_lts.action_index.size() << " action labels" << std::endl;
      for (const auto& action : m_lts.action_index)
      {
          mCRL2log(log::debug) << m_lts.action_index.index(action) << ": " << action << std::endl;
      }

      m_lts.states = visited;
      return visited;
    }

    const std::vector<boost::dynamic_bitset<>>& read_write_patterns() const
    {
      return m_summand_patterns;
    }

    const std::vector<boost::dynamic_bitset<>>& read_write_group_patterns() const
    {
      return m_group_patterns;
    }
    
    std::vector<symbolic::data_expression_index>& data_index()
    {
      return m_lts.data_index;
    }

    utilities::indexed_set<lps::multi_action>& action_index()
    {
      return m_lts.action_index;
    }

    const symbolic_lts& get_symbolic_lts()
    {
      return m_lts;
    }
};

} // namespace mcrl2::lps

#endif // MCRL2_ENABLE_SYLVAN

#endif // MCRL2_LPS_LPSREACH_H
