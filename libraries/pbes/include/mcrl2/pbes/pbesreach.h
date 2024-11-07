// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_PBES_PBESREACH_H
#define MCRL2_PBES_PBESREACH_H

#ifdef MCRL2_ENABLE_SYLVAN

#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/utilities/detail/container_utility.h"
#include "mcrl2/data/merge_data_specifications.h"
#include "mcrl2/pbes/detail/instantiate_global_variables.h"
#include "mcrl2/pbes/detail/pbes_io.h"
#include "mcrl2/pbes/normalize.h"
#include "mcrl2/pbes/pbes_summand_group.h"
#include "mcrl2/pbes/replace_constants_by_variables.h"
#include "mcrl2/pbes/resolve_name_clashes.h"
#include "mcrl2/pbes/rewriters/one_point_rule_rewriter.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/unify_parameters.h"
#include "mcrl2/symbolic/symbolic_reachability.h"
#include "mcrl2/symbolic/print.h"
#include "mcrl2/utilities/stopwatch.h"
#include "mcrl2/utilities/text_utility.h"

namespace mcrl2::pbes_system {

// Returns a data specification containing a structured sort with the names of the propositional variables
// in the PBES as elements.
inline
data::data_specification construct_propositional_variable_data_specification(const pbes_system::srf_pbes& pbesspec, const std::string& sort_name)
{
  // TODO: it should be possible to add a structured sort to pbesspec.data() directly, but I don't know how
  std::vector<std::string> names;
  for (const auto& equation: pbesspec.equations())
  {
    names.push_back(equation.variable().name());
  }
  std::string text = "sort " + sort_name + " = struct " +   utilities::string_join(names, " | ") + ";";
  return data::parse_data_specification(text);
}

struct symbolic_reachability_options: public symbolic::symbolic_reachability_options
{
  bool make_total = false;
  bool reset_parameters = false;
  bool aggressive = false;
  std::size_t solve_strategy = 0;
  std::size_t split_conditions = 0;
  std::string srf;
};

inline
std::ostream& operator<<(std::ostream& out, const symbolic_reachability_options& options)
{
  out << static_cast<symbolic::symbolic_reachability_options>(options);
  out << "solve_strategy = " << options.solve_strategy << std::endl;
  out << "split_conditions = " << options.split_conditions << std::endl;
  out << "total = " << std::boolalpha << options.make_total << std::endl;
  return out;
}

pbes_system::srf_pbes split_conditions(const pbes_system::srf_pbes& pbes, std::size_t granularity)
{
  mCRL2log(log::debug) << "splitting conditions" << std::endl;

  // Find existing identifiers.
  data::set_identifier_generator id_generator;
  for (const srf_equation& equation : pbes.equations())
  {
    id_generator.add_identifier(equation.variable().name());
  }

  // Determine the Xtrue equation.
  pbes_system::propositional_variable Xtrue = pbes.equations()[pbes.equations().size()-2].variable();
  pbes_system::propositional_variable Xfalse = pbes.equations()[pbes.equations().size()-1].variable();

  pbes_system::srf_pbes result = pbes;
  std::vector<srf_equation> added_equations; // These equations are added at the end of the pbes.
  for (srf_equation& equation : result.equations())
  {
    std::vector<srf_summand> split_summands; // The updated summands.
    for (const srf_summand& summand : equation.summands())
    {
      mCRL2log(log::debug) << "splitting summand " << summand << std::endl;

      // Heuristics to determine when to split conjunctive conditions.
      bool should_split = summand.parameters().empty() && granularity > 1; // && find_free_variables(summand.condition()).size() >= 4;

      if (data::sort_bool::is_or_application(summand.condition()))
      {
        // For disjunctive conditions we can introduce one summand per clause.
        for (const data::data_expression& clause : split_disjunction(summand.condition()))
        {
          split_summands.emplace_back(summand.parameters(), clause, summand.variable());
          mCRL2log(log::debug) << "Added summand " << split_summands.back() << std::endl;
        }
      }
      else if (should_split && data::sort_bool::is_and_application(summand.condition()))
      {
        // The summand is simple if there is only a dependency on true or false.
        bool simple = granularity == 3 || summand.variable().name() == Xtrue.name() || summand.variable().name() == Xfalse.name();

        std::vector<srf_summand> split_summands_inner; // The summands for the added equation.
        for (const data::data_expression& clause : split_conjunction(summand.condition()))
        {
          if (simple)
          {
            // For conjunctive equations add !condition => Xfalse, and !condition && Xtrue otherwise.
            split_summands_inner.emplace_back(data::variable_list(),
                                              data::lazy::not_(clause),
                                              !equation.is_conjunctive()
                                              ? propositional_variable_instantiation(Xtrue.name(), {}) :
                                                propositional_variable_instantiation(Xfalse.name(), {})
                                              );
          }
          else
          {
            // Add a new equation per clause.
            const propositional_variable& Y = equation.variable();
            propositional_variable Y1(id_generator(Y.name()), Y.parameters());

            split_summands_inner.emplace_back(data::variable_list(), true_(), propositional_variable_instantiation(Y1.name(), data::make_data_expression_list(Y1.parameters())));
            std::vector<srf_summand> summands;
            summands.emplace_back(data::variable_list(), clause, summand.variable());
            added_equations.emplace_back(equation.symbol(), Y1, summands, !equation.is_conjunctive());
            mCRL2log(log::debug) << "Added equation " << added_equations.back() << std::endl;
          }
          mCRL2log(log::debug) << "Added summand " << split_summands_inner.back() << std::endl;
        }

        if (simple)
        {
          split_summands_inner.emplace_back(data::variable_list(), true_(), summand.variable());
        }

        if (equation.summands().size() == 1)
        {
          // Change the current equation.
          split_summands = split_summands_inner;
          equation.is_conjunctive() = !equation.is_conjunctive();
          mCRL2log(log::debug) << "Changed equation type (conjunctive or disjunctive)" << std::endl;
        }
        else
        {
          // Add a new equation.
          const propositional_variable& Y = equation.variable();
          propositional_variable Y1(id_generator(Y.name()), Y.parameters());

          split_summands.emplace_back(data::variable_list(), true_(), propositional_variable_instantiation(Y1.name(), data::make_data_expression_list(Y1.parameters())));
          added_equations.emplace_back(equation.symbol(), Y1, split_summands_inner, !equation.is_conjunctive());
          mCRL2log(log::debug) << "Added equation " << added_equations.back() << std::endl;
        }
      }
      else
      {
        // Do nothing.
        split_summands.emplace_back(summand);
      }
    }

    equation.summands() = split_summands;
  }

  // The last two equations must be Xfalse and Xtrue.
  result.equations().insert(result.equations().end()-2, added_equations.begin(), added_equations.end());

  return result;
}

// Store information per lace worker.
struct per_worker_information
{
  data::mutable_indexed_substitution<> m_sigma;
  data::rewriter m_rewr;
};

class pbesreach_algorithm
{
    using enumerator_element = data::enumerator_list_element_with_substitution<>;

    template <typename Context, bool ActionLabel>
    friend void symbolic::learn_successors_callback(WorkerP*, Task*, std::uint32_t* v, std::size_t n, void* context);

  protected:
    using ldd = sylvan::ldds::ldd;
    const symbolic_reachability_options& m_options;
    pbes_system::srf_pbes m_pbes;
    data::rewriter m_rewr;
    data::mutable_indexed_substitution<> m_sigma;
    data::enumerator_identifier_generator m_id_generator;
    data::enumerator_algorithm<> m_enumerator;
    data::variable_list m_process_parameters;
    std::size_t m_n;
    std::vector<symbolic::data_expression_index> m_data_index;
    std::vector<pbes_summand_group> m_summand_groups;
    data::data_expression_list m_initial_state;
    std::vector<boost::dynamic_bitset<>> m_summand_patterns;
    std::vector<boost::dynamic_bitset<>> m_group_patterns;
    std::vector<std::size_t> m_variable_order;

    ldd m_visited;
    ldd m_todo;
    ldd m_deadlocks;
    ldd m_initial_vertex;

    /// \brief Updates R.L := R.L U {(x,y) in R | x in X}
    void learn_successors(std::size_t i, pbes_summand_group& R, const ldd& X)
    {
      mCRL2log(log::trace) << "learn successors of summand group " << i << " for X = " << print_states(m_data_index, X, R.read) << std::endl;

      using namespace sylvan::ldds;
      std::pair<pbesreach_algorithm&, pbes_summand_group&> context{*this, R};
      sat_all_nopar(X, symbolic::learn_successors_callback<std::pair<pbesreach_algorithm&, pbes_summand_group&>, false>, &context);
    }

    pbes_system::srf_pbes preprocess(pbes_system::pbes pbesspec, bool make_total)
    {
      pbes_system::detail::instantiate_global_variables(pbesspec);
      normalize(pbesspec);

      if (m_options.one_point_rule_rewrite)
      {
        pbes_system::one_point_rule_rewriter R;
        pbes_system::replace_pbes_expressions(pbesspec, R, false);
      }

      if (m_options.replace_constants_by_variables)
      {
        pbes_system::replace_constants_by_variables(pbesspec, m_rewr, m_sigma);
      }

      pbes_system::srf_pbes result = pbes2srf(pbesspec);
      if (m_options.split_conditions > 0)
      {
        result = split_conditions(result, m_options.split_conditions);
      }

      if (make_total)
      {
        result.make_total();
      }

      unify_parameters(result, m_options.reset_parameters);
      pbes_system::resolve_summand_variable_name_clashes(result, result.equations().front().variable().parameters()); // N.B. This is a required preprocessing step.

      // add a sort for the propositional variable names
      data::data_specification propvar_dataspec = construct_propositional_variable_data_specification(result, "PropositionalVariable");
      result.data() = data::merge_data_specifications(result.data(), propvar_dataspec);
      mCRL2log(log::trace) << "--- srf pbes ---\n" << result.to_pbes() << std::endl;

      return result;
    }

    std::string print_size(const sylvan::ldds::ldd& L)
    {
      return symbolic::print_size(L, m_options.print_exact, m_options.print_nodesize);
    }

  public:
    pbesreach_algorithm(const pbes_system::pbes& pbesspec, const symbolic_reachability_options& options_)
      : m_options(options_),
        m_pbes(preprocess(pbesspec, options_.make_total)),
        m_rewr(symbolic::construct_rewriter(m_pbes.data(), m_options.rewrite_strategy, pbes_system::find_function_symbols(pbesspec), m_options.remove_unused_rewrite_rules)),
        m_enumerator(m_rewr, m_pbes.data(), m_rewr, m_id_generator, false)
    {
      if (!m_options.srf.empty())
      {
        detail::save_pbes(m_pbes.to_pbes(), m_options.srf);
      }

      data::basic_sort propvar_sort("PropositionalVariable"); // todo: choose a unique name
      std::unordered_map<core::identifier_string, data::data_expression> propvar_map;
      for (const auto& equation: m_pbes.equations())
      {
        propvar_map[equation.variable().name()] = data::function_symbol(equation.variable().name(), propvar_sort);
      }

      m_process_parameters = m_pbes.equations().front().variable().parameters();
      m_process_parameters.push_front(data::variable("propvar", propvar_sort)); // todo: choose a unique name
      m_n = m_process_parameters.size();

      // Rewrite the initial expressions to normal form,
      std::vector<data::data_expression> initial_values;
      for (const data::data_expression& expression : make_state(m_pbes.initial_state(), propvar_map))
      {
        initial_values.push_back(m_rewr(expression));
      }

      m_initial_state = data::data_expression_list(initial_values.begin(), initial_values.end());

      m_summand_patterns = compute_read_write_patterns(m_pbes, m_process_parameters);
      symbolic::adjust_read_write_patterns(m_summand_patterns, m_options);

      m_variable_order = symbolic::compute_variable_order(m_options.variable_order, m_process_parameters.size(), m_summand_patterns, true);
      assert(m_variable_order[0] == 0); // It is required that the propositional variable name stays up front
      mCRL2log(log::debug) << "variable order = " << core::detail::print_list(m_variable_order) << std::endl;
      m_summand_patterns = symbolic::reorder_read_write_patterns(m_summand_patterns, m_variable_order);
      mCRL2log(log::debug) << symbolic::print_read_write_patterns(m_summand_patterns);

      m_process_parameters = symbolic::permute_copy(m_process_parameters, m_variable_order);
      m_initial_state = symbolic::permute_copy(m_initial_state, m_variable_order);
      mCRL2log(log::debug) << "process parameters = " << core::detail::print_list(m_process_parameters) << std::endl;

      std::vector<std::set<std::size_t>> groups = symbolic::compute_summand_groups(m_options.summand_groups, m_summand_patterns);
      for (const auto& group: groups)
      {
        mCRL2log(log::debug) << "group " << core::detail::print_set(group) << std::endl;
      }
      m_group_patterns = symbolic::compute_summand_group_patterns(m_summand_patterns, groups);
      for (std::size_t j = 0; j < m_group_patterns.size(); j++)
      {
        m_summand_groups.emplace_back(m_pbes, m_process_parameters, propvar_map, groups[j], m_group_patterns[j], m_summand_patterns, m_variable_order);
      }

      for (std::size_t i = 0; i < m_summand_groups.size(); i++)
      {
        mCRL2log(log::debug) << "=== summand group " << i << " ===\n" << m_summand_groups[i] << std::endl;
      }

      for (const data::variable& param: m_process_parameters)
      {
        m_data_index.push_back(symbolic::data_expression_index(param.sort()));
      }
    }

    virtual ~pbesreach_algorithm() {}

    ldd initial_state()
    {
      return symbolic::state2ldd(m_initial_state, m_data_index);
    }

    /// \returns The set of deadlock states.
    ldd deadlocks()
    {
      return m_deadlocks;
    }

    /// \brief Computes relprod(U, group).
    ldd relprod_impl(const ldd& U, const pbes_summand_group& group, std::size_t i)
    {
      if (m_options.no_relprod)
      {
        ldd z = symbolic::alternative_relprod(U, group);
        mCRL2log(log::trace) << "relprod(" << i << ", todo) = " << print_states(m_data_index, z) << std::endl;
        return z;
      }
      else
      {
        ldd z = relprod(U, group.L, group.Ir);
        mCRL2log(log::trace) << "relprod(" << i << ", todo) = " << print_states(m_data_index, z) << std::endl;
        return z;
      }
    }

    /// \brief Perform a single breadth first step.
    /// \returns The tuple <visited, todo, deadlocks>
    std::tuple<ldd, ldd, ldd> step(const ldd& visited, const ldd& todo, bool learn_transitions = true, bool detect_deadlocks = false)
    {
      using namespace sylvan::ldds;
      auto& R = m_summand_groups;

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

            mCRL2log(log::trace) << "L =\n" << print_relation(m_data_index, R[i].L, R[i].read, R[i].write) << std::endl;
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

            mCRL2log(log::trace) << "L =\n" << print_relation(m_data_index, R[i].L, R[i].read, R[i].write) << std::endl;
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

    ldd run(bool report_states = false)
    {
      using namespace sylvan::ldds;
      auto& R = m_summand_groups;
      std::size_t iteration_count = 0;

      mCRL2log(log::trace) << "initial state = " << core::detail::print_list(m_initial_state) << std::endl;

      stopwatch timer;
      m_initial_vertex = initial_state();
      m_visited = empty_set();
      m_todo = m_initial_vertex;
      m_deadlocks = empty_set();

      while (m_todo != empty_set() && !solution_found() && (m_options.max_iterations == 0 || iteration_count < m_options.max_iterations))
      {
        stopwatch loop_start;
        iteration_count++;
        mCRL2log(log::trace) << "--- iteration " << iteration_count << " ---" << std::endl;
        mCRL2log(log::trace) << "todo = " << print_states(m_data_index, m_todo) << std::endl;
        ldd deadlocks = empty_set();

        std::tie(m_visited, m_todo, deadlocks) = step(m_visited, m_todo, true, m_options.detect_deadlocks);

        if (m_options.detect_deadlocks)
        {
          m_deadlocks = union_(m_deadlocks, deadlocks);
        }

        mCRL2log(log::verbose) << "generated " << std::setw(12) << print_size(union_(m_visited, m_todo)) << " BES equations after "
                               << std::setw(3) << iteration_count << " iterations (time = " << std::setprecision(2)
                               << std::fixed << loop_start.seconds() << "s)" << std::endl;

        if (m_options.detect_deadlocks)
        {
          mCRL2log(log::verbose) << "found " << std::setw(12) << print_size(m_deadlocks) << " deadlocks" << std::endl;
        }

        on_end_while_loop();
        sylvan::sylvan_stats_report(stderr);
      }

      if (report_states)
      {
        std::cout << "number of BES equations = " << print_size(m_visited) << " (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)" << std::endl;
      }
      else
      {
        mCRL2log(log::verbose) << "number of BES equations = " << print_size(m_visited) << " (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)" << std::endl;
      }

      mCRL2log(log::verbose) << "used variable order = " << core::detail::print_list(m_variable_order) << std::endl;

      double total_time = 0.0;
      for (std::size_t i = 0; i < R.size(); i++)
      {
        mCRL2log(log::verbose) << "group " << std::setw(4) << i << " contains " << std::setw(7) << print_size(R[i].L) << " transitions (learn time = "
                               << std::setw(5) << std::setprecision(2) << std::fixed << R[i].learn_time << "s with " << std::setw(9) << R[i].learn_calls 
                               << " calls, cached " << print_size(R[i].Ldomain) << " values" 
                               << std::endl;

        total_time += R[i].learn_time;
      }
      mCRL2log(log::verbose) << "learning transitions took " << total_time << "s" << std::endl;

      std::size_t i = 0;
      for (const auto& param : m_process_parameters)
      {
        auto& table = m_data_index[i];

        mCRL2log(log::verbose) << "Parameter " << i << " (" << param << ")" << " has " << table.size() << " values."<< std::endl;
        for (const auto& data : table)
        {
          mCRL2log(log::debug) << table.index(data) << ": " << data << std::endl;
        }

        ++i;
      }

      return m_visited;
    }

    /// \brief This function is called right after the while loop is finished.
    virtual void on_end_while_loop()
    { }

    /// \returns True iff the solution for the initial state is true.
    virtual bool solution_found() const
    {
      return false;
    }

    /// \returns True iff the vertex is won by even and nothing if the solution has not been determined.
    virtual sylvan::ldds::ldd W0() const
    {
      return sylvan::ldds::empty_set();
    }

    virtual sylvan::ldds::ldd W1() const
    {
      return sylvan::ldds::empty_set();
    }

    std::vector<symbolic::summand_group> summand_groups() const
    {
      std::vector<symbolic::summand_group> result;

      for (const auto& group : m_summand_groups)
      {
        result.push_back(group);
      }

      return result;
    }

    const srf_pbes& pbes() const
    {
      return m_pbes;
    }

    const data::variable_list& process_parameters() const
    {
      return m_process_parameters;
    }

    const std::vector<symbolic::data_expression_index>& data_index() const
    {
      return m_data_index;
    }

    std::vector<symbolic::data_expression_index>& data_index()
    {
      return m_data_index;
    }

    const std::vector<boost::dynamic_bitset<>>& read_write_patterns() const
    {
      return m_summand_patterns;
    }

    const std::vector<boost::dynamic_bitset<>>& read_write_group_patterns() const
    {
      return m_group_patterns;
    }
};

} // namespace mcrl2::pbes_system

#endif // MCRL2_ENABLE_SYLVAN

#endif // MCRL2_PBES_PBESREACH_H