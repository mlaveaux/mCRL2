// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesreach.h

#ifndef MCRL2_PBES_PBESREACH_H
#define MCRL2_PBES_PBESREACH_H

#include <chrono>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>
#include <sylvan_ldd.hpp>
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/data/standard_utility.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/utilities/detail/container_utility.h"
#include "mcrl2/data/consistency.h"
#include "mcrl2/data/enumerator.h"
#include "mcrl2/data/merge_data_specifications.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/lps/symbolic_reachability.h"
#include "mcrl2/pbes/detail/instantiate_global_variables.h"
#include "mcrl2/pbes/detail/pbes_io.h"
#include "mcrl2/pbes/normalize.h"
#include "mcrl2/pbes/replace_constants_by_variables.h"
#include "mcrl2/pbes/resolve_name_clashes.h"
#include "mcrl2/pbes/rewriters/one_point_rule_rewriter.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/unify_parameters.h"
#include "mcrl2/utilities/stopwatch.h"
#include "mcrl2/utilities/text_utility.h"

namespace mcrl2 {

namespace pbes_system {

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

struct symbolic_reachability_options: public lps::symbolic_reachability_options
{
  bool make_total = false;
  std::size_t solve_strategy = 0;
  std::size_t split_conditions = 0;
  std::string srf;
};

inline
std::ostream& operator<<(std::ostream& out, const symbolic_reachability_options& options)
{
  out << static_cast<lps::symbolic_reachability_options>(options);
  out << "solve_strategy = " << options.solve_strategy << std::endl;
  out << "split_conditions = " << options.split_conditions << std::endl;
  out << "total = " << std::boolalpha << options.make_total << std::endl;
  return out;
}

inline
std::pair<std::set<data::variable>, std::set<data::variable>> read_write_parameters(
  const pbes_system::srf_equation& equation,
  const pbes_system::srf_summand& summand,
  const data::variable_list& process_parameters)
{
  using utilities::detail::set_union;
  using utilities::detail::set_intersection;
  using utilities::detail::as_set;

  std::set<data::variable> read_parameters = data::find_free_variables(summand.condition());
  std::set<data::variable> write_parameters;

  // We need special handling for the first parameter, since the clause 'propvar == X' is not in the condition yet
  read_parameters.insert(process_parameters.front());
  if (equation.variable().name() != summand.variable().name())
  {
    write_parameters.insert(process_parameters.front());
  }

  const data::data_expression_list& expressions = summand.variable().parameters();
  auto pi = ++process_parameters.begin(); // skip the first parameter
  auto ei = expressions.begin();
  for (; ei != expressions.end(); ++pi, ++ei)
  {
    if (*pi != *ei)
    {
      write_parameters.insert(*pi);
      data::find_free_variables(*ei, std::inserter(read_parameters, read_parameters.end()));
    }
  }

  auto process_parameter_set = as_set(process_parameters);
  return { set_intersection(read_parameters, process_parameter_set), set_intersection(write_parameters, process_parameter_set) };
}

inline
std::map<data::variable, std::size_t> process_parameter_index(const data::variable_list& process_parameters)
{
  std::map<data::variable, std::size_t> result;
  std::size_t i = 0;
  for (const data::variable& v: process_parameters)
  {
    result[v] = i++;
  }
  return result;
}

inline
std::vector<boost::dynamic_bitset<>> compute_read_write_patterns(const pbes_system::srf_pbes& pbesspec, const data::variable_list& process_parameters)
{
  std::vector<boost::dynamic_bitset<>> result;

  std::size_t n = process_parameters.size();
  std::map<data::variable, std::size_t> index = process_parameter_index(process_parameters);

  for (const auto& equation: pbesspec.equations())
  {
    for (const auto& summand: equation.summands())
    {
      auto [read_parameters, write_parameters] = read_write_parameters(equation, summand, process_parameters);
      auto read = lps::parameter_indices(read_parameters, index);
      auto write = lps::parameter_indices(write_parameters, index);
      boost::dynamic_bitset<> rw(2*n);
      for (std::size_t j: read)
      {
        rw[2*j] = true;
      }
      for (std::size_t j: write)
      {
        rw[2*j + 1] = true;
      }
      result.push_back(rw);
    }
  }

  return result;
}

inline
data::data_expression_list make_state(const pbes_system::propositional_variable_instantiation& x, const std::unordered_map<core::identifier_string, data::data_expression>& propvar_map)
{
  data::data_expression_list result = x.parameters();
  result.push_front(propvar_map.at(x.name()));
  return result;
}

struct summand_group: public lps::summand_group
{
  summand_group(
    const pbes_system::srf_pbes& pbesspec,
    const data::variable_list& process_parameters, // the reordered process parameters
    const std::unordered_map<core::identifier_string, data::data_expression>& propvar_map,
    const std::set<std::size_t>& summand_group_indices,
    const boost::dynamic_bitset<>& read_write_pattern,
    const std::vector<boost::dynamic_bitset<>>& read_write_patterns,
    const std::vector<std::size_t> variable_order // a permutation of [0 .. |process_parameters| - 1]
  )
    : lps::summand_group(process_parameters, read_write_pattern)
  {
    using lps::project;
    using utilities::detail::as_vector;
    using utilities::detail::as_set;
    using utilities::detail::set_union;
    using utilities::detail::contains;

    std::set<std::size_t> used;
    for (std::size_t j: read)
    {
      used.insert(2*j);
    }
    for (std::size_t j: write)
    {
      used.insert(2*j + 1);
    }

    const auto& equations = pbesspec.equations();
    std::size_t k = 0;
    for (std::size_t i = 0; i < equations.size(); i++)
    {
      const core::identifier_string& X_i = equations[i].variable().name();
      const auto& equation_summands = equations[i].summands();
      for (std::size_t j = 0; j < equation_summands.size(); j++, k++)
      {
        if (contains(summand_group_indices, k))
        {
          std::vector<int> copy;
          for (std::size_t q: used)
          {
            bool b = read_write_patterns[k][q];
            copy.push_back(b ? 0 : 1);
          }
          const pbes_system::srf_summand& smd = equation_summands[j];
          auto next_state = make_state(smd.variable(), propvar_map);
          next_state = lps::permute_copy(next_state, variable_order);
          summands.emplace_back(data::and_(data::equal_to(process_parameters.front(), propvar_map.at(X_i)), smd.condition()), smd.parameters(), project(as_vector(next_state), write), copy);
        }
      }
    }
  }
};

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
};

class pbesreach_algorithm
{
    using enumerator_element = data::enumerator_list_element_with_substitution<>;

    template <typename Context>
    friend void lps::learn_successors_callback(WorkerP*, Task*, std::uint32_t* v, std::size_t n, void* context);

  protected:
    using ldd = sylvan::ldds::ldd;
    const symbolic_reachability_options& m_options;
    std::vector<per_worker_information> m_workers;
    pbes_system::srf_pbes m_pbes;
    data::rewriter m_rewr;
    data::enumerator_identifier_generator m_id_generator;
    data::enumerator_algorithm<> m_enumerator;
    data::variable_list m_process_parameters;
    std::size_t m_n;
    std::vector<lps::data_expression_index> m_data_index;
    std::vector<summand_group> m_summand_groups;
    data::data_expression_list m_initial_state;
    std::vector<boost::dynamic_bitset<>> m_summand_patterns;
    std::vector<boost::dynamic_bitset<>> m_group_patterns;
    std::vector<std::size_t> m_variable_order;

    ldd m_visited;
    ldd m_todo;
    ldd m_deadlocks;

    ldd state2ldd(const data::data_expression_list& x)
    {
      MCRL2_DECLARE_STACK_ARRAY(v, std::uint32_t, x.size());

      auto vi = v.begin();
      auto di = m_data_index.begin();
      auto xi = x.begin();
      for (; di != m_data_index.end(); ++vi, ++di, ++xi)
      {
        *vi = di->index(*xi);
      }

      return sylvan::ldds::cube(v.data(), x.size());
    };

    /// \brief Updates R.L := R.L U {(x,y) in R | x in X}
    void learn_successors(std::size_t i, summand_group& R, const ldd& X)
    {
      mCRL2log(log::debug) << "learn successors of summand group " << i << " for X = " << print_states(m_data_index, X, R.read) << std::endl;

      using namespace sylvan::ldds;
      std::pair<pbesreach_algorithm&, summand_group&> context{*this, R};
      sat_all(X, lps::learn_successors_callback<std::pair<pbesreach_algorithm&, summand_group&>>, &context);
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
        pbes_system::replace_constants_by_variables(pbesspec, m_rewr, m_workers[0].m_sigma);
      }

      // Add the initial subsitution for every worker.
      for (std::size_t i = 1; i < m_workers.size(); ++i)
      {
        m_workers[i].m_sigma = m_workers[0].m_sigma;
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

      unify_parameters(result);
      pbes_system::resolve_summand_variable_name_clashes(result, result.equations().front().variable().parameters()); // N.B. This is a required preprocessing step.

      // add a sort for the propositional variable names
      data::data_specification propvar_dataspec = construct_propositional_variable_data_specification(result, "PropositionalVariable");
      result.data() = data::merge_data_specifications(result.data(), propvar_dataspec);
      mCRL2log(log::debug) << "--- srf pbes ---\n" << result.to_pbes() << std::endl;

      return result;
    }

  public:
    pbesreach_algorithm(const pbes_system::pbes& pbesspec, const symbolic_reachability_options& options_, std::size_t num_lace_workers)
      : m_options(options_),
        m_workers(num_lace_workers),
        m_pbes(preprocess(pbesspec, options_.make_total)),
        m_rewr(lps::construct_rewriter(m_pbes.data(), m_options.rewrite_strategy, pbes_system::find_function_symbols(pbesspec), m_options.remove_unused_rewrite_rules)),
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
      lps::adjust_read_write_patterns(m_summand_patterns, m_options);

      m_variable_order = lps::compute_variable_order(m_options.variable_order, m_summand_patterns, true);
      assert(m_variable_order[0] == 0); // It is required that the propositional variable name stays up front
      mCRL2log(log::debug) << "variable order = " << core::detail::print_list(m_variable_order) << std::endl;
      m_summand_patterns = lps::reorder_read_write_patterns(m_summand_patterns, m_variable_order);
      mCRL2log(log::debug) << lps::print_read_write_patterns(m_summand_patterns);

      m_process_parameters = lps::permute_copy(m_process_parameters, m_variable_order);
      m_initial_state = lps::permute_copy(m_initial_state, m_variable_order);
      mCRL2log(log::debug) << "process parameters = " << core::detail::print_list(m_process_parameters) << std::endl;
      mCRL2log(log::debug) << "initial state = " << core::detail::print_list(m_initial_state) << std::endl;

      std::vector<std::set<std::size_t>> groups = lps::compute_summand_groups(m_options.summand_groups, m_summand_patterns);
      for (const auto& group: groups)
      {
        mCRL2log(log::debug) << "group " << core::detail::print_set(group) << std::endl;
      }
      m_group_patterns = lps::compute_summand_group_patterns(m_summand_patterns, groups);
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
        m_data_index.push_back(lps::data_expression_index(param.sort()));
      }
    }

    virtual ~pbesreach_algorithm() {};

    ldd initial_state()
    {
      return state2ldd(m_initial_state);
    }

    /// \returns The set of deadlock states.
    ldd deadlocks()
    {
      return m_deadlocks;
    }

    ldd run(bool report_states = false)
    {
      using namespace sylvan::ldds;
      auto& R = m_summand_groups;
      std::size_t iteration_count = 0;

      mCRL2log(log::debug) << "initial state = " << core::detail::print_list(m_initial_state) << std::endl;

      stopwatch timer;
      ldd initial_state = state2ldd(m_initial_state);
      m_visited = empty_set();
      m_todo = initial_state;
      m_deadlocks = empty_set();

      while (m_todo != empty_set() && !solution_found(initial_state))
      {
        stopwatch loop_start;
        iteration_count++;
        mCRL2log(log::debug) << "--- iteration " << iteration_count << " ---" << std::endl;
        mCRL2log(log::debug) << "todo = " << print_states(m_data_index, m_todo) << std::endl;

        ldd todo1 = m_options.chaining ? m_todo : empty_set();
        ldd potential_deadlocks = m_todo;

        for (std::size_t i = 0; i < R.size(); i++)
        {          
          ldd proj = project(m_options.chaining ? todo1 : m_todo, R[i].Ip);
          
          stopwatch learn_start;
          learn_successors(i, R[i], m_options.cached ? minus(proj, R[i].Ldomain) : proj);
          R[i].learn_time += learn_start.seconds();

          mCRL2log(log::debug) << "L =\n" << print_relation(m_data_index, R[i].L, R[i].read, R[i].write) << std::endl;

          if (m_options.no_relprod)
          {
            ldd z = lps::alternative_relprod(m_options.chaining ? todo1 : m_todo, R[i]);
            mCRL2log(log::debug) << "relprod(" << i << ", todo) = " << print_states(m_data_index, z) << std::endl;
            todo1 = union_(z, todo1);
          }
          else
          {
            mCRL2log(log::debug) << "relprod(" << i << ", todo) = " << print_states(m_data_index, relprod(m_todo, R[i].L, R[i].Ir)) << std::endl;
            todo1 = relprod_union(m_options.chaining ? todo1 : m_todo, R[i].L, R[i].Ir, todo1);
          }

          if (m_options.detect_deadlocks)
          {
            potential_deadlocks = minus(potential_deadlocks, relprev(todo1, R[i].L, R[i].Ir, potential_deadlocks));
          }
        }

        m_visited = union_(m_visited, m_todo);
        m_todo = minus(todo1, m_visited);

        // after all transition groups are applied the remaining potential deadlocks are actual deadlocks.
        if (m_options.detect_deadlocks)
        {
          m_deadlocks = union_(m_deadlocks, potential_deadlocks);
        }

        mCRL2log(log::verbose) << "found " << std::setw(12) << satcount(m_visited) << " states after "
                               << std::setw(3) << iteration_count << " iterations (time = " << std::setprecision(2)
                               << std::fixed << loop_start.seconds() << "s)" << std::endl;

        if (m_options.detect_deadlocks)
        {
          mCRL2log(log::verbose) << "found " << std::setw(12) << satcount(m_deadlocks) << " deadlocks" << std::endl;
        }

        on_end_while_loop();
      }

      if (report_states)
      {
        std::cout << "number of states = " << satcount(m_visited) << " (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)" << std::endl;
      }
      else
      {
        mCRL2log(log::verbose) << "number of states = " << satcount(m_visited) << " (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)" << std::endl;
      }

      mCRL2log(log::verbose) << "visited LDD size = " << nodecount(m_visited) << std::endl;
      mCRL2log(log::verbose) << "used variable order = " << core::detail::print_list(m_variable_order) << std::endl;

      for (std::size_t i = 0; i < R.size(); i++)
      {
        mCRL2log(log::verbose) << "group " << std::setw(4) << i << " contains " << std::setw(7) << satcount(R[i].L) << " transitions (learn time = "
                               << std::setw(5) << std::setprecision(2) << std::fixed << R[i].learn_time << "s)" <<std::endl;
        mCRL2log(log::verbose) << "cached " << satcount(R[i].Ldomain) << " values" << std::endl;
      }

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

    /// \returns True iff the solution for the given vertex is known.
    virtual bool solution_found(const sylvan::ldds::ldd&) const
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

    const std::vector<summand_group>& summand_groups() const
    {
      return m_summand_groups;
    }

    const srf_pbes& pbes() const
    {
      return m_pbes;
    }

    const data::variable_list& process_parameters() const
    {
      return m_process_parameters;
    }

    const std::vector<lps::data_expression_index>& data_index() const
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

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_PBESREACH_H
