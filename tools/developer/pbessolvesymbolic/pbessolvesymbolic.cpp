// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbessolvesymbolic.cpp

#include <array>
#include <chrono>
#include <iomanip>
#include <sylvan_ldd.hpp>
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/pbes/detail/pbes_io.h"
#include "mcrl2/pbes/pbesreach.h"
#include "mcrl2/pbes/symbolic_pbessolve.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/utilities/detail/container_utility.h"

using namespace mcrl2;
using data::tools::rewriter_tool;
using utilities::tools::input_output_tool;


/// \brief maps proposition variable ldd values to (rank, is_disjunctive)
std::map<std::size_t, std::pair<std::size_t, bool>> compute_equation_info(const pbes_system::srf_pbes& pbes, const std::vector<lps::data_expression_index>& data_index)
{
  pbes_system::pbes_equation_index equation_index(pbes);

  // map propositional variable names to the corresponding ldd value
  std::map<core::identifier_string, std::uint32_t> propvar_index;
  for (const data::data_expression& X: data_index[0])
  {
    const auto& X_ = atermpp::down_cast<data::function_symbol>(X);
    std::uint32_t i = propvar_index.size();
    propvar_index[X_.name()] = i;
  }

  // maps ldd values to (rank, is_disjunctive)
  std::map<std::size_t, std::pair<std::size_t, bool>> equation_info;
  for (const auto& equation: pbes.equations())
  {
    const core::identifier_string& name = equation.variable().name();
    std::size_t rank = equation_index.rank(name);
    bool is_disjunctive = !equation.is_conjunctive();
    auto i = propvar_index.find(name);
    if (i != propvar_index.end())
    {
      std::uint32_t ldd_value = i->second;
      equation_info[ldd_value] = { rank, is_disjunctive };
    }
  }

  return equation_info;
}

namespace mcrl2::pbes_system {

class pbesreach_algorithm_partial : public pbes_system::pbesreach_algorithm
{
public:

  pbesreach_algorithm_partial(const pbes_system::pbes& pbesspec, const symbolic_reachability_options& options_) :
    pbes_system::pbesreach_algorithm(pbesspec, options_)
  {
    m_Vwon[0] = sylvan::ldds::empty_set();
    m_Vwon[1] = sylvan::ldds::empty_set();
  }

  void on_end_while_loop() override
  {
    ++iteration_count;
    if (iteration_count % 10 == 0)
    {
      auto equation_info = compute_equation_info(pbes(), data_index());
      pbes_system::symbolic_pbessolve_algorithm solver(m_visited, summand_groups(), equation_info, m_options.no_relprod, false, data_index());

      if (m_options.solve_strategy == 1)
      {
        mCRL2log(log::verbose) << "start cycle detection\n";

        auto result = solver.detect_cycles(m_visited, m_deadlocks);
        m_Vwon[0] = result.first;
        m_Vwon[1] = result.second;

        mCRL2log(log::verbose) << "cycle detection found solution for" << std::setw(12) << satcount(m_Vwon[0]) + satcount(m_Vwon[1]) << " states" << std::endl;

      }
      else if (m_options.solve_strategy == 2)
      {
        mCRL2log(log::verbose) << "start partial solving" << std::endl;

        // Solve assuming that even wins all todo nodes.
        m_Vwon[1] = solver.solve_impl(m_visited, initial_state(), m_deadlocks, union_(m_todo, m_Vwon[0]), m_Vwon[1]).second;

        // Solve assuming that odd wins all todo nodes.
        m_Vwon[0] = solver.solve_impl(m_visited, initial_state(), m_deadlocks, m_Vwon[0], union_(m_todo, m_Vwon[1])).first;
      }
      mCRL2log(log::verbose) << "partial solving found solution for" << std::setw(12) << satcount(m_Vwon[0]) + satcount(m_Vwon[1]) << " states" << std::endl;

      if (m_options.prune_todo_list)
      {
        mCRL2log(log::verbose) << "start pruning todo list" << std::endl;

        using namespace sylvan::ldds;
        auto& R = m_summand_groups;

        stopwatch timer;
        ldd initial_state = state2ldd(m_initial_state);
        ldd visited = union_(m_Vwon[0], m_Vwon[1]);
        ldd todo = initial_state;

        for (std::size_t iter = 1; iter <= iteration_count; ++iter)
        {
          stopwatch loop_start;
          mCRL2log(log::debug) << "--- iteration " << iter << " ---" << std::endl;
          mCRL2log(log::debug) << "todo = " << print_states(m_data_index, todo) << std::endl;

          ldd todo1 = m_options.chaining ? todo : empty_set();

          for (std::size_t i = 0; i < R.size(); i++)
          {
            if (m_options.no_relprod)
            {
              ldd z = lps::alternative_relprod(m_options.chaining ? todo1 : todo, R[i]);
              mCRL2log(log::debug) << "relprod(" << i << ", todo) = " << print_states(m_data_index, z) << std::endl;
              todo1 = union_(z, todo1);
            }
            else
            {
              mCRL2log(log::debug) << "relprod(" << i << ", todo) = " << print_states(m_data_index, relprod(todo, R[i].L, R[i].Ir)) << std::endl;
              todo1 = relprod_union(m_options.chaining ? todo1 : todo, R[i].L, R[i].Ir, todo1);
            }
          }

          visited = union_(visited, todo);
          todo = minus(todo1, visited);

          mCRL2log(log::verbose) << "found " << std::setw(12) << satcount(visited) << " states after "
                                 << std::setw(3) << iter << " iterations (time = " << std::setprecision(2)
                                 << std::fixed << loop_start.seconds() << "s)" << std::endl;
          mCRL2log(log::verbose) << "todo contains " << std::setw(12) << satcount(m_todo) << " states" << std::endl;
          mCRL2log(log::verbose) << "LDD size visited = " << nodecount(m_visited) << ", todo = " << nodecount(m_todo) << " and deadlock = " << nodecount(m_deadlocks) << std::endl;
        }

        mCRL2log(log::verbose) << "pruned todo list to " << satcount(todo) << " states from " << satcount(m_todo) << " states" << std::endl;
        m_todo = todo;
      }
    }
  }

  bool solution_found(const sylvan::ldds::ldd& initial_state) override
  {
    if (includes(m_Vwon[0], initial_state))
    {
      m_result = true;
      return true;
    }
    else if (includes(m_Vwon[1], initial_state))
    {
      m_result = false;
      return true;
    }

    return false;
  }

  bool solution() const override
  {
    return m_result;
  }

private:
  // States for which winners have already been determined.
  std::array<sylvan::ldds::ldd, 2> m_Vwon;
  bool m_result = false;
  std::size_t iteration_count = 0;
};

} // namespace mcrl2::pbes_system

class pbessolvesymbolic_tool: public rewriter_tool<input_output_tool>
{
  typedef rewriter_tool<input_output_tool> super;

  protected:
    pbes_system::symbolic_reachability_options options;

    // Lace options
    std::size_t lace_n_workers = 1;
    std::size_t lace_dqsize = 1024*1024*4; // set large default
    std::size_t lace_stacksize = 0; // use default

    // Sylvan options
    std::size_t min_tablesize = 22;
    std::size_t max_tablesize = 25;
    std::size_t min_cachesize = 22;
    std::size_t max_cachesize = 25;

    void add_options(utilities::interface_description& desc) override
    {
      super::add_options(desc);
      desc.add_option("lace-workers", utilities::make_optional_argument("NAME", "1"), "set number of Lace workers (threads for parallelization), (0=autodetect, default 1)");
      desc.add_option("lace-dqsize", utilities::make_optional_argument("NAME", "4194304"), "set length of Lace task queue (default 1024*1024*4)");
      desc.add_option("lace-stacksize", utilities::make_optional_argument("NAME", "0"), "set size of program stack in kilo bytes (0=default stack size)");
      desc.add_option("min-table-size", utilities::make_optional_argument("NAME", "22"), "minimum Sylvan table size (21-40, default 22)");
      desc.add_option("max-table-size", utilities::make_optional_argument("NAME", "25"), "maximum Sylvan table size (21-40, default 25)");
      desc.add_option("min-cache-size", utilities::make_optional_argument("NAME", "22"), "minimum Sylvan cache size (21-40, default 22)");
      desc.add_option("max-cache-size", utilities::make_optional_argument("NAME", "25"), "maximum Sylvan cache size (21-40, default 25)");
      desc.add_option("cached", "use transition group caching to speed up state space exploration");
      desc.add_option("chaining", "apply the transition groups as a series");
      desc.add_option("groups", utilities::make_optional_argument("GROUPS", "none"),
                      "'none' (default) no summand groups\n"
                      "'simple' summands with the same read/write variables are joined\n"
                      "a list of summand groups separated by semicolons, e.g. '0; 1 3 4; 2 5'");
      desc.add_option("prune-todo-list", "Prune the todo list periodically.");
      desc.add_option("reorder", utilities::make_optional_argument("ORDER", "none"),
                      "'none' (default) no variable reordering\n"
                      "'random' variables are put in a random order\n"
                      "'<order>' a user defined permutation e.g. '1 3 2 0 4'"
      );
      desc.add_option("info", "print read/write information of the summands");
      desc.add_option("solve-strategy",
                      utilities::make_enum_argument<int>("NUM")
                        .add_value_desc(0, "No on-the-fly solving is applied", true)
                        .add_value_desc(1, "Detect winning loops.")
                        .add_value_desc(2, "Solve subgames using the solver."),
                      "Use solve strategy NUM. Strategy 1 periodically applies on-the-fly solving, which may lead to early termination.",
                      's');
      desc.add_option("split-conditions",
                      utilities::make_optional_argument("NUM", "0"),
                      "split conditions to obtain possibly smaller transition groups\n"
                      "0 (default) no splitting performed.\n"
                      "1 only split disjunctive conditions.\n"
                      "2 also split conjunctive conditions into multiple equations which often yield more reachable states.\n"
                      "3 alternative split for conjunctive conditions where even more states can become reachable.",
                      'c');
      desc.add_option("total", "make the SRF PBES total", 't');
      desc.add_hidden_option("no-remove-unused-rewrite-rules", "do not remove unused rewrite rules. ", 'u');
      desc.add_hidden_option("no-one-point-rule-rewrite", "do not apply the one point rule rewriter");
      desc.add_hidden_option("no-discard", "do not discard any parameters");
      desc.add_hidden_option("no-read", "do not discard only-read parameters");
      desc.add_hidden_option("no-write", "do not discard only-write parameters");
      desc.add_hidden_option("no-relprod", "use an inefficient alternative version of relprod (for debugging)");
      desc.add_hidden_option("srf", utilities::make_optional_argument("FILE", ""), "save the preprocessed PBES in SRF format");
      desc.add_hidden_option("dot", utilities::make_optional_argument("FILE", ""), "print the LDD of the parity game in dot format");
    }

    void parse_options(const utilities::command_line_parser& parser) override
    {
      super::parse_options(parser);
      options.cached                                = parser.has_option("cached");
      options.chaining                              = parser.has_option("chaining");
      options.one_point_rule_rewrite                = !parser.has_option("no-one-point-rule-rewrite");
      options.remove_unused_rewrite_rules           = !parser.has_option("no-remove-unused-rewrite-rules");
      options.replace_constants_by_variables        = false; // This option doesn't work in the current implementation
      options.no_discard                            = parser.has_option("no-discard");
      options.no_discard_read                       = parser.has_option("no-read");
      options.no_discard_write                      = parser.has_option("no-write");
      options.no_relprod                            = parser.has_option("no-relprod");
      options.prune_todo_list                       = parser.has_option("prune-todo-list");
      options.info                                  = parser.has_option("info");
      options.summand_groups                        = parser.option_argument("groups");
      options.variable_order                        = parser.option_argument("reorder");
      options.make_total                            = parser.has_option("total");
      if (!options.make_total)
      {
        options.detect_deadlocks                    = true; // This is a required setting if the pbes is not total.
      }
      options.srf                                   = parser.option_argument("srf");
      options.rewrite_strategy                      = rewrite_strategy();
      options.dot_file                              = parser.option_argument("dot");
      if (parser.has_option("lace-workers"))
      {
        lace_n_workers = parser.option_argument_as<int>("lace-workers");
      }
      if (parser.has_option("lace-dqsize"))
      {
        lace_dqsize = parser.option_argument_as<int>("lace-dqsize");
      }
      if (parser.has_option("lace-stacksize"))
      {
        lace_stacksize = parser.option_argument_as<int>("lace-stacksize");
      }
      if (parser.has_option("min-table-size"))
      {
        min_tablesize = parser.option_argument_as<std::size_t>("min-table-size");
      }
      if (parser.has_option("max-table-size"))
      {
        max_tablesize = parser.option_argument_as<std::size_t>("max-table-size");
      }
      if (parser.has_option("min-cache-size"))
      {
        min_cachesize = parser.option_argument_as<std::size_t>("min-cache-size");
      }
      if (parser.has_option("max-cache-size"))
      {
        max_cachesize = parser.option_argument_as<std::size_t>("max-cache-size");
      }
      if (parser.has_option("split-conditions"))
      {
        options.split_conditions = parser.option_argument_as<std::size_t>("split-conditions");
      }

      options.solve_strategy =  parser.option_argument_as<int>("solve-strategy");
      if (options.solve_strategy < 0 || options.solve_strategy > 2)
      {
        throw mcrl2::runtime_error("Invalid strategy " + std::to_string(options.solve_strategy));
      }
    }

  public:
    pbessolvesymbolic_tool()
      : super("pbessolvesymbolic",
              "Wieger Wesselink",
              "Solves a PBES using symbolic data structures",
              "Solves PBES from INFILE. "
              "If INFILE is not present, stdin is used. "
              "The PBES is first instantiated into a parity game, "
              "which is then solved using Zielonka's algorithm. "
             )
    {}

    void solve(pbes_system::pbesreach_algorithm& reach)
    {
      using namespace sylvan::ldds;

      if (options.info)
      {
        std::cout << lps::print_read_write_patterns(reach.read_write_group_patterns());
      }
      else
      {
        timer().start("instantiation");
        ldd V = reach.run();
        timer().finish("instantiation");

        if (reach.solution_found(reach.initial_state()))
        {
          std::cout << (reach.solution() ? "true" : "false") << std::endl;
        }
        else
        {
          auto equation_info = compute_equation_info(reach.pbes(), reach.data_index());
          pbes_system::symbolic_pbessolve_algorithm solver(V, reach.summand_groups(), equation_info, options.no_relprod, options.chaining, reach.data_index());
          mCRL2log(log::debug) << pbes_system::print_pbes_info(reach.pbes()) << std::endl;
          timer().start("solving");
          bool result = solver.solve(V, reach.initial_state(), reach.deadlocks());
          timer().finish("solving");

          std::cout << (result ? "true" : "false") << std::endl;
        }

        if (!options.dot_file.empty())
        {
          print_dot(options.dot_file, V);
        }
      }
    }

    bool run() override
    {
      lace_init(lace_n_workers, lace_dqsize);
      lace_startup(lace_stacksize, nullptr, nullptr);
      sylvan::sylvan_set_sizes(1ULL<<min_tablesize, 1ULL<<max_tablesize, 1ULL<<min_cachesize, 1ULL<<max_cachesize);
      sylvan::sylvan_init_package();
      sylvan::sylvan_init_ldd();

      mCRL2log(log::verbose) << options << std::endl;

      pbes_system::pbes pbesspec = pbes_system::detail::load_pbes(input_filename());

      if (pbesspec.initial_state().empty())
      {
        throw mcrl2::runtime_error("PBESses without parameters are not supported");
      }

      else
      {
        if (options.solve_strategy == 0)
        {
          pbes_system::pbesreach_algorithm reach(pbesspec, options);
          solve(reach);
        }
        else
        {
          pbes_system::pbesreach_algorithm_partial reach(pbesspec, options);
          solve(reach);
        }
      }

      sylvan::sylvan_quit();
      lace_exit();
      return true;
    }
};

int main(int argc, char* argv[])
{
  return pbessolvesymbolic_tool().execute(argc, argv);
}
