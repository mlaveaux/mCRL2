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
#include "mcrl2/utilities/parallel_tool.h"
#include "mcrl2/utilities/power_of_two.h"

using namespace mcrl2;
using data::tools::rewriter_tool;
using utilities::tools::input_output_tool;
using utilities::tools::parallel_tool;

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
    time_exploring += explore_timer.seconds();
    ++iteration_count;

    //if (iteration_count % 10 == 0 || m_options.aggressive)
    if (time_solving * 10 < (time_solving + time_exploring) || m_options.aggressive)
    {
      mCRL2log(log::verbose) << "start partial solving\n";
      stopwatch timer;

      // Store the set of won states to keep track of whether new states have been solved.
      std::array<sylvan::ldds::ldd, 2> Vwon = m_Vwon;

      ldd V = union_(m_visited, m_todo);
      pbes_system::symbolic_parity_game G(pbes(), summand_groups(), data_index(), V, m_options.no_relprod, m_options.chaining);
      G.print_information();
      pbes_system::symbolic_pbessolve_algorithm solver(G);

      if (m_options.solve_strategy == 1)
      {
        std::tie(m_Vwon[0], m_Vwon[1]) = solver.detect_solitair_cycles(m_initial_vertex, V, m_todo, false, m_deadlocks, m_Vwon[0], m_Vwon[1]);
      }      
      else if (m_options.solve_strategy == 2)
      {
        std::tie(m_Vwon[0], m_Vwon[1]) = solver.detect_solitair_cycles(m_initial_vertex, V, m_todo, true, m_deadlocks, m_Vwon[0], m_Vwon[1]);
      }      
      else if (m_options.solve_strategy == 3)
      {
        std::tie(m_Vwon[0], m_Vwon[1]) = solver.detect_forced_cycles(m_initial_vertex, V, m_todo, false, m_deadlocks, m_Vwon[0], m_Vwon[1]);
      }
      else if (m_options.solve_strategy == 4)
      {
        std::tie(m_Vwon[0], m_Vwon[1]) = solver.detect_forced_cycles(m_initial_vertex, V, m_todo, true, m_deadlocks, m_Vwon[0], m_Vwon[1]);
      }
      else if (m_options.solve_strategy == 5)
      {
        std::tie(m_Vwon[0], m_Vwon[1]) = solver.detect_fatal_attractors(m_initial_vertex, V, m_todo, false, m_deadlocks, m_Vwon[0], m_Vwon[1]);
      }
      else if (m_options.solve_strategy == 6)
      {
        std::tie(m_Vwon[0], m_Vwon[1]) = solver.detect_fatal_attractors(m_initial_vertex, V, m_todo, true, m_deadlocks, m_Vwon[0], m_Vwon[1]);
      }
      else if (m_options.solve_strategy == 7)
      {
        std::tie(m_Vwon[0], m_Vwon[1]) = solver.partial_solve(m_initial_vertex, V, m_todo, m_deadlocks, m_Vwon[0], m_Vwon[1]);
      }

      mCRL2log(log::verbose) << "found solution for" << std::setw(12) << satcount(m_Vwon[0]) + satcount(m_Vwon[1]) << " BES equations" << std::endl;
      mCRL2log(log::verbose) << "finished partial solving (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";

      time_solving += timer.seconds();
    }

    explore_timer.reset();
  }

  bool solution_found() const override
  {
    if (includes(m_Vwon[0], m_initial_vertex))
    {
      return true;
    }
    else if (includes(m_Vwon[1], m_initial_vertex))
    {
      return true;
    }

    return false;
  }

  ldd W0() const override
  {
    return m_Vwon[0];
  }

  ldd W1() const override
  {
    return m_Vwon[1];
  }

private:
  // States for which winners have already been determined.
  std::array<sylvan::ldds::ldd, 2> m_Vwon;
  std::size_t iteration_count = 0;

  double time_solving = 0.0;
  double time_exploring = 0.0;
  stopwatch explore_timer;
};

} // namespace mcrl2::pbes_system

class pbessolvesymbolic_tool: public parallel_tool<rewriter_tool<input_output_tool>>
{
  using super = parallel_tool<rewriter_tool<input_output_tool>>;

  protected:
    pbes_system::symbolic_reachability_options options;

    // Lace options
    std::size_t lace_n_workers = 1;
    std::size_t lace_dqsize = 1024*1024*4; // set large default
    std::size_t lace_stacksize = 0; // use default

    // Sylvan options
    std::size_t memory_limit = 3;
    std::size_t initial_ratio = 16;
    std::size_t table_ratio = 1;

    void add_options(utilities::interface_description& desc) override
    {
      super::add_options(desc);
      desc.add_option("lace-dqsize", utilities::make_optional_argument("NUM", "4194304"), "set length of Lace task queue (default 1024*1024*4)");
      desc.add_option("lace-stacksize", utilities::make_optional_argument("NUM", "0"), "set size of program stack in kilobytes (0=default stack size)");
      desc.add_option("memory-limit", utilities::make_optional_argument("NUM", "3"), "Sylvan memory limit in gigabytes (default 3)", 'm');

      desc.add_option("cached", "use transition group caching to speed up state space exploration");
      desc.add_option("chaining", "reduce the amount of breadth-first iterations by applying the transition groups consecutively");
      desc.add_option("groups", utilities::make_optional_argument("GROUPS", "none"),
                      "'none' (default) no summand groups\n"
                      "'used' summands with the same variables are joined\n"
                      "'simple' summands with the same read/write variables are joined\n"
                      "a user defined list of summand groups separated by semicolons, e.g. '0; 1 3 4; 2 5'");
      desc.add_option("reorder", utilities::make_optional_argument("ORDER", "none"),
                      "'none' (default) no variable reordering\n"
                      "'random' variables are put in a random order\n"
                      "'weighted' variables are put in an order defined by their connectivity weight\n"
                      "'a user defined permutation e.g. '1 3 2 0 4'"
      );
      desc.add_option("info", "print read/write information of the summands");
      desc.add_option("max-iterations", utilities::make_optional_argument("NUM", "0"), "limit number of breadth-first iterations to NUM");
      desc.add_option("print-exact", "prints the sizes of LDDs exactly when within the representable range, and in scientific notation otherwise");
      desc.add_option("print-nodesize", "print the number of LDD nodes in addition to the number of elements represented as 'elements[nodes]'");
      desc.add_option("saturation", "reduce the amount of breadth-first iterations by applying the transition groups until fixed point");
      desc.add_option("solve-strategy",
                      utilities::make_enum_argument<int>("NUM")
                        .add_value_desc(0, "No on-the-fly solving is applied", true)
                        .add_value_desc(1, "Detect solitair winning cycles.")
                        .add_value_desc(2, "Detect solitair winning cycles with safe attractors.")
                        .add_value_desc(3, "Detect forced winning cycles.")
                        .add_value_desc(4, "Detect forced winning cycles with safe attractors.")
                        .add_value_desc(5, "Detect fatal attractors.")
                        .add_value_desc(6, "Detect fatal attractors with safe attractors.")
                        .add_value_desc(7, "Solve subgames using a Zielonka solver."),
                      "Use solve strategy NUM. All strategies except 0 periodically apply on-the-fly solving, which may lead to early termination.",
                      's');
      desc.add_option("split-conditions",
                      "split disjunctive conditions to obtain more summands with potentially less dependencies",
                      'c');
      desc.add_option("total", "make the SRF PBES total", 't');
      desc.add_option("reset", "set constant values when introducing parameters");
      desc.add_hidden_option("aggressive", "apply on-the-fly solving after every iteration to detect bugs");
      desc.add_hidden_option("no-remove-unused-rewrite-rules", "do not remove unused rewrite rules. ", 'u');
      desc.add_hidden_option("no-one-point-rule-rewrite", "do not apply the one point rule rewriter");
      desc.add_hidden_option("no-discard", "do not discard any parameters");
      desc.add_hidden_option("no-read", "do not discard only-read parameters");
      desc.add_hidden_option("no-write", "do not discard only-write parameters");
      desc.add_hidden_option("no-relprod", "use an inefficient alternative version of relprod (for debugging)");
      desc.add_hidden_option("initial-ratio", utilities::make_optional_argument("NUM", "16"), "power-of-two ratio of initial and maximum table size (default 16)");
      desc.add_hidden_option("table-ratio", utilities::make_optional_argument("NUM", "16"), "power-of-two ratio of node table and cache table (default 1)");
      desc.add_hidden_option("srf", utilities::make_optional_argument("FILE", ""), "save the preprocessed PBES in SRF format");
      desc.add_hidden_option("dot", utilities::make_optional_argument("FILE", ""), "print the LDD of the parity game in dot format");
      desc.add_hidden_option("split-conditions-unsafe",
                      utilities::make_optional_argument("NUM", "0"),
                      "split conditions to obtain more summands (and equations) with potentially less dependencies\n"
                      "0 (default) no splitting performed.\n"
                      "1 only split disjunctive conditions, same as --split-conditions.\n"
                      "2 also split conjunctive conditions into multiple equations which weakens guards and introduces more reachable BES equations. Note that splitting conditions can lead to expressions that cannot be rewritten if the equations are not sufficiently complete.\n"
                      "3 alternative split for conjunctive conditions where even more states can become reachable.");
    }

    void parse_options(const utilities::command_line_parser& parser) override
    {
      super::parse_options(parser);
      options.aggressive                            = parser.has_option("aggressive");
      options.cached                                = parser.has_option("cached");
      options.chaining                              = parser.has_option("chaining");
      options.one_point_rule_rewrite                = !parser.has_option("no-one-point-rule-rewrite");
      options.print_exact                           = parser.has_option("print-exact");
      options.print_nodesize                        = parser.has_option("print-nodesize");
      options.remove_unused_rewrite_rules           = !parser.has_option("no-remove-unused-rewrite-rules");
      options.replace_constants_by_variables        = false; // This option doesn't work in the current implementation
      options.saturation                            = parser.has_option("saturation");
      options.no_discard                            = parser.has_option("no-discard");
      options.no_discard_read                       = parser.has_option("no-read");
      options.no_discard_write                      = parser.has_option("no-write");
      options.no_relprod                            = parser.has_option("no-relprod");
      options.info                                  = parser.has_option("info");
      options.summand_groups                        = parser.option_argument("groups");
      options.variable_order                        = parser.option_argument("reorder");
      options.make_total                            = parser.has_option("total");
      options.reset_parameters                      = parser.has_option("reset");
      if (!options.make_total)
      {
        options.detect_deadlocks                    = true; // This is a required setting if the pbes is not total.
      }
      options.srf                                   = parser.option_argument("srf");
      options.rewrite_strategy                      = rewrite_strategy();
      options.dot_file                              = parser.option_argument("dot");
      lace_n_workers = number_of_threads();
      if (parser.has_option("lace-dqsize"))
      {
        lace_dqsize = parser.option_argument_as<int>("lace-dqsize");
      }
      if (parser.has_option("lace-stacksize"))
      {
        lace_stacksize = parser.option_argument_as<int>("lace-stacksize");
      }
      if (parser.has_option("memory-limit"))
      {
        memory_limit = parser.option_argument_as<std::size_t>("memory-limit");
      }
      if (parser.has_option("initial-ratio"))
      {
        initial_ratio = parser.option_argument_as<std::size_t>("initial-ratio");
        if (!utilities::is_power_of_two(initial_ratio))
        {
          throw mcrl2::runtime_error("The initial-ratio should be a power of two.");
        }
      }
      if (parser.has_option("table-ratio"))
      {
        table_ratio = parser.option_argument_as<std::size_t>("table-ratio");
        if (!utilities::is_power_of_two(table_ratio))
        {
          throw mcrl2::runtime_error("The table-ratio should be a power of two.");
        }
      }

      if (parser.has_option("split-conditions"))
      {
        options.split_conditions = 1;
      }

      if (parser.has_option("split-conditions-unsafe"))
      {
        options.split_conditions = parser.option_argument_as<std::size_t>("split-conditions-unsafe");
      }

      options.solve_strategy =  parser.option_argument_as<int>("solve-strategy");
      if (options.solve_strategy < 0 || options.solve_strategy > 7)
      {
        throw mcrl2::runtime_error("Invalid strategy " + std::to_string(options.solve_strategy));
      }
      
      if (parser.has_option("max-iterations"))
      {
        options.max_iterations = parser.option_argument_as<std::size_t>("max-iterations");
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
        std::cout << symbolic::print_read_write_patterns(reach.read_write_group_patterns());
      }
      else
      {
        timer().start("instantiation");
        ldd V = reach.run();
        timer().finish("instantiation");

        if (reach.solution_found())
        {
          std::cout << (includes(reach.W0(), (reach.initial_state())) ? "true" : "false") << std::endl;
        }
        else
        {
          if (options.max_iterations == 0)
          {
            pbes_system::symbolic_parity_game G(reach.pbes(), reach.summand_groups(), reach.data_index(), V, options.no_relprod, options.chaining);
            G.print_information();
            pbes_system::symbolic_pbessolve_algorithm solver(G);

            mCRL2log(log::debug) << pbes_system::detail::print_pbes_info(reach.pbes()) << std::endl;
            timer().start("solving");
            bool result = solver.solve(reach.initial_state(), V, reach.deadlocks(), reach.W0(), reach.W1());
            timer().finish("solving");

            std::cout << (result ? "true" : "false") << std::endl;
          }
          else
          {
            // TODO: We could actually try to solve the incomplete parity game.
            std::cout << "Skipped solving since exploration was limited to max-iterations" << std::endl;
          }
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
      sylvan::sylvan_set_limits(memory_limit * 1024 * 1024 * 1024, std::log2(table_ratio), std::log2(initial_ratio));
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
