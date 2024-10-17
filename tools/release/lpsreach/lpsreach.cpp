// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lpsreach.cpp

#include <sylvan_ldd.hpp>
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/lps/lpsreach.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/symbolic_lts_io.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/utilities/parallel_tool.h"
#include "mcrl2/symbolic/ordering.h"

using namespace mcrl2;
using data::tools::rewriter_tool;
using utilities::tools::input_output_tool;
using utilities::tools::parallel_tool;

class lpsreach_tool: public parallel_tool<rewriter_tool<input_output_tool>>
{
  using super = parallel_tool<rewriter_tool<input_output_tool>>;

  protected:
    symbolic::symbolic_reachability_options options;

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
      desc.add_option("deadlock", "report the number of deadlocks (i.e. states with no outgoing transitions).");
      desc.add_option("info", "print read/write information of the summands");
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
      desc.add_option("max-iterations", utilities::make_optional_argument("NUM", "0"), "limit number of breadth-first iterations to NUM");
      desc.add_option("print-exact", "prints the sizes of LDDs exactly when within the representable range, and in scientific notation otherwise");
      desc.add_option("print-nodesize", "print the number of LDD nodes in addition to the number of elements represented as 'elements[nodes]'");
      desc.add_option("saturation", "reduce the amount of breadth-first iterations required by applying the transition groups until fixed point is reached");
      desc.add_hidden_option("no-discard", "do not discard any parameters");
      desc.add_hidden_option("no-read", "do not discard only-read parameters");
      desc.add_hidden_option("no-write", "do not discard only-write parameters");
      desc.add_hidden_option("no-relprod", "use an inefficient alternative version of relprod (for debugging)");
      desc.add_hidden_option("no-remove-unused-rewrite-rules", "do not remove unused rewrite rules. ", 'u');
      desc.add_hidden_option("no-one-point-rule-rewrite", "do not apply the one point rule rewriter");
      desc.add_option("initial-ratio", utilities::make_optional_argument("NUM", "16"), "power-of-two ratio of initial and maximum table size (default 16)");
      desc.add_option("table-ratio", utilities::make_optional_argument("NUM", "16"), "power-of-two ratio of node table and cache table (default 1)");
      desc.add_hidden_option("dot", utilities::make_optional_argument("FILE", ""), "print the LDD of the parity game in dot format");
    }

    void parse_options(const utilities::command_line_parser& parser) override
    {
      super::parse_options(parser);
      options.cached                                = parser.has_option("cached");
      options.chaining                              = parser.has_option("chaining");
      options.detect_deadlocks                      = parser.has_option("deadlock");
      options.one_point_rule_rewrite                = !parser.has_option("no-one-point-rule-rewrite");
      options.print_exact                           = parser.has_option("print-exact");
      options.print_nodesize                        = parser.has_option("print-nodesize");
      options.remove_unused_rewrite_rules           = !parser.has_option("no-remove-unused-rewrite-rules");
      options.replace_constants_by_variables        = false; // This option cannot be used in the symbolic algorithm
      options.saturation                            = parser.has_option("saturation");
      options.no_discard                            = parser.has_option("no-discard");
      options.no_discard_read                       = parser.has_option("no-read");
      options.no_discard_write                      = parser.has_option("no-write");
      options.no_relprod                            = parser.has_option("no-relprod");
      options.info                                  = parser.has_option("info");
      options.summand_groups                        = parser.option_argument("groups");
      options.variable_order                        = parser.option_argument("reorder");
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
      if (parser.has_option("max-iterations"))
      {
        options.max_iterations = parser.option_argument_as<std::size_t>("max-iterations");
      }
    }

  public:
    lpsreach_tool()
      : super("lpsreach",
              "Wieger Wesselink",
              "Applies a symbolic reachability algorithm to an LPS",
              "Read an LPS from INFILE and write a symbolic labelled transition system to OUTFILE. "
              "If OUTFILE is not present, no symbolic LTS is written. If INFILE is not present, stdin is used."
             )
    {}

    bool run() override
    {
      using namespace sylvan::ldds;

      lace_init(lace_n_workers, lace_dqsize);
      lace_startup(lace_stacksize, nullptr, nullptr);
      sylvan::sylvan_set_limits(memory_limit * 1024 * 1024 * 1024, std::log2(table_ratio), std::log2(initial_ratio));
      sylvan::sylvan_init_package();
      sylvan::sylvan_init_ldd();

      mCRL2log(log::verbose) << options << std::endl;

      lps::stochastic_specification stochastic_lpsspec;
      lps::load_lps(stochastic_lpsspec, input_filename());
      lps::specification lpsspec = lps::remove_stochastic_operators(stochastic_lpsspec);

      lps::lpsreach_algorithm algorithm(lpsspec, options);

      if (options.info)
      {
        std::cout << symbolic::print_read_write_patterns(algorithm.read_write_group_patterns());
      }
      else
      {
        ldd V = algorithm.run();
        if (!options.dot_file.empty())
        {
          print_dot(options.dot_file, V);
        }

        if (!output_filename().empty())
        {
          std::ofstream to(output_filename(), std::ofstream::out | std::ofstream::binary);
          if (!to.good())
          {
            throw mcrl2::runtime_error("Could not write to filename " + output_filename());
          }

          to << algorithm.get_symbolic_lts();
        }
      }

      sylvan::sylvan_quit();
      lace_exit();
      return true;
    }
};

int main(int argc, char* argv[])
{
  return lpsreach_tool().execute(argc, argv);
}
