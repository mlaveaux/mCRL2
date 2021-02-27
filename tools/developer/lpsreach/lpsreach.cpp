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
#include "mcrl2/utilities/input_output_tool.h"

using namespace mcrl2;
using data::tools::rewriter_tool;
using utilities::tools::input_output_tool;

class lpsreach_tool: public rewriter_tool<input_output_tool>
{
  typedef rewriter_tool<input_output_tool> super;

  protected:
    lps::symbolic_reachability_options options;

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
      desc.add_option("deadlock", "report the number of deadlocks (i.e. states with no outgoing transitions).");
      desc.add_option("no-discard", "do not discard any parameters");
      desc.add_option("no-read", "do not discard only-read parameters");
      desc.add_option("no-write", "do not discard only-write parameters");
      desc.add_option("no-relprod", "use an inefficient alternative version of relprod (for debugging)");
      desc.add_option("info", "print read/write information of the summands");
      desc.add_option("groups", utilities::make_optional_argument("GROUPS", "none"),
                      "'none' (default) no summand groups\n"
                      "'simple' summands with the same read/write variables are joined\n"
                      "'<groups>' a user defined list of summand groups separated by semicolons, e.g. '0; 1 3 4; 2 5'");
      desc.add_option("reorder", utilities::make_optional_argument("ORDER", "none"),
                      "'none' (default) no variable reordering\n"
                      "'random' variables are put in a random order\n"
                      "'<order>' a user defined permutation e.g. '1 3 2 0 4'"
                      );
      desc.add_hidden_option("no-remove-unused-rewrite-rules", "do not remove unused rewrite rules. ", 'u');
      desc.add_hidden_option("no-one-point-rule-rewrite", "do not apply the one point rule rewriter");
      desc.add_hidden_option("dot", utilities::make_optional_argument("FILE", ""), "print the LDD of the parity game in dot format");
    }

    void parse_options(const utilities::command_line_parser& parser) override
    {
      super::parse_options(parser);
      options.cached                                = parser.has_option("cached");
      options.chaining                              = parser.has_option("chaining");
      options.detect_deadlocks                      = parser.has_option("deadlock");
      options.one_point_rule_rewrite                = !parser.has_option("no-one-point-rule-rewrite");
      options.remove_unused_rewrite_rules           = !parser.has_option("no-remove-unused-rewrite-rules");
      options.replace_constants_by_variables        = false; // This option cannot be used in the symbolic algorithm
      options.no_discard                            = parser.has_option("no-discard");
      options.no_discard_read                       = parser.has_option("no-read");
      options.no_discard_write                      = parser.has_option("no-write");
      options.no_relprod                            = parser.has_option("no-relprod");
      options.info                                  = parser.has_option("info");
      options.summand_groups                        = parser.option_argument("groups");
      options.variable_order                        = parser.option_argument("reorder");
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
    }

  public:
    lpsreach_tool()
      : super("lpsreach",
              "Wieger Wesselink",
              "applies a symbolic reachability algorithm to an LPS",
              "read an LPS from INFILE and write output to OUTFILE. If OUTFILE "
              "is not present, stdout is used. If INFILE is not present, stdin is used."
             )
    {}

    bool run() override
    {
      using namespace sylvan::ldds;

      lace_init(lace_n_workers, lace_dqsize);
      lace_startup(lace_stacksize, nullptr, nullptr);
      sylvan::sylvan_set_sizes(1ULL<<min_tablesize, 1ULL<<max_tablesize, 1ULL<<min_cachesize, 1ULL<<max_cachesize);
      sylvan::sylvan_init_package();
      sylvan::sylvan_init_ldd();

      mCRL2log(log::verbose) << options << std::endl;

      lps::stochastic_specification stochastic_lpsspec;
      lps::load_lps(stochastic_lpsspec, input_filename());
      lps::specification lpsspec = lps::remove_stochastic_operators(stochastic_lpsspec);

      if (lpsspec.process().process_parameters().empty())
      {
        throw mcrl2::runtime_error("Processes without parameters are not supported");
      }

      lps::lpsreach_algorithm algorithm(lpsspec, options, lace_workers());

      if (options.info)
      {
        std::cout << lps::print_read_write_patterns(algorithm.read_write_group_patterns());
      }
      else
      {
        ldd V = algorithm.run();
        if (!options.dot_file.empty())
        {
          print_dot(options.dot_file, V);
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
