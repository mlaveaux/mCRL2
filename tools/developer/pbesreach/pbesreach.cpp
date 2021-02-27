// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbesreach.cpp

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

class pbesreach_tool: public rewriter_tool<input_output_tool>
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
      desc.add_option("no-discard", "do not discard any parameters");
      desc.add_option("no-read", "do not discard only-read parameters");
      desc.add_option("no-write", "do not discard only-write parameters");
      desc.add_option("no-relprod", "use an inefficient alternative version of relprod (for debugging)");
      desc.add_option("info", "print read/write information of the summands");
      desc.add_option("groups", utilities::make_optional_argument("GROUPS", "none"),
                      "'none' (default) no summand groups\n"
                      "'simple' summands with the same read/write variables are joined\n"
                      "a list of summand groups separated by semicolons, e.g. '0; 1 3 4; 2 5'");
      desc.add_option("reorder", utilities::make_optional_argument("ORDER", "none"),
                      "'none' (default) no variable reordering\n"
                      "'random' variables are put in a random order\n"
                      "'<order>' a user defined permutation e.g. '1 3 2 0 4'"
      );
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
      options.replace_constants_by_variables        = false; // This option cannot be used in the symbolic algorithm
      options.no_discard                            = parser.has_option("no-discard");
      options.no_discard_read                       = parser.has_option("no-read");
      options.no_discard_write                      = parser.has_option("no-write");
      options.no_relprod                            = parser.has_option("no-relprod");
      options.info                                  = parser.has_option("info");
      options.summand_groups                        = parser.option_argument("groups");
      options.variable_order                        = parser.option_argument("reorder");
      options.make_total                            = parser.has_option("total");
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
    }

  public:
    pbesreach_tool()
      : super("pbesreach",
              "Wieger Wesselink",
              "applies a symbolic reachability algorithm to a PBES",
              "read a PBES from INFILE and write output to OUTFILE. If OUTFILE "
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

      pbes_system::pbes pbesspec = pbes_system::detail::load_pbes(input_filename());

      if (pbesspec.initial_state().empty())
      {
        throw mcrl2::runtime_error("PBESses without parameters are not supported");
      }

      pbes_system::pbesreach_algorithm algorithm(pbesspec, options, lace_workers());

      if (options.info)
      {
        std::cout << lps::print_read_write_patterns(algorithm.read_write_group_patterns());
      }
      else
      {
        ldd V = algorithm.run(true);
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
  return pbesreach_tool().execute(argc, argv);
}
