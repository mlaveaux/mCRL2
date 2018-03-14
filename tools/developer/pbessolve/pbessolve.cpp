// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbessolve.cpp

#include <iostream>

#include "mcrl2/bes/pbes_input_tool.h"
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/lps/detail/instantiate_global_variables.h"
#include "mcrl2/lts/lts_lts.h"
#include "mcrl2/pbes/pbesinst_structure_graph.h"
#include "mcrl2/pbes/solve_structure_graph.h"
#include "mcrl2/utilities/input_output_tool.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system;
using namespace mcrl2::utilities::tools;
using mcrl2::bes::tools::pbes_input_tool;
using mcrl2::data::tools::rewriter_tool;
using mcrl2::utilities::tools::input_tool;

/// \brief Loads an lps from input_filename, or from stdin if filename equals "".
inline
lps::specification load_lps(const std::string& input_filename)
{
  lps::specification result;
  if (input_filename.empty())
  {
    result.load(std::cin);
  }
  else
  {
    std::ifstream from(input_filename, std::ifstream::in | std::ifstream::binary);
    result.load(from);
  }
  return result;
}

/// \brief Saves an LPS to output_filename, or to stdout if filename equals "".
inline
void save_lps(const lps::specification& lpsspec, const std::string& output_filename)
{
  if (output_filename.empty())
  {
    lpsspec.save(std::cout);
  }
  else
  {
    std::ofstream to(output_filename, std::ofstream::out | std::ofstream::binary);
    lpsspec.save(to);
  }
}

/// \brief Loads a PBES from input_filename, or from stdin if filename equals "".
inline
pbes load_pbes(const std::string& input_filename)
{
    pbes result;
    if (input_filename.empty())
    {
        result.load(std::cin);
    }
    else
    {
        std::ifstream from(input_filename, std::ifstream::in | std::ifstream::binary);
        result.load(from);
    }
    return result;
}

class pbessolve_tool: public rewriter_tool<pbes_input_tool<input_tool>>
{
  protected:
    typedef rewriter_tool<pbes_input_tool<input_tool>> super;

    transformation_strategy m_transformation_strategy; // The strategy to substitute the value of variables with
                                                       // a trivial rhs (aka true or false) in other equations when generating a BES.
    search_strategy m_search_strategy;                 // The search strategy (depth first/breadth first)
    bool check_strategy = false;
    std::string lpsfile;
    std::string ltsfile;
    lts::lts_lts_t ltsspec;

    void parse_options(const utilities::command_line_parser& parser) override
    {
      super::parse_options(parser);
      check_strategy = parser.options.count("check-strategy") > 0;
      if (parser.options.count("lpsfile") > 0 && parser.options.count("ltsfile") > 0)
      {
        throw mcrl2::runtime_error("It is not allowed to use both options --lpsfile and --ltsfile");
      }
      if (parser.options.count("lpsfile") > 0)
      {
        lpsfile = parser.option_argument("lpsfile");
      }
      if (parser.options.count("ltsfile") > 0)
      {
        ltsfile = parser.option_argument("ltsfile");
        ltsspec.load(ltsfile);
      }
      m_transformation_strategy = parser.option_argument_as<mcrl2::pbes_system::transformation_strategy>("strategy");
      m_search_strategy         = parser.option_argument_as<mcrl2::pbes_system::search_strategy>("search");
    }

    void add_options(utilities::interface_description& desc) override
    {
      super::add_options(desc);
      desc.add_hidden_option("check-strategy", "do a sanity check on the computed strategy", 'c');
      desc.add_option("strategy",
                 utilities::make_enum_argument<transformation_strategy>("STRAT")
                   .add_value(lazy, true)
                   .add_value(optimize)
                   .add_value(on_the_fly)
                   .add_value(on_the_fly_with_fixed_points),
                 "use substitution strategy STRAT:",
                 's');
      desc.add_option("search",
                 utilities::make_enum_argument<search_strategy>("SEARCH")
                   .add_value(breadth_first, true)
                   .add_value(depth_first)
                   .add_value(breadth_first_short)
                   .add_value(depth_first_short),
                 "use search strategy SEARCH:",
                 'z');
      desc.add_option("lpsfile",
                 utilities::make_optional_argument("NAME", "name"),
                 "The file containing the LPS that was used to generate the PBES. If this option is set, a counter example LPS will be generated.",
                 'f');
      desc.add_option("ltsfile",
                 utilities::make_optional_argument("NAME", "name"),
                 "The file containing the LTS that was used to generate the PBES. If this option is set, a counter example LTS will be generated.",
                 'g');
    }

  public:
    pbessolve_tool()
      : super("pbessolve",
              "Wieger Wesselink",
              "Generate a BES from a PBES and solve it. ",
              "Solves (P)BES from INFILE. "
              "If INFILE is not present, stdin is used. "
             ),
      m_transformation_strategy(lazy),
      m_search_strategy(breadth_first)
    {}

    bool run() override
    {
      pbes_system::pbes pbesspec = load_pbes(input_filename());
      pbes_system::algorithms::normalize(pbesspec);
      structure_graph G;

      pbesinst_structure_graph_algorithm algorithm(pbesspec, G, rewrite_strategy(), m_search_strategy, m_transformation_strategy);
      mCRL2log(log::verbose) << "Generating parity game..." << std::endl;
      timer().start("instantiation");
      algorithm.run();
      timer().finish("instantiation");

      if (!lpsfile.empty())
      {
        lps::specification lpsspec = load_lps(lpsfile);
        lps::detail::instantiate_global_variables(lpsspec); // N.B. This is necessary, because the global variables might not be valid for the evidence.
        bool result;
        lps::specification evidence;
        timer().start("solving");
        std::tie(result, evidence) = solve_structure_graph_with_counter_example(G, lpsspec, pbesspec, algorithm.equation_index());
        timer().finish("solving");
        std::cout << (result ? "true" : "false") << std::endl;
        std::string output_filename = input_filename() + ".evidence.lps";
        save_lps(evidence, output_filename);
        mCRL2log(log::verbose) << "Saved " << (result ? "witness" : "counter example") << " in " << output_filename << std::endl;
      }
      else if (!ltsfile.empty())
      {
        lts::lts_lts_t evidence;
        timer().start("solving");
        bool result = solve_structure_graph_with_counter_example(G, ltsspec, pbesspec, algorithm.equation_index());
        timer().finish("solving");
        std::cout << (result ? "true" : "false") << std::endl;
        std::string output_filename = input_filename() + ".evidence.lts";
        ltsspec.save(output_filename);
        mCRL2log(log::verbose) << "Saved " << (result ? "witness" : "counter example") << " in " << output_filename << std::endl;
      }
      else
      {
        timer().start("solving");
        bool result = solve_structure_graph(G, check_strategy);
        timer().finish("solving");
        std::cout << (result ? "true" : "false") << std::endl;
      }
      return true;
    }
};

int main(int argc, char* argv[])
{
  return pbessolve_tool().execute(argc, argv);
}