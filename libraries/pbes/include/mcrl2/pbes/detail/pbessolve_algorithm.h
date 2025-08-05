// Author(s): Wieger Wesselink, Maurice Laveaux; Threads are added by Jan Friso Groote 
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_PBES_DETAIL_PBESSOLVE_ALGORITHM_H
#define MCRL2_PBES_DETAIL_PBESSOLVE_ALGORITHM_H

#include "mcrl2/pbes/pbessolve_options.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/utilities/logger.h"
#include "mcrl2/utilities/parallel_tool.h"
#include "mcrl2/utilities/file_utility.h"
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/lps/detail/instantiate_global_variables.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/pbes/pbes_input_tool.h"
#include "mcrl2/pbes/detail/pbes_io.h"
#include "mcrl2/pbes/detail/pbes_remove_counterexample_info.h"
#include "mcrl2/pbes/pbesinst_lazy_counter_example.h"
#include "mcrl2/pbes/pbesinst_structure_graph2.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system::detail {

using namespace mcrl2;
using namespace mcrl2::pbes_system;
using namespace mcrl2::utilities::tools;
using mcrl2::pbes_system::tools::pbes_input_tool;
using mcrl2::data::tools::rewriter_tool;
using mcrl2::utilities::tools::input_tool;
using utilities::tools::parallel_tool;

inline
bool run_solve(const pbes_system::pbes& pbesspec, 
  const data::mutable_map_substitution<>& sigma,
  structure_graph& G,
  const pbes_equation_index& equation_index,  
  pbessolve_options options,
  const std::string& input_filename,
  const std::string& lpsfile,
  const std::string& ltsfile,
  std::string evidence_file,
  mcrl2::utilities::execution_timer& timer)
{  
  bool result;
  if (!lpsfile.empty())
  {
    lps::specification lpsspec;
    lps::load_lps(lpsspec, lpsfile);
    lps::detail::replace_global_variables(lpsspec, sigma);
          
    lps::specification evidence;
    timer.start("solving");
    std::tie(result, evidence) = solve_structure_graph_with_counter_example(
        G, lpsspec, pbesspec, equation_index);
    timer.finish("solving");

    std::cout << (result ? "true" : "false") << std::endl;
    if (evidence_file.empty())
    {
      evidence_file = input_filename + ".evidence.lps";
    }
    lps::save_lps(evidence, evidence_file);
    mCRL2log(log::verbose)
        << "Saved " << (result ? "witness" : "counter example") << " in "
        << evidence_file << std::endl;
  }
  else if (!ltsfile.empty())
  {
    lts::lts_lts_t ltsspec;
    ltsspec.load(ltsfile);

    lts::lts_lts_t evidence;
    timer.start("solving");
    result = solve_structure_graph_with_counter_example(G, ltsspec);
    timer.finish("solving");
    std::cout << (result ? "true" : "false") << std::endl;
    if (evidence_file.empty())
    {
      evidence_file = input_filename + ".evidence.lts";
    }
    ltsspec.save(evidence_file);
    mCRL2log(log::verbose)
        << "Saved " << (result ? "witness" : "counter example") << " in "
        << evidence_file << std::endl;
  }
  else
  {
    timer.start("solving");
    result = solve_structure_graph(G, options.check_strategy);
    timer.finish("solving");
    std::cout << (result ? "true" : "false") << std::endl;
  }

  return result;
}

class pbessolve_tool
    : public parallel_tool<rewriter_tool<pbes_input_tool<input_tool>>>
{
  protected:
    using super = parallel_tool<rewriter_tool<pbes_input_tool<input_tool>>>;

    pbessolve_options options;
    int m_short_strategy = 0;
    partial_solve_strategy m_long_strategy = partial_solve_strategy::no_optimisation;
    std::string lpsfile;
    std::string ltsfile;
    std::string evidence_file;

    void add_options(utilities::interface_description& desc) override
    {
      super::add_options(desc);
      desc.add_hidden_option("check-strategy", "do a sanity check on the computed strategy", 'c');
      desc.add_option("search-strategy",
          utilities::make_enum_argument<search_strategy>("NAME")
              .add_value_desc(breadth_first, "Leads to smaller counter examples", true)
              .add_value_desc(depth_first, ""),
          "Use search strategy NAME:",
          'z');
      desc.add_option("file",
          utilities::make_file_argument("NAME"),
          "The file containing the LPS or LTS that was used to "
          "generate the PBES using lps2pbes -c. If this "
          "option is set, a counter example or witness for the "
          "encoded property will be generated. The "
          "extension of the file should be .lps in case of an LPS "
          "file, in all other cases it is assumed to "
          "be an LTS.",
          'f');
      desc.add_option("prune-todo-list", "Prune the todo list periodically.");
      desc.add_hidden_option("naive-counter-example-instantiation",
          "run the naive instantiation algorithm for pbes with counter example information");
      desc.add_hidden_option("no-remove-unused-rewrite-rules", "do not remove unused rewrite rules. ", 'u');
      desc.add_option("evidence-file",
          utilities::make_file_argument("NAME"),
          "The file to which the evidence is written. If not set, a "
          "default name will be chosen.");
      desc.add_option("solve-strategy",
          utilities::make_enum_argument<int>("NAME")
              .add_value_desc(0, "No on-the-fly solving is applied", true)
              .add_value_desc(1, "Propagate solved equations using an attractor.")
              .add_value_desc(2, "Detect winning loops.")
              .add_value_desc(3, "Solve subgames using a fatal attractor.")
              .add_value_desc(4, "Solve subgames using the solver."),
          "Use solve strategy NAME. Strategies 1-4 periodically apply on-the-fly "
          "solving, which may lead to early termination.",
          's');
      desc.add_hidden_option("long-strategy",
          utilities::make_enum_argument<partial_solve_strategy>("STRATEGY")
              .add_value_desc(partial_solve_strategy::no_optimisation, "Do not apply any optimizations.")
              .add_value_desc(partial_solve_strategy::remove_self_loops, "Remove self loops.")
              .add_value_desc(partial_solve_strategy::propagate_solved_equations_using_substitution,
                  "Propagate solved equations using substitution.")
              .add_value_desc(partial_solve_strategy::propagate_solved_equations_using_attractor,
                  "Propagate solved equations using an attractor.")
              .add_value_desc(partial_solve_strategy::detect_winning_loops_using_fatal_attractor,
                  "Detect winning loops using a fatal attractor.")
              .add_value_desc(partial_solve_strategy::solve_subgames_using_fatal_attractor_local,
                  "Solve subgames using a fatal attractor (local version).")
              .add_value_desc(partial_solve_strategy::solve_subgames_using_fatal_attractor_original,
                  "Solve subgames using a fatal attractor (original version).")
              .add_value_desc(partial_solve_strategy::solve_subgames_using_solver, "Solve subgames using the solver.")
              .add_value_desc(partial_solve_strategy::detect_winning_loops_original,
                  "Detect winning loops (original version)."
                  " N.B. This optimization does not work "
                  "correctly in combination with counter examples."
                  " It may also cause stack overflow."),
          "use strategy STRATEGY (N.B. This is a developer option that overrides "
          "--strategy)",
          'l');
      desc.add_hidden_option("no-replace-constants-by-variables",
          "Do not move constant expressions to a substitution.");
      desc.add_hidden_option("aggressive", "Apply optimizations 4 and 5 at every iteration.");
      desc.add_hidden_option("prune-todo-alternative", "Use a variation of todo list pruning.");
  }

  void parse_options(const utilities::command_line_parser& parser) override
  {
    super::parse_options(parser);

    options.check_strategy = parser.has_option("check-strategy");
    options.replace_constants_by_variables =
        !parser.has_option("no-replace-constants-by-variables");
    options.remove_unused_rewrite_rules =
        !parser.has_option("no-remove-unused-rewrite-rules");
    options.aggressive = parser.has_option("aggressive");
    options.prune_todo_list = parser.has_option("prune-todo-list");
    options.prune_todo_alternative =
        parser.has_option("prune-todo-alternative");
    options.exploration_strategy =
        parser.option_argument_as<mcrl2::pbes_system::search_strategy>(
            "search-strategy");
    options.rewrite_strategy = rewrite_strategy();
    options.number_of_threads = number_of_threads();
    options.naive_counter_example_instantiation = parser.has_option("naive-counter-example-instantiation");

    if (parser.has_option("file"))
    {
      std::string filename = parser.option_argument("file");
      if (mcrl2::utilities::file_extension(filename) == "lps")
      {
        lpsfile = filename;
      }
      else
      {
        ltsfile = filename;
      }
    }

    if (parser.has_option("evidence-file"))
    {
      if (!parser.has_option("file"))
      {
        throw mcrl2::runtime_error(
            "Option --evidence-file cannot be used without option --file");
      }
      evidence_file = parser.option_argument("evidence-file");
    }

    if (parser.has_option("long-strategy"))
    {
      m_long_strategy = parser.option_argument_as<partial_solve_strategy>("long-strategy");
    }
    else
    {
      m_short_strategy = parser.option_argument_as<int>("solve-strategy");
    }
  }

  std::set<utilities::file_format> available_input_formats() const override
  {
    return {pbes_system::pbes_format_internal()};
  }

  public:
  pbessolve_tool(const std::string& toolname)
      : super(toolname, "Wieger Wesselink",
              "Generate a BES from a PBES and solve it. ",
              "Solves (P)BES from INFILE. "
              "If INFILE is not present, stdin is used. "
              "The PBES is first instantiated into a parity game, "
              "which is then solved using Zielonka's algorithm. "
              "It supports the generation of a witness or counter "
              "example for the property encoded by the PBES.")
  {
  }

  template <typename PbesInstAlgorithm, typename PbesInstAlgorithmCE>
  void run_algorithm(pbes_system::pbes& pbesspec,
    const data::mutable_map_substitution<>& sigma)
  {
    bool has_counter_example = detail::has_counter_example_information(pbesspec);
    if (has_counter_example)
    {
      if (lpsfile.empty() && ltsfile.empty())
      {
        mCRL2log(log::warning)
            << "Warning: the PBES has counter example information, but no witness will be generated due to lack of --file"
            << std::endl;
      }
    }
    else if ((!lpsfile.empty() || !ltsfile.empty()))
    {
      mCRL2log(log::warning)
          << "Warning: the PBES has no counter example information. Did you "
            "use the --counter-example option when generating the PBES?"
          << std::endl;
    }

    // When the original has counter example information we remove it and store the provided pbes.
    if (!has_counter_example || options.naive_counter_example_instantiation)
    {      
      mCRL2log(log::verbose) << "Generating parity game..." << std::endl;
      structure_graph G;
      PbesInstAlgorithm instantiate(options, pbesspec, G);

      timer().start("instantiation");
      instantiate.run();
      timer().finish("instantiation");

      detail::run_solve(pbesspec, sigma, G, instantiate.equation_index(), options, input_filename(), lpsfile, ltsfile, evidence_file, timer());
    }
    else
    {
      // Remove the counter example information, but store it for the later step.
      mCRL2log(log::verbose) << "Removing counter example information for first pass." << std::endl;
      pbes_system::pbes pbesspec_without_counterexample = detail::remove_counterexample_info(pbesspec);
      mCRL2log(log::trace) << pbesspec_without_counterexample;

      mCRL2log(log::verbose) << "Generating parity game..." << std::endl;
      structure_graph initial_G;
      PbesInstAlgorithm first_instantiate(options, pbesspec_without_counterexample, initial_G);

      timer().start("first-instantiation");
      first_instantiate.run();
      timer().finish("first-instantiation");

      mCRL2log(log::verbose) << "Number of vertices in the structure graph: "
                             << initial_G.all_vertices().size() << std::endl;      

      // Solve the initial pbes and obtain the strategies in G.
      timer().start("first-solving");
      auto [result, mapping] = solve_structure_graph_winning_mapping(initial_G, true);
      timer().finish("first-solving");
      mCRL2log(log::log_level_t::verbose) << (result ? "true" : "false") << std::endl;

      // Based on the result remove the unnecessary equations related to counter example information. 
      mCRL2log(log::verbose) << "Removing unnecessary example information for other player." << std::endl;
      pbesspec = detail::remove_counterexample_info(pbesspec, !result, result); 
      mCRL2log(log::trace) << pbesspec << std::endl;
      
      structure_graph G;
      PbesInstAlgorithmCE second_instantiate(options, pbesspec, initial_G, !result, mapping, G, first_instantiate.data_rewriter());
      
      // Perform the second instantiation given the proof graph.      
      timer().start("second-instantiation");
      second_instantiate.run();
      timer().finish("second-instantiation");

      mCRL2log(log::verbose) << "Number of vertices in the structure graph: "
                             << G.all_vertices().size() << std::endl;
      
      bool final_result = detail::run_solve(pbesspec, sigma, G, second_instantiate.equation_index(), options, input_filename(), lpsfile, ltsfile, evidence_file, timer());
      if(result != final_result) {
        throw mcrl2::runtime_error("The result of the second instantiation does not match the first instantiation. This is a bug in the tool!");
      }
    }
  }

  bool run() override
  {
    pbes_system::pbes pbesspec =
        pbes_system::detail::load_pbes(input_filename());
    pbes_system::algorithms::normalize(pbesspec);
    data::mutable_map_substitution<> sigma;

    if (!lpsfile.empty())
    {
      // Make sure that the global variables of the LPS and the PBES get the
      // same values
      sigma = pbes_system::detail::instantiate_global_variables(pbesspec);
      pbes_system::detail::replace_global_variables(pbesspec, sigma);
    }

    // Handle tool options here because now we know whether the PBES has counter example information.
    if (m_long_strategy > partial_solve_strategy::no_optimisation)
    {
      options.optimization = m_long_strategy;
    }
    else
    {
      if (m_short_strategy == 0)
      {
        options.optimization = partial_solve_strategy::propagate_solved_equations_using_substitution;
      }
      else if (m_short_strategy == 1)
      {
        options.optimization = partial_solve_strategy::propagate_solved_equations_using_attractor;
      }
      else if (m_short_strategy == 2)
      {
        options.optimization = partial_solve_strategy::detect_winning_loops_using_fatal_attractor;
      }
      else if (m_short_strategy == 3)
      {
        options.optimization = partial_solve_strategy::solve_subgames_using_fatal_attractor_original;
      }
      else if (m_short_strategy == 4)
      {
        options.optimization = partial_solve_strategy::solve_subgames_using_solver;
      }
    }
    
    bool has_counter_example = detail::has_counter_example_information(pbesspec);
    if (has_counter_example)
    {      
        mCRL2log(mcrl2::log::warning) << "Warning: Cannot use partial solving with PBES that has counter example information, using strategy 0 instead." << std::endl;
        options.optimization = partial_solve_strategy::no_optimisation;
    }    

    if (options.optimization < partial_solve_strategy::no_optimisation || options.optimization >  partial_solve_strategy::detect_winning_loops_original)
    {
      throw mcrl2::runtime_error("Invalid strategy " +
                                 std::to_string(static_cast<int>(options.optimization)));
    }
    if (options.prune_todo_list && options.optimization < partial_solve_strategy::propagate_solved_equations_using_substitution)
    {
      mCRL2log(log::warning) << "Option --prune-todo-list has no effect for "
                                "strategies less than 2."
                             << std::endl;
    }
    if (options.optimization == partial_solve_strategy::detect_winning_loops_original && has_counter_example)
    {
      throw mcrl2::runtime_error("optimisation 8 cannot be used with a PBES that has counter example information");
    }
    if (options.optimization == partial_solve_strategy::detect_winning_loops_original && options.number_of_threads > 1)
    {
      throw mcrl2::runtime_error("optimisation 8 does not work correctly with multiple threads, using 1 thread instead.");
    }

    mCRL2log(log::log_level_t::verbose) << "Using optimisation " << options.optimization << "\n";

    if (options.optimization <= partial_solve_strategy::remove_self_loops)
    {
      run_algorithm<pbesinst_structure_graph_algorithm, pbesinst_counter_example_structure_graph_algorithm>(pbesspec, sigma);
    }
    else
    {
      run_algorithm<pbesinst_structure_graph_algorithm2, pbesinst_counter_example_structure_graph_algorithm2>(pbesspec, sigma);
    }
    return true;
  }
};

inline
bool pbessolve(const pbes& p)
{
  pbessolve_options options;
  pbes pbesspec = p;
  pbes_system::algorithms::normalize(pbesspec);
  structure_graph G;
  pbesinst_structure_graph_algorithm algorithm(options, pbesspec, G);
  algorithm.run();
  return solve_structure_graph(G);
}

} // mcrl2::pbes_system::detail

#endif // MCRL2_PBES_DETAIL_PBESSOLVE_ALGORITHM_H