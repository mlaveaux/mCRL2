// Author(s): Jeroen Keiren, Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbessolvesymbolic.cpp

#include <array>
#include <cassert>
#include <cstddef>
#include <iomanip>
#include <string>
#include <sylvan_ldd.hpp>

#include "mcrl2/core/identifier_string.h"
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/pbes/detail/pbes_io.h"
#include "mcrl2/pbes/detail/pbes_remove_counterexample_info.h"
#include "mcrl2/pbes/detail/pbessolve_algorithm.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/pbes_expression.h"
#include "mcrl2/pbes/pbesinst_structure_graph.h"
#include "mcrl2/pbes/pbesreach.h"
#include "mcrl2/pbes/pbesreach_partial.h"
#include "mcrl2/pbes/rewriters/data_rewriter.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/symbolic_pbessolve.h"
#include "mcrl2/pbes/unify_parameters.h"
#include "mcrl2/utilities/exception.h"
#include "mcrl2/utilities/execution_timer.h"
#include "mcrl2/utilities/file_utility.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/utilities/logger.h"
#include "mcrl2/utilities/parallel_tool.h"
#include "mcrl2/utilities/power_of_two.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system;

using data::tools::rewriter_tool;
using utilities::tools::input_output_tool;
using utilities::tools::parallel_tool;

/// \brief The arguments for the lace task.
struct arguments
{
  pbes_system::symbolic_reachability_options options;
  std::string input_filename;
  std::string output_filename;
  std::string evidence_filename;
  std::string lpsfile;
  std::string ltsfile;
  mcrl2::utilities::execution_timer& timer;
};

TASK_DECL_1(bool, pbessolvesymbolic_task, arguments*); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
#define pbessolvesymbolic_task(a) RUN(pbessolvesymbolic_task, a) // NOLINT(cppcoreguidelines-macro-usage)

namespace mcrl2::pbes_system
{

namespace detail
{

/// \brief Computes the LDD cube for the vertex X, i.e. [index(name), index(param_0), ...].
/// \returns false if X is unknown to the symbolic exploration, because its name is not in
///          propvar_map or one of its parameter values is not in data_index; the contents of cube
///          are then unspecified. Callers must never treat this as grounds for pruning: index()
///          returns npos for unknown values, which truncates to 0xFFFFFFFF in the cube and merely
///          makes member_cube return false.
inline bool vertex_cube(const propositional_variable_instantiation& X,
  const std::vector<symbolic::data_expression_index>& data_index,
  const std::unordered_map<core::identifier_string, data::data_expression>& propvar_map,
  std::vector<std::uint32_t>& cube)
{
  cube.clear();

  std::unordered_map<core::identifier_string, data::data_expression>::const_iterator propvar_it
    = propvar_map.find(X.name());
  if (propvar_it == propvar_map.end())
  {
    mCRL2log(log::debug) << "vertex_cube: " << X.name() << " is unknown to the symbolic exploration" << std::endl;
    return false;
  }

  std::size_t name_index = data_index[0].index(propvar_it->second);
  if (name_index == symbolic::data_expression_index::npos)
  {
    mCRL2log(log::debug) << "vertex_cube: " << X.name() << " is unknown to the symbolic exploration" << std::endl;
    return false;
  }
  cube.push_back(static_cast<std::uint32_t>(name_index));

  std::size_t i = 1;
  for (const data::data_expression& param: X.parameters())
  {
    std::size_t param_index = data_index[i].index(param);
    if (param_index == symbolic::data_expression_index::npos)
    {
      mCRL2log(log::debug) << "vertex_cube: " << X << " has an unknown parameter value" << std::endl;
      return false;
    }
    cube.push_back(static_cast<std::uint32_t>(param_index));
    ++i;
  }
  return true;
}

} // namespace detail

class pbesinst_symbolic_counter_example_structure_graph_algorithm : public pbesinst_structure_graph_algorithm
{
public:
  pbesinst_symbolic_counter_example_structure_graph_algorithm(structure_graph& G,
    const pbessolve_options& options,
    const pbes& p,
    bool _alpha,
    const std::unordered_map<core::identifier_string, data::data_expression>& _propvar_map,
    const std::vector<symbolic::data_expression_index>& _data_index,
    const sylvan::ldds::ldd& Valpha_,
    const sylvan::ldds::ldd& S,
    std::optional<data::rewriter> rewriter = std::nullopt)
    : pbesinst_structure_graph_algorithm(options, p, G, rewriter),
      alpha(_alpha),
      strategy(S),
      Valpha(Valpha_),
      data_index(_data_index),
      propvar_map(_propvar_map),
      X_false(p.equations()[p.equations().size() - 2].variable().name()),
      X_true(p.equations()[p.equations().size() - 1].variable().name())
  {}

  std::function<pbes_expression(const propositional_variable_instantiation&)> phi_substitution(
    const std::size_t thread_index,
    const fixpoint_symbol& symbol,
    const propositional_variable_instantiation& X,
    const pbes_expression& phi) override
  {
    return compose_substitutions(X_false_X_true_substitution(X_false, X_true),
      compose_substitutions(rewrite_star_substitution(data_index, propvar_map, strategy, Valpha, X, alpha),
        pbesinst_structure_graph_algorithm::phi_substitution(thread_index, symbol, X, phi)));
  }

private:
  bool alpha;
  sylvan::ldds::ldd strategy;
  sylvan::ldds::ldd Valpha;
  const std::vector<symbolic::data_expression_index>& data_index;
  const std::unordered_map<core::identifier_string, data::data_expression>& propvar_map;
  const core::identifier_string& X_false;
  const core::identifier_string& X_true;

  struct X_false_X_true_substitution
  {
    const core::identifier_string& X_false;
    const core::identifier_string& X_true;

    X_false_X_true_substitution(const core::identifier_string& X_false, const core::identifier_string& X_true)
      : X_false(X_false),
        X_true(X_true)
    {}

    pbes_expression operator()(const propositional_variable_instantiation& x)
    {
      if (x.name() == X_false)
      {
        return false_();
      }
      else if (x.name() == X_true)
      {
        return true_();
      }
      else
      {
        return x;
      }
    }
  };

  /// Removes PBES expressions that are irrelevant w.r.t the given strategy
  struct rewrite_star_substitution
  {
    mutable std::vector<std::uint32_t> singleton;

    const std::vector<symbolic::data_expression_index>& data_index;
    const std::unordered_map<core::identifier_string, data::data_expression>& propvar_map;
    const sylvan::ldds::ldd& strategy;
    const sylvan::ldds::ldd& Valpha;
    const propositional_variable_instantiation& X;
    const bool alpha;

    // Everything below depends only on X, so it is computed once in the constructor instead of
    // once per successor Y (phi_substitution, and hence this substitution, is constructed once
    // per equation, but operator() below is called once per successor in the right-hand side).
    std::vector<std::uint32_t> m_X_cube;
    bool m_X_is_counter_example; // X is a counter example equation, i.e. in L
    bool m_X_known;              // X was resolved to a vertex of the symbolic game
    bool m_X_is_alpha;

    rewrite_star_substitution(const std::vector<symbolic::data_expression_index>& data_index,
      const std::unordered_map<core::identifier_string, data::data_expression>& propvar_map,
      const sylvan::ldds::ldd& strategy,
      const sylvan::ldds::ldd& Valpha,
      const propositional_variable_instantiation& X,
      bool alpha)
      : data_index(data_index),
        propvar_map(propvar_map),
        strategy(strategy),
        Valpha(Valpha),
        X(X),
        alpha(alpha)
    {
      // TODO: This depends on the encoding used in pbesreach.
      // The counter example equations (those in L) are not part of the symbolic game, so they are
      // knowingly absent from propvar_map. Recognising them here keeps a failure of vertex_cube a
      // reliable indication that the two instantiations have gone out of step.
      m_X_is_counter_example = mcrl2::pbes_system::detail::is_counter_example_name(X.name());
      m_X_known = !m_X_is_counter_example && detail::vertex_cube(X, data_index, propvar_map, m_X_cube);
      m_X_is_alpha = m_X_known && sylvan::ldds::member_cube(Valpha, m_X_cube);
    }

    pbes_expression operator()(const propositional_variable_instantiation& Y) const
    {
      // the rewrite_star substitution is only applicable to closed PVIs.
      if (!find_free_variables(Y).empty())
      {
        mCRL2log(log::trace) << "rewrite_star " << Y << " contains free variables, not applying substitution\n";
        return Y;
      }

      if (mcrl2::pbes_system::detail::is_counter_example_name(Y.name()))
      {
        // If Y in L return Y
        mCRL2log(log::debug) << "rewrite_star " << Y << " is counter example equation (in L)" << std::endl;
        return Y;
      }

      if (!m_X_known)
      {
        // X could not be resolved to a known vertex in the symbolic exploration. This is a second,
        // independent source of blow-up if it ever happens; preserve the pre-existing behaviour of
        // treating X as not belonging to alpha, i.e. do not prune. For counter example equations
        // this is expected and not worth reporting.
        if (!m_X_is_counter_example)
        {
          mCRL2log(log::debug) << "rewrite_star " << X << " could not be resolved to a known vertex, not pruning"
                               << std::endl;
        }
        return Y;
      }

      if (m_X_is_alpha)
      {
        // Determine whether (X, Y) is in the strategy.

        // Add the propositional variables.
        singleton.clear();
        singleton.emplace_back(m_X_cube[0]);
        singleton.emplace_back(data_index[0].index(propvar_map.at(Y.name())));

        // Add the interleaved data expressions.
        std::size_t i = 1;
        auto param_Y_it = Y.parameters().begin();

        for (std::size_t k = 1; k < m_X_cube.size(); ++k)
        {
          singleton.emplace_back(m_X_cube[k]);
          singleton.emplace_back(data_index[i].index(*param_Y_it));

          ++param_Y_it;
          ++i;
        }

        if (sylvan::ldds::member_cube(strategy, singleton))
        {
          // If Y in E0
          mCRL2log(log::debug) << "rewrite_star " << Y << " is reachable" << std::endl;
          return Y;
        }
        else
        {
          if (alpha == 0)
          {
            // If Y is not reachable, replace it by false
            mCRL2log(log::debug) << "rewrite_star " << Y << " is not reachable, becomes false" << std::endl;
            return false_();
          }
          else
          {
            // If Y is not reachable, replace it by true
            mCRL2log(log::debug) << "rewrite_star " << Y << " is not reachable, becomes true" << std::endl;
            return true_();
          }
        }
      }
      else
      {
        mCRL2log(log::debug) << "rewrite_star " << Y << " is reachable" << std::endl;
        return Y;
      }
    }
  };
};

} // namespace mcrl2::pbes_system

class pbessolvesymbolic_tool : public parallel_tool<rewriter_tool<input_output_tool>>
{
  using super = parallel_tool<rewriter_tool<input_output_tool>>;

protected:
  pbes_system::symbolic_reachability_options options;

  // Lace options
  std::size_t lace_dqsize = static_cast<std::size_t>(1024 * 1024 * 4); // set large default
  std::size_t lace_stacksize = 0; // use default

  // Sylvan options
  std::size_t memory_limit = 3;
  std::size_t initial_ratio = 16;
  std::size_t table_ratio = 1;

  // Counter example options.
  std::string lpsfile;
  std::string ltsfile;
  std::string evidence_file;

  void add_options(utilities::interface_description& desc) override
  {
    super::add_options(desc);
    desc.add_option("lace-dqsize",
      utilities::make_optional_argument("NUM", "4194304"),
      "set the length of Lace task queue (default 1024*1024*4)");
    desc.add_option("lace-stacksize",
      utilities::make_optional_argument("NUM", "0"),
      "set the size of Sylvan program stack in gigabytes (0=default stack size). "
      "This is the main stack for all calculations. If it is too small a bus error occurs. ");
    desc.add_option("memory-limit",
      utilities::make_optional_argument("NUM", "3"),
      "Sylvan memory limit in gigabytes (default 3)",
      'm');

    desc.add_option("cached", "use transition group caching to speed up state space exploration");
    desc.add_option("chaining",
      "reduce the amount of breadth-first iterations by applying the transition groups consecutively");
    desc.add_option("groups",
      utilities::make_optional_argument("GROUPS", "none"),
      "'none' (default) no summand groups\n"
      "'used' summands with the same variables are joined\n"
      "'simple' summands with the same read/write variables are joined\n"
      "a user defined list of summand groups separated by semicolons, e.g. '0; 1 3 4; 2 5'");
    desc.add_option("reorder",
      utilities::make_optional_argument("ORDER", "none"),
      "'none' (default) no variable reordering\n"
      "'random' variables are put in a random order\n"
      "'weighted' variables are put in an order defined by their connectivity weight\n"
      "'a user defined permutation e.g. '1 3 2 0 4'");
    desc.add_option("info", "print read/write information of the summands");
    desc.add_option("max-iterations",
      utilities::make_optional_argument("NUM", "0"),
      "limit number of breadth-first iterations to NUM");
    desc.add_option("print-exact",
      "prints the sizes of LDDs exactly when within the representable range, and in scientific notation otherwise");
    desc.add_option("print-nodesize",
      "print the number of LDD nodes in addition to the number of elements represented as 'elements[nodes]'");
    desc.add_option("saturation",
      "reduce the amount of breadth-first iterations by applying the transition groups until fixed point");
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
      "Use solve strategy NUM. All strategies except 0 periodically apply on-the-fly solving, which may lead to early "
      "termination.",
      's');
    desc.add_option("split-conditions",
      "split disjunctive conditions to obtain more summands with potentially less dependencies",
      'c');
    desc.add_option("total", "make the SRF PBES total", 't');
    desc.add_option("reset", "set constant values when introducing parameters");

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
    desc.add_option("evidence-file",
      utilities::make_file_argument("NAME"),
      "The file to which the evidence is written. If not set, a "
      "default name will be chosen.");

    desc.add_hidden_option("aggressive", "apply on-the-fly solving after every iteration to detect bugs");
    desc.add_hidden_option("check-strategy", "do a sanity check on the computed strategy", 'y');
    desc.add_hidden_option("no-remove-unused-rewrite-rules", "do not remove unused rewrite rules. ", 'u');
    desc.add_hidden_option("no-one-point-rule-rewrite", "do not apply the one point rule rewriter");
    desc.add_hidden_option("no-discard", "do not discard any parameters");
    desc.add_hidden_option("no-read", "do not discard only-read parameters");
    desc.add_hidden_option("no-write", "do not discard only-write parameters");
    desc.add_hidden_option("no-relprod", "use an inefficient alternative version of relprod (for debugging)");
    desc.add_hidden_option("initial-ratio",
      utilities::make_optional_argument("NUM", "16"),
      "power-of-two ratio of initial and maximum table size (default 16)");
    desc.add_hidden_option("table-ratio",
      utilities::make_optional_argument("NUM", "16"),
      "power-of-two ratio of node table and cache table (default 1)");
    desc.add_hidden_option("srf",
      utilities::make_optional_argument("FILE", ""),
      "save the preprocessed PBES in SRF format");
    desc.add_hidden_option("dot",
      utilities::make_optional_argument("FILE", ""),
      "print the LDD of the parity game in dot format");
    desc.add_hidden_option("split-conditions-unsafe",
      utilities::make_optional_argument("NUM", "0"),
      "split conditions to obtain more summands (and equations) with potentially less dependencies\n"
      "0 (default) no splitting performed.\n"
      "1 only split disjunctive conditions, same as --split-conditions.\n"
      "2 also split conjunctive conditions into multiple equations which weakens guards and introduces more reachable "
      "BES equations. Note that splitting conditions can lead to expressions that cannot be rewritten if the equations "
      "are not sufficiently complete.\n"
      "3 alternative split for conjunctive conditions where even more states can become reachable.");
    desc.add_hidden_option("naive-counter-example-instantiation",
      "run the naive instantiation algorithm for pbes with counter example information");
  }

  void parse_options(const utilities::command_line_parser& parser) override
  {
    super::parse_options(parser);
    options.aggressive = parser.has_option("aggressive");
    options.cached = parser.has_option("cached");
    options.chaining = parser.has_option("chaining");
    options.check_strategy = parser.has_option("check-strategy");
    options.one_point_rule_rewrite = !parser.has_option("no-one-point-rule-rewrite");
    options.print_exact = parser.has_option("print-exact");
    options.print_nodesize = parser.has_option("print-nodesize");
    options.remove_unused_rewrite_rules = !parser.has_option("no-remove-unused-rewrite-rules");
    options.replace_constants_by_variables = false; // This option doesn't work in the current implementation
    options.saturation = parser.has_option("saturation");
    options.no_discard = parser.has_option("no-discard");
    options.no_discard_read = parser.has_option("no-read");
    options.no_discard_write = parser.has_option("no-write");
    options.no_relprod = parser.has_option("no-relprod");
    options.info = parser.has_option("info");
    options.summand_groups = parser.option_argument("groups");
    options.variable_order = parser.option_argument("reorder");
    options.make_total = parser.has_option("total");
    options.reset_parameters = parser.has_option("reset");
    options.naive_counter_example_instantiation = parser.has_option("naive-counter-example-instantiation");
    if (!options.make_total)
    {
      options.detect_deadlocks = true; // This is a required setting if the pbes is not total.
    }
    options.srf = parser.option_argument("srf");
    options.rewrite_strategy = rewrite_strategy();
    options.dot_file = parser.option_argument("dot");
    options.max_workers = number_of_threads();
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

    options.solve_strategy = parser.option_argument_as<int>("solve-strategy");
    if (options.solve_strategy > 7)
    {
      throw mcrl2::runtime_error("Invalid strategy " + std::to_string(options.solve_strategy));
    }

    if (parser.has_option("max-iterations"))
    {
      options.max_iterations = parser.option_argument_as<std::size_t>("max-iterations");
    }

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
        throw mcrl2::runtime_error("Option --evidence-file cannot be used without option --file");
      }
      evidence_file = parser.option_argument("evidence-file");
    }

    if (options.check_strategy)
    {
      if (options.summand_groups.compare("none") != 0)
      {
        throw mcrl2::runtime_error("Cannot check strategy for merged summand groups");
      }
      options.compute_strategy = true;
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
        "which is then solved using Zielonka's algorithm. ")
  {}

  bool run() override
  {
    lace_set_stacksize(lace_stacksize*1024*1024*1024);
    lace_start(number_of_threads(), lace_dqsize);
    sylvan::sylvan_set_limits(memory_limit * 1024 * 1024 * 1024, static_cast<int>(std::log2(table_ratio)), static_cast<int>(std::log2(initial_ratio)));
    sylvan::sylvan_init_package();
    sylvan::sylvan_init_ldd();
    sylvan::ldds::initialise();

    auto args = arguments{.options = options,
      .input_filename = input_filename(),
      .output_filename = output_filename(),
      .evidence_filename = evidence_file,
      .lpsfile = lpsfile,
      .ltsfile = ltsfile,
      .timer = timer()};
    pbessolvesymbolic_task(&args);

    sylvan::sylvan_quit();
    lace_stop();
    return true;
  }
};

template<typename PbesReachAlgorithm, typename PbesInstAlgorithm>
void solve(pbes_system::pbes pbesspec,
  symbolic_reachability_options& options_,
  const std::string& input_filename,
  const std::string& evidence_filename,
  const std::string& lpsfile,
  const std::string& ltsfile,
  mcrl2::utilities::execution_timer& timer)
{
  using namespace sylvan::ldds;

  bool has_counter_example = mcrl2::pbes_system::detail::has_counter_example_information(pbesspec);
  if (has_counter_example && (options_.solve_strategy == 5 || options_.solve_strategy == 6))
  {
    // TODO: Cannot use the partial solvers.
    mCRL2log(mcrl2::log::warning)
      << "Warning: Cannot use partial solving using fatal attractor solving (solve strategies 5 and 6) with PBES that "
         "has counter example information, using solving strategy 0 instead."
      << std::endl;
    options_.solve_strategy = 0;
  }

  if (has_counter_example)
  {
    if (lpsfile.empty() && ltsfile.empty())
    {
      mCRL2log(log::warning)
        << "Warning: the PBES has counter example information, but no witness will be generated due to lack of --file"
        << std::endl;
    }
  }
  else if (!lpsfile.empty() || !ltsfile.empty())
  {
    mCRL2log(log::warning) << "Warning: the PBES has no counter example information. Did you "
                              "use the"
                              " --counter-example option when generating the PBES?"
                           << std::endl;
  }

  // This has to be done consistently with the LPS for the counter examples.
  data::mutable_map_substitution<> sigma = pbes_system::detail::instantiate_global_variables(pbesspec);
  pbes_system::detail::replace_global_variables(pbesspec, sigma);
  pbes_system::srf_pbes_with_ce pre_srf_pbes = preprocess(pbesspec, options_);

  mCRL2log(log::trace) << "============== Pre-SRF PBES ==============" << std::endl;
  mCRL2log(log::trace) << pre_srf_pbes.to_pbes() << std::endl;

  pbes_system::srf_pbes srf_pbes = pre_srf2srfpbes(pre_srf_pbes);

  pbesspec = pre_srf_pbes.to_pbes();
  if (options_.info)
  {
    PbesReachAlgorithm reach(srf_pbes, options_);
    std::cout << symbolic::print_read_write_patterns(reach.read_write_group_patterns());
  }
  else
  {
    // If you provide a file, but the PBES has no counter example information, then use the two pass instantiation. This
    // will be useless, but at least the file will be written.
    if ((!has_counter_example && lpsfile.empty() && ltsfile.empty()) || options_.naive_counter_example_instantiation)
    {
      PbesReachAlgorithm reach(srf_pbes, options_);
      mCRL2log(log::debug) << pbes_system::detail::print_pbes_info(reach.pbes()) << std::endl;

      timer.start("instantiation");
      reach.run();
      timer.finish("instantiation");
      if (!options_.dot_file.empty())
      {
        print_dot(options_.dot_file, reach.V());
      }

      if (reach.solution_found())
      {
        std::cout << (includes(reach.partial_solution().winning[0], (reach.initial_state())) ? "true" : "false")
                  << std::endl;
      }
      else
      {
        if (options_.max_iterations == 0)
        {
          pbes_system::symbolic_parity_game G(reach.pbes(),
            reach.summand_groups(),
            reach.data_index(),
            reach.V(),
            options_.no_relprod,
            options_.chaining,
            options_.compute_strategy);
          G.print_information();
          pbes_system::symbolic_pbessolve_algorithm solver(G, options_.check_strategy, options_.compute_strategy);

          mCRL2log(log::debug) << pbes_system::detail::print_pbes_info(reach.pbes()) << std::endl;
          timer.start("solving");
          auto [result, solution]
            = solver.solve(reach.initial_state(), reach.V(), reach.deadlocks(), reach.partial_solution());
          timer.finish("solving");

          std::cout << (result ? "true" : "false") << std::endl;
        }
        else
        {
          // TODO: We could actually try to solve the incomplete parity game.
          std::cout << "Skipped solving since exploration was limited to max-iterations" << std::endl;
        }
      }
    }
    else
    {
      // We are generating a counterexample, so the strategy must be computed, irregardless of
      // whether we check the strategy afterwards.
      options_.compute_strategy = true;

      PbesReachAlgorithm reach(srf_pbes, options_);

      timer.start("first-instantiation");
      reach.run();
      timer.finish("first-instantiation");

      // Instantiation may lead to a partially explored parity game;
      // V represents all seen vertices, reach.V() are the explored, and reach.I() the unexplored
      // vertices
      ldd V = union_(reach.V(), reach.I());

      if (!options_.dot_file.empty())
      {
        print_dot(options_.dot_file, V);
      }

      pbes_system::symbolic_parity_game G(reach.pbes(),
        reach.summand_groups(),
        reach.data_index(),
        V,
        options_.no_relprod,
        options_.chaining,
        options_.compute_strategy);
      G.print_information();
      pbes_system::symbolic_pbessolve_algorithm solver(G, options_.check_strategy, options_.compute_strategy);

      timer.start("first-solving");
      // Solve the remainder of the partially solved game.
      bool solution_found = false;
      bool result;
      symbolic_solution_t solution(true);
      if (reach.I() == empty_set())
      {
        std::tie(result, solution)
          = solver.solve(reach.initial_state(), V, reach.deadlocks(), reach.partial_solution());
        solution_found = true;
      }
      else
      {
        solution
          = solver.partial_solve(reach.initial_state(), V, reach.I(), reach.deadlocks(), reach.partial_solution());

        if (includes(solution.winning[0], reach.initial_state()))
        {
          solution_found = true;
          result = true;
        }
        else if (includes(solution.winning[1], reach.initial_state()))
        {
          solution_found = true;
          result = false;
        }
      }
      timer.finish("first-solving");

      if (!solution_found)
      {
        std::cout << "Exploration was limited to max-iterations, and the partially explored parity game does not "
                     "contain enough information to compute the solution."
                  << std::endl;
      }
      else
      {
        mCRL2log(log::log_level_t::verbose) << (result ? "true" : "false") << std::endl;

        // Based on the result remove the unnecessary equations related to counter example information.
        mCRL2log(log::verbose) << "Removing unnecessary counter example information for other player." << std::endl;
        pbes_system::pbes pbesspec_simplified
          = mcrl2::pbes_system::detail::remove_counterexample_info(pbesspec, !result, result);
        mCRL2log(log::trace) << pbesspec_simplified << std::endl;

        structure_graph SG;

        // Set some options for the second instantiation.
        pbessolve_options pbessolve_options;
        // All optimizations disabled for the second run. They are not needed due to the
        // availability of a winning strategy
        pbessolve_options.rewrite_strategy = options_.rewrite_strategy;
        pbessolve_options.remove_unused_rewrite_rules = options_.remove_unused_rewrite_rules;
        pbessolve_options.check_strategy = options_.check_strategy;
        pbessolve_options.number_of_threads
          = 1; // If we spawn multiple threads here, the threads of Sylvan and the explicit exploration will interfere

        // At this point compute_strategy=true, so both strategies must have been populated by solve().
        if (!solution.strategy[0].has_value() || !solution.strategy[1].has_value())
        {
          throw mcrl2::runtime_error("Expected strategies to be computed, but they were not.");
        }

        PbesInstAlgorithm second_instantiate(SG,
          pbessolve_options,
          pbesspec_simplified,
          !result,
          reach.propvar_map(),
          reach.data_index(),
          G.players(V)[result ? 0 : 1],
          result ? *solution.strategy[0] : *solution.strategy[1], // NOLINT(bugprone-unchecked-optional-access)
          reach.rewriter());

        // Perform the second instantiation given the proof graph.
        timer.start("second-instantiation");
        second_instantiate.run();
        timer.finish("second-instantiation");

        mCRL2log(log::verbose) << "Number of vertices in the structure graph: " << SG.all_vertices().size()
                               << std::endl;
        [[maybe_unused]]
        bool final_result
          = pbes_system::detail::run_solve(pbesspec,
            sigma,
            SG,
            second_instantiate.equation_index(),
            pbessolve_options,
            input_filename,
            lpsfile,
            ltsfile,
            evidence_filename,
            timer);
        if (result != final_result)
        {
          throw mcrl2::runtime_error(
            "The result of the first and second instantiations do not match, this is a bug in the tool!");
        }
      }
    }
  }
}

TASK_IMPL_1(bool, pbessolvesymbolic_task, arguments*, args) // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
{
  mCRL2log(log::verbose) << args->options << std::endl;

  pbes_system::pbes pbesspec = pbes_system::detail::load_pbes(args->input_filename);

  if (pbesspec.initial_state().empty())
  {
    throw mcrl2::runtime_error("PBESses without parameters are not supported");
  }
  else
  {
    if (args->options.solve_strategy == 0)
    {
      solve<pbesreach_algorithm, pbesinst_symbolic_counter_example_structure_graph_algorithm>(pbesspec,
        args->options,
        args->input_filename,
        args->evidence_filename,
        args->lpsfile,
        args->ltsfile,
        args->timer);
    }
    else
    {
      solve<pbesreach_algorithm_partial, pbesinst_symbolic_counter_example_structure_graph_algorithm>(pbesspec,
        args->options,
        args->input_filename,
        args->evidence_filename,
        args->lpsfile,
        args->ltsfile,
        args->timer);
    }
  }

  return true;
}

int main(int argc, char* argv[])
{
  return pbessolvesymbolic_tool().execute(argc, argv);
}
