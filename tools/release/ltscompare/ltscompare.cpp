// Author(s): Muck van Weerdenburg
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file ltscompare.cpp

#define NAME "ltscompare"
#define AUTHOR "Muck van Weerdenburg"

#include "mcrl2/utilities/input_tool.h"

#include "mcrl2/lts/lts_algorithm.h"
#include "mcrl2/lts/lts_io.h"

using namespace mcrl2::lts;
using namespace mcrl2::lts::detail;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::utilities;
using namespace mcrl2::core;
using namespace mcrl2::log;

struct t_tool_options
{
  std::string     name_for_first  = "";
  std::string     name_for_second = "";
  lts_type        format_for_first = lts_none;
  lts_type        format_for_second = lts_none;
  lts_equivalence equivalence = lts_eq_none;
  lts_preorder    preorder = lts_preorder::lts_pre_none;
  mcrl2::lps::exploration_strategy strategy = mcrl2::lps::es_breadth;
  std::vector<std::string> tau_actions;   // Actions with these labels must be considered equal to tau.
  bool generate_counter_examples = false;
  std::string counter_example_file = "";
  bool structured_output = false;
  bool enable_preprocessing      = true;
};

typedef  input_tool ltscompare_base;
class ltscompare_tool : public ltscompare_base
{
  private:
    t_tool_options tool_options;

    // Check whether preconditions w.r.t. equivalence and pre-order are satisfied.
    // This is needed to make sure wrappers can call this tool without the --help,
    // --equivalence or --preorder options
    void check_preconditions()
    {
      if (tool_options.equivalence != lts_eq_none && tool_options.preorder != lts_preorder::lts_pre_none)
      {
        throw mcrl2::runtime_error("options -e/--equivalence and -p/--preorder cannot be used simultaneously");
      }

      if (tool_options.equivalence == lts_eq_none && tool_options.preorder == lts_preorder::lts_pre_none)
      {
        throw mcrl2::runtime_error("one of the options -e/--equivalence and -p/--preorder must be used");
      }

      if (tool_options.name_for_first.empty() && tool_options.name_for_second.empty())
      {
        throw mcrl2::runtime_error("too few file arguments");
      }
    }

  public:
    ltscompare_tool() :
      ltscompare_base(NAME,AUTHOR,
                      "compare two LTSs",
                      "Determine whether or not the labelled transition systems (LTSs) in INFILE1 and INFILE2 are related by some equivalence or preorder. "
                      "If INFILE1 is not supplied, stdin is used. "
                      "If INFILE1 and/or INFILE2 is '-', stdin is used. "
                      "Reading two LTSs via stdin is only supported for the 'aut' format, these LTSs must be separated by an EOT character (\\x04).\n"
                      "\n"
                      "The input formats are determined by the contents of INFILE1 and INFILE2. "
                      "Options --in1 and --in2 can be used to force the input format of INFILE1 and INFILE2, respectively. "
                      "The supported formats are:\n"
                      + mcrl2::lts::detail::supported_lts_formats_text()
                     )
    {
    }

  private:

    template <class LTS_TYPE>
    bool lts_compare(void)
    {
      LTS_TYPE l1,l2;
      l1.load(tool_options.name_for_first);
      l2.load(tool_options.name_for_second);

      l1.record_hidden_actions(tool_options.tau_actions);
      l2.record_hidden_actions(tool_options.tau_actions);

      bool result = true;
      if (tool_options.equivalence != lts_eq_none)
      {
        mCRL2log(verbose) << "comparing LTSs using " <<
                     tool_options.equivalence << "..." << std::endl;

        result = destructive_compare(l1, l2, tool_options.equivalence, tool_options.generate_counter_examples, tool_options.counter_example_file, tool_options.structured_output);

        mCRL2log(info) << "LTSs are " << ((result) ? "" : "not ")
                       << "equal ("
                       << description(tool_options.equivalence) << ")\n";
      }

      if (tool_options.preorder != lts_preorder::lts_pre_none)
      {
        mCRL2log(verbose) << "comparing LTSs for " <<
                     description(tool_options.preorder) << "..."
                     " using the " << print_exploration_strategy(tool_options.strategy) << " strategy.\n";

        result = destructive_compare(l1, l2, tool_options.preorder, tool_options.generate_counter_examples, tool_options.counter_example_file, tool_options.structured_output, tool_options.strategy, tool_options.enable_preprocessing);

        if (!tool_options.structured_output)
        {
          mCRL2log(info) << "The LTS in " << tool_options.name_for_first
                         << " is " << ((result) ? "" : "not ")
                         << "included in"
                         << " the LTS in " << tool_options.name_for_second
                         << " (using " << description(tool_options.preorder)
                         << ")." << std::endl;
        }
      }

      std::cout << (tool_options.structured_output ? "result: " : "") << std::boolalpha << result << std::endl;
      return true; // The tool terminates in a correct way.
    }

  public:
    bool run() override
    {
      check_preconditions();

      if (tool_options.format_for_first==lts_none)
      {
        tool_options.format_for_first = guess_format(tool_options.name_for_first);
      }

      if (tool_options.format_for_second==lts_none)
      {
        tool_options.format_for_second = guess_format(tool_options.name_for_second);
      }

      if (tool_options.format_for_first!=tool_options.format_for_second)
      {
        throw mcrl2::runtime_error("The input labelled transition systems have different types");
      }

      switch (tool_options.format_for_first)
      {
        case lts_lts:
        case lts_lts_probabilistic:
        {
          return lts_compare<lts_lts_t>();
        }
        case lts_none:
          mCRL2log(mcrl2::log::warning) << "No input format is specified. Assuming .aut format.\n";
          [[fallthrough]];
        case lts_aut:
        case lts_aut_probabilistic:
        {
          return lts_compare<lts_aut_t>();
        }
        case lts_fsm:
        case lts_fsm_probabilistic:
        {
          return lts_compare<lts_fsm_t>();
        }
        case lts_dot:
        {
          throw mcrl2::runtime_error("Reading the .dot format is not supported anymore.");
        }
      }

      return true;
    }

  protected:
    std::string synopsis() const override
    {
      return "[OPTION]... [INFILE1] INFILE2";
    }

    void check_positional_options(const command_line_parser& parser) override
    {
      if (2 < parser.arguments.size())
      {
        parser.error("too many file arguments");
      }
    }

    void set_tau_actions(std::vector <std::string>& tau_actions, std::string const& act_names)
    {
      std::string::size_type lastpos = 0, pos;
      while ((pos = act_names.find(',',lastpos)) != std::string::npos)
      {
        tau_actions.push_back(act_names.substr(lastpos,pos-lastpos));
        lastpos = pos+1;
      }
      tau_actions.push_back(act_names.substr(lastpos));
    }

    void add_options(interface_description& desc) override
    {
      ltscompare_base::add_options(desc);

      desc.
      add_option("in1", make_mandatory_argument("FORMAT"),
                 "use FORMAT as the format for INFILE1 (or stdin)", 'i').
      add_option("in2", make_mandatory_argument("FORMAT"),
                 "use FORMAT as the format for INFILE2", 'j').
      add_option("equivalence", make_enum_argument<lts_equivalence>("NAME)")
                 .add_value(lts_eq_none, true)
                 .add_value(lts_eq_bisim)
                 .add_value(lts_eq_bisim_gv)
                 .add_value(lts_eq_bisim_gjkw)
#ifdef BRANCH_BIS_EXPERIMENT_JFG
                 .add_value(lts_eq_bisim_gj)
#endif    
                 .add_value(lts_eq_branching_bisim)
                 .add_value(lts_eq_branching_bisim_gv)
                 .add_value(lts_eq_branching_bisim_gjkw)
#ifdef BRANCH_BIS_EXPERIMENT_JFG
                 .add_value(lts_eq_branching_bisim_gj)
#endif    
                 .add_value(lts_eq_divergence_preserving_branching_bisim)
                 .add_value(lts_eq_divergence_preserving_branching_bisim_gv)
                 .add_value(lts_eq_divergence_preserving_branching_bisim_gjkw)
#ifdef BRANCH_BIS_EXPERIMENT_JFG
                 .add_value(lts_eq_divergence_preserving_branching_bisim_gj)
#endif    
                 .add_value(lts_eq_weak_bisim)
                 .add_value(lts_eq_divergence_preserving_weak_bisim)
                 .add_value(lts_eq_sim)
                 .add_value(lts_eq_ready_sim)
                 .add_value(lts_eq_trace)
                 .add_value(lts_eq_weak_trace)
                 .add_value(lts_eq_coupled_sim),
                 "use equivalence NAME (not allowed in combination with -p/--preorder):", 'e').
      add_option("preorder", make_enum_argument<lts_preorder>("NAME")
                 .add_value(lts_preorder::lts_pre_none, true)
                 .add_value(lts_preorder::lts_pre_sim)
                 .add_value(lts_preorder::lts_pre_ready_sim)
                 .add_value(lts_preorder::lts_pre_trace)
                 .add_value(lts_preorder::lts_pre_weak_trace)
                 .add_value(lts_preorder::lts_pre_trace_anti_chain)
                 .add_value(lts_preorder::lts_pre_weak_trace_anti_chain)
                 .add_value(lts_preorder::lts_pre_failures_refinement)
                 .add_value(lts_preorder::lts_pre_weak_failures_refinement)
                 .add_value(lts_preorder::lts_pre_failures_divergence_refinement)
                 .add_value(lts_preorder::lts_pre_impossible_futures),
                 "use preorder NAME (not allowed in combination with -e/--equivalence):", 'p').
      add_option("strategy", make_enum_argument<mcrl2::lps::exploration_strategy>("NAME")
                 .add_value_short(mcrl2::lps::es_breadth, "b", true)
                 .add_value_short(mcrl2::lps::es_depth, "d")
                 , "explore the state space using strategy NAME (only for antichain based algorithms; includes all failures refinements) :"
                 , 's').
      add_option("tau", make_mandatory_argument("ACTNAMES"),
                 "consider actions with a name in the comma separated list ACTNAMES to "
                 "be internal (tau) actions in addition to those defined as such by "
                 "the input").
      add_option("counter-example",
                 "generate counter example if the input lts's are not equivalent",'c').
      add_option("counter-example-file", mcrl2::utilities::make_file_argument("NAME"),
                 "the file to which the counterexample should be written");
      desc.add_hidden_option("structured-output",
                 "generate counter examples on stdout");
      desc.add_hidden_option("no-preprocessing",
                        "disable preprocessing applied to the input LTSs for refinement checking",'\0');
    }

    void parse_options(const command_line_parser& parser) override
    {
      ltscompare_base::parse_options(parser);

      tool_options.equivalence = parser.option_argument_as<lts_equivalence>("equivalence");
      tool_options.preorder = parser.option_argument_as<lts_preorder>("preorder");

      if (parser.has_option("counter-example") && parser.has_option("preorder"))
      {
        if (tool_options.preorder == lts_preorder::lts_pre_sim)
        {
          parser.error("counter examples cannot be used with simulation pre-order");
        }
        if (tool_options.preorder == lts_preorder::lts_pre_ready_sim)
        {
          parser.error("counter examples cannot be used with ready simulation pre-order");
        }
        if (tool_options.preorder == lts_preorder::lts_pre_trace)
        {
          parser.error("counter examples cannot be used with the plain trace pre-order (use trace-ac instead)");
        }
        if (tool_options.preorder == lts_preorder::lts_pre_weak_trace)
        {
          parser.error("counter examples cannot be used with the plain weak trace pre-order (use weak-trace-ac instead");
        }
      }

      if (parser.has_option("no-preprocessing"))
      {
        tool_options.enable_preprocessing = false;
      }

      if (parser.has_option("tau"))
      {
        set_tau_actions(tool_options.tau_actions, parser.option_argument("tau"));
      }

      tool_options.generate_counter_examples = parser.has_option("counter-example");

      if (parser.has_option("counter-example-file"))
      {
        tool_options.counter_example_file = parser.option_argument("counter-example-file");
      }

      tool_options.structured_output = parser.has_option("structured-output");

      if (parser.arguments.size() == 1)
      {
        tool_options.name_for_second = parser.arguments[0];
      }
      else if (parser.arguments.size() == 2)
      {
        tool_options.name_for_first  = parser.arguments[0];
        tool_options.name_for_second = parser.arguments[1];
      } // else something strange going on, caught in check_preconditions()

      if (parser.has_option("strategy"))
      {
        tool_options.strategy = mcrl2::lps::parse_exploration_strategy(parser.option_argument("strategy"));

        if (tool_options.strategy != mcrl2::lps::es_breadth && tool_options.generate_counter_examples)
        {
          mCRL2log(mcrl2::log::warning) << "Generated counter example might not be the shortest with the " << print_exploration_strategy(tool_options.strategy) << " strategy.\n";
        }

        if (tool_options.preorder != lts_preorder::lts_pre_trace_anti_chain
            && tool_options.preorder != lts_preorder::lts_pre_weak_trace_anti_chain
            && tool_options.preorder != lts_preorder::lts_pre_failures_refinement
            && tool_options.preorder != lts_preorder::lts_pre_weak_failures_refinement
            && tool_options.preorder != lts_preorder::lts_pre_failures_divergence_refinement)
        {
          parser.error("strategy can only be chosen for antichain based algorithms.");
        }
      }

      if (parser.has_option("in1"))
      {
        tool_options.format_for_first = mcrl2::lts::detail::parse_format(parser.option_argument("in1"));

        if (tool_options.format_for_first == lts_none)
        {
          mCRL2log(warning) << "format '" << parser.option_argument("in1") <<
                    "' is not recognised; option ignored" << std::endl;
        }
      }
      else if (!tool_options.name_for_first.empty())
      {
        tool_options.format_for_first = mcrl2::lts::detail::guess_format(tool_options.name_for_first);
      }
      else
      {
        mCRL2log(warning) << "cannot detect format from stdin and no input format specified; assuming aut format" << std::endl;
        tool_options.format_for_first = lts_aut;
      }
      if (parser.has_option("in2"))
      {
        tool_options.format_for_second = mcrl2::lts::detail::parse_format(parser.option_argument("in2"));

        if (tool_options.format_for_second == lts_none)
        {
          mCRL2log(warning) << "format '" << parser.option_argument("in2") <<
                    "' is not recognised; option ignored" << std::endl;
        }
      }
      else if (!tool_options.name_for_first.empty())
      {
        tool_options.format_for_second = mcrl2::lts::detail::guess_format(tool_options.name_for_second);
      }
    }

};

int main(int argc, char** argv)
{
  return ltscompare_tool().execute(argc,argv);
}
