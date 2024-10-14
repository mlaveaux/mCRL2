// Author(s): Muck van Weerdenburg, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file ltsconvert.cpp

#define NAME "ltsconvert"
#define AUTHOR "Muck van Weerdenburg, Jan Friso Groote"

#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/lts/lts_io.h"
#include "mcrl2/lts/lts_algorithm.h"

using namespace mcrl2::lts;
using namespace mcrl2::lts::detail;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::utilities;
using namespace mcrl2::core;
using namespace mcrl2::log;

class t_tool_options
{
  public:
    std::string     infilename;
    std::string     outfilename;
    std::string     lpsfile;
    lts_type        intype;
    lts_type        outtype;
    lts_equivalence equivalence;
    std::vector<std::string> tau_actions;   // Actions with these labels must be considered equal to tau.
    bool            remove_state_information=false;
    bool            determinise=false;
    bool            check_reach=true;
    bool            add_state_as_state_label=false;

    inline t_tool_options() 
     : intype(lts_none), 
       outtype(lts_none), 
       equivalence(lts_eq_none)
    {
    }

    inline std::string source_string() const
    {
      return (infilename.empty()) ? std::string("standard input") :
             std::string("'" + infilename + "'");
    }

    inline std::string target_string() const
    {
      return (outfilename.empty()) ? std::string("standard output") :
             std::string("'" + outfilename + "'");
    }

    void set_source(std::string const& filename)
    {
      infilename = filename;

      if (intype==lts_none)
      {
        intype = mcrl2::lts::detail::guess_format(infilename,false);
      }

      if (intype == lts_none)
      {
        intype = lts_aut;
        mCRL2log(warning) << "Cannot determine type of input. Assuming .aut.\n";
      }
      // When there is no equivalence and determinisation is not applied, the input lts can be probabilistic. 
      if (equivalence == lts_eq_none && !determinise)
      {
        if (intype == lts_lts)
        {
          intype = lts_lts_probabilistic;
        }
        else if (intype == lts_aut)
        {
          intype = lts_aut_probabilistic;
        }
        else if (intype == lts_fsm)
        {
          intype = lts_fsm_probabilistic;
        }
      }
    }

    void set_target(std::string const& filename)
    {
      outfilename = filename;

      if (outtype == lts_none)
      {
        mCRL2log(verbose) << "Trying to detect output format by extension..." << std::endl;

        outtype = mcrl2::lts::detail::guess_format(outfilename,true);

        if (outtype == lts_none)
        {
          if (!lpsfile.empty())
          {
            mCRL2log(warning) << "no output format set; using fsm because --lps was used" << std::endl;
            outtype = lts_fsm;
          }
          else
          {
            mCRL2log(warning) << "no output format set or detected; using default (mcrl2)" << std::endl;
            outtype = lts_lts;
          }
        }
        // When there is no equivalence and determinisation is not applied, the output lts can be probabilistic. 
        if (equivalence == lts_eq_none && !determinise)
        {
          if (outtype == lts_lts)
          {
            outtype = lts_lts_probabilistic;
          }
          else if (outtype == lts_aut)
          {
            outtype = lts_aut_probabilistic;
          }
          else if (outtype == lts_fsm)
          {
            outtype = lts_fsm_probabilistic;
          }
        }
      }
    }

};

class ltsconvert_tool : public input_output_tool
{
  private:
    t_tool_options tool_options;

  public:
    ltsconvert_tool() :
      input_output_tool(NAME,AUTHOR,
                      "convert and optionally minimise an LTS",
                      "Convert the labelled transition system (LTS) from INFILE to OUTFILE in the\n"
                      "requested format after applying the selected minimisation method (default is\n"
                      "none). If OUTFILE is not supplied, stdout is used. If INFILE is not supplied,\n"
                      "stdin is used.\n"
                      "\n"
                      "The output format is determined by the extension of OUTFILE, whereas the input\n"
                      "format is determined by the content of INFILE. Options --in and --out can be\n"
                      "used to force the input and output formats. The supported formats are:\n"
                      + mcrl2::lts::detail::supported_lts_formats_text(lts_lts)
                     )
    {
    }

  private:

    template < class LTS_TYPE >
    bool load_convert_and_save()
    {
      using namespace mcrl2::lts;
      using namespace mcrl2::lts::detail;

      LTS_TYPE l;
      l.load(tool_options.infilename);
      l.apply_hidden_actions(tool_options.tau_actions);

      if (tool_options.check_reach)
      {
        timer().start("reachability check");
        reachability_check(l, true); // Remove unreachable states from the input transition system.
        timer().finish("reachability check");
      }

      if (tool_options.remove_state_information)
      {
        l.clear_state_labels();
      }

      if (tool_options.add_state_as_state_label)
      {
        // Add the state numbers as the labels of the states. 
        if constexpr (std::is_same<LTS_TYPE,probabilistic_lts_fsm_t>::value ||
                      std::is_same<LTS_TYPE,lts_fsm_t>::value)
        {
          l.clear_process_parameters();
          l.add_process_parameter("state_number","Nat");
          for(std::size_t i=0; i<l.num_states(); ++i)
          {
            l.add_state_element_value(0,std::to_string(i));
          }
        }
        if constexpr (std::is_same<LTS_TYPE,probabilistic_lts_lts_t>::value ||
                      std::is_same<LTS_TYPE,lts_lts_t>::value)
        {
          // l.set_process_parameters(atermpp::term_list<mcrl2::data::variable>(mcrl2::data::variable("state_number",mcrl2::data::sort_nat::nat())));
          l.set_process_parameters({mcrl2::data::variable("state_number",mcrl2::data::sort_nat::nat())});
        }
        l.add_state_number_as_state_information();
      }

      if (tool_options.equivalence != lts_eq_none)
      {
        if constexpr (LTS_TYPE::is_probabilistic_lts)
        {
           throw mcrl2::runtime_error("The LTS is probabilistic. The reduction cannot be applied. Consider the tool ltspbisim instead.");
        }
        else 
        {
          mCRL2log(verbose) << "Reducing LTS (modulo " <<  description(tool_options.equivalence) << ")..." << std::endl;
          mCRL2log(verbose) << "Before reduction: " << l.num_states() << " states and " << l.num_transitions() << " transitions." << std::endl;
          reduce(l,tool_options.equivalence);
          mCRL2log(verbose) << "After reduction: " << l.num_states() << " states and " << l.num_transitions() << " transitions." << std::endl;
        }
      }

      if (tool_options.determinise)
      {
        if constexpr (LTS_TYPE::is_probabilistic_lts)
        {
           throw mcrl2::runtime_error("The LTS is probabilistic. It cannot be determinised.");
        }
        else 
        {
          mCRL2log(verbose) << "determinising LTS..." << std::endl;
          mCRL2log(verbose) << "before determinisation: " << l.num_states() << " states and " << l.num_transitions() << " transitions" << std::endl;
          timer().start("determinisation");
          determinise(l);
          timer().finish("determinisation");
          mCRL2log(verbose) << "after determinisation: " << l.num_states() << " states and " << l.num_transitions() << " transitions" << std::endl;
        }
      }

      mcrl2::lps::stochastic_specification spec;

      if (!tool_options.lpsfile.empty())
      {
        // No lpsfile is given. Only straightforward translations are possible.
        load_lps(spec, tool_options.lpsfile);
      }

      switch (tool_options.outtype)
      {
        case lts_lts:
        {
          lts_lts_t l_out;
          lts_convert(l,l_out,spec.data(),spec.action_labels(),spec.process().process_parameters(),!tool_options.lpsfile.empty());
          l_out.save(tool_options.outfilename);
          return true;
        }
        case lts_lts_probabilistic:
        {
          probabilistic_lts_lts_t l_out;
          lts_convert(l,l_out,spec.data(),spec.action_labels(),spec.process().process_parameters(),!tool_options.lpsfile.empty());
          l_out.save(tool_options.outfilename);
          return true;
        }
        case lts_none:
          mCRL2log(warning) << "Cannot determine type of output. Assuming .aut.\n";
          [[fallthrough]];
        case lts_aut:
        {
          lts_aut_t l_out;
          lts_convert(l,l_out,spec.data(),spec.action_labels(),spec.process().process_parameters(),!tool_options.lpsfile.empty());
          l_out.save(tool_options.outfilename);
          return true;
        }
        case lts_aut_probabilistic:
        {
          probabilistic_lts_aut_t l_out;
          lts_convert(l,l_out,spec.data(),spec.action_labels(),spec.process().process_parameters(),!tool_options.lpsfile.empty());
          l_out.save(tool_options.outfilename);
          return true;
        }
        case lts_fsm:
        {
          lts_fsm_t l_out;
          lts_convert(l,l_out,spec.data(),spec.action_labels(),spec.process().process_parameters(),!tool_options.lpsfile.empty());
          l_out.save(tool_options.outfilename);
          return true;
        }
        case lts_fsm_probabilistic:
        {
          probabilistic_lts_fsm_t l_out;
          lts_convert(l,l_out,spec.data(),spec.action_labels(),spec.process().process_parameters(),!tool_options.lpsfile.empty());
          l_out.save(tool_options.outfilename);
          return true;
        }
        case lts_dot:
        {
          lts_dot_t l_out;
          lts_convert(l,l_out,spec.data(),spec.action_labels(),spec.process().process_parameters(),!tool_options.lpsfile.empty());
          l_out.save(tool_options.outfilename);
          return true;
        }
      }
      return true;
    }

  public:
    bool run()
    {
      switch (tool_options.intype)
      {
        case lts_lts:
        {
          return load_convert_and_save<lts_lts_t>();
        }
        case lts_lts_probabilistic:
        {
          return load_convert_and_save<probabilistic_lts_lts_t>();
        }
        case lts_none:
        case lts_aut:
        {
          return load_convert_and_save<lts_aut_t>();
        }
        case lts_aut_probabilistic:
        {
          return load_convert_and_save<probabilistic_lts_aut_t>();
        }
        case lts_fsm:
        {
          return load_convert_and_save<lts_fsm_t>();
        }
        case lts_fsm_probabilistic:
        {
          return load_convert_and_save<probabilistic_lts_fsm_t>();
        }
        case lts_dot:
        {
          throw mcrl2::runtime_error("Cannot read a .dot file anymore.");
        }
      }
      return true;
    }

  protected:
    void add_options(interface_description& desc)
    {
      input_output_tool::add_options(desc);

      desc.add_option("no-reach",
                      "do not perform a reachability check on the input LTS.");
      desc.add_option("no-state",
                      "remove the state information. This can be useful when state labels are huge.", 'n');
      desc.add_option("determinise", "determinise LTS", 'D');
      desc.add_option("lps", make_file_argument("FILE"),
                      "use FILE as the LPS from which the input LTS was generated; this might "
                      "be needed to store the correct parameter names of states when saving "
                      "in fsm format and to convert non-mCRL2 LTSs to a mCRL2 LTS.", 'l');
      desc.add_option("in", make_mandatory_argument("FORMAT"),
                      "use FORMAT as the input format.", 'i').
      add_option("out", make_mandatory_argument("FORMAT"),
                 "use FORMAT as the output format.", 'o');
      desc.add_option("equivalence",make_enum_argument<lts_equivalence>("NAME")
                      .add_value(lts_eq_none, true)
                      .add_value(lts_eq_bisim)
                      .add_value(lts_eq_bisim_gv)
                      .add_value(lts_eq_bisim_gjkw)
#ifdef BRANCH_BIS_EXPERIMENT_JFG
                      .add_value(lts_eq_bisim_gj)
#endif
                      .add_value(lts_eq_bisim_sigref)
                      .add_value(lts_eq_branching_bisim)
                      .add_value(lts_eq_branching_bisim_gv)
                      .add_value(lts_eq_branching_bisim_gjkw)
#ifdef BRANCH_BIS_EXPERIMENT_JFG
                      .add_value(lts_eq_branching_bisim_gj)
#endif
                      .add_value(lts_eq_branching_bisim_martens)
                      .add_value(lts_eq_branching_bisim_sigref)
                      .add_value(lts_eq_divergence_preserving_branching_bisim)
                      .add_value(lts_eq_divergence_preserving_branching_bisim_gv)
                      .add_value(lts_eq_divergence_preserving_branching_bisim_gjkw)
#ifdef BRANCH_BIS_EXPERIMENT_JFG
                      .add_value(lts_eq_divergence_preserving_branching_bisim_gj)
#endif
                      .add_value(lts_eq_divergence_preserving_branching_bisim_sigref)
                      .add_value(lts_eq_weak_bisim)
                      .add_value(lts_eq_divergence_preserving_weak_bisim)
                      .add_value(lts_eq_sim)
                      .add_value(lts_eq_ready_sim)		      
                      .add_value(lts_eq_trace)
                      .add_value(lts_eq_weak_trace)
                      .add_value(lts_red_tau_star),
                      "generate an equivalent LTS, preserving equivalence NAME:"
                      , 'e');
      desc.add_option("tau", make_mandatory_argument("ACTNAMES"),
                      "consider actions with a name in the comma separated list ACTNAMES to "
                      "be internal (tau) actions in addition to those defined as such by "
                      "the input.");
      desc.add_hidden_option("add-state-as-state-label",
                             "add the state number as the label of the states in the input file, "
                             "and remove other state labels if they exist");
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

    void parse_options(const command_line_parser& parser)
    {
      input_output_tool::parse_options(parser);

      if (parser.options.count("lps"))
      {
        if (1 < parser.options.count("lps"))
        {
          mCRL2log(warning) << "multiple LPS files specified; can only use one\n";
        }

        tool_options.lpsfile = parser.option_argument("lps");
      }
      if (parser.options.count("in"))
      {
        if (1 < parser.options.count("in"))
        {
          mCRL2log(warning) << "multiple input formats specified; can only use one\n";
        }

        tool_options.intype = mcrl2::lts::detail::parse_format(parser.option_argument("in"));

        if (tool_options.intype == lts_none)
        {
          mCRL2log(warning) << "format '" << parser.option_argument("in") <<
                    "' is not recognised; option ignored" << std::endl;
        }
      }
      if (parser.options.count("out"))
      {
        if (1 < parser.options.count("out"))
        {
          mCRL2log(warning) << "multiple output formats specified; can only use one\n";
        }

        tool_options.outtype = mcrl2::lts::detail::parse_format(parser.option_argument("out"));

        if (tool_options.outtype == lts_none)
        {
          mCRL2log(warning) << "format '" << parser.option_argument("out") <<
                    "' is not recognised; option ignored" << std::endl;
        }
      }

      tool_options.equivalence = parser.option_argument_as<lts_equivalence>("equivalence");

      if (parser.options.count("tau"))
      {
        set_tau_actions(tool_options.tau_actions, parser.option_argument("tau"));
      }

      if (parser.options.count("add-state-as-state-label"))
      {
        if (tool_options.intype == lts_aut)
        {
          parser.error("cannot use --add-state-as-state-label on a .aut input file, as .aut files do not have state labels\n");
        }
        tool_options.add_state_as_state_label=true;
      }

      tool_options.determinise                       = 0 < parser.options.count("determinise");
      tool_options.check_reach                       = parser.options.count("no-reach") == 0;
      tool_options.remove_state_information          = parser.options.count("no-state") != 0;

      if (tool_options.determinise && (tool_options.equivalence != lts_eq_none))
      {
        parser.error("cannot use option -D/--determinise together with LTS reduction options\n");
      }

      if (2 < parser.arguments.size())
      {
        parser.error("too many file arguments");
      }
      else
      {
        if (0 < parser.arguments.size())
        {
          tool_options.set_source(parser.arguments[0]);
        }
        else
        {
          if (tool_options.intype == lts_none)
          {
            mCRL2log(warning) << "cannot detect format from stdin and no input format specified; assuming aut format" << std::endl;
            tool_options.intype = lts_aut;
          }
        }
        if (1 < parser.arguments.size())
        {
          tool_options.set_target(parser.arguments[1]);
        }
        else
        {
          if (tool_options.outtype == lts_none)
          {
            if (!tool_options.lpsfile.empty())
            {
              mCRL2log(warning) << "no output format set; using fsm because --lps was used" << std::endl;
              tool_options.outtype = lts_fsm;
            }
            else
            {
              mCRL2log(warning) << "no output format set or detected; using default (aut)" << std::endl;
              tool_options.outtype = lts_aut;
            }
          }
        }
      }
    }

};

int main(int argc, char** argv)
{
  return ltsconvert_tool().execute(argc,argv);
}
