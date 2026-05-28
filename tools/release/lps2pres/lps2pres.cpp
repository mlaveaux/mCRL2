// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lps2pres.cpp
/// \brief This tool translates an .lps and a real modal formula to
///        a parameterised real equation system. 

#include "mcrl2/lps/io.h"
#include "mcrl2/modal_formula/parse.h"
#include "mcrl2/pres/io.h"
#include "mcrl2/pres/lps2pres.h"
#include "mcrl2/pres/pres_output_tool.h"
#include "mcrl2/utilities/input_output_tool.h"

using namespace mcrl2;
using namespace mcrl2::pres_system;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::log;
using pres_system::tools::pres_output_tool;

namespace mcrl2::pres_system::detail
{
/// \brief Prints a warning if formula contains an action that is not used in lpsspec.
inline void check_lps2pres_actions(const state_formulas::state_formula& formula, const lps::stochastic_specification& lpsspec)
{
  std::set<process::action_label> used_lps_actions = lps::find_action_labels(lpsspec.process());
  std::set<process::action_label> used_state_formula_actions = state_formulas::find_action_labels(formula);
  std::set<process::action_label> diff = utilities::detail::set_difference(used_state_formula_actions, used_lps_actions);
  if (!diff.empty())
  {
    mCRL2log(log::warning) << "Warning: the modal formula contains actions "
                           << core::detail::print_list(diff)
                           << " that are in the data specification, but does not appear in the LPS!" << std::endl;
  }
}

} // namespace mcrl2::pres_system::detail

class lps2pres_tool : public pres_output_tool<input_output_tool>
{
  using super = pres_output_tool<input_output_tool>;

protected:
  std::string formula_filename;
  bool timed = false;
  bool unoptimized = false;
  bool preprocess_modal_operators = false;
  bool check_only = false;

  std::string synopsis() const override { return "[OPTION]... --formula=FILE [INFILE [OUTFILE]]\n"; }

  void add_options(interface_description& desc) override
  {
    super::add_options(desc);
    desc.add_option("formula", make_file_argument("FILE"), "use the state formula from FILE", 'f');
    desc.add_option("preprocess-modal-operators",
        "insert dummy fixpoints in modal operators, which may lead to smaller PRESs",
        'm');
    desc.add_option("timed", "use the timed version of the algorithm, even for untimed LPS's", 't');
    desc.add_option("unoptimized", "do not simplify boolean expressions", 'u');
    desc.add_hidden_option("check-only", "check syntax and semantics of state formula; do not generate PRES", 'e');
    }

    void parse_options(const command_line_parser& parser) override
    {
      super::parse_options(parser);
      if (parser.options.count("formula"))
      {
        formula_filename = parser.option_argument("formula");
      }
      preprocess_modal_operators = parser.options.count("preprocess-modal-operators") > 0;
      timed       = parser.options.count("timed") > 0;
      unoptimized = parser.options.count("unoptimized") > 0;
      check_only = parser.options.count("check-only") > 0;
    }

  public:
    lps2pres_tool() : super(
        "lps2pres",
        "Jan Friso Groote. Based on the tool lps2pbes",
        "generate a PRES from a stochastic LPS and a modal formula over reals",
        "Convert the state formula in FILE and the LPS in INFILE to a parameterised "
        "real equation system (PRES) and save it to OUTFILE. If OUTFILE is not "
        "present, stdout is used. If INFILE is not present, stdin is used."
      )
    {}

    bool run() override
    {
      if (formula_filename.empty())
      {
        throw mcrl2::runtime_error("option -f is not specified");
      }
      if (input_filename().empty())
      {
        mCRL2log(log::verbose) << "reading LPS from stdin..." << std::endl;
      }
      else
      {
        mCRL2log(log::verbose) << "reading LPS from file '" <<  input_filename() << "'..." << std::endl;
      }
      lps::stochastic_specification lpsspec;
      load_lps(lpsspec, input_filename());
      mCRL2log(log::verbose) << "reading input from file '" <<  formula_filename << "'..." << std::endl;
      std::ifstream from(formula_filename.c_str(), std::ifstream::in | std::ifstream::binary);
      if (!from)
      {
        throw mcrl2::runtime_error("cannot open state formula file: " + formula_filename);
      }
      std::string text = utilities::read_text(from);
      const bool formula_is_quantitative = true;
      state_formulas::state_formula_specification formspec = state_formulas::algorithms::parse_state_formula_specification(text, lpsspec, formula_is_quantitative);
      pres_system::detail::check_lps2pres_actions(formspec.formula(), lpsspec);
      mCRL2log(log::verbose) << "converting state formula and LPS to a PRES..." << std::endl;
      const bool structured=false;
      pres_system::pres result = pres_system::lps2pres(lpsspec, formspec, timed, structured, unoptimized, preprocess_modal_operators, check_only);

      if (check_only)
      {
        mCRL2log(mcrl2::log::info)
          << "the file '" << formula_filename
          << "' contains a well-formed state formula" << std::endl;
        return true;
      }

      if (output_filename().empty())
      {
        mCRL2log(log::verbose) << "writing PRES to stdout..." << std::endl;
      }
      else
      {
        mCRL2log(log::verbose) << "writing PRES to file '" <<  output_filename() << "'..." << std::endl;
      }
      save_pres(result, output_filename(), m_pres_output_format);
      return true;
    }

};

int main(int argc, char** argv)
{
  return lps2pres_tool().execute(argc, argv);
}
