// Author(s): Jan Friso Groote. Based on lts2pbes by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts2pres.cpp
/// \brief This tool reads an LTS and a quantitative formula,  and translates them into a PRES.

#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/lts/detail/lts_load.h"
#include "mcrl2/modal_formula/parse.h"
#include "mcrl2/pres/lts2pres.h"
#include "mcrl2/pres/pres_output_tool.h"

using namespace mcrl2;
using namespace mcrl2::utilities;
using pres_system::tools::pres_output_tool;
using utilities::tools::input_output_tool;

inline
void check_lts(const lts::probabilistic_lts_lts_t& ltsspec)
{
  for (const lts::transition& tr: ltsspec.get_transitions())
  {
    std::size_t label = tr.label();
    if (label >= ltsspec.action_labels().size())
    {
      throw mcrl2::runtime_error("Invalid LTS detected: there are not enough action labels");
    }
  }
}

inline
std::set<process::action_label> find_action_labels(const lts::probabilistic_lts_lts_t& ltsspec)
{
  std::set<std::size_t> used_labels;
  for (const lts::transition& tr: ltsspec.get_transitions())
  {
    used_labels.insert(tr.label());
  }
  std::set<process::action_label> result;
  const std::vector<lts::action_label_lts>& action_labels = ltsspec.action_labels();
  for (std::size_t index: used_labels)
  {
    for (const process::action& a: action_labels[index].actions())
    {
      result.insert(a.label());
    }
  }
  return result;
}

/// \brief Prints a warning if formula contains an action that is not used in ltsspec.
inline void check_lts2pres_actions(const state_formulas::state_formula& formula, const lts::probabilistic_lts_lts_t& ltsspec)
{
  std::set<process::action_label> used_lts_actions = find_action_labels(ltsspec);
  std::set<process::action_label> used_state_formula_actions = state_formulas::find_action_labels(formula);
  std::set<process::action_label> diff = utilities::detail::set_difference(used_state_formula_actions, used_lts_actions);
  if (!diff.empty())
  {
    mCRL2log(log::warning) << "Warning: the modal formula contains actions "
                           << core::detail::print_list(diff)
                           << " that are in the data specification, but does not appear in the LTS!" << std::endl;
  }
}

class lts2pres_tool : public pres_output_tool<input_output_tool>
{
  private:
    using super = pres_output_tool<input_output_tool>;

  protected:

    std::string formfilename;
    bool preprocess_modal_operators = false;
    lts::probabilistic_lts_lts_t ltsspec;

    void add_options(interface_description& desc) override
    {
      super::add_options(desc);

      desc.add_option("formula", make_file_argument("FILE"),
                      "use the state formula from FILE", 'f');

      desc.add_option("preprocess-modal-operators",
                      "insert dummy fixpoints in modal operators, which may lead to smaller PRESs", 'p');

      lts::detail::add_options(desc);
    }

    void parse_options(const command_line_parser& parser) override
    {
      super::parse_options(parser);

      if (parser.options.count("formula"))
      {
        formfilename = parser.option_argument("formula");
      }

      preprocess_modal_operators = parser.options.count("preprocess-modal-operators") > 0;
      lts::detail::load_lts(parser, input_filename(), ltsspec);
      check_lts(ltsspec);
    }

  public:

    lts2pres_tool() : super(
        "lts2pres",
        "Jan Friso Groote",
        "translates an LTS into an PRES",
        "Translates an LTS in INFILE and writes the resulting PRES to "
        "OUTFILE. If OUTFILE is not present, standard output is used. If INFILE is not "
        "present, standard input is used.")
    {}

  public:
    bool run() override
    {
      lps::stochastic_specification lpsspec = lts::detail::extract_specification(ltsspec);
      std::ifstream from(formfilename.c_str());
      if (!from)
      {
        throw mcrl2::runtime_error("Cannot open state formula file: " + formfilename + ".");
      }
      const bool formula_is_quantitative = true;
      state_formulas::state_formula_specification formspec = state_formulas::parse_state_formula_specification(from, lpsspec, formula_is_quantitative);
      check_lts2pres_actions(formspec.formula(), ltsspec);
      from.close();
      ltsspec.set_data(data::merge_data_specifications(ltsspec.data(), formspec.data()));
      if (!formspec.action_labels().empty())
      {
        mCRL2log(log::warning) << "The modal formula contains action declarations. These are ignored.\n";
      }
      pres_system::pres result = pres_system::lts2pres(ltsspec, formspec, preprocess_modal_operators);

      //save the result
      if (output_filename().empty())
      {
        mCRL2log(log::verbose) << "Writing PRES to stdout..." << std::endl;
      }
      else
      {
        mCRL2log(log::verbose) << "Writing PRES to file '" <<  output_filename() << "'..." << std::endl;
      }
      save_pres(result, output_filename(), pres_output_format());
      return true;
    }
};

int main(int argc, char** argv)
{
  return lts2pres_tool().execute(argc, argv);
}
