// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lps2pbes.cpp
/// \brief Add your file description here.

#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/pbes/tools.h"
#include "mcrl2/pbes/pbes_output_tool.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/modal_formula/parse.h"
#include "mcrl2/pbes/lps2pbes.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::log;
using pbes_system::tools::pbes_output_tool;

class lps2pbes_tool : public pbes_output_tool<input_output_tool>
{
    typedef pbes_output_tool<input_output_tool> super;

  protected:
    std::string formula_filename;
    bool timed = false;
    bool structured = false;
    bool unoptimized = false;
    bool preprocess_modal_operators = false;
    bool generate_counter_example = false;
    bool check_only = false;
    bool print_ast = false;

    std::string synopsis() const override
    {
      return "[OPTION]... --formula=FILE [INFILE [OUTFILE]]\n";
    }

    void add_options(interface_description& desc) override
    {
      super::add_options(desc);
      desc.add_option("formula", make_file_argument("FILE"),
                      "use the state formula from FILE", 'f');
      desc.add_option("preprocess-modal-operators",
                      "insert dummy fixpoints in modal operators, which may lead to smaller PBESs", 'm');
      desc.add_option("timed",
                      "use the timed version of the algorithm, even for untimed LPS's", 't');
      desc.add_option("structured",
                      "generate equations such that no mixed conjunctions and disjunctions occur", 's');
      desc.add_option("unoptimized",
                      "do not simplify boolean expressions", 'u');
      desc.add_option("counter-example",
                      "add counter example equations to the generated PBES", 'c');
      desc.add_option("check-only",
                      "check syntax and semantics of state formula; do not generate PBES", 'e');
      desc.add_hidden_option("print-ast",
                      "prints the abstract syntax tree of the state formula; do not generate PBES", 'A');
    }

    void parse_options(const command_line_parser& parser) override
    {
      super::parse_options(parser);
      if (parser.options.count("formula"))
      {
        formula_filename = parser.option_argument("formula");
      }
      preprocess_modal_operators = parser.options.count("preprocess-modal-operators") > 0;
      structured  = parser.options.count("structured") > 0;
      timed       = parser.options.count("timed") > 0;
      unoptimized = parser.options.count("unoptimized") > 0;
      generate_counter_example = parser.options.count("counter-example") > 0;
      check_only = parser.options.count("check-only") > 0;
      print_ast = parser.options.count("print-ast") > 0;
    }

  public:
    lps2pbes_tool() : super(
        "lps2pbes",
        "Wieger Wesselink; Tim Willemse",
        "generate a PBES from an LPS and a state formula",
        "Convert the state formula in FILE and the LPS in INFILE to a parameterised "
        "boolean equation system (PBES) and save it to OUTFILE. If OUTFILE is not "
        "present, stdout is used. If INFILE is not present, stdin is used."
      )
    {}

    bool run() override
    {
      if (print_ast) 
      {
        lps::specification plain_lpsspec;
        load_lps(plain_lpsspec, input_filename());  // Read as a non stochastic lps, because lps2pbes cannot handle stochastic lps's.
        lps::stochastic_specification lpsspec(plain_lpsspec);
        mCRL2log(log::verbose) << "Reading input from file '" <<  formula_filename << "'..." << std::endl;
        std::ifstream from(formula_filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (!from)
        {
          throw mcrl2::runtime_error("Cannot open state formula file: " + formula_filename);
        }
        std::string text = utilities::read_text(from);
        mcrl2::state_formulas::parse_state_formula_options options;
        options.type_check = false;
        options.translate_regular_formulas = false;
        options.translate_user_notation = false;
        options.resolve_name_clashes = false;
        options.check_monotonicity = false;
        state_formulas::state_formula_specification formspec = mcrl2::state_formulas::parse_state_formula_specification(text, lpsspec, false, options);
        mCRL2log(log::verbose) << "Converting state formula and LPS to a PBES..." << std::endl;
        std::cout << pp(formspec, false) << std::endl;
        return true;
      }

      lps2pbes(input_filename(),
               output_filename(),
               pbes_output_format(),
               formula_filename,
               timed,
               structured,
               unoptimized,
               preprocess_modal_operators,
               generate_counter_example,
               check_only
             );
      return true;
    }

};

int main(int argc, char** argv)
{
  return lps2pbes_tool().execute(argc, argv);
}
