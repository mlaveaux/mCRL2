// Author(s): Menno Bartels and Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbestransform.cpp

#include "mcrl2/data/rewriter.h"
#include "mcrl2/pbes/detail/stategraph_pbes.h"
#include "mcrl2/pbes/pbes_symmetry.h"
#include "mcrl2/utilities/detail/transform_tool.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/pbes/pbes_input_tool.h"
#include "mcrl2/pbes/detail/pbes_io.h"
#include "mcrl2/pbes/detail/stategraph_influence.h"

using namespace mcrl2;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::data::tools;
using namespace mcrl2::pbes_system;
using namespace mcrl2::pbes_system::tools;

class pbessymmetry_tool: public rewriter_tool<pbes_input_tool<input_tool>>
{
  using super = rewriter_tool<pbes_input_tool>;

public:
  pbessymmetry_tool()
      : super("pbessymmetry",
            "Menno Bartels and Maurice Laveaux",
            "Determines symmetries within a given PBES",
            "Detects symmetries within the PBES in INFILE and write the result to STDOUT. If INFILE is not present, stdin is used.")
  {}
  
  void parse_options(const command_line_parser& parser) override
  {
      super::parse_options(parser);
      if (parser.has_option("permutation"))
      {
        m_permutation = pbes_system::detail::permutation(parser.option_argument("permutation"));
      }
  }

  void add_options(interface_description& desc) override
  {
      desc.add_option("permutation",
        utilities::make_mandatory_argument("PERMUTATION"),
        "Checks whether a permutation is a symmetry for the PBES.",
        'y');
      super::add_options(desc);
  }

  bool run() override
  {   
    pbes input;
    mcrl2::pbes_system::load_pbes(input, input_filename(), pbes_input_format());
    
    pbes_symmetry algorithm(input);
    if (m_permutation.mapping().size() > 0)
    {
      if (algorithm.check_permutation(m_permutation))
      {
        std::cout << "true" << std::endl;
      }
      else
      {
        std::cout << "false" << std::endl;
      }
      return true;
    }
    else
    {
      algorithm.run();
    }
    
    return true;
  }

private:
  pbes_system::detail::permutation m_permutation;
};

int main(int argc, char* argv[])
{
  return pbessymmetry_tool().execute(argc, argv);
}
