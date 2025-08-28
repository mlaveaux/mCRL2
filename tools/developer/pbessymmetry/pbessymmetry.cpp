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
  }

  void add_options(interface_description& desc) override
  {
      super::add_options(desc);
  }

  bool run() override
  {   
    // TODO: where does the input format go?
    pbes input;
    mcrl2::pbes_system::load_pbes(input, input_filename(), pbes_input_format());
    
    mcrl2::data::rewriter rewr = create_rewriter();
    pbes_symmetry algorithm(input, rewr);

    return true;
  }
};

int main(int argc, char* argv[])
{
  return pbessymmetry_tool().execute(argc, argv);
}
