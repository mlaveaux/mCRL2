// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/input_input_output_tool.h
/// \brief add your file description here.

#ifndef MCRL2_UTILITIES_INPUT_INPUT_OUTPUT_TOOL_H
#define MCRL2_UTILITIES_INPUT_INPUT_OUTPUT_TOOL_H

#include "mcrl2/utilities/input_input_tool.h"





namespace mcrl2::utilities::tools
{

/// \brief Base class for tools that take two files as input,
/// and write the results to a file.
class input_input_output_tool: public input_input_tool
{
  protected:
    /// The output file name
    std::string m_output_filename;

    /// \brief Checks if the number of positional options is OK.
    /// \param parser A command line parser
    void check_positional_options(const command_line_parser& parser) override
    {
      if (3 < parser.arguments.size())
      {
        parser.error("too many file arguments");
      }
    }

    /// \brief Returns the synopsis of the tool.
    std::string synopsis() const override { return "[OPTION]... INFILE1 [INFILE2 [OUTFILE]]\n"; }

    /// \brief Parse non-standard options
    /// \param parser A command line parser
    void parse_options(const command_line_parser& parser) override
    {
      input_input_tool::parse_options(parser);
      if (2 < parser.arguments.size())
      {
        m_output_filename = parser.arguments[2];
      }
    }

    /// \brief Returns a message about the output filename
    std::string output_file_message() const
    {
      std::ostringstream out;
      out << "Output written to " << ((m_output_filename.empty())? "standard output" : ("'" + m_output_filename + "'"));
      return out.str();
    }

    /// \brief Adds a message about input and output files to the given description.
    std::string make_tool_description(const std::string& description) const
    {
      return description + " If INFILE2 is not present, standard input is used. If OUTFILE is not present, standard output is used.";
    }

  public:
    /// \brief Constructor.
    input_input_output_tool(const std::string& name,
                            const std::string& author,
                            const std::string& what_is,
                            const std::string& tool_description,
                            std::string known_issues = ""
                           )
      : input_input_tool(name, author, what_is, tool_description, known_issues)
    {
    }

    /// \brief Returns a const reference to the output filename.
    const std::string& output_filename() const
    {
      return m_output_filename;
    }

    /// \brief Returns a reference to the output filename.
    std::string& output_filename()
    {
      return m_output_filename;
    }
};

} // namespace mcrl2::utilities::tools





#endif // MCRL2_UTILITIES_INPUT_INPUT_OUTPUT_TOOL_H
