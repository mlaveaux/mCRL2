// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/parallel_tool.h
/// \brief Class that transforms a tool into a parallel tool with a --threads flag. 

#ifndef MCRL2_UTILITIES_PARALLEL_TOOL_H
#define MCRL2_UTILITIES_PARALLEL_TOOL_H

#include "mcrl2/utilities/command_line_interface.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/configuration.cxx"
#else
  import utilities;
#endif

namespace mcrl2::utilities::tools
{

/// \brief Base class for tools that are using multiple threads.
template <typename Tool>
class parallel_tool: public Tool
{
  protected:
    /// The number of threads in the parallel tool.
    std::size_t m_number_of_threads=1;

    /// \brief Add the --threads option to an interface description. 
    /// \param desc An interface description
    void add_options(utilities::interface_description& desc) override
    {
      Tool::add_options(desc);
      if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
      {
        desc.add_option("threads", utilities::make_mandatory_argument("NUM"), "run with NUM threads (default=1). With multiple threads the stack size on a Mac is limited which can lead to bus errors. ");
      }
      else
      {
        desc.add_option("threads", utilities::make_mandatory_argument("NUM"), "run with NUM threads (default=1). This option is not available as this tool is compiled as a sequential tool.");
      }
    }


    /// \brief Parse non-standard options
    /// \param parser A command line parser
    void parse_options(const utilities::command_line_parser& parser) override
    {
      Tool::parse_options(parser);

      if (parser.has_option("threads"))
      {
        m_number_of_threads = parser.option_argument_as<std::size_t>("threads");
        if (m_number_of_threads<1)
        {
          throw mcrl2::runtime_error("The number of threads should at least be 1.");
        }
        if (!mcrl2::utilities::detail::GlobalThreadSafe && number_of_threads()!=1)
        {
          throw mcrl2::runtime_error("This tool is compiled for sequential use. The number of threads (now: " + 
                                     std::to_string(number_of_threads()) +
                                     ") can only be 1.");
        }
      }
    }

  public:

    /// \brief Constructor.
    parallel_tool(const std::string& name,
                  const std::string& author,
                  const std::string& what_is,
                  const std::string& tool_description,
                  std::string known_issues = ""
                 )
      : Tool(name, author, what_is, tool_description, known_issues)
    {}

    /// \brief Returns the number of threads in this tool.
    /// \return The number of threads.
    std::size_t number_of_threads() const
    {
      return m_number_of_threads;
    }
};

} // namespace mcrl2::utilities::tools

#endif // MCRL2_UTILITIES_PARALLEL_TOOL_H
