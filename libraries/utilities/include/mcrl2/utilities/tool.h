// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/tool.h
/// \brief add your file description here.

#ifndef MCRL2_UTILITIES_TOOL_H
#define MCRL2_UTILITIES_TOOL_H

#include "mcrl2/utilities/logger.h"

#include "mcrl2/utilities/command_line_interface.h"
#include "mcrl2/utilities/platform.h"

#ifdef MCRL2_PLATFORM_WINDOWS
  #include <io.h>
  #include <fcntl.h>
#endif

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/execution_timer.cxx"
#else  
  import utilities;
#endif

/// \brief The namespace for command line tool classes
namespace mcrl2::utilities::tools
{

/// \brief Base class for command line tools.
/// as result.
class tool
{
  protected:
    /// The name of the tool
    std::string m_name;

    /// The name of the developer(s)
    std::string m_author;

    /// One-line "what is" description of the tool
    std::string m_what_is;

    /// The description of the tool
    std::string m_tool_description;

    /// Known issues of the tool
    std::string m_known_issues;

    /// The filename to which timings must be written
    std::string m_timing_filename;

    /// The timer which can be used by the tools
    execution_timer m_timer;

    /// Determines whether timing output should be written
    bool m_timing_enabled = false;

    /// \brief Add options to an interface description.
    /// \param desc An interface description
    virtual void add_options(interface_description& desc)
    {
      desc.add_option("timings", make_optional_argument<std::string>("FILE", ""),
                      "append timing measurements to FILE. Measurements are written to "
                      "standard error if no FILE is provided");
    }

    /// \brief Parse non-standard options
    /// \param parser A command line parser
    virtual void parse_options(const command_line_parser& parser)
    {
      if (parser.options.count("timings") > 0)
      {
        m_timing_enabled = true;
        log::logger::set_report_time_info();
        m_timing_filename = parser.option_argument("timings");
      }
    }

    /// \brief Executed only if run would be executed and invoked before run.
    /// \return Whether run should still be executed
    virtual bool pre_run(int& /*argc*/, char** /*argv*/)
    {
      return true;
    }

    /// \brief Parse standard options
    /// \param parser A command line parser
    virtual void check_standard_options(const command_line_parser& parser)
    {
      if (parser.options.count("quiet"))
      {
        if (parser.options.count("debug"))
        {
          parser.error("options -q/--quiet and -d/--debug cannot be used together\n");
        }
        if (parser.options.count("verbose"))
        {
          parser.error("options -q/--quiet and -v/--verbose cannot be used together\n");
        }
      }
#ifndef MCRL2_TOOL_CLASSES_NO_CORE
      if (parser.options.count("quiet"))
      {
        log::logger::set_reporting_level(log::quiet);
      }
      if (parser.options.count("verbose"))
      {
        log::logger::set_reporting_level(log::verbose);
      }
      if (parser.options.count("debug"))
      {
        log::logger::set_reporting_level(log::debug);
      }
      if (parser.options.count("log-level"))
      {
        log::logger::set_reporting_level(log::log_level_from_string(parser.option_argument("log-level")));
      }
#endif
    }

    /// \brief Checks if the number of positional options is OK.
    /// By default this function handles standard options: -v, -d and -q
    /// Furthermore, it checks that all options occur at most once
    /// \param parser A command line parser
    virtual void check_positional_options(const command_line_parser& parser)
    {
      parser.check_no_duplicate_arguments();
    }

    /// \brief Returns the synopsis of the tool
    /// \return The string "[OPTION]...\n"
    virtual std::string synopsis() const
    {
      return "[OPTION]...\n";
    }

  public:
    /// \brief Constructor.
    tool(const std::string& name,
         const std::string& author,
         const std::string& what_is,
         const std::string& tool_description,
         std::string known_issues = ""
        )
      : m_name(name),
        m_author(author),
        m_what_is(what_is),
        m_tool_description(tool_description),
        m_known_issues(known_issues),
        m_timing_filename(""),
        m_timer(name)
    {}

    /// \brief Destructor.
    virtual ~tool() = default;

    /// \brief Run the tool. The options must be set manually.
    /// \return True if the tool execution was successful.
    virtual bool run() = 0;

    /// \brief Return the filename in which timings must be saved.
    const std::string& timing_filename() const
    {
      return m_timing_filename;
    }

    /// \brief Return reference to the timer that can be used.
    execution_timer& timer()
    {
      return m_timer;
    }

    /// \brief Run the tool with the given command line options.
    /// \param argc Number of command line arguments
    /// \param argv Command line arguments
    /// \return The execution result
    /// \post If timing was enabled, timer().report() has been called
    int execute(int argc, char* argv[])
    {
#ifdef MCRL2_PLATFORM_WINDOWS
      // All tools expect their std::cin to be binary streams. In Windows, 
      // this is the way to guarantee this. In applications that do not have
      // a standard input attached (Windows-mode applications for instance),
      // _fileno(stdin) returns a value less than 0.
      if (_fileno(stdin) >= 0 && _setmode(_fileno(stdin), _O_BINARY) == -1)
      {
        throw std::runtime_error("Cannot set stdin to binary mode");
      }
#endif
      try
      {
        interface_description clinterface(argv[0], m_name, m_author, m_what_is, synopsis(), m_tool_description, m_known_issues);
        add_options(clinterface);
        command_line_parser parser(clinterface, argc, argv);
        check_standard_options(parser);

        if (parser.continue_execution())
        {
          check_positional_options(parser);
          parse_options(parser);

          // If pre_run succeeds, then do the actual running.
          bool result = pre_run(argc, argv);
          if (result)
          {
            // Create timer, and by default measure running time of run()
            // method.
            m_timer = execution_timer(m_name, timing_filename());

            timer().start("total");
            result = run();
            timer().finish("total");

            if (m_timing_enabled)
            {
              timer().report();
            }
          }

          // Either pre_run or run failed.
          if (!result)
          {
            return EXIT_FAILURE;
          }
        }

        return EXIT_SUCCESS;
      }
      catch (mcrl2::runtime_error& e)
      {
        mCRL2log(mcrl2::log::error) << e.what() << std::endl;
      }
      catch (std::exception& e)
      {
        mCRL2log(mcrl2::log::error) << e.what() << std::endl;
      }
      return EXIT_FAILURE;
    }
};

} // namespace mcrl2::utilities::tools

#endif // MCRL2_UTILITIES_TOOL_H
