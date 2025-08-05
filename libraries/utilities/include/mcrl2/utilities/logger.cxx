// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file logger.h

#ifndef MCRL2_UTILITIES_LOGGER_H
#define MCRL2_UTILITIES_LOGGER_H
MCRL2_MODULE;

#include <array>
#include <cassert>
#include <atomic>
#include <cstdio>
#include <ctime>
#include <set>
#include <stdexcept>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/noncopyable.h"
  #include "mcrl2/utilities/text_utility.h"
#else
  export module utilities:logger;

  import :noncopyable;
  import :text_utility;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::log {

/// \brief Log levels that are supported
/// \note log_debugi with i>=1 automatically indent 2*i spaces.
enum log_level_t
{
  quiet, // No log message should ever be printed to this log level!
  error,
  warning,
  info,
  status, // This is overwritten in the output when printed consecutively, and the highest level that is shown by default.
  verbose,
  debug,
  trace,
};

/// \brief Convert log level to string
/// This string is used to prefix messages in the logging output.
inline
std::string_view log_level_to_string(const log_level_t level)
{
  static std::array<std::string_view,87> buffer = { "quiet", "error", "warning", "info", "status", "verbose", "debug", "trace" };
  if ((unsigned) level >= buffer.size()) {
    return "unknown level";
  }
  return buffer[level];
}

/// \brief Convert string to log level
inline
log_level_t log_level_from_string(const std::string_view s)
{
  if (s == "quiet")
  {
    return quiet;
  }
  else if (s == "error")
  {
    return error;
  }
  else if (s == "warning")
  {
    return warning;
  }
  else if (s == "info")
  {
    return info;
  }
  else if (s == "status")
  {
    return status;
  }
  else if (s == "verbose")
  {
    return verbose;
  }
  else if (s == "debug")
  {
    return debug;
  }
  else if (s == "trace")
  {
    return trace;
  }
  else
  {
    throw std::runtime_error("Unknown log-level " + std::string(s) + " provided.");
  }
}

std::string format_time(const time_t* t);

/// \brief Interface class for output policy.
///
/// Separates the exact way of doing output from the logger class.
class output_policy
{
  public:
    /// \brief Constructor.
    output_policy() = default;

    /// \brief Destructor.
    virtual ~output_policy() = default;

    /// \brief Output message.
    /// \param[in] msg Message that is written to output.
    /// \param[in] timestamp Timestamp to use in the output.
    /// \param[in] level The log level to print the message to.
    /// \param[in] print_time_information An indication whether the time must be printed.
    ///  \details Any implementation must assure that output is written using an atomic action, to prevent
    /// mixing of different lines into one line in the output.
    virtual void output(const log_level_t level, const time_t timestamp, const std::string& msg, const bool print_time_information) = 0;
};

std::set<output_policy*> initialise_output_policies();

/// \brief Class for logging messages
///
/// Based on a description in the article "Logging In C++", Petru Marginean
/// Dr. Dobb's Journal, September 5, 2007
/// (url: http://drdobbs.com/cpp/201804215)
/// Requires that OutputPolicy is a class which as a static member output(const std::string&)
class logger: private utilities::noncopyable
{
  protected:
    /// \brief Stream that is printed to internally
    /// Collects the full debug message that we are currently printing.
    std::ostringstream m_os;

    /// \brief The loglevel of the current message
    log_level_t m_level;

    /// \brief Timestamp of the current message
    time_t m_timestamp = 0L;

    static std::atomic<log_level_t>& log_level()
    {
      static std::atomic<log_level_t> g_log_level(log_level_t::info);
      return g_log_level;
    }

    /// \brief An indication whether time information should be printed.
    static bool& m_print_time_information()
    {
      thread_local bool print_timing_info=false;
      return print_timing_info;
    }

    /// \brief Output policies
    static
    std::set<output_policy*>& output_policies()
    {
      static std::set<output_policy*> m_output_policies = initialise_output_policies();
      return m_output_policies;
    }

  public:
    /// \brief Default constructor
    logger(const log_level_t l)
    {
      assert(quiet != l);
      m_level = l;
      std::time(&m_timestamp);
    }

    /// \brief Destructor; flushes output.
    /// Flushing during destruction is important to confer thread safety to the
    /// logging mechanism. Requires that output performs output in an atomic way.
    ~logger()
    {
      for(output_policy* policy: output_policies())
      {
        policy->output(m_level, m_timestamp, m_os.str(), m_print_time_information());
      }
    }

    /// \brief Register output policy
    static
    void register_output_policy(output_policy& policy)
    {
      output_policies().insert(&policy);
    }

    /// \brief Unregister output policy
    static
    void unregister_output_policy(output_policy& policy)
    {
      std::set<output_policy*>::iterator i = output_policies().find(&policy);
      if(i != output_policies().end())
      {
        output_policies().erase(i);
      }
    }

    /// \brief Clear all output policies
    static
    void clear_output_policies()
    {
      output_policies().clear();
    }

    /// \brief Set reporting level
    /// \param[in] level Log level
    static
    void set_reporting_level(const log_level_t level)
    {
      log_level() = level;
    }

    /// \brief Get reporting level
    static
    log_level_t get_reporting_level()
    {
      return log_level();
    }

    /// \brief Indicate that timing information should be printed.
    static void set_report_time_info()
    {
      m_print_time_information() = true;
    }

    /// \brief Indicate that timing information should not be printed.
    static void clear_report_time_info()
    {
      m_print_time_information() = false;
    }

    /// \brief Get whether timing information is printed.
    /// \return True if time information is printed, otherwise false.
    static bool get_report_time_info() 
    {
      return m_print_time_information();
    }

    /// Get access to the stream provided by the logger.
    /// \param[in] l Log level for the stream
    std::ostringstream& get()
    {
      return m_os;
    }
};

class formatter_interface
{
public:
  /// \brief Format msg,
  /// \param[in] level The log level of the message
  /// \param[in] timestamp The timestamp of the log message
  /// \param[in] msg The message to be formatted
  /// \param[in] print_time_information A boolean that if true indicates that time usage information must be printed. If false this 
  ///                                   information is suppressed. 
  /// \return The formatted message (\a msg)
  static std::string format(const log_level_t level, const time_t timestamp, const std::string& msg, const bool print_time_information)
  {
    /// suppress non used variable warnings.
    (void)level; (void)timestamp; (void)print_time_information;

    assert(quiet != level);
    return msg;
  }
}; 

/// \brief Mixin that takes care of formatting of a message.
///
/// In this case, the formatter
class formatter: public formatter_interface
{
protected:
  /// \brief Records whether the last message that was printed ended with
  ///        a new line.
  static
  std::atomic<bool>& last_message_ended_with_newline()
  {
    static std::atomic<bool> m_last_message_ended_with_newline = true;
    return m_last_message_ended_with_newline;
  }

  static
  std::atomic<bool>& last_message_was_status()
  {
    static std::atomic<bool> m_last_message_was_status = false;
    return m_last_message_was_status;
  }

  static
  std::atomic<std::size_t>& caret_pos()
  {
    static std::atomic<std::size_t> m_caret_pos = 0;
    return m_caret_pos;
  }

  static
  std::atomic<std::size_t>& last_caret_pos()
  {
    static std::atomic<std::size_t> m_last_caret_pos = 0;
    return m_last_caret_pos;
  }

public:
  /// \brief Prefix each line in s with some extra information.
  /// The things that are added are:
  /// - current time
  /// - log level
  /// - indentation
  static std::string format(const log_level_t level, const time_t timestamp, const std::string& msg, const bool print_time_information);
};

/// \brief File output class.
///
/// Provides facilities to output to a file. By default output is sent to stderr.
class file_output: public output_policy
{
  protected:
    /// \brief Obtain the underlying stream used to print to a file.
    static std::atomic<FILE*>& get_stream()
    {
      static std::atomic<FILE*> g_stream(stderr);
      return g_stream;
    }

  public:
    file_output() = default;

    ~file_output() override = default;

    /// \param[in] stream A file handle
    static
    void set_stream(FILE* stream)
    {
      get_stream() = stream;
    }

    /// Output message to stream.
    /// \param[in] level The log level on which to output the message
    /// \param[in] timestamp The timestamp to use for the message
    /// \param[in] msg The message to be printed
    /// \param[in] print_time_information A boolean that if true indicates that time usage information must be printed. If false this 
    ///                                   information is suppressed. 
    ///
    /// \note This uses fprintf (and not e.g. <<) because fprintf is guaranteed to be
    /// atomic.
    void output(const log_level_t level,
        const time_t timestamp,
        const std::string& msg,
        const bool print_time_information) override
    {
      assert(quiet != level);
      FILE* p_stream = get_stream();
      if (!p_stream)
      {
        return;
      }

      // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg) Maurice: waiting for std::format to be supported.
      fprintf(p_stream, "%s", formatter::format(level, timestamp, msg, print_time_information).c_str());
      fflush(p_stream);
    }
};

/// \brief The default output policy used by the logger
inline
output_policy& default_output_policy()
{
  static file_output m_default = file_output();
  return m_default;
}

/// \brief Initialise the output policies. This returns the singleton set
///        containing the default output policy.
inline
std::set<output_policy*> initialise_output_policies()
{
  std::set<output_policy*> result;
  result.insert(&default_output_policy());
  return result;
}

/// \returns True whenever the logging for the given level is enabled.
inline bool mCRL2logEnabled(const log_level_t level)
{
  return level <= mcrl2::log::logger::get_reporting_level();
}

} // namespace mcrl2::log

#endif // MCRL2_UTILITIES_LOGGER_H
