// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/execution_timer.h
/// \brief Class to obtain running times of code.

#ifndef MCRL2_UTILITIES_EXCEPTION_H
#define MCRL2_UTILITIES_EXCEPTION_H

module;

#include <cassert>
#include <sstream>
#include <stdexcept>

#ifdef MCRL2_ENABLE_MODULES
  export module utilities;
#endif

namespace mcrl2
{

/**
 * \brief Standard exception class for reporting runtime errors.
 **/
MCRL2_EXPORT class runtime_error : public std::runtime_error
{
public:
  /// \brief Constructor
  /// \param[in] message the exception message
  runtime_error(const std::string& message) : std::runtime_error(message)
  {}
};

/**
 * \brief Exception class for errors raised by the command-line parser.
 **/
MCRL2_EXPORT class command_line_error : public runtime_error
{
private:
  std::string m_msg;

  /// \returns A string that contains "<name>: <message>" followed by a message pointing to --help on the next line.
  static std::string format(const std::string& name, const std::string& message)
  {
    std::stringstream s;
    s << name << ": " << message << "\n"
      << "Try '" << name << " --help' for more information.";
    return s.str();
  }

public:
  command_line_error(const std::string& name, const std::string& message) noexcept
      : runtime_error(format(name, message))
  {}

  ~command_line_error() noexcept override = default;
};

} // namespace mcrl2

#endif
