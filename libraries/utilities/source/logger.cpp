// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file logger.cpp
MCRL2_MODULE;

#include <cassert>
#include <string>
#include <iostream>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/logger.cxx"
#else
  module utilities:logger_impl;

  import :logger;
#endif


namespace mcrl2::log {

  std::string format_time(const time_t* t)
  {
    struct tm* r;
    r = localtime(t);

    char buffer[11];
    std::size_t res = strftime(buffer, sizeof(buffer), "%H:%M:%S", r);
    if(res == 0)
    {
      std::clog << "Could not write time to buffer" << std::endl;
      return std::string();
    }
    return buffer;
  }

  std::string formatter::format(const log_level_t level, 
                                const time_t timestamp, 
                                const std::string& msg,
                                const bool print_time_information)
  {
    assert(quiet != level);
    // Construct the message header, with or without time and debug level info.
    const bool print_log_level= level >= debug;
    std::stringstream start_of_line;
    if (print_time_information || print_log_level)
    {
      start_of_line << "[";
    }
    if (print_time_information)
    {
      start_of_line << format_time(&timestamp);
    }
    if (print_time_information && print_log_level)
    {
      start_of_line << " ";
    }
    if (print_log_level)
    {
      start_of_line << log_level_to_string(level);
    }
    if (print_time_information || print_log_level)
    {
      start_of_line << "] ";
    }
    if (print_log_level)
    {
      start_of_line << std::string(7 - log_level_to_string(level).size(), ' ');
    }

    // Check if message is finished. We
    bool msg_ends_with_newline = false;
    if (msg.size() > 0)
    {
      msg_ends_with_newline = (msg[msg.size()-1] == '\n');
    }

    std::stringstream result;

    // Determine whether the last line should be overwritten (by inserting a
    // carriage return).
    bool overwrite = false;
    if (last_message_was_status())
    {
      if (level == status)
      {
        if (last_message_ended_with_newline())
        {
          if (msg != "\n")
          {
            result << "\r" << start_of_line.str();
            overwrite = true;
          }
          else
          {
            result << "\n";
          }
        }
      }
      else
      {
        result << "\n" << start_of_line.str();
      }
    }
    else
    if (last_message_ended_with_newline())
    {
      result << start_of_line.str();
    }

    // Copy the message to result, keeping track of the caret position and
    // inserting message headers at newlines. A trailing newline is ignored.
    for (std::string::const_iterator it = msg.begin(); it != msg.end(); )
    {
      if (*it != '\n')
      {
        result << *it++;
        ++caret_pos();
      }
      else
      {
        if (++it != msg.end())
        {
          result << "\n" << start_of_line.str();
          caret_pos() = 0;
          last_caret_pos() = 0;
        }
      }
    }

    // Pad message with spaces when overwriting a previous line
    if (msg_ends_with_newline && overwrite && caret_pos() < last_caret_pos())
    {
      for (std::size_t i = 0; i < last_caret_pos() - caret_pos(); ++i)
      {
        result << ' ';
      }
    }

    // If this is a finished status message, remember the caret position
    // so the next line can overwrite properly. If this is a finished message
    // of another type, restore the trailing newline.
    if(msg_ends_with_newline)
    {
      if (level == status)
      {
        last_caret_pos() = caret_pos().load();
      }
      else
      {
        result << "\n";
      }
      caret_pos() = 0;
    }

    // Store information about this message so next time we can determine
    // whether we should overwrite the last line or not.
    last_message_ended_with_newline() = msg_ends_with_newline;
    if (level == status)
    {
      last_message_was_status() = msg != "\n";
    }
    else
    {
      last_message_was_status() = false;
    }

    return result.str();
  }

} // namespace mcrl2::log

