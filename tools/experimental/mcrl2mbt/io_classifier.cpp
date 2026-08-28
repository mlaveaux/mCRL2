// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "io_classifier.h"

#include <fstream>
#include <stdexcept>

#include "mcrl2/utilities/exception.h"

namespace mcrl2
{

void io_classifier::load(const std::string& filename)
{
  std::ifstream file(filename);
  if (!file)
  {
    throw mcrl2::runtime_error("cannot open action classification file: " + filename);
  }

  read(file);
}

void io_classifier::read(std::istream& input)
{
  enum class section
  {
    none,
    input,
    output
  };
  section current = section::none;
  std::string line;

  while (std::getline(input, line))
  {
    // Strip comments and trailing whitespace
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos)
    {
      line.erase(hash);
    }

    // Trim leading/trailing whitespace
    const std::size_t first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos)
    {
      continue; // blank line
    }
    const std::size_t last = line.find_last_not_of(" \t\r");
    line = line.substr(first, last - first + 1);

    if (line == "input")
    {
      current = section::input;
    }
    else if (line == "output")
    {
      current = section::output;
    }
    else
    {
      switch (current)
      {
      case section::input:
        m_inputs.insert(line);
        break;
      case section::output:
        m_outputs.insert(line);
        break;
      case section::none:
        throw mcrl2::runtime_error("action '" + line + "' appears before any section header");
      }
    }
  }
}

io_classifier::action_type io_classifier::classify(const std::string& action_name) const
{
  if (m_inputs.contains(action_name))
  {
    return action_type::input;
  }

  if (m_outputs.contains(action_name))
  {
    return action_type::output;
  }
  
  return action_type::internal;
}

} // namespace mcrl2
