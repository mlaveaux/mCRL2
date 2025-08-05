// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/detail/io.h
/// \brief add your file description here.

#ifndef MCRL2_UTILITIES_DETAIL_IO_H
#define MCRL2_UTILITIES_DETAIL_IO_H

#include "mcrl2/utilities/text_utility.h"

#include <fstream>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::utilities::detail
{

/// \brief Reads text from the file filename, or from stdin if filename equals "-".
inline
std::string read_text(const std::string& filename)
{
  if (filename.empty())
  {
    return utilities::read_text(std::cin);
  }
  else
  {
    std::ifstream from(filename.c_str());
    if (!from.good())
    {
      throw mcrl2::runtime_error("Could not read from filename " + filename);
    }
    return utilities::read_text(from);
  }
}

/// \brief Saves text to the file filename, or to stdout if filename equals "-".
inline
void write_text(const std::string& filename, const std::string& text)
{
  if (filename.empty())
  {
    std::cout << text;
  }
  else
  {
    std::ofstream to(filename);
    if (!to.good())
    {
      throw mcrl2::runtime_error("Could not write to filename " + filename);
    }
    to << text;
  }
}

} // namespace mcrl2::utilities::detail

#endif // MCRL2_UTILITIES_DETAIL_IO_H
