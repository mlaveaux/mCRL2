// Author(s): Jeroen Keiren.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file solution_strategy.h
/// \brief Strategies to solve a BES.

#ifndef MCRL2_PBES_SOLUTION_STRATEGY_H
#define MCRL2_PBES_SOLUTION_STRATEGY_H

#include <fstream>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system
{

using solution_strategy_t = enum
{
  gauss,
  small_progr_measures
};

static
std::string solution_strategy_to_string(const solution_strategy_t s)
{
  switch (s)
  {
    case gauss:
      return "gauss";
      break;
    case small_progr_measures:
      return "spm";
      break;
  }
  throw mcrl2::runtime_error("unknown solution strategy");
}

static
std::ostream& operator<<(std::ostream& os, const solution_strategy_t s)
{
  os << solution_strategy_to_string(s);
  return os;
}

static
solution_strategy_t parse_solution_strategy(const std::string& s)
{
  if (s == "gauss")
  {
    return gauss;
  }
  else if (s == "spm")
  {
    return small_progr_measures;
  }
  else
  {
    throw mcrl2::runtime_error("unsupported solution strategy '" + s + "'");
  }
}

static
std::istream& operator>>(std::istream& is, solution_strategy_t& s)
{
  try
  {
    std::string str;
    is >> str;
    s = parse_solution_strategy(str);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

static
std::string description(const solution_strategy_t s)
{
  switch (s)
  {
    case gauss:
      return "Gauss elimination (inefficient; plain implementation)";
      break;
    case small_progr_measures:
      return "Small progress measures";
      break;
  }
  throw mcrl2::runtime_error("unknown solution strategy");
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_SOLUTION_STRATEGY_H
