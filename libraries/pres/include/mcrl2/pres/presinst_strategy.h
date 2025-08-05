// Author(s): Jan Friso Groote. Based on mcrl2/pres/presinst_strategy.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/presinst_strategy.h
/// \brief add your file description here.

#ifndef MCRL2_PRES_PRESINST_STRATEGY_H
#define MCRL2_PRES_PRESINST_STRATEGY_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pres_system {

/// \brief presinst transformation strategies
enum presinst_strategy
{
  presinst_lazy_strategy,
  presinst_finite_strategy
};

/// \brief Parse a presinst transformation strategy.
inline
presinst_strategy parse_presinst_strategy(const std::string& s)
{
  if (s == "finite")
  {
    return presinst_finite_strategy;
  }
  else if (s == "lazy")
  {
    return presinst_lazy_strategy;
  }
  else
  {
    throw mcrl2::runtime_error("unknown presinst strategy specified (got `" + s + "')");
  }
}

inline
std::istream& operator>>(std::istream& is, presinst_strategy& s)
{
  try
  {
    std::string str;
    is >> str;
    s = parse_presinst_strategy(str);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

/// \brief Returns a string representation of a presinst transformation strategy.
inline
std::string print_presinst_strategy(const presinst_strategy strategy)
{
  if (strategy == presinst_finite_strategy)
  {
    return "finite";
  }
  else if (strategy == presinst_lazy_strategy)
  {
    return "lazy";
  }
  throw mcrl2::runtime_error("unknown presinst strategy");
}

inline
std::ostream& operator<<(std::ostream& os, const presinst_strategy strategy)
{
  os << print_presinst_strategy(strategy);
  return os;
}

/// \brief Returns a string representation of a presinst transformation strategy.
inline
std::string description(const presinst_strategy strategy)
{
  if (strategy == presinst_finite_strategy)
  {
    return "for computing all possible boolean equations";
  }
  else if (strategy == presinst_lazy_strategy)
  {
    return "for computing only boolean equations which can be reached from the initial state";
  }
  throw mcrl2::runtime_error("unknown presinst strategy");
}

} // namespace mcrl2::pres_system



#endif // MCRL2_PRES_PRESINST_STRATEGY_H
