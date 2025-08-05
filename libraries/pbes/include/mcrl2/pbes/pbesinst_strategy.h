// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesinst_strategy.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PBESINST_STRATEGY_H
#define MCRL2_PBES_PBESINST_STRATEGY_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system {

/// \brief pbesinst transformation strategies
enum pbesinst_strategy
{
  pbesinst_lazy_strategy,
  pbesinst_alternative_lazy_strategy,
  pbesinst_finite_strategy
};

/// \brief Parse a pbesinst transformation strategy.
inline
pbesinst_strategy parse_pbesinst_strategy(const std::string& s)
{
  if (s == "finite")
  {
    return pbesinst_finite_strategy;
  }
  else if (s == "alternative-lazy")
  {
    return pbesinst_alternative_lazy_strategy;
  }
  else if (s == "lazy")
  {
    return pbesinst_lazy_strategy;
  }
  else
  {
    throw mcrl2::runtime_error("unknown pbesinst strategy specified (got `" + s + "')");
  }
}

inline
std::istream& operator>>(std::istream& is, pbesinst_strategy& s)
{
  try
  {
    std::string str;
    is >> str;
    s = parse_pbesinst_strategy(str);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

/// \brief Returns a string representation of a pbesinst transformation strategy.
inline
std::string print_pbesinst_strategy(const pbesinst_strategy strategy)
{
  if (strategy == pbesinst_finite_strategy)
  {
    return "finite";
  }
  else if (strategy == pbesinst_alternative_lazy_strategy)
  {
    return "alternative-lazy";
  }
  else if (strategy == pbesinst_lazy_strategy)
  {
    return "lazy";
  }
  throw mcrl2::runtime_error("unknown pbesinst strategy");
}

inline
std::ostream& operator<<(std::ostream& os, const pbesinst_strategy strategy)
{
  os << print_pbesinst_strategy(strategy);
  return os;
}

/// \brief Returns a string representation of a pbesinst transformation strategy.
inline
std::string description(const pbesinst_strategy strategy)
{
  if (strategy == pbesinst_finite_strategy)
  {
    return "for computing all possible boolean equations";
  }
  else if (strategy == pbesinst_alternative_lazy_strategy)
  {
    return "an alternative version of the lazy strategy that supports more options";
  }
  else if (strategy == pbesinst_lazy_strategy)
  {
    return "for computing only boolean equations which can be reached from the initial state";
  }
  throw mcrl2::runtime_error("unknown pbesinst strategy");
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_PBESINST_STRATEGY_H
