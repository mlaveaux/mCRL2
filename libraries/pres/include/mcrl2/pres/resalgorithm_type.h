// Author(s): Jan Friso Groote. 
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/resalgorithm_type.h
/// \brief Contains the various types of algorithms to solve res's.

#ifndef MCRL2_PRES_RESALGORITHM_TYPE_H
#define MCRL2_PRES_RESALGORITHM_TYPE_H

#include <iomanip>
#include "mcrl2/core/detail/print_utility.h"



namespace mcrl2::pres_system {

enum solution_algorithm { gauss_elimination, numerical, numerical_directed };

inline
std::string print_algorithm(const solution_algorithm alg)
{
  switch(alg)
  {
    case gauss_elimination: return "gauss";
    case numerical: return "numerical";
    case numerical_directed: return "numerical_directed";
    default: throw mcrl2::runtime_error("unknown res algorithm");
  }
}


inline
std::string description(const solution_algorithm a)
{
  switch(a)
  {
    case gauss_elimination: return "solve the res using gauss elimination; this is guaranteed to terminate but may require an excessive amount of time.";
    case numerical: return "solve the res by a numerical recursive algorithm; this is not guaranteed to terminate.";
    case numerical_directed: return "solve the res by a numerical recursive algorithm with directed propagation; this is not guaranteed to terminate.";
    default: throw mcrl2::runtime_error("unknown algorithm");
  }
}

/// \brief Parse an algorithm 
/// \param[in] s A string
/// \return The algorithm represented by \a s
inline
solution_algorithm parse_algorithm(const std::string& s)
{
  if(s == "g" || s == "gauss")
  {
    return gauss_elimination;
  }
  else if (s == "n" || s == "numerical")
  {
    return numerical;
  }
  else if (s == "m" || s == "numerical_directed")
  {
    return numerical_directed;
  }
  else
  {
    throw mcrl2::runtime_error("Unknown algorithm " + s);
  }
}


// \overload
inline
std::istream& operator>>(std::istream& is, solution_algorithm& l)
{
  try {
    std::stringbuf buffer;
    is >> &buffer;
    l = parse_algorithm(buffer.str());
  }
  catch (mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }

  return is;
}

// \overload
inline
std::ostream& operator<<(std::ostream& os, const solution_algorithm l)
{
  os << print_algorithm(l);
  return os;
}

} // namespace mcrl2::pres_system



#endif // MCRL2_PRES_RESALGORITHM_TYPE_H
