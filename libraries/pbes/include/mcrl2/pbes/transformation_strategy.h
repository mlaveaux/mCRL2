// Author(s): Jan Friso Groote, XIAO Qi
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file transformation_strategy.h
/// \brief Strategies for the generation of a BES from a PBES.

#ifndef MCRL2_PBES_TRANSFORMATION_STRATEGY_H
#define MCRL2_PBES_TRANSFORMATION_STRATEGY_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system
{

/// \brief Strategies for the generation of a BES from a PBES
enum transformation_strategy
{
  lazy,          // generate equations but do not optimize on the fly
  optimize,      // optimize by substituting true and false for already
  // investigated variables in a rhs, while generating this rhs.
  on_the_fly,    // make a distinction between variables that occur somewhere,
  // and variables that do not occur somewhere. When generating
  // a rhs, optimize this rhs as in "optimize". If the rhs is
  // equal to T or F, substitute this value throughout the
  // equation system, and maintain which variables become unused
  // by doing so, as these do not have to be investigated further.
  // E.g. if a rhs is  X1 && X2, X1 does not occur elsewhere and
  // X2 turns out to be equal to false, then X1 is moved to the
  // set of irrelevant variables, and not investigated further.
  on_the_fly_with_fixed_points
  // Do the same as with on the fly, but for each generated variable
  // in the rhs, investigate whether this variable lies on a loop
  // such that depending on its fixed point, it can be set to true
  // or false. Due to the breadth first nature of the main algorithm
  // the existence of such loops must be investigated separately
  // for each variable, which can take a lot of time.
};

inline
transformation_strategy parse_transformation_strategy(const std::string& s)
{
  if (s == "0")
  {
    return lazy;
  }
  else if (s == "1")
  {
    return optimize;
  }
  else if (s == "2")
  {
    return on_the_fly;
  }
  else if (s == "3")
  {
    return on_the_fly_with_fixed_points;
  }
  else
  {
    throw mcrl2::runtime_error("unknown transformation strategy " + s);
  }
}

inline
std::string print_transformation_strategy(const transformation_strategy s)
{
  switch(s)
  {
    case lazy: return "0";
    case optimize: return "1";
    case on_the_fly: return "2";
    case on_the_fly_with_fixed_points: return "3";
  }
  throw mcrl2::runtime_error("unknown transformation strategy");
}

inline
std::istream& operator>>(std::istream& is, transformation_strategy& strategy)
{
  try
  {
    std::string s;
    is >> s;
    strategy = parse_transformation_strategy(s);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

inline
std::ostream& operator<<(std::ostream& os, const transformation_strategy s)
{
  os << print_transformation_strategy(s);
  return os;
}

inline
std::string description(const transformation_strategy s)
{
  switch(s)
  {
    case lazy: return "Compute all boolean equations which can be reached"
        " from the initial state, without optimization."
        " This is is the most data efficient"
        " option per generated equation.";
    case optimize: return "Optimize by immediately substituting the right"
        " hand sides for already investigated variables"
        " that are true or false when generating an"
        " expression. This is as memory efficient as 0.";
    case on_the_fly: return "In addition to 1, also substitute variables that"
        " are true or false into an already generated right"
        " hand side. This can mean that certain variables"
        " become unreachable (e.g. X0 in X0 and X1, when X1"
        " becomes false, assuming X0 does not occur"
        " elsewhere. It will be maintained which variables"
        " have become unreachable as these do not have to be"
        " investigated. Depending on the PBES, this can"
        " reduce the size of the generated BES substantially"
        " but requires a larger memory footprint.";
    case on_the_fly_with_fixed_points: return "In addition to 2, investigate for generated"
        " variables whether they occur on a loop, such that"
        " they can be set to true or false, depending on the"
        " fixed point symbol. This can increase the time"
        " needed to generate an equation substantially.";
  }
  throw mcrl2::runtime_error("unknown transformation strategy");
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_TRANSFORMATION_STRATEGY_H
