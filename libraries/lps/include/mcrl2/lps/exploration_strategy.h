// Author(s): Muck van Weerdenburg
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef MCRL2_LPS_EXPLORATION_STRATEGY_H
#define MCRL2_LPS_EXPLORATION_STRATEGY_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::lps
{

enum exploration_strategy { es_none,
                            es_breadth,
                            es_depth,
                            es_random,
                            es_value_prioritized,
                            es_value_random_prioritized,
                            es_highway
                          };

inline
exploration_strategy parse_exploration_strategy(const std::string& s)
{
  if (s=="b" || s=="breadth")
  {
    return es_breadth;
  }
  if (s=="d" || s=="depth")
  {
    return es_depth;
  }
  if (s=="r" || s=="random")
  {
    return es_random;
  }
  if (s=="p" || s=="prioritized")
  {
    return es_value_prioritized;
  }
  if (s=="q" || s=="rprioritized")
  {
    return es_value_random_prioritized;
  }
  if (s=="h" || s == "highway")
  {
    return es_highway;
  }
  return es_none;
}

inline
std::string print_exploration_strategy(const exploration_strategy es)
{
  switch (es)
  {
    case es_breadth:
      return "breadth";
    case es_depth:
      return "depth";
    case es_random:
      return "random";
    case es_value_prioritized:
      return "prioritized";
    case es_value_random_prioritized:
      return "rprioritized";
    case es_highway:
      return "highway";
    default:
      throw mcrl2::runtime_error("unknown exploration strategy");
  }
}

inline
std::istream& operator>>(std::istream& is, exploration_strategy& strat)
{
  try
  {
    std::string s;
    is >> s;
    strat = parse_exploration_strategy(s);
  }
  catch (mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

inline
std::ostream& operator<<(std::ostream& os, const exploration_strategy strat)
{
  os << print_exploration_strategy(strat);
  return os;
}

inline std::string description(const exploration_strategy strat)
{
  switch (strat)
  {
    case es_breadth:
      return "breadth-first search";
    case es_depth:
      return "depth-first search";
    case es_random:
      return "random simulation. Out of all next states one is chosen at random independently of whether this state has already been observed. Consequently, random simulation only terminates when a deadlocked state is encountered.";
    case es_value_prioritized:
      return "prioritize single actions on its first argument being of sort Nat where only those actions with the lowest value for this parameter are selected. E.g. if there are actions a(3) and b(4), a(3) remains and b(4) is skipped. Actions without a first parameter of sort Nat and multactions with more than one action are always chosen (option is experimental)";
    case es_value_random_prioritized:
      return "prioritize actions on its first argument being of sort Nat (see option --prioritized), and randomly select one of these to obtain a prioritized random simulation (option is experimental)";
    case es_highway:
      return "highway search. Only part of the state space is explored, by restricting the size of the todo list. N.B. The implementation deviates slightly from the published version.";
    default:
      throw mcrl2::runtime_error("unknown exploration_strategy");
  }
}

} // namespace mcrl2::lps

#endif // MCRL2_LPS_EXPLORATION_STRATEGY_H
