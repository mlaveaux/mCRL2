// Author(s): Jan Friso Groote. Based on pbes/pbes_rewriter_type.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/pres_rewriter_type.h
/// \brief add your file description here.

#ifndef MCRL2_PRES_PRES_REWRITER_TYPE_H
#define MCRL2_PRES_PRES_REWRITER_TYPE_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pres_system
{

/// \brief An enumerated type for the available pres rewriters
enum pres_rewriter_type
{
  simplify,
  quantifier_all,
  quantifier_finite,
  quantifier_inside,
  quantifier_one_point,

};

/// \brief Parses a pres rewriter type
inline
pres_rewriter_type parse_pres_rewriter_type(const std::string& type)
{
  if (type == "simplify")
  {
    return simplify;
  }
  if (type == "quantifier-all")
  {
    return quantifier_all;
  }
  if (type == "quantifier-finite")
  {
    return quantifier_finite;
  }
  if (type == "quantifier-inside")
  {
    return quantifier_inside;
  }
  if (type == "quantifier-one-point")
  {
    return quantifier_one_point;
  }
  throw mcrl2::runtime_error("unknown pres rewriter option " + type);
}

/// \brief Prints a pres rewriter type
inline
std::string print_pres_rewriter_type(const pres_rewriter_type type)
{
  switch (type)
  {
    case simplify:
      return "simplify";
    case quantifier_all:
      return "quantifier-all";
    case quantifier_finite:
      return "quantifier-finite";
    case quantifier_inside:
      return "quantifier-inside";
    case quantifier_one_point:
      return "quantifier-one-point";
    default:
    return "unknown pres rewriter";
  }
}

/// \brief Returns a description of a pres rewriter
inline
std::string description(const pres_rewriter_type type)
{
  switch (type)
  {
    case simplify          :
      return "for simplification";
    case quantifier_all    :
      return "for eliminating all quantifiers";
    case quantifier_finite :
      return "for eliminating finite quantifier variables";
    case quantifier_inside :
      return "for pushing quantifiers inside";
    case quantifier_one_point :
      return "for one point rule quantifier elimination";

  }
  throw mcrl2::runtime_error("unknown pres rewriter");
}

/// \brief Stream operator for rewriter type
/// \param is An input stream
/// \param t A rewriter type
/// \return The input stream
inline
std::istream& operator>>(std::istream& is, pres_rewriter_type& t)
{
  std::string s;
  is >> s;
  try
  {
    t = parse_pres_rewriter_type(s);
  }
  catch (const mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

inline
std::ostream& operator<<(std::ostream& os, const pres_rewriter_type t)
{
  os << print_pres_rewriter_type(t);
  return os;
}

} // namespace mcrl2::pres_system

#endif // MCRL2_PRES_PRES_REWRITER_TYPE_H
