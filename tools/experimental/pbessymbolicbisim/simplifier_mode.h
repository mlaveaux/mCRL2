// Author(s): Thomas Neele
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_PBESSYMBOLICBISIM_SIMPLIFIER_MODE_H
#define MCRL2_PBESSYMBOLICBISIM_SIMPLIFIER_MODE_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

enum simplifier_mode
{
  simplify_fm,
#ifdef DBM_PACKAGE_AVAILABLE
  simplify_dbm,
#endif
  simplify_finite_domain,
  simplify_identity,
  simplify_auto
};

inline
simplifier_mode parse_mode(const std::string& s)
{
  if(s == "fm")
  {
    return simplify_fm;
  }
#ifdef DBM_PACKAGE_AVAILABLE
  else if(s == "dbm")
  {
    return simplify_dbm;
  }
#endif
  else if(s == "finite" || s == "f")
  {
    return simplify_finite_domain;
  }
  else if(s == "identity" || s == "id")
  {
    return simplify_identity;
  }
  else if(s == "auto" || s == "a")
  {
    return simplify_auto;
  }
  else
  {
    throw mcrl2::runtime_error("Unknown simplifier " + s);
  }
}

// \overload
inline
std::istream& operator>>(std::istream& is, simplifier_mode& eq)
{
  try
  {
    std::string s;
    is >> s;
    eq = parse_mode(s);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

inline std::string print_simplifier_mode(const simplifier_mode r)
{
  switch(r)
  {
    case simplify_fm:
      return "fm";
#ifdef DBM_PACKAGE_AVAILABLE
    case simplify_dbm:
      return "dbm";
#endif
    case simplify_finite_domain:
      return "finite";
    case simplify_identity:
      return "identity";
    case simplify_auto:
      return "auto";
    default:
      throw mcrl2::runtime_error("Unknown simplifier_mode.");
  }
}

// \overload
inline
std::ostream& operator<<(std::ostream& os, const simplifier_mode eq)
{
  os << print_simplifier_mode(eq);
  return os;
}

/** \brief Gives a description of an equivalence.
 * \param[in] eq The equivalence type.
 * \return A string describing the equivalence specified by \a
 * eq. */
inline std::string description(const simplifier_mode r)
{
  switch(r)
  {
    case simplify_fm:
      return "Use functions from the mCRL2 data library to eliminate redundant inequalities";
#ifdef DBM_PACKAGE_AVAILABLE
    case simplify_dbm:
      return "Use the UPPAAL DBM package to simplify inequalities";
#endif
    case simplify_finite_domain:
      return "Only simplify expressions over finite discrete data";
    case simplify_identity:
      return "Do not simplify expressions";
    case simplify_auto:
      return "Automatically select the best simplifier";
    default:
      throw mcrl2::runtime_error("Unknown simplifier_mode.");
  }
}

#endif // MCRL2_PBESSYMBOLICBISIM_SIMPLIFIER_MODE_H
