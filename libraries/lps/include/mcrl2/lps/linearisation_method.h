// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/linearisation_method.h

#ifndef MCRL2_LPS_LINEARISATION_METHOD_H
#define MCRL2_LPS_LINEARISATION_METHOD_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::lps
{

/// \brief The available linearisation methods
using t_lin_method = enum
{
  lmStack,
  lmRegular,
  lmRegular2
};

/// \brief String representation of a linearisation method
/// \param[in] lin_method A linerisation method
/// \return String representation of \a lin_method
inline
std::string print_lin_method(const t_lin_method lin_method)
{
  switch(lin_method)
  {
    case lmStack: return "stack";
    case lmRegular: return "regular";
    case lmRegular2: return "regular2";
    default: throw mcrl2::runtime_error("unknown linearisation method");
  }
}

inline
std::string description(const t_lin_method lin_method)
{
  switch(lin_method)
  {
    case lmStack: return "for using stack data types (useful when 'regular' and 'regular2' do not work)";
    case lmRegular: return "for generating an LPS in regular form (specification should be regular)";
    case lmRegular2: return "for a variant of 'regular' that uses more data variables (useful when 'regular' does not work)";
    default: throw mcrl2::runtime_error("unknown linearisation method");
  }
}

/// \brief Parse a linearisation method
/// \param[in] s A string
/// \return The linearisation method represented by \a s
inline
t_lin_method parse_lin_method(const std::string& s)
{
  if(s == "stack")
  {
    return lmStack;
  }
  else if (s == "regular")
  {
    return lmRegular;
  }
  else if (s == "regular2")
  {
    return lmRegular2;
  }
  else
  {
    throw mcrl2::runtime_error("unknown linearisation strategy " + s);
  }
}

// \overload
inline
std::istream& operator>>(std::istream& is, t_lin_method& l)
{
  try {
    std::stringbuf buffer;
    is >> &buffer;
    l = parse_lin_method(buffer.str());
  }
  catch (mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }

  return is;
}

// \overload
inline
std::ostream& operator<<(std::ostream& os, const t_lin_method l)
{
  os << print_lin_method(l);
  return os;
}

} // namespace mcrl2::lps

#endif // MCRL2_LPS_LINEARISATION_METHOD_H
