// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/identifier_string.h
/// \brief aterm representations of identifier strings.

#ifndef MCRL2_CORE_IDENTIFIER_STRING_H
#define MCRL2_CORE_IDENTIFIER_STRING_H

#include "mcrl2/atermpp/aterm_string.h"
#include "mcrl2/core/detail/function_symbols.h"

namespace mcrl2
{

namespace core
{

/// \brief String type of the LPS library.
/// Identifier strings are represented internally as ATerms.
typedef atermpp::aterm_string identifier_string;

//--- start generated class identifier_string ---//
/// \\brief list of identifier_strings
typedef atermpp::term_list<identifier_string> identifier_string_list;

/// \\brief vector of identifier_strings
typedef std::vector<identifier_string>    identifier_string_vector;

// prototype declaration
std::string pp(const identifier_string& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const identifier_string& x)
{
  return out << core::pp(x);
}

/// \\brief swap overload
inline void swap(identifier_string& t1, identifier_string& t2)
{
  t1.swap(t2);
}
//--- end generated class identifier_string ---//

/// \brief Tests if a term is an identifier string.
/// \param[in] t A term
/// \return Whether t is an identifier string.
inline
bool is_identifier_string(const atermpp::aterm& t)
{
  return t.type_is_appl() && t.size() == 0;
}

/// \brief Provides the empty identifier string.
/// \return The empty identifier string. 
inline
identifier_string empty_identifier_string()
{
  return atermpp::empty_string();
}

} // namespace core

} // namespace mcrl2


namespace std {

/// \brief hash specialization
template<>
struct hash<mcrl2::core::identifier_string>
{
  std::size_t operator()(const mcrl2::core::identifier_string& x) const
  {
    return std::hash<atermpp::aterm>()(x);
  }
};

} // namespace std

#endif // MCRL2_CORE_IDENTIFIER_STRING_H
