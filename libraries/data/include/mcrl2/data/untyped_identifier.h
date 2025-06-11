// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/untyped_identifier.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_UNTYPED_IDENTIFIER_H
#define MCRL2_DATA_UNTYPED_IDENTIFIER_H

#include "mcrl2/data/data_expression.h"

namespace mcrl2 {

namespace data {

//--- start generated class untyped_identifier ---//
/// \\brief An untyped identifier
class untyped_identifier: public data_expression
{
  public:
    /// \\brief Default constructor X3.
    untyped_identifier()
      : data_expression(core::detail::default_values::UntypedIdentifier)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit untyped_identifier(const atermpp::aterm& term)
      : data_expression(term)
    {
      assert(core::detail::check_term_UntypedIdentifier(*this));
    }

    /// \\brief Constructor Z14.
    explicit untyped_identifier(const core::identifier_string& name)
      : data_expression(atermpp::aterm(core::detail::function_symbol_UntypedIdentifier(), name))
    {}

    /// \\brief Constructor Z2.
    untyped_identifier(const std::string& name)
      : data_expression(atermpp::aterm(core::detail::function_symbol_UntypedIdentifier(), core::identifier_string(name)))
    {}

    /// Move semantics
    untyped_identifier(const untyped_identifier&) noexcept = default;
    untyped_identifier(untyped_identifier&&) noexcept = default;
    untyped_identifier& operator=(const untyped_identifier&) noexcept = default;
    untyped_identifier& operator=(untyped_identifier&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }
};

/// \\brief Make_untyped_identifier constructs a new term into a given address.
/// \\ \param t The reference into which the new untyped_identifier is constructed. 
template <class... ARGUMENTS>
inline void make_untyped_identifier(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_UntypedIdentifier(), args...);
}

// prototype declaration
std::string pp(const untyped_identifier& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const untyped_identifier& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(untyped_identifier& t1, untyped_identifier& t2)
{
  t1.swap(t2);
}
//--- end generated class untyped_identifier ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_UNTYPED_IDENTIFIER_H
