// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/untyped_sort_variable.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_UNTYPED_SORT_VARIABLE_H
#define MCRL2_DATA_UNTYPED_SORT_VARIABLE_H

#include "mcrl2/data/sort_expression.h"

namespace mcrl2 {

namespace data {

//--- start generated class untyped_sort_variable ---//
/// \\brief Untyped sort variable
class untyped_sort_variable: public sort_expression
{
  public:
    /// \\brief Default constructor X3.
    untyped_sort_variable()
      : sort_expression(core::detail::default_values::UntypedSortVariable)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit untyped_sort_variable(const atermpp::aterm& term)
      : sort_expression(term)
    {
      assert(core::detail::check_term_UntypedSortVariable(*this));
    }

    /// \\brief Constructor Z14.
    explicit untyped_sort_variable(const atermpp::aterm_int& value)
      : sort_expression(atermpp::aterm(core::detail::function_symbol_UntypedSortVariable(), value))
    {}

    /// Move semantics
    untyped_sort_variable(const untyped_sort_variable&) noexcept = default;
    untyped_sort_variable(untyped_sort_variable&&) noexcept = default;
    untyped_sort_variable& operator=(const untyped_sort_variable&) noexcept = default;
    untyped_sort_variable& operator=(untyped_sort_variable&&) noexcept = default;

    const atermpp::aterm_int& value() const
    {
      return atermpp::down_cast<atermpp::aterm_int>((*this)[0]);
    }
//--- start user section untyped_sort_variable ---//
    /// \brief Constructor.
    untyped_sort_variable(std::size_t value)
      : sort_expression(atermpp::aterm(core::detail::function_symbol_UntypedSortVariable(), atermpp::aterm_int(value)))
    {}
//--- end user section untyped_sort_variable ---//
};

/// \\brief Make_untyped_sort_variable constructs a new term into a given address.
/// \\ \param t The reference into which the new untyped_sort_variable is constructed. 
template <class... ARGUMENTS>
inline void make_untyped_sort_variable(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_UntypedSortVariable(), args...);
}

/// \\brief Test for a untyped_sort_variable expression
/// \\param x A term
/// \\return True if \\a x is a untyped_sort_variable expression
inline
bool is_untyped_sort_variable(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::UntypedSortVariable;
}

// prototype declaration
std::string pp(const untyped_sort_variable& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const untyped_sort_variable& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(untyped_sort_variable& t1, untyped_sort_variable& t2)
{
  t1.swap(t2);
}
//--- end generated class untyped_sort_variable ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_UNTYPED_SORT_VARIABLE_H
