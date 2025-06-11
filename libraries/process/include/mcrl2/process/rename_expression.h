// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/rename_expression.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_RENAME_EXPRESSION_H
#define MCRL2_PROCESS_RENAME_EXPRESSION_H

#include "mcrl2/core/detail/default_values.h"
#include "mcrl2/core/detail/soundness_checks.h"

namespace mcrl2
{

namespace process
{

//--- start generated class rename_expression ---//
/// \\brief A rename expression
class rename_expression: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    rename_expression()
      : atermpp::aterm(core::detail::default_values::RenameExpr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit rename_expression(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_RenameExpr(*this));
    }

    /// \\brief Constructor Z12.
    rename_expression(core::identifier_string& source, core::identifier_string& target)
      : atermpp::aterm(core::detail::function_symbol_RenameExpr(), source, target)
    {}

    /// \\brief Constructor Z1.
    rename_expression(const std::string& source, const std::string& target)
      : atermpp::aterm(core::detail::function_symbol_RenameExpr(), core::identifier_string(source), core::identifier_string(target))
    {}

    /// Move semantics
    rename_expression(const rename_expression&) noexcept = default;
    rename_expression(rename_expression&&) noexcept = default;
    rename_expression& operator=(const rename_expression&) noexcept = default;
    rename_expression& operator=(rename_expression&&) noexcept = default;

    const core::identifier_string& source() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const core::identifier_string& target() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[1]);
    }
};

/// \\brief Make_rename_expression constructs a new term into a given address.
/// \\ \param t The reference into which the new rename_expression is constructed. 
template <class... ARGUMENTS>
inline void make_rename_expression(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_RenameExpr(), args...);
}

/// \\brief list of rename_expressions
typedef atermpp::term_list<rename_expression> rename_expression_list;

/// \\brief vector of rename_expressions
typedef std::vector<rename_expression>    rename_expression_vector;

/// \\brief Test for a rename_expression expression
/// \\param x A term
/// \\return True if \\a x is a rename_expression expression
inline
bool is_rename_expression(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::RenameExpr;
}

// prototype declaration
std::string pp(const rename_expression& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const rename_expression& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(rename_expression& t1, rename_expression& t2)
{
  t1.swap(t2);
}
//--- end generated class rename_expression ---//

} // namespace process

} // namespace mcrl2

#endif // MCRL2_PROCESS_RENAME_EXPRESSION_H
