// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/communication_expression.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_COMMUNICATION_EXPRESSION_H
#define MCRL2_PROCESS_COMMUNICATION_EXPRESSION_H

#include "mcrl2/process/action_name_multiset.h"

namespace mcrl2
{

namespace process
{

//--- start generated class communication_expression ---//
/// \\brief A communication expression
class communication_expression: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    communication_expression()
      : atermpp::aterm(core::detail::default_values::CommExpr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit communication_expression(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_CommExpr(*this));
    }

    /// \\brief Constructor Z12.
    communication_expression(const action_name_multiset& action_name, const core::identifier_string& name)
      : atermpp::aterm(core::detail::function_symbol_CommExpr(), action_name, name)
    {}

    /// \\brief Constructor Z1.
    communication_expression(const action_name_multiset& action_name, const std::string& name)
      : atermpp::aterm(core::detail::function_symbol_CommExpr(), action_name, core::identifier_string(name))
    {}

    /// Move semantics
    communication_expression(const communication_expression&) noexcept = default;
    communication_expression(communication_expression&&) noexcept = default;
    communication_expression& operator=(const communication_expression&) noexcept = default;
    communication_expression& operator=(communication_expression&&) noexcept = default;

    const action_name_multiset& action_name() const
    {
      return atermpp::down_cast<action_name_multiset>((*this)[0]);
    }

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[1]);
    }
};

/// \\brief Make_communication_expression constructs a new term into a given address.
/// \\ \param t The reference into which the new communication_expression is constructed. 
template <class... ARGUMENTS>
inline void make_communication_expression(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_CommExpr(), args...);
}

/// \\brief list of communication_expressions
typedef atermpp::term_list<communication_expression> communication_expression_list;

/// \\brief vector of communication_expressions
typedef std::vector<communication_expression>    communication_expression_vector;

/// \\brief Test for a communication_expression expression
/// \\param x A term
/// \\return True if \\a x is a communication_expression expression
inline
bool is_communication_expression(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::CommExpr;
}

// prototype declaration
std::string pp(const communication_expression& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const communication_expression& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(communication_expression& t1, communication_expression& t2)
{
  t1.swap(t2);
}
//--- end generated class communication_expression ---//

} // namespace process

} // namespace mcrl2

#endif // MCRL2_PROCESS_COMMUNICATION_EXPRESSION_H
