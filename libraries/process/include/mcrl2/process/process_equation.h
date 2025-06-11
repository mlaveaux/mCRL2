// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/process_equation.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_PROCESS_EQUATION_H
#define MCRL2_PROCESS_PROCESS_EQUATION_H

#include "mcrl2/process/process_expression.h"

namespace mcrl2
{

namespace process
{

//--- start generated class process_equation ---//
/// \\brief A process equation
class process_equation: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    process_equation()
      : atermpp::aterm(core::detail::default_values::ProcEqn)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit process_equation(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_ProcEqn(*this));
    }

    /// \\brief Constructor Z12.
    process_equation(const process_identifier& identifier, const data::variable_list& formal_parameters, const process_expression& expression)
      : atermpp::aterm(core::detail::function_symbol_ProcEqn(), identifier, formal_parameters, expression)
    {}

    /// Move semantics
    process_equation(const process_equation&) noexcept = default;
    process_equation(process_equation&&) noexcept = default;
    process_equation& operator=(const process_equation&) noexcept = default;
    process_equation& operator=(process_equation&&) noexcept = default;

    const process_identifier& identifier() const
    {
      return atermpp::down_cast<process_identifier>((*this)[0]);
    }

    const data::variable_list& formal_parameters() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[1]);
    }

    const process_expression& expression() const
    {
      return atermpp::down_cast<process_expression>((*this)[2]);
    }
};

/// \\brief Make_process_equation constructs a new term into a given address.
/// \\ \param t The reference into which the new process_equation is constructed. 
template <class... ARGUMENTS>
inline void make_process_equation(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ProcEqn(), args...);
}

/// \\brief list of process_equations
typedef atermpp::term_list<process_equation> process_equation_list;

/// \\brief vector of process_equations
typedef std::vector<process_equation>    process_equation_vector;

/// \\brief Test for a process_equation expression
/// \\param x A term
/// \\return True if \\a x is a process_equation expression
inline
bool is_process_equation(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ProcEqn;
}

// prototype declaration
std::string pp(const process_equation& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const process_equation& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(process_equation& t1, process_equation& t2)
{
  t1.swap(t2);
}
//--- end generated class process_equation ---//

// template function overloads
std::string pp(const process_equation_list& x);
std::string pp(const process_equation_vector& x);
void normalize_sorts(process_equation_vector& x, const data::sort_specification& sortspec);
std::set<data::sort_expression> find_sort_expressions(const process::process_equation_vector& x);

} // namespace process

} // namespace mcrl2

#endif // MCRL2_PROCESS_PROCESS_EQUATION_H
