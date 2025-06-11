// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/where_clause.h
/// \brief The class where_clause.

#ifndef MCRL2_DATA_WHERE_CLAUSE_H
#define MCRL2_DATA_WHERE_CLAUSE_H

#include "mcrl2/data/assignment.h"

namespace mcrl2
{

namespace data
{

//--- start generated class where_clause ---//
/// \\brief A where expression
class where_clause: public data_expression
{
  public:
    /// \\brief Default constructor X3.
    where_clause()
      : data_expression(core::detail::default_values::Whr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit where_clause(const atermpp::aterm& term)
      : data_expression(term)
    {
      assert(core::detail::check_term_Whr(*this));
    }

    /// \\brief Constructor Z14.
    where_clause(const data_expression& body, const assignment_expression_list& declarations)
      : data_expression(atermpp::aterm(core::detail::function_symbol_Whr(), body, declarations))
    {}

    /// \\brief Constructor Z2.
    template <typename Container>
    where_clause(const data_expression& body, const Container& declarations, typename atermpp::enable_if_container<Container, assignment_expression>::type* = nullptr)
      : data_expression(atermpp::aterm(core::detail::function_symbol_Whr(), body, assignment_expression_list(declarations.begin(), declarations.end())))
    {}

    /// Move semantics
    where_clause(const where_clause&) noexcept = default;
    where_clause(where_clause&&) noexcept = default;
    where_clause& operator=(const where_clause&) noexcept = default;
    where_clause& operator=(where_clause&&) noexcept = default;

    const data_expression& body() const
    {
      return atermpp::down_cast<data_expression>((*this)[0]);
    }

    const assignment_expression_list& declarations() const
    {
      return atermpp::down_cast<assignment_expression_list>((*this)[1]);
    }
//--- start user section where_clause ---//
    const assignment_list& assignments() const
    {
      return atermpp::down_cast<const assignment_list>((*this)[1]);
    }
//--- end user section where_clause ---//
};

/// \\brief Make_where_clause constructs a new term into a given address.
/// \\ \param t The reference into which the new where_clause is constructed. 
template <class... ARGUMENTS>
inline void make_where_clause(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Whr(), args...);
}

// prototype declaration
std::string pp(const where_clause& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const where_clause& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(where_clause& t1, where_clause& t2)
{
  t1.swap(t2);
}
//--- end generated class where_clause ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_WHERE_CLAUSE_H

