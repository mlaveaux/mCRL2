// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/assignment.h
/// \brief The class assignment.

#ifndef MCRL2_DATA_ASSIGNMENT_H
#define MCRL2_DATA_ASSIGNMENT_H

#include "mcrl2/data/undefined.h"
#include "mcrl2/data/untyped_identifier.h"

namespace mcrl2::data
{

//--- start generated classes ---//
/// \\brief Assignment expression
class assignment_expression: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    assignment_expression()
      : atermpp::aterm(core::detail::default_values::WhrDecl)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit assignment_expression(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_WhrDecl(*this));
    }

    /// Move semantics
    assignment_expression(const assignment_expression&) noexcept = default;
    assignment_expression(assignment_expression&&) noexcept = default;
    assignment_expression& operator=(const assignment_expression&) noexcept = default;
    assignment_expression& operator=(assignment_expression&&) noexcept = default;
};

/// \\brief list of assignment_expressions
using assignment_expression_list = atermpp::term_list<assignment_expression>;

/// \\brief vector of assignment_expressions
using assignment_expression_vector = std::vector<assignment_expression>;

// prototypes
inline bool is_assignment(const atermpp::aterm& x);
inline bool is_untyped_identifier_assignment(const atermpp::aterm& x);

/// \\brief Test for a assignment_expression expression
/// \\param x A term
/// \\return True if \\a x is a assignment_expression expression
inline
bool is_assignment_expression(const atermpp::aterm& x)
{
  return data::is_assignment(x) ||
         data::is_untyped_identifier_assignment(x);
}

// prototype declaration
std::string pp(const assignment_expression& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const assignment_expression& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(assignment_expression& t1, assignment_expression& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Assignment of a data expression to a variable
class assignment: public assignment_expression
{
  public:
    /// \\brief Default constructor X3.
    assignment()
      : assignment_expression(core::detail::default_values::DataVarIdInit)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit assignment(const atermpp::aterm& term)
      : assignment_expression(term)
    {
      assert(core::detail::check_term_DataVarIdInit(*this));
    }

    /// \\brief Constructor Z14.
    assignment(const variable& lhs, const data_expression& rhs)
      : assignment_expression(atermpp::aterm(core::detail::function_symbol_DataVarIdInit(), lhs, rhs))
    {}

    /// Move semantics
    assignment(const assignment&) noexcept = default;
    assignment(assignment&&) noexcept = default;
    assignment& operator=(const assignment&) noexcept = default;
    assignment& operator=(assignment&&) noexcept = default;

    const variable& lhs() const
    {
      return atermpp::down_cast<variable>((*this)[0]);
    }

    const data_expression& rhs() const
    {
      return atermpp::down_cast<data_expression>((*this)[1]);
    }
//--- start user section assignment ---//
    /// \brief Applies the assignment to a variable
    /// \param[in] x A variable
    /// \return The value <tt>x[lhs() := rhs()]</tt>.
    const data_expression& operator()(const variable& x) const
    {
      return x == lhs() ? rhs() : x;
    }
//--- end user section assignment ---//
};

/// \\brief Make_assignment constructs a new term into a given address.
/// \\ \param t The reference into which the new assignment is constructed. 
template <class... ARGUMENTS>
inline void make_assignment(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_DataVarIdInit(), args...);
}

/// \\brief list of assignments
using assignment_list = atermpp::term_list<assignment>;

/// \\brief vector of assignments
using assignment_vector = std::vector<assignment>;

/// \\brief Test for a assignment expression
/// \\param x A term
/// \\return True if \\a x is a assignment expression
inline
bool is_assignment(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::DataVarIdInit;
}

// prototype declaration
std::string pp(const assignment& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const assignment& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(assignment& t1, assignment& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Assignment of a data expression to a string
class untyped_identifier_assignment: public assignment_expression
{
  public:
    /// \\brief Default constructor X3.
    untyped_identifier_assignment()
      : assignment_expression(core::detail::default_values::UntypedIdentifierAssignment)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit untyped_identifier_assignment(const atermpp::aterm& term)
      : assignment_expression(term)
    {
      assert(core::detail::check_term_UntypedIdentifierAssignment(*this));
    }

    /// \\brief Constructor Z14.
    untyped_identifier_assignment(const core::identifier_string& lhs, const data_expression& rhs)
      : assignment_expression(atermpp::aterm(core::detail::function_symbol_UntypedIdentifierAssignment(), lhs, rhs))
    {}

    /// \\brief Constructor Z2.
    untyped_identifier_assignment(const std::string& lhs, const data_expression& rhs)
      : assignment_expression(atermpp::aterm(core::detail::function_symbol_UntypedIdentifierAssignment(), core::identifier_string(lhs), rhs))
    {}

    /// Move semantics
    untyped_identifier_assignment(const untyped_identifier_assignment&) noexcept = default;
    untyped_identifier_assignment(untyped_identifier_assignment&&) noexcept = default;
    untyped_identifier_assignment& operator=(const untyped_identifier_assignment&) noexcept = default;
    untyped_identifier_assignment& operator=(untyped_identifier_assignment&&) noexcept = default;

    const core::identifier_string& lhs() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const data_expression& rhs() const
    {
      return atermpp::down_cast<data_expression>((*this)[1]);
    }
//--- start user section untyped_identifier_assignment ---//
    /// \brief Applies the assignment to a variable
    /// \param[in] x An identifier string
    /// \return The value <tt>x[lhs() := rhs()]</tt>.
    data_expression operator()(const untyped_identifier& x) const
    {
      return x == lhs() ? rhs() : data_expression(x);
    }
//--- end user section untyped_identifier_assignment ---//
};

/// \\brief Make_untyped_identifier_assignment constructs a new term into a given address.
/// \\ \param t The reference into which the new untyped_identifier_assignment is constructed. 
template <class... ARGUMENTS>
inline void make_untyped_identifier_assignment(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_UntypedIdentifierAssignment(), args...);
}

/// \\brief list of untyped_identifier_assignments
using untyped_identifier_assignment_list = atermpp::term_list<untyped_identifier_assignment>;

/// \\brief vector of untyped_identifier_assignments
using untyped_identifier_assignment_vector = std::vector<untyped_identifier_assignment>;

/// \\brief Test for a untyped_identifier_assignment expression
/// \\param x A term
/// \\return True if \\a x is a untyped_identifier_assignment expression
inline
bool is_untyped_identifier_assignment(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::UntypedIdentifierAssignment;
}

// prototype declaration
std::string pp(const untyped_identifier_assignment& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const untyped_identifier_assignment& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(untyped_identifier_assignment& t1, untyped_identifier_assignment& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated classes ---//

/// \brief Constructs an assignment_list by pairwise combining a variable and expression
/// \param variables A sequence of data variables
/// \param expressions A sequence of data expressions
/// \return The corresponding assignment list.
template < typename VariableSequence, typename ExpressionSequence >
assignment_vector make_assignment_vector(VariableSequence const& variables, ExpressionSequence const& expressions)
{
  std::vector<assignment> result;
  typename ExpressionSequence::const_iterator j=expressions.begin();
  for(typename VariableSequence::const_iterator i=variables.begin(); i!=variables.end(); ++i, ++j)
  {
    assert(j!=expressions.end());
    result.push_back(assignment(*i,*j));
  }
  assert(j==expressions.end());
  return result;
}

/// \brief Converts an iterator range to data_expression_list
/// \param variables A sequence of variables.
/// \param expressions A sequence of variables.
/// \pre length(variables) == length(expressions)
/// \note This function uses implementation details of the iterator type
/// and hence is sometimes efficient than copying all elements of the list.
template < typename VariableSequence, typename ExpressionSequence >
assignment_list make_assignment_list(const VariableSequence& variables, const ExpressionSequence& expressions)
{
  std::vector<assignment> result(make_assignment_vector(variables,expressions));
  return assignment_list(result.begin(),result.end());
}

/// \brief Returns the left hand sides of an assignment list
/// \param x An assignment list
inline
variable_list left_hand_sides(const assignment_list& x)
{
  return variable_list(x.begin(), x.end(), [](const assignment& a) { return a.lhs(); });
}

/// \brief Returns the right hand sides of an assignment list
/// \param x An assignment list
inline
data_expression_list right_hand_sides(const assignment_list& x)
{
  return data_expression_list(x.begin(), x.end(), [](const assignment& a) { return a.rhs(); });
}

// template function overloads
std::string pp(const assignment_list& x, bool precedence_aware = true);
std::string pp(const assignment_vector& x, bool precedence_aware = true);

} // namespace mcrl2::data



#endif // MCRL2_DATA_ASSIGNMENT_H

