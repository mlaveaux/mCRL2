// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/data_equation.h
/// \brief The class data_equation.

#ifndef MCRL2_DATA_DATA_EQUATION_H
#define MCRL2_DATA_DATA_EQUATION_H

#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_symbol.h"

namespace mcrl2
{

namespace data
{

// predeclare
namespace sort_bool
{
function_symbol const& true_();
} // namespace sort_bool

//--- start generated class data_equation ---//
/// \\brief A data equation
class data_equation: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    data_equation()
      : atermpp::aterm(core::detail::default_values::DataEqn)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit data_equation(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_DataEqn(*this));
    }

    /// \\brief Constructor Z12.
    data_equation(const variable_list& variables, const data_expression& condition, const data_expression& lhs, const data_expression& rhs)
      : atermpp::aterm(core::detail::function_symbol_DataEqn(), variables, condition, lhs, rhs)
    {}

    /// \\brief Constructor Z1.
    template <typename Container>
    data_equation(const Container& variables, const data_expression& condition, const data_expression& lhs, const data_expression& rhs, typename atermpp::enable_if_container<Container, variable>::type* = nullptr)
      : atermpp::aterm(core::detail::function_symbol_DataEqn(), variable_list(variables.begin(), variables.end()), condition, lhs, rhs)
    {}

    /// Move semantics
    data_equation(const data_equation&) noexcept = default;
    data_equation(data_equation&&) noexcept = default;
    data_equation& operator=(const data_equation&) noexcept = default;
    data_equation& operator=(data_equation&&) noexcept = default;

    const variable_list& variables() const
    {
      return atermpp::down_cast<variable_list>((*this)[0]);
    }

    const data_expression& condition() const
    {
      return atermpp::down_cast<data_expression>((*this)[1]);
    }

    const data_expression& lhs() const
    {
      return atermpp::down_cast<data_expression>((*this)[2]);
    }

    const data_expression& rhs() const
    {
      return atermpp::down_cast<data_expression>((*this)[3]);
    }
//--- start user section data_equation ---//
    /// \brief Constructor.
    ///
    /// \param[in] variables The free variables of the data_equation.
    /// \param[in] lhs The left hand side of the data_equation.
    /// \param[in] rhs The right hand side of the data_equation.
    /// \post this is the data equation representing the input, with
    ///       condition true
    template < typename Container >
    data_equation(const Container& variables,
                  const data_expression& lhs,
                  const data_expression& rhs,
                  typename atermpp::enable_if_container< Container, variable >::type* = nullptr)
      : atermpp::aterm(core::detail::function_symbol_DataEqn(), variable_list(variables.begin(),variables.end()), sort_bool::true_(), lhs, rhs)
    {}

    /// \brief Constructor.
    ///
    /// \param[in] lhs The left hand side of the data equation.
    /// \param[in] rhs The right hand side of the data equation.
    /// \post this is the data equations representing the input, without
    ///       variables, and condition true
    data_equation(const data_expression& lhs,
                  const data_expression& rhs)
      : atermpp::aterm(core::detail::function_symbol_DataEqn(), variable_list(), sort_bool::true_(), lhs, rhs)
    {}
//--- end user section data_equation ---//
};

/// \\brief Make_data_equation constructs a new term into a given address.
/// \\ \param t The reference into which the new data_equation is constructed. 
template <class... ARGUMENTS>
inline void make_data_equation(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_DataEqn(), args...);
}

/// \\brief list of data_equations
typedef atermpp::term_list<data_equation> data_equation_list;

/// \\brief vector of data_equations
typedef std::vector<data_equation>    data_equation_vector;

// prototype declaration
std::string pp(const data_equation& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const data_equation& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(data_equation& t1, data_equation& t2)
{
  t1.swap(t2);
}
//--- end generated class data_equation ---//

/// \brief Recognizer function.
/// \param[in] t A aterm of which it is checked whether it is a data_equation.
/// \returns True if the provided argument is a data_equation. 
inline bool is_data_equation(const atermpp::aterm& t)
{
  return t.function()==core::detail::function_symbol_DataEqn();
}

// template function overloads
std::string pp(const data_equation_list& x, bool precedence_aware = true);
std::string pp(const data_equation_vector& x, bool precedence_aware = true);
data::data_equation translate_user_notation(const data::data_equation& x);
std::set<data::sort_expression> find_sort_expressions(const data::data_equation& x);
std::set<data::function_symbol> find_function_symbols(const data::data_equation& x);

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_DATA_EQUATION_H

