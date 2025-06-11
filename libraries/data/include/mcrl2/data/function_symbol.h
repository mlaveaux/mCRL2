// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/function_symbol.h
/// \brief The class function symbol.

#ifndef MCRL2_DATA_FUNCTION_SYMBOL_H
#define MCRL2_DATA_FUNCTION_SYMBOL_H

#include "mcrl2/data/application.h"

namespace mcrl2
{

namespace data
{

typedef std::pair<core::identifier_string, sort_expression> function_symbol_key_type;

//--- start generated class function_symbol ---//
/// \\brief A function symbol
class function_symbol: public data_expression
{
  public:


    /// Move semantics
    function_symbol(const function_symbol&) noexcept = default;
    function_symbol(function_symbol&&) noexcept = default;
    function_symbol& operator=(const function_symbol&) noexcept = default;
    function_symbol& operator=(function_symbol&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const sort_expression& sort() const
    {
      return atermpp::down_cast<sort_expression>((*this)[1]);
    }
//--- start user section function_symbol ---//
    /// \brief Default constructor.
    function_symbol()
      : data_expression(core::detail::default_values::OpId)
    {}

    /// \brief Constructor.
    /// \param term A term
    explicit function_symbol(const atermpp::aterm& term)
      : data_expression(term)
    {
      assert(core::detail::check_term_OpId(*this));
    }

    /// \brief Constructor.
    function_symbol(const core::identifier_string& name, const sort_expression& sort)
    {
      atermpp::make_term_appl_with_index<function_symbol, std::pair<core::identifier_string, sort_expression> >
                      (*this, core::detail::function_symbol_OpId(), name, sort);
    }


    /// \brief Constructor.
    function_symbol(const std::string& name, const sort_expression& sort)
    { 
      atermpp::make_term_appl_with_index<function_symbol, std::pair<core::identifier_string, sort_expression> >
                      (*this, core::detail::function_symbol_OpId(), core::identifier_string(name), sort);
    }

//--- end user section function_symbol ---//
};

/// \\brief Make_function_symbol constructs a new term into a given address.
/// \\ \param t The reference into which the new function_symbol is constructed. 
template <class... ARGUMENTS>
inline void make_function_symbol(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl_with_index<function_symbol,std::pair<core::identifier_string, sort_expression>>(t, core::detail::function_symbol_OpId(), args...);
}

/// \\brief list of function_symbols
typedef atermpp::term_list<function_symbol> function_symbol_list;

/// \\brief vector of function_symbols
typedef std::vector<function_symbol>    function_symbol_vector;

// prototype declaration
std::string pp(const function_symbol& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const function_symbol& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(function_symbol& t1, function_symbol& t2)
{
  t1.swap(t2);
}
//--- end generated class function_symbol ---//

// template function overloads
std::string pp(const function_symbol_list& x, bool precedence_aware = true);
std::string pp(const function_symbol_vector& x, bool precedence_aware = true);
std::set<data::variable> find_all_variables(const data::function_symbol& x);

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_FUNCTION_SYMBOL_H

