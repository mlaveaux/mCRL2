// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/variable.h
/// \brief The class variable.

#ifndef MCRL2_DATA_VARIABLE_H
#define MCRL2_DATA_VARIABLE_H

#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/core/identifier_string.h"
#include "mcrl2/data/data_expression.h"

namespace mcrl2
{

namespace data
{

//--- start generated class variable ---//
/// \\brief A data variable
class variable: public data_expression
{
  public:


    /// Move semantics
    variable(const variable&) noexcept = default;
    variable(variable&&) noexcept = default;
    variable& operator=(const variable&) noexcept = default;
    variable& operator=(variable&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const sort_expression& sort() const
    {
      return atermpp::down_cast<sort_expression>((*this)[1]);
    }
//--- start user section variable ---//
    /// \brief Default constructor.
    variable()
      : data_expression(core::detail::default_values::DataVarId)
    {}

    /// \brief Constructor.
    /// \param term A term
    explicit variable(const atermpp::aterm& term)
      : data_expression(term)
    {
      assert(core::detail::check_term_DataVarId(*this));
    }

    /// \brief Constructor.
    variable(const core::identifier_string& name, const sort_expression& sort)
    {
      atermpp::make_term_appl(*this, core::detail::function_symbol_DataVarId(), name, sort);
    }


    /// \brief Constructor.
    variable(const std::string& name, const sort_expression& sort)
    {
      atermpp::make_term_appl(*this, core::detail::function_symbol_DataVarId(), core::identifier_string(name), sort);
    }

//--- end user section variable ---//
};

/// \\brief Make_variable constructs a new term into a given address.
/// \\ \param t The reference into which the new variable is constructed. 
template <class... ARGUMENTS>
inline void make_variable(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_DataVarId(), args...);
}

/// \\brief list of variables
typedef atermpp::term_list<variable> variable_list;

/// \\brief vector of variables
typedef std::vector<variable>    variable_vector;

// prototype declaration
std::string pp(const variable& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const variable& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(variable& t1, variable& t2)
{
  t1.swap(t2);
}
//--- end generated class variable ---//


// template function overloads
std::string pp(const variable_list& x, bool precedence_aware = true);
std::string pp(const variable_vector& x, bool precedence_aware = true);
std::string pp(const std::set<variable>& x, bool precedence_aware = true);
std::set<data::variable> find_all_variables(const data::variable& x);
std::set<data::variable> find_all_variables(const data::variable_list& x);
std::set<core::identifier_string> find_identifiers(const data::variable_list& x);

} // namespace data

} // namespace mcrl2


namespace std 
{

template<>
struct hash<mcrl2::data::variable>
{
  // Default constructor, required for const qualified hash functions. 
    hash()
    {}
  
    std::size_t operator()(const mcrl2::data::variable& v) const
    {
      const hash<atermpp::aterm> hasher; 
      return hasher(v);
    }
};
  
} // namespace std

#endif // MCRL2_DATA_VARIABLE_H

