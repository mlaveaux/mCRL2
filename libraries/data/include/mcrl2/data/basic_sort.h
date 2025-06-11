// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/basic_sort.h
/// \brief The class basic_sort.

#ifndef MCRL2_DATA_BASIC_SORT_H
#define MCRL2_DATA_BASIC_SORT_H

#include "mcrl2/data/sort_expression.h"

namespace mcrl2
{

namespace data
{

//--- start generated class basic_sort ---//
/// \\brief A basic sort
class basic_sort: public sort_expression
{
  public:
    /// \\brief Default constructor X3.
    basic_sort()
      : sort_expression(core::detail::default_values::SortId)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit basic_sort(const atermpp::aterm& term)
      : sort_expression(term)
    {
      assert(core::detail::check_term_SortId(*this));
    }

    /// \\brief Constructor Z14.
    explicit basic_sort(const core::identifier_string& name)
      : sort_expression(atermpp::aterm(core::detail::function_symbol_SortId(), name))
    {}

    /// \\brief Constructor Z2.
    basic_sort(const std::string& name)
      : sort_expression(atermpp::aterm(core::detail::function_symbol_SortId(), core::identifier_string(name)))
    {}

    /// Move semantics
    basic_sort(const basic_sort&) noexcept = default;
    basic_sort(basic_sort&&) noexcept = default;
    basic_sort& operator=(const basic_sort&) noexcept = default;
    basic_sort& operator=(basic_sort&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }
};

/// \\brief Make_basic_sort constructs a new term into a given address.
/// \\ \param t The reference into which the new basic_sort is constructed. 
template <class... ARGUMENTS>
inline void make_basic_sort(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_SortId(), args...);
}

// prototype declaration
std::string pp(const basic_sort& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const basic_sort& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(basic_sort& t1, basic_sort& t2)
{
  t1.swap(t2);
}
//--- end generated class basic_sort ---//

/// \brief list of basic sorts
typedef atermpp::term_list<basic_sort> basic_sort_list;
/// \brief vector of basic sorts
typedef std::vector<basic_sort> basic_sort_vector;

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_BASIC_SORT_H

