// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/untyped_possible_sorts.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_UNTYPED_POSSIBLE_SORTS_H
#define MCRL2_DATA_UNTYPED_POSSIBLE_SORTS_H

#include "mcrl2/data/sort_expression.h"

namespace mcrl2 {

namespace data {

//--- start generated class untyped_possible_sorts ---//
/// \\brief Multiple possible sorts
class untyped_possible_sorts: public sort_expression
{
  public:
    /// \\brief Default constructor X3.
    untyped_possible_sorts()
      : sort_expression(core::detail::default_values::UntypedSortsPossible)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit untyped_possible_sorts(const atermpp::aterm& term)
      : sort_expression(term)
    {
      assert(core::detail::check_term_UntypedSortsPossible(*this));
    }

    /// \\brief Constructor Z14.
    explicit untyped_possible_sorts(const sort_expression_list& sorts)
      : sort_expression(atermpp::aterm(core::detail::function_symbol_UntypedSortsPossible(), sorts))
    {}

    /// \\brief Constructor Z2.
    template <typename Container>
    untyped_possible_sorts(const Container& sorts, typename atermpp::enable_if_container<Container, sort_expression>::type* = nullptr)
      : sort_expression(atermpp::aterm(core::detail::function_symbol_UntypedSortsPossible(), sort_expression_list(sorts.begin(), sorts.end())))
    {}

    /// Move semantics
    untyped_possible_sorts(const untyped_possible_sorts&) noexcept = default;
    untyped_possible_sorts(untyped_possible_sorts&&) noexcept = default;
    untyped_possible_sorts& operator=(const untyped_possible_sorts&) noexcept = default;
    untyped_possible_sorts& operator=(untyped_possible_sorts&&) noexcept = default;

    const sort_expression_list& sorts() const
    {
      return atermpp::down_cast<sort_expression_list>((*this)[0]);
    }
};

/// \\brief Make_untyped_possible_sorts constructs a new term into a given address.
/// \\ \param t The reference into which the new untyped_possible_sorts is constructed. 
template <class... ARGUMENTS>
inline void make_untyped_possible_sorts(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_UntypedSortsPossible(), args...);
}

// prototype declaration
std::string pp(const untyped_possible_sorts& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const untyped_possible_sorts& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(untyped_possible_sorts& t1, untyped_possible_sorts& t2)
{
  t1.swap(t2);
}
//--- end generated class untyped_possible_sorts ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_UNTYPED_POSSIBLE_SORTS_H
