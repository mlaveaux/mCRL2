// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/untyped_sort.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_UNTYPED_SORT_H
#define MCRL2_DATA_UNTYPED_SORT_H

#include "mcrl2/data/sort_expression.h"

namespace mcrl2 {

namespace data {

//--- start generated class untyped_sort ---//
/// \\brief Unknown sort expression
class untyped_sort: public sort_expression
{
  public:
    /// \\brief Default constructor X3.
    untyped_sort()
      : sort_expression(core::detail::default_values::UntypedSortUnknown)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit untyped_sort(const atermpp::aterm& term)
      : sort_expression(term)
    {
      assert(core::detail::check_term_UntypedSortUnknown(*this));
    }

    /// Move semantics
    untyped_sort(const untyped_sort&) noexcept = default;
    untyped_sort(untyped_sort&&) noexcept = default;
    untyped_sort& operator=(const untyped_sort&) noexcept = default;
    untyped_sort& operator=(untyped_sort&&) noexcept = default;
};

// prototype declaration
std::string pp(const untyped_sort& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const untyped_sort& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(untyped_sort& t1, untyped_sort& t2)
{
  t1.swap(t2);
}
//--- end generated class untyped_sort ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_UNTYPED_SORT_H
