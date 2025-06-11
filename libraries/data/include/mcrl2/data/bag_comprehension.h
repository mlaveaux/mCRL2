// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/bag_comprehension.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_BAG_COMPREHENSION_H
#define MCRL2_DATA_BAG_COMPREHENSION_H

#include "mcrl2/data/abstraction.h"

namespace mcrl2 {

namespace data {

/// \brief universal quantification.
///
class bag_comprehension: public abstraction
{
  public:

    /// Constructor.
    ///
    /// \param[in] d An aterm.
    /// \pre d has the interal structure of an abstraction.
    /// \pre d is a universal quantification.
    explicit bag_comprehension(const aterm& d)
      : abstraction(d)
    {
      assert(is_abstraction(d));
      assert(static_cast<abstraction>(d).binding_operator() == bag_comprehension_binder());
    }

    /// Constructor.
    ///
    /// \param[in] variables A nonempty list of binding variables (objects of type variable).
    /// \param[in] body The body of the bag_comprehension abstraction.
    /// \pre variables is not empty.
    template < typename Container >
    bag_comprehension(const Container& variables,
           const data_expression& body,
           typename atermpp::enable_if_container< Container, variable >::type* = nullptr)
      : abstraction(bag_comprehension_binder(), variables, body)
    {
      assert(!variables.empty());
    }

    /// Move semantics
    bag_comprehension(const bag_comprehension&) noexcept = default;
    bag_comprehension(bag_comprehension&&) noexcept = default;
    bag_comprehension& operator=(const bag_comprehension&) noexcept = default;
    bag_comprehension& operator=(bag_comprehension&&) noexcept = default;

}; // class bag_comprehension

template <class... ARGUMENTS>
void make_bag_comprehension(atermpp::aterm& result, ARGUMENTS... arguments)
{
  make_abstraction(result, bag_comprehension_binder(), arguments...);
}


//--- start generated class bag_comprehension ---//
// prototype declaration
std::string pp(const bag_comprehension& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const bag_comprehension& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(bag_comprehension& t1, bag_comprehension& t2)
{
  t1.swap(t2);
}
//--- end generated class bag_comprehension ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_BAG_COMPREHENSION_H
