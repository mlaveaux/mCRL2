// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/exists.h
/// \brief The class exists.

#ifndef MCRL2_DATA_EXISTS_H
#define MCRL2_DATA_EXISTS_H

#include "mcrl2/data/abstraction.h"

namespace mcrl2
{

namespace data
{

/// \brief existential quantification.
///
class exists: public abstraction
{
  public:

    /// Constructor.
    ///
    /// \param[in] d An aterm.
    /// \pre d has the internal structure of an abstraction.
    /// \pre d is an existential quantification.
    explicit exists(const aterm& d)
      : abstraction(d)
    {
      assert(is_abstraction(d));
      assert(abstraction(d).binding_operator() == exists_binder());
    }

    /// Constructor.
    ///
    /// \param[in] variables A nonempty list of binding variables (objects of type variable).
    /// \param[in] body The body of the exists abstraction.
    /// \pre variables is not empty.
    template < typename Container >
    exists(const Container& variables,
           const data_expression& body,
           typename atermpp::enable_if_container< Container, variable >::type* = nullptr)
      : abstraction(exists_binder(), variables, body)
    {
      assert(!variables.empty());
    }

    /// Move semantics
    exists(const exists&) noexcept = default;
    exists(exists&&) noexcept = default;
    exists& operator=(const exists&) noexcept = default;
    exists& operator=(exists&&) noexcept = default;

}; // class exists

template <class... ARGUMENTS>
void make_exists(atermpp::aterm& result, ARGUMENTS... arguments)
{
  make_abstraction(result, exists_binder(), arguments...);
}


//--- start generated class exists ---//
// prototype declaration
std::string pp(const exists& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const exists& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(exists& t1, exists& t2)
{
  t1.swap(t2);
}
//--- end generated class exists ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_EXISTS_H

