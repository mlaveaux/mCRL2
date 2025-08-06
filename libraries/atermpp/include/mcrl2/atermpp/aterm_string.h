// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_STRING_H
#define MCRL2_ATERMPP_ATERM_STRING_H

#include <cassert>
#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm.cxx"
#else
  import atermpp;
#endif

namespace atermpp
{

/// \brief Term containing a string.
class aterm_string: public aterm
{
public:
  /// \brief Default constructor.
  aterm_string() = default;

  /// \brief Constructor.
  /// \param t A term without arguments of type appl. The string is given by the function symbol.
  explicit aterm_string(const aterm& t)
    : aterm(t)
  {
    assert(size() == 0);
  }

  /// \brief Constructor.
  /// \param f A function symbol with arity 0.
  aterm_string(const function_symbol& f)
    : aterm(f)
  {
    assert(f.arity() == 0);
    assert(size() == 0);
  }

  /// \brief Constructor that allows construction from a string.
  /// \param s A string.
  aterm_string(const std::string& s)
    : aterm(function_symbol(s,0))
  {
    assert(size() == 0);
  }

  /// Copy constructor.
  /// \param t An aterm_string.
  aterm_string(const aterm_string& t) noexcept = default;

  /// Move constructor.
  /// \param t An aterm_string.
  aterm_string(aterm_string&& t) noexcept = default;

  /// Assignment operator.
  /// \param t An aterm_string.
  aterm_string& operator=(const aterm_string& t) noexcept = default;

  /// Move assignment operator.
  /// \param t An aterm_string.
  aterm_string& operator=(aterm_string&& t) noexcept = default;

  /// \brief Conversion operator
  /// \return The term converted to string
  operator const std::string&() const
  {
    return function().name();
  }
};

/// \brief Returns the empty aterm_string
inline
const aterm_string& empty_string()
{
  static atermpp::aterm_string t("");
  return t;
}

} // namespace atermpp


namespace std
{

/// \brief Swaps two aterm_strings.
/// \details This operation is more efficient than exchanging terms by an assignment,
///          as swapping does not require to change the protection of terms.
///          In order to be used in the standard containers, the declaration must
///          be preceded by an empty template declaration.
/// \param t1 The first term
/// \param t2 The second term

template <>
inline void swap(atermpp::aterm_string& t1, atermpp::aterm_string& t2) noexcept
{
  t1.swap(t2);
}
} // namespace std

#endif // MCRL2_ATERMPP_ATERM_STRING_H
