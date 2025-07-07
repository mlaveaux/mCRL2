// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/aterm_int.h
/// \brief Term containing an integer.

#ifndef MCRL2_ATERMPP_ATERM_INT_H
#define MCRL2_ATERMPP_ATERM_INT_H

#include "mcrl2/atermpp/aterm.h"

#include <mcrl3_ffi.h>

namespace atermpp
{

/// \brief An integer term stores a single std::size_t value. It carries
/// 	       no arguments.
/// \details A special function symbol is used to identify integer terms
/// 			 internally.
class aterm_int : public aterm
{
public:
  /// \brief Default constructor.
  aterm_int() noexcept = default;

  /// \brief Constructs an integer term from a value.
  /// \param value The value of the new integer.
  explicit aterm_int(std::size_t value)
    : aterm( mcrl3::ffi::term_create_int(value)) 
  {}

  /// Construct a term from the internal representation of an integer term.
  explicit aterm_int(mcrl3::ffi::aterm_t term)
    : aterm(term)
  {
    assert(type_is_int() || !defined());
  }

  /// \brief Constructs an integer term from an aterm.
  explicit aterm_int(const aterm& t)
      : aterm(t)
  {
    assert(type_is_int() || !defined());
  }

  /// This class has user-declared copy constructor so declare default copy and move operators.
  aterm_int(const aterm_int& other) noexcept = default;
  aterm_int& operator=(const aterm_int& other) noexcept = default;
  aterm_int(aterm_int&& other) noexcept = default;
  aterm_int& operator=(aterm_int&& other) noexcept = default;

  /// \brief Provide the value stored in an aterm.
  /// \returns The value of the integer term.
  [[nodiscard]] std::size_t value() const noexcept { return mcrl3::ffi::term_get_int_value(m_term); }

  /// \brief Swaps two integer terms without changing the protection.
  /// \param t The term that is swapped with the current term.
  void swap(aterm_int& t) noexcept { aterm::swap(t); }
};

/// \brief Constructs an integer term from a value.
/// \param target The term into which the term is constructed.
inline void make_aterm_int(aterm_int& target, std::size_t value)
{
  target.unprotected_assign(aterm_int(mcrl3::ffi::term_create_int(value)));
}

} // namespace atermpp

namespace std
{

/// \brief Swaps two aterm_ints.
/// \details This operation is more efficient than exchanging terms by an assignment,
///          as swapping does not require to change the protection of terms.
///          In order to be used in the standard containers, the declaration must
///          be preceded by an empty template declaration.
/// \param t1 The first term
/// \param t2 The second term
template <>
inline void swap(atermpp::aterm_int& t1, atermpp::aterm_int& t2) noexcept
{
  t1.swap(t2);
}
} // namespace std

#endif // MCRL2_ATERMPP_ATERM_INT_H
