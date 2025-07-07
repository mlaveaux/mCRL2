// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/aterm.h
/// \brief The term_appl class represents function application.

#ifndef MCRL2_ATERMPP_UNPROTECTED_ATERM_H
#define MCRL2_ATERMPP_UNPROTECTED_ATERM_H

#include "mcrl2/atermpp/detail/aterm_appl_iterator.h"
#include "mcrl2/atermpp/function_symbol.h"

#include <mcrl3_ffi.h>

#include <cassert>
#include <limits>
#include <utility>

namespace atermpp
{

/// Forward declaration of the unprotected_aterm class.
class unprotected_aterm;

/// Forward declaration of the aterm class.
class aterm;

namespace detail
{
/// \brief Returns the address of the unprotected term.
mcrl3::ffi::unprotected_aterm_t address(const unprotected_aterm& t) noexcept;
} // namespace detail

/// \brief An unprotected term does not change the reference count of the
///        shared term when it is copied or moved.
class unprotected_aterm
{
  friend mcrl3::ffi::unprotected_aterm_t detail::address(const unprotected_aterm& t) noexcept;

public:
  /// Iterator used to iterate through an term_appl.
  using iterator = term_appl_iterator<aterm>;

  /// Const iterator used to iterate through an term_appl.
  using const_iterator = term_appl_iterator<aterm>;

  /// \brief Default constuctor.
  unprotected_aterm() noexcept = default;

  /// \brief Constructor.
  /// \param term The term from which the new term is constructed.
  explicit unprotected_aterm(mcrl3::ffi::unprotected_aterm_t term) noexcept
    : m_term(term)
  {}

  /// \brief Dynamic check whether the term is an aterm.
  /// \return True iff this term is an term_appl.
  /// \details This function has constant complexity.
  ///          It is defined as !type_is_int() && !type_is_list().
  [[nodiscard]] bool type_is_appl() const noexcept { return !type_is_int() && !type_is_list(); }

  /// \brief Dynamic check whether the term is an aterm_int.
  /// \return True iff this term has internal structure of an aterm_int.
  /// \details This function has constant complexity.
  [[nodiscard]] bool type_is_int() const noexcept { return mcrl3::ffi::term_is_int(m_term); }

  /// \brief Dynamic check whether the term is an aterm_list.
  /// \returns True iff this term has the structure of an term_list
  /// \details This function has constant complexity.
  [[nodiscard]] bool type_is_list() const noexcept { return mcrl3::ffi::term_is_list(m_term); }

  /// \brief Dynamic check whether the term is an empty aterm_list.
  /// \returns True iff this term has the structure of an term_list
  /// \details This function has constant complexity.
  [[nodiscard]] bool type_is_empty_list() const noexcept { return mcrl3::ffi::term_is_empty_list(m_term); }

  /// \brief Comparison operator.
  /// \details Terms are stored in a maximally shared way. This
  ///         means that this equality operator can be calculated
  ///         in constant time.
  /// \return true iff t is equal to the current term.
  bool operator==(const unprotected_aterm& t) const { return m_term.ptr == t.m_term.ptr; }

  /// \brief Comparison operator for two unprotected aterms.
  /// \details This operator is constant time. It compares
  ///         the addresses where terms are stored. That means
  ///         that the outcome of this operator is only stable
  ///         as long as aterms are not garbage collected.
  /// \param t A term to which the current term is compared.
  /// \return True iff the current term is smaller than the argument.
  std::weak_ordering operator<=>(const unprotected_aterm& t) const { return m_term.ptr <=> t.m_term.ptr; }

  /// \brief Returns true if this term is not equal to the term assigned by
  ///        the default constructor of aterms, aterm_appls and aterm_int.
  /// \details The default constructor of a term_list<T> is the empty list, on which
  ///          the operator defined yields true. This operation is more efficient
  ///          than comparing the current term with an aterm(), aterm_list() or an
  ///          aterm_int().
  /// \return A boolean indicating whether this term equals the default constructor.
  [[nodiscard]] bool defined() const { return m_term.ptr != nullptr; }

  /// \brief Swaps this term with its argument.
  /// \details This operation is more efficient than exchanging terms by an assignment,
  ///          as swapping does not require to change the protection of terms.
  /// \param t The term with which this term is swapped.
  void swap(unprotected_aterm& t) noexcept { std::swap(m_term, t.m_term); }

  /// \brief Yields the function symbol in an aterm.
  /// \returns The function symbol of the term, which can also be an AS_EMPTY_LIST,
  ///          AS_INT and AS_LIST.
  /// \details This is for internal use only.
  [[nodiscard]] const function_symbol& function() const
  {
    // return reinterpret_cast<const function_symbol&>(mcrl3::ffi::term_get_function_symbol(m_term).ptr);
  }
  
  /// \brief Returns the number of arguments of this term.
  /// \return The number of arguments of this term.
  [[nodiscard]] std::size_t size() const
  {
    mcrl3::ffi::function_symbol_t func = mcrl3::ffi::term_get_function_symbol(m_term);
    return mcrl3::ffi::function_symbol_get_arity(func);
  }

  /// \brief Returns the largest possible number of arguments.
  /// \return The largest possible number of arguments.
  [[nodiscard]] constexpr std::size_t max_size() const { return std::numeric_limits<std::size_t>::max(); }
  
  /// \brief Returns the i-th argument.
  /// \param i A positive integer.
  /// \return The argument with the given index.
  const unprotected_aterm operator[](const std::size_t i) const
  {
    assert(i < size()); // Check the bounds.
    return unprotected_aterm(mcrl3::ffi::term_get_argument(m_term, i));
  }
  
  /// \brief Returns an iterator pointing to the first argument.
  /// \return An iterator pointing to the first argument.
  [[nodiscard]] const_iterator begin() const
  {
    return const_iterator(static_cast<const aterm*>(operator[](0).m_term.ptr));
  }

  /// \brief Returns a const_iterator pointing past the last argument.
  /// \return A const_iterator pointing past the last argument.
  [[nodiscard]] const_iterator end() const
  {
    return const_iterator(static_cast<const aterm*>(operator[](size()).m_term.ptr));
  }

  /// Marks the term as used during garbage collection.
  void mark() const
  {
    if (defined())
    {
      mcrl3::ffi::term_mark(m_term);
    }
  }

protected:
  mcrl3::ffi::unprotected_aterm_t m_term{nullptr};
};

namespace detail
{
mcrl3::ffi::unprotected_aterm_t address(const unprotected_aterm& t) noexcept
{
  return t.m_term;
}

} // namespace detail

} // namespace atermpp

namespace std
{
/// \brief Standard hash function.
template <>
struct hash<atermpp::unprotected_aterm>
{
  std::size_t operator()(const atermpp::unprotected_aterm& t) const
  {
    return reinterpret_cast<std::size_t>(atermpp::detail::address(t).ptr) >> 4;
  }
};
} // namespace std

#endif // MCRL2_ATERMPP_UNPROTECTED_ATERM_H