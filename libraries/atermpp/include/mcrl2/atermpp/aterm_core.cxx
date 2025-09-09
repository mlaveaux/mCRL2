// Author(s): Jan Friso Groote, Maurice Laveaux, Wieger Wesselink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_H
#define MCRL2_ATERMPP_ATERM_H
MCRL2_MODULE;

#include <algorithm>
#include <cassert>
#include <compare>
#include <sstream>
#include <type_traits>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/aterm_core_data.cxx"
  #include "mcrl2/atermpp/type_traits.cxx"
#else
  export module atermpp:aterm_core;

  import :function_symbol;
  import :type_traits;
  import :detail.aterm_core_data;
#endif


// Forward declaration
MCRL2_MODULE_EXPORT namespace atermpp::detail
{
  class thread_aterm_pool;
  
  detail::_aterm* address(const unprotected_aterm_core& t);
}

/// \brief The main namespace for the aterm++ library.
MCRL2_MODULE_EXPORT namespace atermpp
{

/// \brief An unprotected term does not change the reference count of the
///        shared term when it is copied or moved.
class unprotected_aterm_core
{
  friend detail::_aterm* detail::address(const unprotected_aterm_core& t);

protected:
  const detail::_aterm* m_term;

public:

  /// \brief Default constuctor.
  unprotected_aterm_core() noexcept
   : m_term(nullptr)
  {}

  /// \brief Constructor.
  /// \param term The term from which the new term is constructed.
  unprotected_aterm_core(const detail::_aterm* term) noexcept
   : m_term(term)
  {}

  /// \brief Dynamic check whether the term is an aterm.
  /// \return True iff this term is an term_appl.
  /// \details This function has constant complexity.
  ///          It is defined as !type_is_int() && !type_is_list().
  bool type_is_appl() const noexcept
  {
    return !type_is_int() && !type_is_list();
  }

  /// \brief Dynamic check whether the term is an aterm_int.
  /// \return True iff this term has internal structure of an aterm_int.
  /// \details This function has constant complexity.
  bool type_is_int() const noexcept
  {
    const function_symbol& f=m_term->function();
    return f == detail::g_as_int;
  }

  /// \brief Dynamic check whether the term is an aterm_list.
  /// \returns True iff this term has the structure of an term_list
  /// \details This function has constant complexity.
  bool type_is_list() const noexcept
  {
    const function_symbol& f=m_term->function();
    return f == detail::g_as_list || f == detail::g_as_empty_list;
  }

  /// \brief Comparison operator.
  /// \details Terms are stored in a maximally shared way. This
  ///         means that this equality operator can be calculated
  ///         in constant time.
  /// \return true iff t is equal to the current term.
  bool operator ==(const unprotected_aterm_core& t) const
  {
    return m_term == t.m_term;
  }

  /// \brief Comparison operator for two unprotected aterms.
  /// \details This operator is constant time. It compares
  ///         the addresses where terms are stored. That means
  ///         that the outcome of this operator is only stable
  ///         as long as aterms are not garbage collected.
  /// \param t A term to which the current term is compared.
  /// \return True iff the current term is smaller than the argument.
  std::weak_ordering operator<=>(const unprotected_aterm_core& t) const
  {
    return m_term<=>t.m_term;
  }

  /// \brief Returns true if this term is not equal to the term assigned by
  ///        the default constructor of aterms, aterm_appls and aterm_int.
  /// \details The default constructor of a term_list<T> is the empty list, on which
  ///          the operator defined yields true. This operation is more efficient
  ///          than comparing the current term with an aterm(), aterm_list() or an
  ///          aterm_int().
  /// \return A boolean indicating whether this term equals the default constructor.
  bool defined() const
  {
    return m_term != nullptr;
  }

  /// \brief Swaps this term with its argument.
  /// \details This operation is more efficient than exchanging terms by an assignment,
  ///          as swapping does not require to change the protection of terms.
  /// \param t The term with which this term is swapped.
  void swap(unprotected_aterm_core& t) noexcept
  {
    std::swap(m_term, t.m_term);
  }

  /// \brief Yields the function symbol in an aterm.
  /// \returns The function symbol of the term, which can also be an AS_EMPTY_LIST,
  ///          AS_INT and AS_LIST.
  /// \details This is for internal use only.
  const function_symbol& function() const
  {
    return m_term->function();
  }
};

/// \brief The aterm_core base class that provides protection of the underlying shared terms.
/// \details Terms are protected using one of the following two invariants:
///          (1) A term that can be accessed is a subterm of a term with a reference count
///              larger than 0 (when reference counting is used). Or, 
///          (2) A term that can be accessed if it is a subterm of a term that occurs at
///              an address which exist in the protection set of a process, or which sits
///              in an atermpp container, which automatically is a container protection set.
///              Furthermore, every address in a protection set contains a valid term.
///          During garbage collection or rehashing, this situation is stable in the sense
///          that all terms that are protected remain protected until the end of the 
///          garbage collection or rehashing phase by the same address or term. This means
///          that during garbage collection no terms can be deleted, for instance in an
///          assignment, or in a destruct. 
//               
class aterm_core : public unprotected_aterm_core
{
public:

  /// \brief Default constructor.
  aterm_core() noexcept;

  /// \brief Standard destructor.
  ~aterm_core() noexcept;

  /// \brief Constructor based on an internal term data structure. This is not for public use.
  /// \details Takes ownership of the passed underlying term.
  /// \param t A pointer to an internal aterm_core data structure.
  /// \todo Should be protected, but this cannot yet be done due to a problem
  ///       in the compiling rewriter.
  explicit aterm_core(const detail::_aterm *t) noexcept;

  /// \brief Copy constructor.
  /// \param other The aterm_core that is copied.
  /// \details  This class has a non-trivial destructor so explicitly define the copy and move operators.
  aterm_core(const aterm_core& other) noexcept;

  /// \brief Move constructor.
  /// \param other The aterm_core that is moved into the new term. This term may have changed after this operation.
  /// \details This operation does not employ increments and decrements of reference counts and is therefore more
  ///          efficient than the standard copy construct.
  aterm_core(aterm_core&& other) noexcept;

  /// \brief Assignment operator.
  /// \param other The aterm_core that will be assigned.
  /// \return A reference to the assigned term.
  aterm_core& operator=(const aterm_core& other) noexcept;

  /// \brief Assignment operator, to be used if busy and forbidden flags are explicitly available.
  //  \detail This can be used as an optimisation, because it avoids getting access to thread local variables,
  //          which is as it stands relatively expensive. The effect is equal to the assignment operator =. 
  /// \param other The aterm_core that will be assigned.
  aterm_core& assign(const aterm_core& other,
                detail::thread_aterm_pool& pool) noexcept;

  /// \brief Assignment operator, to be used when the busy flags do not need to be set.
  /// \details This is only safe in the parallel context when the busy flag is already
  ///          known to be set. This is also checked by an assert. This can be used for
  ///          instance in a lambda function that is passed in a make_.... function, as
  ///          this unprotected assign will only be called when a term is constructed. 
  /// \param other The aterm_core that will be assigned.
  template <bool CHECK_BUSY_FLAG=true>
  aterm_core& unprotected_assign(const aterm_core& other) noexcept;

  /// \brief Move assignment operator.
  /// \param other The aterm_core that will be assigned.
  /// \return A reference to the assigned term.
  aterm_core& operator=(aterm_core&& other) noexcept;
};

namespace detail
{
  /// \returns A pointer to the underlying aterm_core.
  inline _aterm* address(const unprotected_aterm_core& t)
  {
    return const_cast<_aterm*>(t.m_term);
  }
}

} // namespace atermpp

#endif // MCRL2_ATERMPP_ATERM_H
