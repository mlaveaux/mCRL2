// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_DETAIL_ATERM_H
#define MCRL2_ATERMPP_DETAIL_ATERM_H

#include "mcrl2/atermpp/function_symbol.h"
#include "mcrl2/utilities/type_traits.h"

#include <limits>

namespace atermpp
{

// Forward declaration.
class unprotected_aterm;
class aterm;

constexpr static std::size_t MarkedReferenceCount = std::numeric_limits<std::size_t>::max();

namespace detail
{

/// \brief Can be used to check whether all elements in the parameter pack are derived from aterms.
template<typename ...Terms>
using are_terms = mcrl2::utilities::forall<std::is_convertible<Terms, unprotected_aterm>...>;

/// \brief This is the class to which an aterm points. Each _aterm consists
///        of a function symbol and a reference count used for garbage
///        collection.
class _aterm : public
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  mcrl2::utilities::shared_reference_counted<_aterm, GlobalThreadSafe>,
#endif
  mcrl2::utilities::noncopyable
{
public:
  /// \brief Create a term from a function symbol.
  _aterm(const function_symbol& symbol) :
    m_function_symbol(symbol)
  {}

  const function_symbol& function() const noexcept
  {
    return m_function_symbol;
  }

  /// \brief Mark this term to be garbage collected.
  /// \details Changes the reference count, so only apply whenever !is_reachable().
  void mark() const
  {
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
    assert(!is_reachable());
    m_reference_count = MarkedReferenceCount;
    increment_reference_count_changes();
#else
    m_function_symbol.m_function_symbol.tag();
#endif
  }

  /// \brief Remove the mark from a term.
  /// \details Changes the reference count, so only apply whenever it was marked.
  void reset() const
  {    
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
    assert(is_marked());
    m_reference_count = 0;
    increment_reference_count_changes();
#else
    m_function_symbol.m_function_symbol.reset();
#endif
  }

  /// \brief During garbage collection a term will be marked whenever it occurs as the argument
  ///        of a reachable term, but is not protected itself.
  /// \returns True whenever this term has been marked.
  bool is_marked() const noexcept
  {
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
    return m_reference_count == MarkedReferenceCount;
#else
    return m_function_symbol.m_function_symbol.tagged();
#endif
  }

  /// \brief A term is reachable in the garbage collection graph if it is protected or whenever it occurs
  ///        as an argument of a reachable term. The latter will be ensured in the marking phase of garbage collection.
  /// \returns True whenever the term is reachable ie either marked or protected.
  bool is_reachable() const noexcept
  {   
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
    return m_reference_count > 0;
#else
    return m_function_symbol.m_function_symbol.tagged();
#endif
  }

  /// \brief Access to the function symbol pointer when the term is tagged.
  const detail::_function_symbol* get_symbol() const
  {
    return pointer(m_function_symbol.m_function_symbol.get());
  }

private:
  function_symbol m_function_symbol;
};

inline _aterm* address(const unprotected_aterm& t);

} // namespace detail
} // namespace atermpp

#endif /* MCRL2_ATERMPP_DETAIL_ATERM_H */
