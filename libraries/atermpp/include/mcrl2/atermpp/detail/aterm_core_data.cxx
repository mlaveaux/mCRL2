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
MCRL2_MODULE;

#include <limits>
#include <type_traits>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/aterm_configuration.cxx"
  #include "mcrl2/utilities/type_traits.cxx"
  #include "mcrl2/utilities/noncopyable.cxx"
  #include "mcrl2/utilities/shared_reference.cxx"
#else
  export module atermpp:detail.aterm_core_data;

  import :detail.aterm_configuration;
  import :function_symbol;

  import utilities;
#endif

// Forward declaration.
MCRL2_MODULE_EXPORT namespace atermpp
{
class unprotected_aterm_core;
class aterm_core;
}

MCRL2_MODULE_EXPORT namespace atermpp
{

namespace detail
{

/// Can be used to check whether all elements in the parameter pack are derived from aterms.
template<typename ...Terms>
using are_terms = mcrl2::utilities::forall<std::is_convertible<Terms, unprotected_aterm_core>...>;

/// Check whether all arguments of a parameter pack are terms, constant functions yielding a term
/// or a function putting a term in a result parameter. 

template <typename Term>
using is_term_or_function = std::disjunction<std::is_convertible<Term, unprotected_aterm_core>,
                            std::disjunction<mcrl2::utilities::is_applicable<Term, unprotected_aterm_core&>,
                                             mcrl2::utilities::is_constant_function_yielding<Term, unprotected_aterm_core> > >;

template<typename ...Terms>
using are_terms_or_functions = mcrl2::utilities::forall<is_term_or_function<Terms>...>;

/// \brief This is the class to which an aterm points.
class _aterm : public mcrl2::utilities::noncopyable
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
  void mark() const
  {
    m_function_symbol.m_function_symbol.tag();
  }

  /// \brief Remove the mark from a term.
  void unmark() const
  {
    m_function_symbol.m_function_symbol.untag();
  }

  /// \brief Check if the term is already marked.
  bool is_marked() const
  {
    return m_function_symbol.m_function_symbol.tagged();
  }

private:
  function_symbol m_function_symbol;
};

} // namespace detail
} // namespace atermpp

#endif /* MCRL2_ATERMPP_DETAIL_ATERM_H */
