// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/detail/term_traits_optimized.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_DETAIL_TERM_TRAITS_OPTIMIZED_H
#define MCRL2_PBES_DETAIL_TERM_TRAITS_OPTIMIZED_H

#include "mcrl2/pbes/pbes_expression.h"



namespace mcrl2::core {

/// \brief Contains type information for terms.
template <typename T>
struct term_traits_optimized
{
};

/// \brief Contains type information for pbes expressions.
template <>
struct term_traits_optimized<pbes_system::pbes_expression>: public core::term_traits<pbes_system::pbes_expression>
{

  using super = core::term_traits<pbes_system::pbes_expression>;

  static inline
  term_type not_(const term_type& x)
  {
    return pbes_system::optimized_not(x);
  } 

  static inline
  void make_not(term_type& result, const term_type& x)
  {
    pbes_system::make_optimized_not(result, x);
  }

  static inline
  term_type and_(const term_type& x, const term_type& y)
  {
    return pbes_system::optimized_and(x, y);
  }

  static inline
  void make_and(term_type& result, const term_type& x, const term_type& y)
  {
    pbes_system::make_optimized_and(result, x, y);
  }

  static inline
  term_type or_(const term_type& x, const term_type& y)
  {
    return pbes_system::optimized_or(x, y);
  } 

  static inline
  void make_or(term_type& result, const term_type& x, const term_type& y)
  {
    pbes_system::make_optimized_or(result, x, y);
  }

  static inline
  term_type imp(const term_type& x, const term_type& y)
  {
    return pbes_system::optimized_imp(x, y);
  } 

  static inline
  void make_imp(term_type& result, const term_type& x, const term_type& y)
  {
    pbes_system::make_optimized_imp(result, x, y);
  }

  static inline
  term_type forall(const variable_sequence_type& d, const term_type& x)
  {
    return pbes_system::optimized_forall(d, x);
  }

  static inline
  void make_forall(term_type& result, const variable_sequence_type& d, const term_type& x)
  {
    pbes_system::make_optimized_forall(result, d, x);
  }

  static inline
  term_type exists(const variable_sequence_type& d, const term_type& x)
  {
    return pbes_system::optimized_exists(d, x);
  }

  static inline
  void make_exists(term_type& result, const variable_sequence_type& d, const term_type& x)
  {
    pbes_system::make_optimized_exists(result, d, x);
  }

  template <typename FwdIt>
  static inline
  term_type join_or(FwdIt first, FwdIt last)
  {
    return utilities::detail::join(first, last, or_, false_());
  }

  template <typename FwdIt>
  static inline
  term_type join_and(FwdIt first, FwdIt last)
  {
    return utilities::detail::join(first, last, and_, true_());
  }
};

} // namespace mcrl2::core



#endif // MCRL2_PBES_DETAIL_TERM_TRAITS_OPTIMIZED_H
