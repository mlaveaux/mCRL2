// Author(s): Jan Friso Groote. Based on pres/detail/term_traits_optimized.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/detail/term_traits_optimized.h
/// \brief add your file description here.

#ifndef MCRL2_PRES_DETAIL_TERM_TRAITS_OPTIMIZED_H
#define MCRL2_PRES_DETAIL_TERM_TRAITS_OPTIMIZED_H

#include "mcrl2/pres/pres_expression.h"



namespace mcrl2::core {

/// \brief Contains type information for terms.
template <typename T>
struct term_traits_optimized
{
};

/// \brief Contains type information for pres expressions.
template <>
struct term_traits_optimized<pres_system::pres_expression>: public core::term_traits<pres_system::pres_expression>
{
  using super = core::term_traits<pres_system::pres_expression>;

  static inline
  term_type minus(const term_type& x)
  {
    return pres_system::optimized_minus(x);
  } 

  static inline
  void make_minus(term_type& result, const term_type& x)
  {
    pres_system::make_optimized_minus(result, x);
  }

  static inline
  term_type and_(const term_type& x, const term_type& y)
  {
    return pres_system::optimized_and(x, y);
  }

  static inline
  void make_and(term_type& result, const term_type& x, const term_type& y)
  {
    pres_system::make_optimized_and(result, x, y);
  }

  static inline
  term_type or_(const term_type& x, const term_type& y)
  {
    return pres_system::optimized_or(x, y);
  } 

  static inline
  void make_or(term_type& result, const term_type& x, const term_type& y)
  {
    pres_system::make_optimized_or(result, x, y);
  }

  static inline
  term_type plus(const term_type& x, const term_type& y)
  {
    return pres_system::optimized_plus(x, y);
  } 

  static inline
  void make_plus(term_type& result, const term_type& x, const term_type& y)
  {
    pres_system::make_optimized_plus(result, x, y);
  }

  static inline
  void make_imp(term_type& result, const term_type& x, const term_type& y)
  {
    pres_system::make_optimized_minus(result, x);
    pres_system::make_optimized_or(result, result, y);
  } 

  static inline
  term_type imp(const term_type& x, const term_type& y)
  {
    pres_system::pres_expression result;
    make_imp(result, x, y);
    return result;
  } 

  static inline
  term_type infimum(const variable_sequence_type& d, const term_type& x)
  {
    return pres_system::optimized_infimum(d, x, true);
  }

  static inline
  void make_infimum(term_type& result, const variable_sequence_type& d, const term_type& x)
  {
    pres_system::make_optimized_infimum(result, d, x, true);
  }

  static inline
  term_type supremum(const variable_sequence_type& d, const term_type& x)
  {
    return pres_system::optimized_supremum(d, x, true);
  }

  static inline
  void make_supremum(term_type& result, const variable_sequence_type& d, const term_type& x)
  {
    pres_system::make_optimized_supremum(result, d, x, true);
  }

  static inline
  term_type sum(const variable_sequence_type& d, const term_type& x)
  {
    return pres_system::optimized_sum(d, x, true);
  }

  static inline
  void make_sum(term_type& result, const variable_sequence_type& d, const term_type& x)
  {
    pres_system::make_optimized_sum(result, d, x, true);
  }

  static inline
  term_type condsm(const term_type& x, const term_type& y, const term_type& z)
  {
    return pres_system::optimized_condsm(x, y, z);
  } 

  static inline
  void make_condsm(term_type& result, const term_type& x, const term_type& y, const term_type& z)
  {
    pres_system::make_optimized_condsm(result, x, y, z);
  }

  static inline
  term_type condeq(const term_type& x, const term_type& y, const term_type& z)
  {
    return pres_system::optimized_condeq(x, y, z);
  } 

  static inline
  void make_condeq(term_type& result, const term_type& x, const term_type& y, const term_type& z)
  {
    pres_system::make_optimized_condeq(result, x, y, z);
  }

  static inline
  term_type eqinf(const term_type& x)
  {
    return pres_system::optimized_eqinf(x);
  } 

  static inline
  void make_eqinf(term_type& result, const term_type& x)
  {
    pres_system::make_optimized_eqinf(result, x);
  }

  static inline
  term_type eqninf(const term_type& x)
  {
    return pres_system::optimized_eqninf(x);
  } 

  static inline
  void make_eqninf(term_type& result, const term_type& x)
  {
    pres_system::make_optimized_eqninf(result, x);
  }

  static inline
  term_type const_multiply(const data::data_expression& d, const term_type& x)
  {
    return pres_system::optimized_const_multiply(d, x);
  } 

  static inline
  void make_const_multiply(term_type& result, const data::data_expression& d, const term_type& x)
  {
    pres_system::make_optimized_const_multiply(result, d, x);
  }

  static inline
  term_type const_multiply_alt(const data::data_expression& d, const term_type& x)
  {
    return pres_system::optimized_const_multiply_alt(d, x);
  } 

  static inline
  void make_const_multiply_alt(term_type& result, const data::data_expression& d, const term_type& x)
  {
    pres_system::make_optimized_const_multiply_alt(result, d, x);
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



#endif // MCRL2_PRES_DETAIL_TERM_TRAITS_OPTIMIZED_H
