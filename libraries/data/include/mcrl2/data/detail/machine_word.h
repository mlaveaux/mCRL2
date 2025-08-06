// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/machine_number.h
/// \brief The class machine_number, which is a subclass of data_expression.

#ifndef MCRL2_DATA_DETAIL_MACHINE_WORD_H
#define MCRL2_DATA_DETAIL_MACHINE_WORD_H

#include <cassert>
#include <cmath>
#include <limits>

#include "mcrl2/data/bool.h"
#include "mcrl2/data/machine_number.h"
#include "mcrl2/data/machine_word.h"

namespace mcrl2::data::sort_machine_word
{

namespace detail
{

/* Below standard operations on machine numbers are implemented.  */

inline const machine_number& zero_word()
{
  static const machine_number zero=machine_number(0);
  return zero;
}

inline const machine_number& one_word()
{
  static const machine_number one=machine_number(1);
  return one;
}

inline const machine_number& two_word()
{
  static const machine_number two=machine_number(2);
  return two;
}

inline const machine_number& three_word()
{
  static const machine_number three=machine_number(3);
  return three;
}

inline const machine_number& four_word()
{
  static const machine_number four=machine_number(4);
  return four;
}

inline const machine_number& max_word()
{
  static const machine_number max=machine_number(std::numeric_limits<std::size_t>::max());
  return max;
}

inline bool equals_zero_word(const std::size_t n)
{
  return n==0;
}

inline bool equals_one_word(const std::size_t n)
{
  return n==1;
}

inline bool equals_max_word(const std::size_t n)
{
  return n==std::numeric_limits<std::size_t>::max();
}

inline std::size_t succ_word(const std::size_t n)
{
  return 1+n;
}

inline bool equal_word(const std::size_t n1, const std::size_t n2)
{
  return n1==n2;
}

inline bool less_word(const std::size_t n1, const std::size_t n2)
{
  return n1<n2;
}

inline bool less_equal_word(const std::size_t n1, const std::size_t n2)
{
  return n1<=n2;
}

inline std::size_t add_word(const std::size_t n1, const std::size_t n2)
{
  return n1+n2;
}

inline std::size_t add_with_carry_word(const std::size_t n1, const std::size_t n2)
{
  return n1+n2+1;
}

inline bool add_overflow_word(const std::size_t n1, const std::size_t n2)
{
  if (n1+n2<n1)
  {
    return true;  // In this case there is an overflow.
  }
 return false; // No overflow. 
}

inline bool add_with_carry_overflow_word(const std::size_t n1, const std::size_t n2)
{
  if (n1+n2+1<n1)
  {
    return true;  // In this case there is an overflow.
  }
  return false; // No overflow. 
}

inline std::size_t times_word(const std::size_t n1, const std::size_t n2)
{
  return n1*n2;
}

inline std::size_t times_with_carry_word(const std::size_t n1, const std::size_t n2, const std::size_t n3)
{
  return n1*n2+n3;
}

std::size_t times_overflow_word(const std::size_t n1, const std::size_t n2);

std::size_t times_with_carry_overflow_word(const std::size_t n1, const std::size_t n2, std::size_t n3);

/// Performs minus modulo the largest std::size_t.
inline std::size_t minus_word(const std::size_t n1, const std::size_t n2)
{
  // In C++ underflow of unsigned integer operators is defined behaviour, and std::size_t is unsigned.
  return n1-n2;
}

inline std::size_t monus_word(const std::size_t n1, const std::size_t n2)
{
  if (n1>n2) 
  { 
    return n1-n2;
  }
  return 0;
}

inline std::size_t div_word(const std::size_t n1, const std::size_t n2)
{
  return n1/n2;
}

inline std::size_t mod_word(const std::size_t n1, const std::size_t n2)
{
  return n1 % n2;
}

std::size_t div_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3);

std::size_t mod_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3);

std::size_t div_double_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4);

std::size_t mod_double_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4);

std::size_t div_triple_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4, const std::size_t n5);

inline std::size_t sqrt_word(const std::size_t n)
{
  return sqrt(n);
}

std::size_t sqrt_doubleword(const std::size_t n1, const std::size_t n2);

std::size_t sqrt_tripleword(const std::size_t n1, const std::size_t n2, const std::size_t n3);

std::size_t sqrt_tripleword_overflow(const std::size_t n1, const std::size_t n2, const std::size_t n3);

std::size_t sqrt_quadrupleword(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4);

std::size_t sqrt_quadrupleword_overflow(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4);

inline std::size_t pred_word(const std::size_t n)
{
  return n-1;
}

inline void rightmost_bit(data_expression& result, const std::size_t n)
{
  if ((n & 1) == 1)
  {
    result=sort_bool::true_();
    return;
  }
  result=sort_bool::false_();
}

inline std::size_t shift_right(const data_expression& e1, const std::size_t n)
{
  std::size_t m = n>>1;
  if (e1==sort_bool::true_())
  {
    return m | static_cast<std::size_t>(1)<<((8*sizeof(std::size_t))-1);
  }
  return m;
}

} // end namespace detail 

/// \brief The machine number representing 0.
/// \param result The data_expression into which the created object is stored.
/// \return The machine number 0.
inline void zero_word_manual_implementation(data_expression& result)
{
  result=detail::zero_word();
}

/// \brief The machine number representing 1.
/// \param result The data_expression into which the created object is stored.
/// \return The machine number 1.
inline void one_word_manual_implementation(data_expression& result)
{
  result=detail::one_word();
}

/// \brief The machine number representing 2.
/// \param result The data_expression into which the created object is stored.
/// \return The machine number 2.
inline void two_word_manual_implementation(data_expression& result)
{
  result=detail::two_word();
}

/// \brief The machine number representing 3.
/// \param result The data_expression into which the created object is stored.
/// \return The machine number 3.
inline void three_word_manual_implementation(data_expression& result)
{
  result=detail::three_word();
}

/// \brief The machine number representing 4.
/// \param result The data_expression into which the created object is stored.
/// \return The machine number 4.
inline void four_word_manual_implementation(data_expression& result)
{
  result=detail::four_word();
}

/// \brief The largest representable machine number.
/// \param result The data_expression into which the created object is stored.
/// \return The largest number a machine word can hold. 
inline void max_word_manual_implementation(data_expression& result)
{
  result=detail::max_word();
}

/// \brief Checks whether the argument is equal to 0.
/// \param e 
/// \param result The data_expression into which the created object is stored.
/// \return True if e equals 0, otherwise false.
inline void equals_zero_word_manual_implementation(data_expression& result, const data_expression& e)
{
  assert(is_machine_number(e));
  if (detail::equals_zero_word(atermpp::down_cast<machine_number>(e).value()))
  {
    result=sort_bool::true_();
    return;
  }
  result=sort_bool::false_();
  return;
}

/// \brief Checks whether the argument is not equal to 0.
/// \param e 
/// \param result The data_expression into which the created object is stored.
/// \return True if e equals 0, otherwise false.
inline void not_equals_zero_word_manual_implementation(data_expression& result, const data_expression& e)
{
  assert(is_machine_number(e));
  if (detail::equals_zero_word(atermpp::down_cast<machine_number>(e).value()))
  {
    result=sort_bool::false_();
    return;
  }
  result=sort_bool::true_();
  return;
}

/// \brief Checks whether the argument is equal to 1.
/// \param e 
/// \param result The data_expression into which the created object is stored.
/// \return True if e equals 0, otherwise false.
inline void equals_one_word_manual_implementation(data_expression& result, const data_expression& e)
{
  assert(is_machine_number(e));
  if (detail::equals_one_word(atermpp::down_cast<machine_number>(e).value()))
  {
    result=sort_bool::true_();
    return;
  }
  result=sort_bool::false_();
  return;
}

/// \brief Checks whether the argument is equal to the largest 64 bit number.
/// \param e 
/// \param result The data_expression into which the created object is stored.
/// \return True if e equals the largest 64 number, otherwise false.
inline void equals_max_word_manual_implementation(data_expression& result, const data_expression& e)
{
  assert(is_machine_number(e));
  if (detail::equals_max_word(atermpp::down_cast<machine_number>(e).value()))
  {
    result=sort_bool::true_();
    return;
  }
  result=sort_bool::false_();
  return;
}

/// \brief The successor function on a machine numbers, that wraps around.
/// \param result The data_expression into which the created object is stored.
/// \param e 
/// \return e+1, or zero if n is the maximum number.
inline void succ_word_manual_implementation(data_expression& result, const data_expression& e)
{
  assert(is_machine_number(e));
  make_machine_number(result, detail::succ_word(atermpp::down_cast<machine_number>(e).value()));
}

/// \brief The equality function on two machine words. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1==e2.
inline void equal_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  const bool b=detail::equal_word(
                            atermpp::down_cast<machine_number>(e1).value(),
                            atermpp::down_cast<machine_number>(e2).value());
  result=(b? sort_bool::true_(): sort_bool::false_());
}

/// \brief The non equality function on two machine words. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1!=e2.
inline void not_equal_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  const bool b=detail::equal_word(
                          atermpp::down_cast<machine_number>(e1).value(),
                          atermpp::down_cast<machine_number>(e2).value());
  result=(b? sort_bool::false_(): sort_bool::true_());
}

/// \brief The less than function on two machine words. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1<e2.
inline void less_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  const bool b=detail::less_word(atermpp::down_cast<machine_number>(e1).value(),atermpp::down_cast<machine_number>(e2).value());
  result=(b? sort_bool::true_(): sort_bool::false_());
}

/// \brief The less than or equal function on two machine words. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1<=e2.
inline void less_equal_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  const bool b=detail::less_equal_word(atermpp::down_cast<machine_number>(e1).value(),atermpp::down_cast<machine_number>(e2).value());
  result=(b? sort_bool::true_(): sort_bool::false_());
}

/// \brief The greater than function on two machine words. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1>e2.
inline void greater_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  const bool b=detail::less_word(atermpp::down_cast<machine_number>(e2).value(),atermpp::down_cast<machine_number>(e1).value());
  result=(b? sort_bool::true_(): sort_bool::false_());
}

/// \brief The greater than or equal function on two machine words. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1>=e2.
inline void greater_equal_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  const bool b=detail::less_equal_word(atermpp::down_cast<machine_number>(e2).value(),atermpp::down_cast<machine_number>(e1).value());
  result=(b? sort_bool::true_(): sort_bool::false_());
}

/// \brief The result of adding two words modulo the maximal representable machine word plus 1.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1+e2 modulo the machine word. 
inline void add_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::add_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief The result of adding two words plus 1 modulo the maximal representable machine word plus 1.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1+e2+1 modulo the machine word. 
inline void add_with_carry_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::add_with_carry_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief An indication whether an overflow occurs when e1 and e2 are added. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return If e1+e2 is larger than a machine word, then true, else false. 
inline void add_overflow_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  if (detail::add_overflow_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()))
  {
    result=sort_bool::true_();
    return;
  }
  result=sort_bool::false_();
}

/// \brief An indication whether an overflow occurs when e1 and e2 are added. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return If e1+e2+1 is larger than a machine word, then true, else false. 
inline void add_with_carry_overflow_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  if (detail::add_with_carry_overflow_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()))
  {
    result=sort_bool::true_();
    return;
  }
  result=sort_bool::false_();
}

/// \brief The result of multiplying two words modulo the maximal representable machine word plus 1.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1*e2 modulo the machine word. 
inline void times_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::times_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief The result of multiplying two words and adding the third modulo the maximal representable machine word plus 1.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1*e2+e3 modulo the machine word. 
inline void times_with_carry_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2, const data_expression& e3)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3));
  make_machine_number(result, detail::times_with_carry_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value(),
                              atermpp::down_cast<machine_number>(e3).value()));
}

/// \brief The result of multiplying two words divided by the maximal representable machine word plus 1.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1*e2 div the maximal machine word+1. 
inline void times_overflow_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::times_overflow_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief The result of multiplying two words and adding a third divided by the maximal representable machine word plus 1.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1*e2 div the maximal machine word+1. 
inline void times_with_carry_overflow_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2, const data_expression& e3)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3));
  make_machine_number(result, detail::times_with_carry_overflow_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value(),
                              atermpp::down_cast<machine_number>(e3).value()));
}

/// \brief The result of subtracting two words modulo the maximal representable machine word plus 1.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1-e2 modulo the machine word. 
inline void minus_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::minus_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief The result of subtracting two words. If the result is negative 0 is returned. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return max(0,e1-e2).
inline void monus_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::monus_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief Calculates the division of the first word by the second. 
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1/e2.
inline void div_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::div_word(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief Calculates e1 modulo e2.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return e1 modulo e2. 
inline void mod_word_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::mod_word(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief The square root of e, rounded down to a machine word. 
/// \param result The data_expression into which the created object is stored.
/// \param e The argument.
/// \return The square root of e rounded down. 
inline void sqrt_word_manual_implementation(data_expression& result, const data_expression& e)
{
  assert(is_machine_number(e));
  make_machine_number(result, detail::sqrt_word(atermpp::down_cast<machine_number>(e).value()));
}

/// \brief Calculates (base*e1 + e2) div e3.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \return (base*e1 + e2) div e3
inline void div_doubleword_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2, const data_expression& e3)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3));
  make_machine_number(result, detail::div_doubleword(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value(), 
                              atermpp::down_cast<machine_number>(e3).value()));
}

/// \brief Calculates (base*e1 + e2) div (base*e3 + e4).
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \param e4 The fourth argument. 
/// \return (base*e1 + e2) div (base*e3 + e4)
inline void div_double_doubleword_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2, const data_expression& e3, const data_expression& e4)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3) && is_machine_number(e4));
  make_machine_number(result, detail::div_double_doubleword(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value(), 
                              atermpp::down_cast<machine_number>(e3).value(), 
                              atermpp::down_cast<machine_number>(e4).value()));
}

/// \brief Calculates (base*(base*e1 + e2)+e3) div (base*e4 + e5).
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \param e4 The fourth argument. 
/// \param e5 The fifth argument. 
/// \return (base*(base*e1 + e2)+e3) div (base*e4 + e5)
inline void div_triple_doubleword_manual_implementation(
                        data_expression& result, 
                        const data_expression& e1, 
                        const data_expression& e2, 
                        const data_expression& e3, 
                        const data_expression& e4, 
                        const data_expression& e5)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3) && is_machine_number(e4) && is_machine_number(e5));
  make_machine_number(result, detail::div_triple_doubleword(
                                 atermpp::down_cast<machine_number>(e1).value(), 
                                 atermpp::down_cast<machine_number>(e2).value(), 
                                 atermpp::down_cast<machine_number>(e3).value(), 
                                 atermpp::down_cast<machine_number>(e4).value(),
                                 atermpp::down_cast<machine_number>(e5).value()));
}

/// \brief Calculates (base*e1 + e2) mod (base*e3 + e4).
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \param e4 The fourth argument. 
/// \return (base*e1 + e2) mod (base*e3 + e4)
inline void mod_double_doubleword_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2, const data_expression& e3, const data_expression& e4)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3) && is_machine_number(e4));
  make_machine_number(result, detail::mod_double_doubleword(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value(), 
                              atermpp::down_cast<machine_number>(e3).value(), 
                              atermpp::down_cast<machine_number>(e4).value()));
}

/// \brief The square root of base*e1+e2 rounded down.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \return The square root of base*e1+e2 rounded down.
inline void sqrt_doubleword_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(is_machine_number(e1) && is_machine_number(e2));
  make_machine_number(result, detail::sqrt_doubleword(
                              atermpp::down_cast<machine_number>(e1).value(),
                              atermpp::down_cast<machine_number>(e2).value()));
}

/// \brief Calculates (base*e1 + e2) mod e3. The result fits in one word.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \return (base*e1 + e2) mod e3
inline void mod_doubleword_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2, const data_expression& e3)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3));
  make_machine_number(result, detail::mod_doubleword(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value(), 
                              atermpp::down_cast<machine_number>(e3).value()));
}

/// \brief Calculates the least significant word of the square root of base*(base*e1+e2)+e3.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \return the least significant word of the square root of base*(base*e1+e2)+e3.
inline void sqrt_tripleword_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2, const data_expression& e3)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3));
  make_machine_number(result, detail::sqrt_tripleword(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value(), 
                              atermpp::down_cast<machine_number>(e3).value()));
}

/// \brief Calculates the most significant word of the square root of base*(base*e1+e2)+e3.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \return the most significant word of the square root of base*(base*e1+e2)+e3.
inline void sqrt_tripleword_overflow_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2, const data_expression& e3)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3));
  make_machine_number(result, detail::sqrt_tripleword_overflow(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value(), 
                              atermpp::down_cast<machine_number>(e3).value()));
}

/// \brief Calculates the least significant word of the square root of base*(base*(base*e1+e2)+e3)+e4.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \param e4 The fourth argument. 
/// \return The least significant word of the square root of base*(base*(base*e1+e2)+e3)+e4.
inline void sqrt_quadrupleword_manual_implementation(
                                      data_expression& result, 
                                      const data_expression& e1, 
                                      const data_expression& e2, 
                                      const data_expression& e3, 
                                      const data_expression& e4)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3) && is_machine_number(e4));
  make_machine_number(result, detail::sqrt_quadrupleword(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value(), 
                              atermpp::down_cast<machine_number>(e3).value(), 
                              atermpp::down_cast<machine_number>(e4).value()));
}

/// \brief Calculates the most significant word of the square root of base*(base*(base*e1+e2)+e3)+e4.
/// \param result The data_expression into which the created object is stored.
/// \param e1 The first argument.
/// \param e2 The second argument. 
/// \param e3 The third argument. 
/// \param e4 The fourth argument. 
/// \return The most significant word of the square root of base*(base*(base*e1+e2)+e3)+e4.
inline void sqrt_quadrupleword_overflow_manual_implementation(
                                      data_expression& result, 
                                      const data_expression& e1, 
                                      const data_expression& e2, 
                                      const data_expression& e3, 
                                      const data_expression& e4)
{
  assert(is_machine_number(e1) && is_machine_number(e2) && is_machine_number(e3) && is_machine_number(e4));
  make_machine_number(result, detail::sqrt_quadrupleword_overflow(
                              atermpp::down_cast<machine_number>(e1).value(), 
                              atermpp::down_cast<machine_number>(e2).value(), 
                              atermpp::down_cast<machine_number>(e3).value(), 
                              atermpp::down_cast<machine_number>(e4).value()));
}

/// \brief The predecessor function on a machine numbers, that wraps around.
/// \param result The data_expression into which the created object is stored.
/// \param e 
/// \return e-1, or maxword if e is zero.
inline void pred_word_manual_implementation(data_expression& result, const data_expression& e)
{
  make_machine_number(result, detail::pred_word(atermpp::down_cast<machine_number>(e).value())); 
}

/// \brief The right most bit of a machine number.
/// \param result The data_expression into which the created object is stored.
/// \param e 
/// \return true if the rightmost bit is 1.
inline void rightmost_bit_manual_implementation(data_expression& result, const data_expression& e)
{
  detail::rightmost_bit(result, atermpp::down_cast<machine_number>(e).value()); 
}

/// \brief The machine word shifted one position to the right.
/// \param result The data_expression into which the created object is stored.
/// \param e1 A boolean indicating what the left most bit must be.  
/// \param e2 The value shifted to the right. 
/// \return The machine number e2 divided by 2 prepended with a bit 1 if e1 is true. 
inline void shift_right_manual_implementation(data_expression& result, const data_expression& e1, const data_expression& e2)
{
  assert(e1==sort_bool::true_() || e1==sort_bool::false_());
  make_machine_number(result, detail::shift_right(e1, atermpp::down_cast<machine_number>(e2).value())); 
}


} // namespace mcrl2::data::sort_machine_word





#endif // MCRL2_DATA_MACHINE_NUMBER_H

