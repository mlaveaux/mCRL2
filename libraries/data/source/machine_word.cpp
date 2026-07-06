// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/source/machine_word.cpp
/// \brief This file implements all functions that rely on the boost
///        multiprecision library, such that this library does not need
///        to be in the include path.

#include <boost/multiprecision/cpp_int.hpp>
#include "mcrl2/data/detail/machine_word.h"







namespace mcrl2::data::sort_machine_word::detail
{

/* Below standard operations on machine numbers are implemented.  */

std::size_t times_overflow_word(const std::size_t n1, const std::size_t n2)
{
  boost::multiprecision::uint128_t m1=n1;
  boost::multiprecision::uint128_t m2=n2;
  return static_cast<std::size_t>((m1*m2) >> std::numeric_limits<std::size_t>::digits);
}

std::size_t times_with_carry_overflow_word(const std::size_t n1, const std::size_t n2, const std::size_t n3)
{
  boost::multiprecision::uint128_t m1=n1;
  boost::multiprecision::uint128_t m2=n2;
  boost::multiprecision::uint128_t m3=n3;
  return static_cast<std::size_t>((m1*m2+m3) >> std::numeric_limits<std::size_t>::digits);
}

std::size_t div_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3)
{
  boost::multiprecision::uint128_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  boost::multiprecision::uint128_t m2=n3;
  return static_cast<std::size_t>(m1 / m2);
}

std::size_t mod_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3)
{
  boost::multiprecision::uint128_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  boost::multiprecision::uint128_t m2=n3;
  return static_cast<std::size_t>(m1 % m2);
}

std::size_t div_double_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4)
{
  boost::multiprecision::uint128_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  boost::multiprecision::uint128_t m2=n3;
  m2 = (m2 << std::numeric_limits<std::size_t>::digits)+n4;
  return static_cast<std::size_t>(m1 / m2);
}

std::size_t mod_double_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4)
{
  boost::multiprecision::uint128_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  boost::multiprecision::uint128_t m2=n3;
  m2 = (m2 << std::numeric_limits<std::size_t>::digits)+n4;
  return static_cast<std::size_t>(m1 % m2);
}

std::size_t div_triple_doubleword(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4, const std::size_t n5)
{
  boost::multiprecision::uint256_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n3;
  boost::multiprecision::uint256_t m2=n4;
  m2 = (m2 << std::numeric_limits<std::size_t>::digits)+n5;
  return static_cast<std::size_t>(m1 / m2);
}

std::size_t sqrt_doubleword(const std::size_t n1, const std::size_t n2)
{
  boost::multiprecision::uint128_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  return static_cast<std::size_t>(sqrt(m1));
}

std::size_t sqrt_tripleword(const std::size_t n1, const std::size_t n2, const std::size_t n3)
{
  boost::multiprecision::uint256_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n3;
  return static_cast<std::size_t>(sqrt(m1));
}

std::size_t sqrt_tripleword_overflow(const std::size_t n1, const std::size_t n2, const std::size_t n3)
{
  boost::multiprecision::uint256_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n3;
  return static_cast<std::size_t>(sqrt(m1)>>std::numeric_limits<std::size_t>::digits);
}

std::size_t sqrt_quadrupleword(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4)
{
  boost::multiprecision::uint256_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n3;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n4;
  return static_cast<std::size_t>(sqrt(m1));
}

std::size_t sqrt_quadrupleword_overflow(const std::size_t n1, const std::size_t n2, const std::size_t n3, const std::size_t n4)
{
  boost::multiprecision::uint256_t m1=n1;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n2;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n3;
  m1 = (m1 << std::numeric_limits<std::size_t>::digits)+n4;
  return static_cast<std::size_t>(sqrt(m1)>>std::numeric_limits<std::size_t>::digits);
}


} // end namespace detail
