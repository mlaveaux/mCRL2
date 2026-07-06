// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/math.h
/// \brief Some general mathematical functions not provided by standard libraries.

#ifndef MCRL2_UTILITIES_MATH_H
#define MCRL2_UTILITIES_MATH_H

#include <cassert>
#include <cmath>
#include <cstdlib>

namespace mcrl2::utilities
{

// Compute base 2 logarithm of n, by checking which is the leftmost
// bit that has been set.

inline
std::size_t ceil_log2(std::size_t n)
{
  assert(n>0);
  std::size_t result = 0;
  while(n != 0)
  {
    n = n >> 1;
    ++result;
  }
  return result;
}

// Calculate n^m for numbers n,m of type std::size_t
inline 
std::size_t power_size_t(const std::size_t n_in, const std::size_t m_in)
{
  std::size_t result=1;
  std::size_t n=n_in;
  std::size_t m=m_in;
  while (m>0)  // Invariant: result*n^m=n_in^m_in;
  {
    if (m % 2==1)
    {
      result=result*n;
    }
    n=n*n;
    m=m/2;
  }
  return result;
}

template <typename T>
bool is_odd(T t)
{
  return t % 2 != 0;
}

template <typename T>
bool is_even(T t)
{
  return t % 2 == 0;
}

} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_MATH_H

