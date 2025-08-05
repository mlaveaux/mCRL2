// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/parse_numbers.h
/// \brief add your file description here.

#ifndef MCRL2_UTILITIES_PARSE_NUMBERS_H
#define MCRL2_UTILITIES_PARSE_NUMBERS_H

#include <cctype>
#include <cstddef>
#include <string>
#include <vector>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::utilities
{

namespace detail {

// Reads the next integer from the range [first, last), and the spaces behind it
// Returns the position in the range after the next integer
// Precondition: (first != last) && !std::isspace(*first)
template <typename Iterator>
Iterator parse_next_natural_number(Iterator first, Iterator last, std::size_t& result)
{
  assert((first != last) && !std::isspace(*first));

  Iterator i = first;
  result = 0;

  for (;;)
  {
    if (*i < '0' || *i > '9')
    {
      throw mcrl2::runtime_error("could not read an integer from " + std::string(first, last));
    }
    result *= 10;
    result += *i - '0';
    ++i;
    if (i == last)
    {
      break;
    }
    if (std::isspace(*i))
    {
      ++i;
      while (i != last && std::isspace(*i))
      {
        ++i;
      }
      break;
    }
  }
  return i;
}

} // namespace detail

/// \brief Parses a natural number from a string
inline
std::size_t parse_natural_number(const std::string& text)
{
  auto first = text.begin();
  auto last = text.end();

  // skip leading spaces
  while (first != last && std::isspace(*first))
  {
    ++first;
  }

  if (first == last)
  {
    throw mcrl2::runtime_error("could not read an integer from " + text);
  }

  std::size_t value;
  first = detail::parse_next_natural_number(first, last, value);

  if (first != last)
  {
    throw mcrl2::runtime_error("could not read an integer from " + text);
  }

  return value;
}

/// \brief Parses a sequence of natural numbers (separated by spaces) from a string
inline
std::vector<std::size_t> parse_natural_number_sequence(const std::string& text)
{
  std::vector<std::size_t> result;

  auto first = text.begin();
  auto last = text.end();

  // skip leading spaces
  while (first != last && std::isspace(*first))
  {
    ++first;
  }

  while (first != last)
  {
    std::size_t value;
    first = detail::parse_next_natural_number(first, last, value);
    result.push_back(value);
  }

  return result;
}

} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_PARSE_NUMBERS_H
