// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/text_utility.h
/// \brief String manipulation functions.

#ifndef MCRL2_UTILITIES_TEXT_UTILITY_H
#define MCRL2_UTILITIES_TEXT_UTILITY_H
MCRL2_MODULE;

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifdef MCRL2_ENABLE_MODULES
  export module utilities:text_utility;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities
{

/// \brief Transform parameter into string.
/// \param x Some expression
/// \pre type T has operator <<
/// \return The string representation of x.
template <typename T>
inline
std::string to_string(const T& x)
{
  std::stringstream ss;
  ss << x;
  return ss.str();
}

/// \brief Split a string into paragraphs.
/// \param text A string
/// \return The paragraphs of <tt>text</tt>
std::vector<std::string> split_paragraphs(const std::string& text);

/// \brief Split the text.
/// \param line A string
/// \param separators A string
/// \return The splitted text
std::vector<std::string> split(const std::string& line, const std::string& separators);

/// \brief Read text from a file.
/// \param filename A string
/// \param warn If true, a warning is printed to standard error if the file is not found
/// \return The contents of the file
std::string read_text(const std::string& filename, bool warn=false);

/// \brief Read text from a stream.
/// \param in An input stream
/// \return The text read from the stream
inline
std::string read_text(std::istream& in)
{
  in.unsetf(std::ios::skipws); //  Turn of white space skipping on the stream
  std::string s;
  std::copy(
    std::istream_iterator<char>(in),
    std::istream_iterator<char>(),
    std::back_inserter(s)
  );
  return s;
}

/// \brief Remove comments from a text (everything from '%' until end of line).
/// \param text A string
/// \return The removal result
std::string remove_comments(const std::string& text);

/// \brief Removes whitespace from a string.
/// \param text A string
/// \return The removal result
std::string remove_whitespace(const std::string& text);

/// \brief Regular expression replacement in a string.
/// \param src A string
/// \param dest A string
/// \param text A string
/// \return The transformed string
std::string regex_replace(const std::string& src, const std::string& dest, const std::string& text);

/// \brief Split a string using a regular expression separator.
/// \param text A string
/// \param sep A string
/// \return The splitted string
std::vector<std::string> regex_split(const std::string& text, const std::string& sep);

/// \brief Remove all trailing and leading spaces from the input.
/// \param text A string
/// \return The trimmed string
std::string trim_copy(const std::string& text);

/// \brief Remove all trailing and leading spaces from the input.
/// \param text A string
void trim(std::string& text);

/// \brief Joins a sequence of strings. This is a replacement for
/// boost::algorithm::join, since it gives stack overflow errors with
/// Visual C++ express 9.0 under some circumstances.
template <typename Container>
std::string string_join(const Container& c, const std::string& separator)
{
  std::ostringstream out;
  for (typename Container::const_iterator i = c.begin(); i != c.end(); ++i)
  {
    if (i != c.begin())
    {
      out << separator;
    }
    out << *i;
  }
  return out.str();
}

/// \brief Apply word wrapping to a text.
/// \param text A string of text.
/// \param max_line_length The maximum line length.
/// \return The wrapped text.
std::string word_wrap_text(const std::string& text, unsigned int max_line_length = 78);

/// \brief Test if a string is a number.
/// \param s A string of text.
/// \return True if s is of the form "0 | -? [1-9][0-9]*", false otherwise
bool is_numeric_string(const std::string& s);

/// \brief Convert a number to a string in the buffer starting at position start_position.
/// \param number The number to be converted.
/// \param buffer A buffer in which the string will be stored that is sufficiently large.
/// \param start_position The first position where a number is written. 
inline void number2string(std::size_t number, std::string& buffer, std::size_t start_position)
{
  // First calculate the number of digital digits of number;
  std::size_t number_of_digits=0;
  if (number==0)
  {
    number_of_digits=1;
  }
  else
  {
    for(std::size_t copy=number ; copy!=0; ++number_of_digits, copy=copy/10)
    {}
  }

  // Put the number in the buffer at the right position.
  std::size_t position=start_position+number_of_digits;
  buffer.resize(position);

  while (position>start_position)
  {
    --position;
    buffer[position] = '0' + number % 10;
    number = number/10;
  }
}

/// \brief Convert a number to string. 
/// \details This function is much faster than std::to_string and 
///          its use is therefore preferred in those cases were performance counts. 
/// \param number A number to be transformed. 
inline std::string number2string(std::size_t number) 
{ 
  std::string buffer;
  buffer.reserve(std::numeric_limits<std::size_t>::digits10 + 1); 
  number2string(number, buffer, 0); 
  return std::string(buffer); 
}  

} // namespace mcrl2::utilities



#endif // MCRL2_UTILITIES_TEXT_UTILITY_H
