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
MCRL2_MODULE;

#include <fstream>
#include <regex>

#include <boost/algorithm/string.hpp>  // for the functions trim, split, trim_copy, is_any_of, trim_right, trim_right_copy.

#include "mcrl2/utilities/logger.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
  #include "mcrl2/utilities/text_utility.cxx"
#else
  module utilities;

  import :exception;
  import :logger;
#endif

namespace mcrl2::utilities
{

/// \brief Split a string into paragraphs.
/// \param text A string
/// \return The paragraphs of <tt>text</tt>
std::vector<std::string> split_paragraphs(const std::string& text)
{
  std::vector<std::string> result;

  // find multiple line endings
  std::regex paragraph_split{R"(\n\s*\n)"};

  // the -1 below directs the token iterator to display the parts of
  // the string that did NOT match the regular expression.
  std::regex_token_iterator cur(text.begin(), text.end(), paragraph_split, -1);
  std::regex_token_iterator<std::string::const_iterator> end;

  for (; cur != end; ++cur)
  {
    std::string paragraph = *cur;
    boost::trim(paragraph);
    if (paragraph.size() > 0)
    {
      result.push_back(paragraph);
    }
  }
  return result;
}

/// \brief Split the text.
/// \param line A string
/// \param separators A string
/// \return The splitted text
std::vector<std::string> split(const std::string& line, const std::string& separators)
{
  std::vector<std::string> result;
  boost::algorithm::split(result, line, boost::algorithm::is_any_of(separators));
  return result;
}

/// \brief Read text from a file.
/// \param filename A string
/// \param warn If true, a warning is printed to standard error if the file is not found
/// \return The contents of the file
std::string read_text(const std::string& filename, bool warn)
{
  std::ifstream in(filename.c_str());
  if (!in)
  {
    if (warn)
    {
      mCRL2log(log::warning) << "Could not open input file: " << filename << std::endl;
      return "";
    }
    else
    {
      throw mcrl2::runtime_error("Could not open input file: " + filename);
    }
  }
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
std::string remove_comments(const std::string& text)
{
  // matches everything from '%' until end of line
  std::regex src {"%[^\\n]*\\n"};

  std::string dest("\n");
  return std::regex_replace(text, src, dest);
}

/// \brief Removes whitespace from a string.
/// \param text A string
/// \return The removal result
std::string remove_whitespace(const std::string& text)
{
  std::regex src {"\\s"};
  std::string dest("");
  return std::regex_replace(text, src, dest);
}

/// \brief Regular expression replacement in a string.
/// \param src A string
/// \param dest A string
/// \param text A string
/// \return The transformed string
std::string regex_replace(const std::string& src, const std::string& dest, const std::string& text)
{
  return std::regex_replace(text, std::regex(src), dest);
}

/// \brief Split a string using a regular expression separator.
/// \param text A string
/// \param sep A string
/// \return The splitted string
std::vector<std::string> regex_split(const std::string& text, const std::string& sep)
{
  std::vector<std::string> result;
  // find multiple line endings
  std::regex paragraph_split { sep };
  // the -1 below directs the token iterator to display the parts of
  // the string that did NOT match the regular expression.
  std::sregex_token_iterator cur(text.begin(), text.end(), paragraph_split, -1);
  std::sregex_token_iterator end;
  for (; cur != end; ++cur)
  {
    std::string word = *cur;
    boost::trim(word);
    if (word.size() > 0)
    {
      result.push_back(word);
    }
  }
  return result;
}

/// \brief Apply word wrapping to a text that doesn't contain newlines.
/// \param line A string of text.
/// \param max_line_length The maximum line length.
/// \return The wrapped text.
static
std::vector<std::string> word_wrap_line(const std::string& line, unsigned int max_line_length)
{
  std::vector<std::string> result;
  std::string text = line;

  for (;;)
  {
    if (text.size() <= max_line_length)
    {
      result.push_back(boost::trim_right_copy(text));
      break;
    }
    std::string::size_type i = text.find_last_of(" \t", max_line_length);
    if (i == std::string::npos)
    {
      result.push_back(text.substr(0, max_line_length));
      text = text.substr(max_line_length);
    }
    else
    {
      result.push_back(text.substr(0, i));
      text = text.substr(i + 1);
    }
  }
  return result;
}

/// \brief Apply word wrapping to a text.
/// \param text A string of text.
/// \param max_line_length The maximum line length.
/// \return The wrapped text.
std::string word_wrap_text(const std::string& text, unsigned int max_line_length)
{
  std::vector<std::string> result;

  // split the lines and remove trailing white space
  std::vector<std::string> lines = split(text, "\n");
  for (auto & line : lines)
  {
    boost::trim_right(line);
  }

  // word wrap each of the lines
  for (auto & line : lines)
  {
    std::vector<std::string> v = word_wrap_line(line, max_line_length);
    result.insert(result.end(), v.begin(), v.end());
  }

  return string_join(result, "\n");
}

/// \brief Test if a string is a number.
/// \param s A string of text.
/// \return True if s is of the form "0 | -? [1-9][0-9]*", false otherwise
bool is_numeric_string(const std::string& s)
{
  // The static below prevents the regular expression recognizer to be compiled
  // each time a string is matched, which is far too time consuming.
  static std::regex re {"0|(-?[1-9][0-9]*)"};
  return std::regex_match(s, re);
}

std::string trim_copy(const std::string& text)
{
  return boost::trim_copy(text);
}

void trim(std::string& text)
{
  boost::trim(text);
}

} // namespace mcrl2::utilities


