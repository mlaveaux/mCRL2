// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_IO_TEXT_H
#define MCRL2_ATERMPP_ATERM_IO_TEXT_H

#include "mcrl2/atermpp/aterm_io.h"

#include <deque>

namespace atermpp
{

/// \brief Writes terms in textual format to an output stream.
class text_aterm_ostream final : public aterm_ostream
{
public:
  /// \param os The output stream.
  /// \param newline When true each term is written on a new line.
  text_aterm_ostream(std::ostream& os, bool newline = false);

  void put(const aterm& term) override;

private:
  /// \brief Writes a term in textual format on the same line.
  void write_term_line(const aterm& term);

  std::ostream& m_stream;

  bool m_newline = false; ///< Indicates that terms are separated by a newline.
};

/// \brief Reads terms in textual format from an input stream.
class text_aterm_istream final : public aterm_istream
{
public:
  text_aterm_istream(std::istream& os);

  void get(aterm& t) override;

private:
  /// \brief Parse a term from the input stream and return it.
  aterm parse_aterm(int& character);

  /// \brief Parses an "f"(t0, ..., tn) application as an aterm.
  aterm parse_aterm_appl(const std::string& function_name, int& character);

  /// \brief Parses an std::size_t as an aterm_int.
  aterm_int parse_aterm_int(int& character);

  /// \brief Parses a list of arguments [...] as terms.
  aterm_list parse_aterm_list(int& character, char begin, char end);

  /// \brief Reads a quoted string from the stream.
  /// \returns The parsed string and the first character after this string that is not whitespace.
  std::string parse_quoted_string(int& character);

  /// \brief Reads an unquoted string from the stream.
  std::string parse_unquoted_string(int& character);

  /// \returns A string indicating the parse error position.
  std::string print_parse_error_position();

  /// \returns The first character that is not whitespace or end-of-file (EOF).
  /// \param skip_whitespace, returns the next non space character.
  /// \param required Throw error when the next character is EOL.
  int next_char(bool skip_whitespace = true, bool required = false);

  std::istream& m_stream;

  std::size_t m_line = 0; ///< The line number of the current character.
  std::size_t m_column = 0; ///< The column of the current character.

  std::size_t m_history_limit = 64; ///< Determines the maximum number of characters that are stored.
  std::deque<char> m_history; ///< Stores the characters that have been read so-far.

  int character; ///< The last character that was read.
};

} // namespace atermpp

#endif // ATERM_IO_TEXT_H
