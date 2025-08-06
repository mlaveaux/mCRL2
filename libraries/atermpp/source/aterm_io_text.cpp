// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/atermpp/aterm_io_text.h"

#include <fstream>

namespace atermpp
{

// utility functions

static void write_string_with_escape_symbols(const std::string& s, std::ostream& os)
{
  // Check whether the string starts with a - or a number, or contains the symbols
  // \, ", (, ), [, ], comma, space, \t, \n or \r. If yes, the string will be
  // surrounded by quotes, and the symbols \, ", \t, \n, \r
  // will be preceded by an escape symbol.

  char c = s[0];
  bool contains_special_symbols = ((c =='-') || isdigit(c));

  for(std::string::const_iterator i=s.begin(); !contains_special_symbols && i!=s.end(); ++i)
  {
    if (*i=='\\' || *i=='"' || *i=='(' || *i==')' || *i=='[' || *i==']' || *i==',' || *i==' ' || *i=='\n' || *i=='\t' || *i=='\r')
    {
      contains_special_symbols=true;
    }
  }

  if (contains_special_symbols)
  {
    // This function symbol needs quotes.
    os << "\"";
    for (char i: s)
    {
      // We need to escape special characters.
      switch (i)
      {
        case '\\':
        case '"':
          os << "\\" << i;
          break;
        case '\n':
          os << "\\n";
          break;
        case '\t':
          os << "\\t";
          break;
        case '\r':
          os << "\\r";
          break;
        default:
          os << i;
          break;
      }
    }
    os << "\"";
  }
  else
  {
    os << s;
  }
}

// Public functions

text_aterm_ostream::text_aterm_ostream(std::ostream& os, bool newline)
  : m_stream(os),
    m_newline(newline)
{}

void text_aterm_ostream::put(const aterm& term)
{
  write_term_line(term);

  if (m_newline)
  {
    m_stream << "\n";
  }
}

text_aterm_istream::text_aterm_istream(std::istream& is)
  : m_stream(is)
{
  character = next_char();
}

void text_aterm_istream::get(aterm& term)
{
  try
  {
    if (character != EOF)
    {
      term = parse_aterm(character);
    }
  }
  catch (std::runtime_error& e)
  {
    throw std::runtime_error(e.what() + std::string("\n") + print_parse_error_position());
  }

  // Reset the parsing error buffers.
  m_column = 0;
  m_history.clear();

  return;
}

// Private functions

void text_aterm_ostream::write_term_line(const aterm& t)
{
  if (t.type_is_int())
  {
    // Write a single integer as is.
    m_stream << atermpp::down_cast<aterm_int>(t).value();
  }
  else if (t.type_is_list())
  {
    // A list l0...ln is formatted as [l0, ..., ln].
    m_stream << "[";
    const aterm_list& list = down_cast<aterm_list>(t);
    for (aterm_list::const_iterator it = list.begin(); it != list.end(); ++it)
    {
      if (it!=list.begin())
      {
        m_stream << ",";
      }
      write_term_line(*it);
    }

    m_stream << "]";
  }
  else
  {
    // An aterm f(t0, ..., tn) is written as f(t0, ..., tn)
    assert(t.type_is_appl());

    aterm appl = m_transformer(t);

    write_string_with_escape_symbols(appl.function().name(), m_stream);

    if (appl.function().arity() > 0)
    {
      m_stream << "(";
      write_term_line(appl[0]);
      for (std::size_t i = 1; i < appl.function().arity(); i++)
      {
        m_stream << ",";
        write_term_line(appl[i]);
      }
      m_stream << ")";
    }
  }
}

aterm text_aterm_istream::parse_aterm(int& character)
{
  // Parse the term.
  switch (character)
  {
    case '"':
    {
      std::string function_name = parse_quoted_string(character);
      return parse_aterm_appl(function_name, character);
    }
    case '[':
    {
      return parse_aterm_list(character, '[', ']');
    }
    default:
    {
      if (isdigit(character) || character == '-')
      {
        return parse_aterm_int(character);
      }

      std::string function_symbol = parse_unquoted_string(character);
      return parse_aterm_appl(function_symbol, character);
    }
  }
}

aterm text_aterm_istream::parse_aterm_appl(const std::string& function_name, int& character)
{
  // Parse the arguments.
  aterm_list arguments = parse_aterm_list(character, '(', ')');

  // Wrap up this function application.
  function_symbol symbol(function_name, arguments.size());
  return m_transformer(aterm(symbol, arguments.begin(), arguments.end()));
}

aterm_int text_aterm_istream::parse_aterm_int(int& character)
{
  std::array<char, 32> number{};
  auto it = number.begin();

  if (character == '-')
  {
    *it = static_cast<char>(character);
    ++it;
    character = next_char(true, true);
  }

  while (isdigit(character) && it != number.end())
  {
    *it = static_cast<char>(character);
    ++it;
    character = next_char();
  }

  *it = '\0';
  return aterm_int(static_cast<std::size_t>(atol(number.data())));
}

aterm_list text_aterm_istream::parse_aterm_list(int& character, char begin, char end)
{
  aterm_list list;

  // A list is [t0, ..., tn] or surrounded by ().
  if (character == begin)
  {
    character = next_char(true, true);
    if (character != end)
    {
      list.push_front(parse_aterm(character));

      while (character == ',')
      {
        character = next_char(true, true);
        list.push_front(parse_aterm(character));
      }

      if (character != end)
      {
        throw std::runtime_error(std::string("Missing ") + end + " while parsing a list term");
      }
    }

    character = next_char(true);
  }

  return reverse(list);
}


std::string text_aterm_istream::print_parse_error_position()
{
  std::stringstream s;
  s << "Error occurred at line " << m_line << ", col " << m_column << " near: ";
  for(const auto& element : m_history)
  {
    s << element;
  }
  return s.str();
}

int text_aterm_istream::next_char(bool skip_whitespace, bool required)
{
  character = EOF;

  do
  {
    try
    {
      // In liblts_lts the exception bit is set, so we need to use exception to handle EOF.
      character = m_stream.get();
    }
    catch (std::ios::failure&)
    {
      return EOF;
    }

    if (character != EOF)
    {
      if (character == '\n')
      {
        m_line++;
        m_column = 0;
      }
      else
      {
        m_column++;
      }

      if (m_history.size() >= m_history_limit)
      {
        // If the history is full the first element must be removed.
        m_history.erase(m_history.begin());
      }

      m_history.emplace_back(character);
    }
    else if (required)
    {
      throw std::runtime_error("Premature end of file while parsing.");
    }
  }
  while (isspace(character) && skip_whitespace);

  // The stream also returns a newline for the last symbol.
  return character == '\n' ? EOF : character;
}

std::string text_aterm_istream::parse_quoted_string(int& character)
{
  // We need a buffer for printing and parsing.
  std::string string;

  assert(character == '"');

  // First obtain the first symbol after ".
  character = next_char();

  while (character != '"')
  {
    switch (character)
    {
      case '\\':
        character = next_char(false, true);
        switch (character)
        {
          case 'n':
            string += '\n';
            break;
          case 'r':
            string += '\r';
            break;
          case 't':
            string += '\t';
            break;
          default:
            string += static_cast<char>(character);
            break;
        }
        break;
      default:
        string += static_cast<char>(character);
        break;
    }
    character = next_char(false, true);
  }

  character = next_char(true, false);
  return string;
}

std::string text_aterm_istream::parse_unquoted_string(int& character)
{
  std::string string;

  if (character != '(')
  {
    // First parse the identifier
    while (character != '"' && character != '(' && character != ')' && character != ']'
      && character != ']' && character != ',' && character != ' ' && character != '\n'
      && character != '\t' && character != '\r' && character != EOF)
    {
      string += static_cast<char>(character);
      character = next_char(false);
    }
  }

  return string;
}

void write_term_to_text_stream(const aterm& term, std::ostream& os)
{
  text_aterm_ostream(os) << term;
}

aterm read_term_from_string(const std::string& s)
{
  std::stringstream ss(s);
  aterm t;
  read_term_from_text_stream(ss, t);
  return t;
}

void read_term_from_text_stream(std::istream& is, aterm& t)
{
  text_aterm_istream(is).get(t);
}

} // namespace atermpp
