// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/parse.h
///
/// \brief Parse mCRL2 specifications and expressions.

#ifndef MCRL2_CORE_PARSE_H
#define MCRL2_CORE_PARSE_H

#include <iostream>

#include "mcrl2/atermpp/aterm_io.h"
#include "mcrl2/core/detail/dparser_functions.h"
#include "mcrl2/core/identifier_string.h"

struct D_ParserTables; // prototype

extern "C"
{
  extern D_ParserTables parser_tables_mcrl2;
  extern D_ParserTables parser_tables_fsm;
  extern D_ParserTables parser_tables_dot;
}

namespace mcrl2::core
{

class parse_node_exception : public mcrl2::runtime_error
{
protected:
  static std::string get_error_message(const parse_node& node, const std::string& user_message)
  {
    try
    {
      return node.add_context(user_message);
    }
    catch (...)
    {
      return user_message;
    }
  }

  parse_node_exception(const std::string& message)
    : mcrl2::runtime_error(message)
  { }

public:
  parse_node_exception(const parse_node& node, const std::string& message)
    : mcrl2::runtime_error(get_error_message(node, message))
  { }
};

class parse_node_unexpected_exception : public parse_node_exception
{
private:
  static std::string get_error_message(const parser& p, const parse_node& node)
  {
    std::string inherited = parse_node_exception::get_error_message(node, "unexpected parse node!");
    try
    {
      std::ostringstream out;
      out << inherited << std::endl;
      p.print_node(out, node);
      return out.str();
    }
    catch (...)
    {
      return inherited;
    }
  }
public:
  parse_node_unexpected_exception(const parser& p, const parse_node& node)
    : parse_node_exception(get_error_message(p, node))
  { }
};

struct parser_actions
{
  const parser& m_parser;

  parser_actions(const parser& parser_)
    : m_parser(parser_)
  {}

  // starts a traversal in node, and calls the function f to each subnode of the given type
  template <typename Function>
  void traverse(const parse_node& node, const Function& f) const
  {
    if (!node)
    {
      return;
    }
    if (!f(node))
    {
      for (int i = 0; i < node.child_count(); i++)
      {
        traverse(node.child(i), f);
      }
    }
  }

  // callback function that applies a function to nodes of a given type
  template <typename Function>
  struct visitor
  {
    const parser_table& table;
    const std::string& type;
    const Function& f;

    visitor(const parser_table& table_, const std::string& type_, const Function& f_)
      : table(table_),
        type(type_),
        f(f_)
    {}

    bool operator()(const parse_node& node) const
    {
      if (table.symbol_name(node) == type)
      {
        f(node);
        return true;
      }
      return false;
    }
  };

  template <typename Function>
  visitor<Function> make_visitor(const parser_table& table, const std::string& type, const Function& f) const
  {
    return visitor<Function>(table, type, f);
  }

  // callback function that applies a function to a node, and adds the result to a container
  template <typename Container, typename Function>
  struct collector
  {
    const parser_table& table;
    const std::string& type;
    Container& container;
    const Function& f;

    collector(const parser_table& table_, const std::string& type_, Container& container_, const Function& f_)
      : table(table_),
        type(type_),
        container(container_),
        f(f_)
    {}

    bool operator()(const parse_node& node) const
    {
      if (table.symbol_name(node) == type)
      {
        container.push_back(f(node));
        return true;
      }
      return false;
    }
  };

  template <typename Container, typename Function>
  collector<Container, Function> make_collector(const parser_table& table, const std::string& type, Container& container, const Function& f) const
  {
    return collector<Container, Function>(table, type, container, f);
  }

  // callback function that applies a function to a node, and adds the result to a set container
  template <typename SetContainer, typename Function>
  struct set_collector
  {
    const parser_table& table;
    const std::string& type;
    SetContainer& container;
    const Function& f;

    set_collector(const parser_table& table_, const std::string& type_, SetContainer& container_, const Function& f_)
      : table(table_),
        type(type_),
        container(container_),
        f(f_)
    {}

    bool operator()(const parse_node& node) const
    {
      if (table.symbol_name(node) == type)
      {
        container.insert(f(node));
        return true;
      }
      return false;
    }
  };

  template <typename SetContainer, typename Function>
  set_collector<SetContainer, Function> make_set_collector(const parser_table& table, const std::string& type, SetContainer& container, const Function& f) const
  {
    return set_collector<SetContainer, Function>(table, type, container, f);
  }

  std::string symbol_name(const parse_node& node) const
  {
    return m_parser.symbol_table().symbol_name(node.symbol());
  }
};

struct default_parser_actions: public parser_actions
{
  default_parser_actions(const parser& parser_)
    : parser_actions(parser_)
  {}

  template <typename T, typename Function>
  atermpp::term_list<T> parse_list(const parse_node& node, const std::string& type, const Function& f) const
  {
    std::vector<T> result;
    traverse(node, make_collector(m_parser.symbol_table(), type, result, f));
    return atermpp::term_list<T>(result.begin(), result.end());
  }

  template <typename T, typename Function>
  std::vector<T> parse_vector(const parse_node& node, const std::string& type, const Function& f) const
  {
    std::vector<T> result;
    traverse(node, make_collector(m_parser.symbol_table(), type, result, f));
    return result;
  }

  core::identifier_string parse_Id(const parse_node& node) const
  {
    return core::identifier_string(node.string());
  }

  core::identifier_string parse_Number(const parse_node& node) const
  {
    return core::identifier_string(node.string());
  }

  core::identifier_string_list parse_IdList(const parse_node& node) const
  {
    return parse_list<core::identifier_string>(node, "Id", [&](const core::parse_node& node) { return parse_Id(node); });
  }
};

template <typename T>
void print_aterm(const T&)
{
}

template<>
inline
void print_aterm(const atermpp::aterm& x)
{
  std::clog << "aterm: " << x << std::endl;
}

/// \brief Parse an identifier.
inline
identifier_string parse_identifier(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("Id");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  identifier_string result = default_parser_actions(p).parse_Id(node);
  return result;
}

inline
bool is_user_identifier(std::string const& s)
{
  try
  {
    parse_identifier(s);
  }
  catch (...)
  {
    return false;
  }
  return true;
}

} // namespace mcrl2::core

#endif // MCRL2_CORE_PARSE_H
