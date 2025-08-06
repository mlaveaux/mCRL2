// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/parser_utility.h
/// \brief add your file description here.

#ifndef MCRL2_CORE_PARSER_UTILITY_H
#define MCRL2_CORE_PARSER_UTILITY_H

#include <iostream>

#include "mcrl2/core/parse.h"

namespace mcrl2::core
{

namespace detail {

/// \brief Calls the function f on each node in the parse tree with x as root
template <typename Function>
void foreach_parse_node(const parse_node& x, Function f)
{
  if (x)
  {
    f(x);
    for (int i = 0; i < x.child_count(); i++)
    {
      foreach_parse_node(x.child(i), f);
    }
  }
}

/// \brief Checks if a node is the binary operation op
struct is_binary_operator_node
{
  const parser_table& table;
  std::string op;

  is_binary_operator_node(const parser_table& table_, const std::string& op_)
    : table(table_),
      op(op_)
  {}

  bool operator()(const parse_node& x)
  {
    if (!x)
    {
      return false;
    }
    bool result = true;
    if (x.child_count() != 3)
    {
      result = false;
    }
    else
    {
      if (table.symbol_name(x.child(1).symbol()) != op)
      {
        result = false;
      }
    }
    return result;
  }
};

/// \brief Checks if a node is the binary operation '&&'
struct is_and_node: public is_binary_operator_node
{
  is_and_node(const parser_table& table_)
    : is_binary_operator_node(table_, "&&")
  {}
};

/// \brief Checks if a node is the binary operation '||'
struct is_or_node: public is_binary_operator_node
{
  is_or_node(const parser_table& table_)
    : is_binary_operator_node(table_, "||")
  {}
};

/// \brief Checks if a node is the merge operation '||'
struct is_merge_node: public is_binary_operator_node
{
  is_merge_node(const parser_table& table_)
    : is_binary_operator_node(table_, "||")
  {}
};

/// \brief Checks if a node is the left merge operation '||_'
struct is_left_merge_node: public is_binary_operator_node
{
  is_left_merge_node(const parser_table& table_)
    : is_binary_operator_node(table_, "||_")
  {}
};

/// \brief Checks if a node is of type 'x && (y || z)'
struct is_and_or_node
{
  const parser_table& table;

  is_and_or_node(const parser_table& table_)
    : table(table_)
  {}

  bool operator()(const parse_node& x)
  {
    bool result = true;
    if (!is_and_node(table)(x))
    {
      result = false;
    }
    else
    {
      // check if the types match
      parse_node left = x.child(0);
      parse_node right = x.child(2);
      if (x.symbol() != left.symbol())
      {
        result = false;
      }
      if (x.symbol() != right.symbol())
      {
        result = false;
      }
      if (!is_or_node(table)(right))
      {
        result = false;
      }
    }
    return result;
  }
};

/// \brief Checks if a node is of type 'x ||_ (y || z)'
struct is_left_merge_merge
{
  const parser_table& table;

  is_left_merge_merge(const parser_table& table_)
    : table(table_)
  {}

  bool operator()(const parse_node& x)
  {
    bool result = true;
    if (!is_left_merge_node(table)(x))
    {
      result = false;
    }
    else
    {
      // check if the types match
      parse_node left = x.child(0);
      parse_node right = x.child(2);
      if (x.symbol() != left.symbol())
      {
        result = false;
      }
      if (x.symbol() != right.symbol())
      {
        result = false;
      }
      if (!is_merge_node(table)(right))
      {
        result = false;
      }
    }
    return result;
  }
};

struct find_and_or
{
  const parser_table& table;

  find_and_or(const parser_table& table_)
    : table(table_)
  {}

  void operator()(const parse_node& x)
  {
    if (is_and_or_node(table)(x))
    {
      std::cout << "Warning: the expression of the form 'x && y || z' on location " << x.line() << ":" << x.column() << " may be parsed differently than before" << std::endl;
    }
  }
};

struct find_left_merge_merge
{
  const parser_table& table;

  find_left_merge_merge(const parser_table& table_)
    : table(table_)
  {}

  void operator()(const parse_node& x)
  {
    if (is_left_merge_merge(table)(x))
    {
      std::cout << "Warning: the expression of the form 'x ||_ y || z' on location " << x.line() << ":" << x.column() << " may be parsed differently than before" << std::endl;
    }
  }
};

} // namespace detail

/// \brief Prints a warning for each occurrence of 'x && y || z' in the parse tree.
inline
void warn_and_or(const parse_node& /* node */)
{
  // core::detail::foreach_parse_node(node, core::detail::find_and_or(parser_tables_mcrl2));
}

/// \brief Prints a warning for each occurrence of 'x ||_ y || z' in the parse tree.
inline
void warn_left_merge_merge(const parse_node& /* node */)
{
  // core::detail::foreach_parse_node(node, core::detail::find_left_merge_merge(parser_tables_mcrl2));
}

} // namespace mcrl2::core

#endif // MCRL2_CORE_PARSER_UTILITY_H
