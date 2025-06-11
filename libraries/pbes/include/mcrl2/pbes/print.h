// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/print.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PRINT_H
#define MCRL2_PBES_PRINT_H

#include "mcrl2/data/print.h"
#include "mcrl2/pbes/traverser.h"

namespace mcrl2 {

namespace pbes_system {

constexpr inline int precedence(const forall&) { return 21; }
constexpr inline int precedence(const exists&) { return 21; }
constexpr inline int precedence(const imp&)    { return 22; }
constexpr inline int precedence(const or_&)    { return 23; }
constexpr inline int precedence(const and_&)   { return 24; }
constexpr inline int precedence(const not_&)   { return 25; }
inline int precedence(const pbes_expression& x)
{
  if      (is_forall(x)) { return precedence(atermpp::down_cast<forall>(x)); }
  else if (is_exists(x)) { return precedence(atermpp::down_cast<exists>(x)); }
  else if (is_imp(x))    { return precedence(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))     { return precedence(atermpp::down_cast<or_>(x)); }
  else if (is_and(x))    { return precedence(atermpp::down_cast<and_>(x)); }
  else if (is_not(x))    { return precedence(atermpp::down_cast<not_>(x)); }
  return core::detail::max_precedence;
}

// only defined for binary operators
inline bool is_left_associative(const imp&)  { return false; }
inline bool is_left_associative(const or_&)  { return true; }
inline bool is_left_associative(const and_&) { return true; }
inline bool is_left_associative(const pbes_expression& x)
{
  if (is_imp(x))      { return is_left_associative(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))  { return is_left_associative(atermpp::down_cast<or_>(x)); }
  else if (is_and(x)) { return is_left_associative(atermpp::down_cast<and_>(x)); }
  return false;
}

inline bool is_right_associative(const imp&)  { return true; }
inline bool is_right_associative(const or_&)  { return true; }
inline bool is_right_associative(const and_&) { return true; }
inline bool is_right_associative(const pbes_expression& x)
{
  if (is_imp(x))      { return is_right_associative(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))  { return is_right_associative(atermpp::down_cast<or_>(x)); }
  else if (is_and(x)) { return is_right_associative(atermpp::down_cast<and_>(x)); }
  return false;
}

namespace detail {

template <typename Derived>
struct printer: public pbes_system::add_traverser_sort_expressions<data::detail::printer, Derived>
{
  typedef pbes_system::add_traverser_sort_expressions<data::detail::printer, Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::derived;
  using super::print_abstraction;
  using super::print_list;
  using super::print_variables;

  // N.B. We need a special version due to the "val" operator that needs to be
  // put around data expressions.
  template <typename T>
  void print_pbes_expression(const T& x, bool needs_parentheses = false)
  {
    bool is_data_expr = is_data(x);
    if (is_data_expr)
    {
      derived().print("val");
    }
    if (needs_parentheses || is_data_expr)
    {
      derived().print("(");
    }
    derived().apply(x);
    if (needs_parentheses || is_data_expr)
    {
      derived().print(")");
    }
  }

  void print_pbes_unary_operand(const pbes_expression& x, const pbes_expression& operand)
  {
    print_pbes_expression(operand, precedence(operand) < precedence(x));
  }

  // N.B. We need a special version due to the "val" operator that needs to be
  // put around data expressions.
  template <typename T>
  void print_pbes_unary_left_operation(const T& x, const std::string& op)
  {
    derived().print(op);
    print_pbes_unary_operand(x, x.operand());
  }

  // N.B. We need a special version due to the "val" operator that needs to be
  // put around data expressions.
  template <typename T>
  void print_pbes_binary_operation(const T& x, const std::string& op)
  {
    const auto& x1 = x.left();
    const auto& x2 = x.right();
    auto p = precedence(x);
    auto p1 = precedence(x1);
    auto p2 = precedence(x2);
    print_pbes_expression(x1, (p1 < p) || (p1 == p && !is_left_associative(x)));
    derived().print(op);
    print_pbes_expression(x2, (p2 < p) || (p2 == p && !is_right_associative(x)));
  }

  // N.B. We need a special version due to the "val" operator that needs to be
  // put around data expressions.
  template <typename Abstraction>
  void print_pbes_abstraction(const Abstraction& x, const std::string& op)
  {
    derived().enter(x);
    derived().print(op + " ");
    print_variables(x.variables(), true, true, false, "", "", ", ");
    derived().print(". ");
    print_pbes_unary_operand(x, x.body());
    derived().leave(x);
  }

  void apply(const pbes_system::propositional_variable& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    print_variables(x.parameters());
    derived().leave(x);
  }

  void apply(const pbes_system::fixpoint_symbol& x)
  {
    derived().enter(x);
    derived().print(x.is_mu() ? "mu" : "nu");
    derived().leave(x);
  }

  void apply(const pbes_system::pbes_equation& x)
  {
    derived().enter(x);
    derived().apply(x.symbol());
    derived().print(" ");
    derived().apply(x.variable());
    // TODO: change the weird convention of putting the rhs of an equation on a new line
    derived().print(" =\n       ");
    bool print_val = data::is_data_expression(x.formula());
    if (print_val)
    {
      derived().print("val(");
    }
    derived().apply(x.formula());
    if (print_val)
    {
      derived().print(")");
    }
    derived().print(";");
    derived().leave(x);
  }

  void apply(const pbes_system::pbes& x)
  {
    derived().enter(x);
    derived().apply(x.data());
    print_variables(x.global_variables(), true, true, true, "glob ", ";\n\n", ";\n     ");
    print_list(x.equations(), "pbes ", "\n\n", "\n     ");
    derived().print("init ");
    print_pbes_expression(x.initial_state());
    derived().print(";\n");
    derived().leave(x);
  }

  void apply(const pbes_system::propositional_variable_instantiation& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    print_list(x.parameters(), "(", ")", ", ", false);
    derived().leave(x);
  }

  void apply(const pbes_system::not_& x)
  {
    derived().enter(x);
    print_pbes_unary_left_operation(x, "!");
    derived().leave(x);
  }

  void apply(const pbes_system::and_& x)
  {
    derived().enter(x);
    print_pbes_binary_operation(x, " && ");
    derived().leave(x);
  }

  void apply(const pbes_system::or_& x)
  {
    derived().enter(x);
    print_pbes_binary_operation(x, " || ");
    derived().leave(x);
  }

  void apply(const pbes_system::imp& x)
  {
    derived().enter(x);
    print_pbes_binary_operation(x, " => ");
    derived().leave(x);
  }

  void apply(const pbes_system::forall& x)
  {
    derived().enter(x);
    print_pbes_abstraction(x, "forall");
    derived().leave(x);
  }

  void apply(const pbes_system::exists& x)
  {
    derived().enter(x);
    print_pbes_abstraction(x, "exists");
    derived().leave(x);
  }
};

} // namespace detail

/// \brief Prints the object x to a stream.
struct stream_printer
{
  template <typename T>
  void operator()(const T& x, std::ostream& out, bool precedence_aware)
  {
    core::detail::apply_printer<pbes_system::detail::printer> printer(out, precedence_aware);
    printer.apply(x);
  }
};

/// \brief Returns a string representation of the object x.
template <typename T>
std::string pp(const T& x, bool precedence_aware = true)
{
  std::ostringstream out;
  stream_printer()(x, out, precedence_aware);
  return out.str();
}

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_PRINT_H
