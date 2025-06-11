// Author(s): Jan Friso Groote. Based on pbes/print.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/print.h
/// \brief add your file description here.

#ifndef MCRL2_PRES_PRINT_H
#define MCRL2_PRES_PRINT_H

#include "mcrl2/data/print.h"
#include "mcrl2/pres/traverser.h"

namespace mcrl2 {

namespace pres_system {

constexpr inline int precedence(const infimum&)  { return 21; }
constexpr inline int precedence(const supremum&) { return 21; }
constexpr inline int precedence(const sum&)      { return 21; }
constexpr inline int precedence(const plus&)     { return 22; }
constexpr inline int precedence(const imp&)      { return 23; }
constexpr inline int precedence(const or_&)      { return 24; }
constexpr inline int precedence(const and_&)     { return 25; }
constexpr inline int precedence(const minus&)    { return 25; }
constexpr inline int precedence(const const_multiply&) { return 26; }
constexpr inline int precedence(const const_multiply_alt&) { return 26; }
inline int precedence(const pres_expression& x)
{
  if      (is_infimum(x))           { return precedence(atermpp::down_cast<infimum>(x)); }
  else if (is_supremum(x))          { return precedence(atermpp::down_cast<supremum>(x)); }
  else if (is_sum(x))               { return precedence(atermpp::down_cast<sum>(x)); }
  else if (is_imp(x))               { return precedence(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))                { return precedence(atermpp::down_cast<or_>(x)); }
  else if (is_and(x))               { return precedence(atermpp::down_cast<and_>(x)); }
  else if (is_plus(x))              { return precedence(atermpp::down_cast<plus>(x)); }
  else if (is_const_multiply(x))    { return precedence(atermpp::down_cast<const_multiply>(x)); }
  else if (is_const_multiply_alt(x)){ return precedence(atermpp::down_cast<const_multiply_alt>(x)); }
  else if (is_minus(x))             { return precedence(atermpp::down_cast<minus>(x)); }
  return core::detail::max_precedence;
}

// only defined for binary operators
inline bool is_left_associative(const imp&)  { return false; }
inline bool is_left_associative(const or_&)  { return true; }
inline bool is_left_associative(const and_&) { return true; }
inline bool is_left_associative(const plus&) { return true; }
inline bool is_left_associative(const const_multiply&) { return true; }
inline bool is_left_associative(const const_multiply_alt&) { return true; }
inline bool is_left_associative(const pres_expression& x)
{
  if (is_imp(x))       { return is_left_associative(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))   { return is_left_associative(atermpp::down_cast<or_>(x)); }
  else if (is_and(x))  { return is_left_associative(atermpp::down_cast<and_>(x)); }
  else if (is_plus(x)) { return is_left_associative(atermpp::down_cast<plus>(x)); }
  else if (is_const_multiply(x)) { return is_left_associative(atermpp::down_cast<const_multiply>(x)); }
  else if (is_const_multiply_alt(x)) { return is_left_associative(atermpp::down_cast<const_multiply_alt>(x)); }
  return false;
}

inline bool is_right_associative(const imp&)  { return true; }
inline bool is_right_associative(const or_&)  { return true; }
inline bool is_right_associative(const and_&) { return true; }
inline bool is_right_associative(const plus&) { return true; }
inline bool is_right_associative(const const_multiply&) { return true; }
inline bool is_right_associative(const const_multiply_alt&) { return true; }
inline bool is_right_associative(const pres_expression& x)
{
  if (is_imp(x))       { return is_right_associative(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))   { return is_right_associative(atermpp::down_cast<or_>(x)); }
  else if (is_and(x))  { return is_right_associative(atermpp::down_cast<and_>(x)); }
  else if (is_plus(x)) { return is_right_associative(atermpp::down_cast<plus>(x)); }
  else if (is_const_multiply(x)) { return is_right_associative(atermpp::down_cast<const_multiply>(x)); }
  else if (is_const_multiply_alt(x)) { return is_right_associative(atermpp::down_cast<const_multiply_alt>(x)); }
  return false;
}

namespace detail {

template <typename Derived>
struct printer: public pres_system::add_traverser_sort_expressions<data::detail::printer, Derived>
{
  typedef pres_system::add_traverser_sort_expressions<data::detail::printer, Derived> super;

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
  void print_pres_expression(const T& x, bool needs_parentheses = false)
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

  void print_pres_unary_operand(const pres_expression& x, const pres_expression& operand)
  {
    print_pres_expression(operand, precedence(operand) < precedence(x));
  }

  // N.B. We need a special version due to the "val" operator that needs to be
  // put around data expressions.
  template <typename T>
  void print_pres_unary_left_operation(const T& x, const std::string& op)
  {
    derived().print(op);
    print_pres_unary_operand(x, x.operand());
  }

  // N.B. We need a special version due to the "val" operator that needs to be
  // put around data expressions.
  template <typename T>
  void print_pres_binary_operation(const T& x, const std::string& op)
  {
    const auto& x1 = x.left();
    const auto& x2 = x.right();
    auto p = precedence(x);
    auto p1 = precedence(x1);
    auto p2 = precedence(x2);
    print_pres_expression(x1, (p1 < p) || (p1 == p && !is_left_associative(x)));
    derived().print(op);
    print_pres_expression(x2, (p2 < p) || (p2 == p && !is_right_associative(x)));
  }

  // N.B. We need a special version due to the "val" operator that needs to be
  // put around data expressions.
  template <typename Abstraction>
  void print_pres_abstraction(const Abstraction& x, const std::string& op)
  {
    derived().enter(x);
    derived().print(op + " ");
    print_variables(x.variables(), true, true, false, "", "", ", ");
    derived().print(". ");
    print_pres_unary_operand(x, x.body());
    derived().leave(x);
  }

  void apply(const pres_system::propositional_variable& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    print_variables(x.parameters());
    derived().leave(x);
  }

  /* void apply(const pres_system::constant_mult& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    print_variables(x.parameters());
    derived().leave(x);
  } */

  void apply(const pres_system::fixpoint_symbol& x)
  {
    derived().enter(x);
    derived().print(x.is_mu() ? "mu" : "nu");
    derived().leave(x);
  }

  void apply(const pres_system::pres_equation& x)
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

  void apply(const pres_system::pres& x)
  {
    derived().enter(x);
    derived().apply(x.data());
    print_variables(x.global_variables(), true, true, true, "glob ", ";\n\n", ";\n     ");
    print_list(x.equations(), "pres ", "\n\n", "\n     ");
    derived().print("init ");
    print_pres_expression(x.initial_state());
    derived().print(";\n");
    derived().leave(x);
  }

  void apply(const pres_system::propositional_variable_instantiation& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    print_list(x.parameters(), "(", ")", ", ", false);
    derived().leave(x);
  }

  void apply(const pres_system::minus& x)
  {
    derived().enter(x);
    print_pres_unary_left_operation(x, " -");
    derived().leave(x);
  }

  void apply(const pres_system::and_& x)
  {
    derived().enter(x);
    print_pres_binary_operation(x, " && ");
    derived().leave(x);
  }

  void apply(const pres_system::or_& x)
  {
    derived().enter(x);
    print_pres_binary_operation(x, " || ");
    derived().leave(x);
  }

  void apply(const pres_system::imp& x)
  {
    derived().enter(x);
    print_pres_binary_operation(x, " => ");
    derived().leave(x);
  }

  void apply(const pres_system::plus& x)
  {
    derived().enter(x);
    print_pres_binary_operation(x, " + ");
    derived().leave(x);
  }

  void apply(const pres_system::const_multiply& x)
  {
    derived().enter(x);
    print_pres_binary_operation(x, "*");
    derived().leave(x);
  }

  void apply(const pres_system::const_multiply_alt& x)
  {
    derived().enter(x);
    print_pres_binary_operation(x, "*");
    derived().leave(x);
  }

  void apply(const pres_system::infimum& x)
  {
    derived().enter(x);
    print_pres_abstraction(x, "inf");
    derived().leave(x);
  }

  void apply(const pres_system::supremum& x)
  {
    derived().enter(x);
    print_pres_abstraction(x, "sup");
    derived().leave(x);
  }

  void apply(const pres_system::sum& x)
  {
    derived().enter(x);
    print_pres_abstraction(x, "sum");
    derived().leave(x);
  }

  void apply(const pres_system::condsm& x)
  {
    derived().enter(x);
    derived().print("(");
    apply(x.arg1());
    derived().print(" -> ");
    apply(x.arg2());
    derived().print(" <> ");
    apply(x.arg3());
    derived().print(")");
    derived().leave(x);
  }

  void apply(const pres_system::condeq& x)
  {
    derived().enter(x);
    derived().print("(");
    apply(x.arg1());
    derived().print(" => ");
    apply(x.arg2());
    derived().print(" <> ");
    apply(x.arg3());
    derived().print(")");
    derived().leave(x);
  }

  void apply(const pres_system::eqinf& x)
  {
    derived().enter(x);
    derived().print("eqinf(");
    apply(x.operand());
    derived().print(")");
    derived().leave(x);
  }

  void apply(const pres_system::eqninf& x)
  {
    derived().enter(x);
    derived().print("eqninf(");
    apply(x.operand());
    derived().print(")");
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
    core::detail::apply_printer<pres_system::detail::printer> printer(out, precedence_aware);
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

} // namespace pres_system

} // namespace mcrl2

#endif // MCRL2_PRES_PRINT_H
