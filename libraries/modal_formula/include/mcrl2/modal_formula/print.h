// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/print.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_PRINT_H
#define MCRL2_MODAL_FORMULA_PRINT_H

#include "mcrl2/lps/print.h"
#include "mcrl2/modal_formula/traverser.h"

namespace mcrl2 {

namespace action_formulas {

constexpr inline int precedence(const forall&) { return 21; }
constexpr inline int precedence(const exists&) { return 21; }
constexpr inline int precedence(const imp&)    { return 22; }
constexpr inline int precedence(const or_&)    { return 23; }
constexpr inline int precedence(const and_&)   { return 24; }
constexpr inline int precedence(const at&)     { return 25; }
constexpr inline int precedence(const not_&)   { return 26; }
inline int precedence(const action_formula& x)
{
  if (is_forall(x))      { return precedence(atermpp::down_cast<forall>(x)); }
  else if (is_exists(x)) { return precedence(atermpp::down_cast<exists>(x)); }
  else if (is_imp(x))    { return precedence(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))     { return precedence(atermpp::down_cast<or_>(x)); }
  else if (is_and(x))    { return precedence(atermpp::down_cast<and_>(x)); }
  else if (is_at(x))     { return precedence(atermpp::down_cast<at>(x)); }
  else if (is_not(x))    { return precedence(atermpp::down_cast<not_>(x)); }
  return core::detail::max_precedence;
}

// only defined for binary operators
inline bool is_left_associative(const imp&)  { return false; }
inline bool is_left_associative(const or_&)  { return true; }
inline bool is_left_associative(const and_&) { return true; }
inline bool is_left_associative(const action_formula& x)
{
  if (is_imp(x))      { return is_left_associative(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))  { return is_left_associative(atermpp::down_cast<or_>(x)); }
  else if (is_and(x)) { return is_left_associative(atermpp::down_cast<and_>(x)); }
  return false;
}

inline bool is_right_associative(const imp&)  { return true; }
inline bool is_right_associative(const or_&)  { return true; }
inline bool is_right_associative(const and_&) { return true; }
inline bool is_right_associative(const action_formula& x)
{
  if (is_imp(x))      { return is_right_associative(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))  { return is_right_associative(atermpp::down_cast<or_>(x)); }
  else if (is_and(x)) { return is_right_associative(atermpp::down_cast<and_>(x)); }
  return false;
}

namespace detail
{

template <typename Derived>
struct printer: public action_formulas::add_traverser_sort_expressions<lps::detail::printer, Derived>
{
  typedef action_formulas::add_traverser_sort_expressions<lps::detail::printer, Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::derived;
  using super::print_abstraction;
  using super::print_list;
  using super::print_action_declarations;
  using super::print_expression;
  using super::print_unary_left_operation;
  using super::print_binary_operation;

  void apply(const action_formulas::true_& x)
  {
    derived().enter(x);
    derived().print("true");
    derived().leave(x);
  }

  void apply(const action_formulas::false_& x)
  {
    derived().enter(x);
    derived().print("false");
    derived().leave(x);
  }

  void apply(const action_formulas::not_& x)
  {
    derived().enter(x);
    print_unary_left_operation(x, "!");
    derived().leave(x);
  }

  void apply(const action_formulas::and_& x)
  {
    derived().enter(x);
    print_binary_operation(x, " && ");
    derived().leave(x);
  }

  void apply(const action_formulas::or_& x)
  {
    derived().enter(x);
    print_binary_operation(x, " || ");
    derived().leave(x);
  }

  void apply(const action_formulas::imp& x)
  {
    derived().enter(x);
    print_binary_operation(x, " => ");
    derived().leave(x);
  }

  void apply(const action_formulas::forall& x)
  {
    derived().enter(x);
    print_abstraction(x, "forall");
    derived().leave(x);
  }

  void apply(const action_formulas::exists& x)
  {
    derived().enter(x);
    print_abstraction(x, "exists");
    derived().leave(x);
  }

  void apply(const action_formulas::at& x)
  {
    derived().enter(x);
    derived().apply(x.operand());
    derived().print(" @ ");
    print_expression(x.time_stamp(), precedence(x.time_stamp()) < core::detail::max_precedence);
    derived().leave(x);
  }

  void apply(const action_formulas::multi_action& x)
  {
    derived().enter(x);
    if (x.actions().empty())
    {
      derived().print("tau");
    }
    else
    {
      print_list(x.actions(), "", "", "|");
    }
    derived().leave(x);
  }
};

} // namespace detail

/// \brief Prints the object t to a stream.
template <typename T>
void pp(const T& t, std::ostream& out, bool precendence_aware)
{
  core::detail::apply_printer<action_formulas::detail::printer> printer(out, precendence_aware);
  printer.apply(t);
}

/// \brief Returns a string representation of the object t.
template <typename T>
std::string pp(const T& t, bool precendence_aware = true)
{
  std::ostringstream out;
  action_formulas::pp(t, out, precendence_aware);
  return out.str();
}

} // namespace action_formulas

namespace regular_formulas {

constexpr inline int precedence(const seq&)          { return 31; }
constexpr inline int precedence(const alt&)          { return 32; }
constexpr inline int precedence(const trans&)        { return 33; }
constexpr inline int precedence(const trans_or_nil&) { return 33; }
inline int precedence(const regular_formula& x)
{
  if (is_seq(x))               { return precedence(atermpp::down_cast<seq>(x)); }
  else if (is_alt(x))          { return precedence(atermpp::down_cast<alt>(x)); }
  else if (is_trans(x))        { return precedence(atermpp::down_cast<trans>(x)); }
  else if (is_trans_or_nil(x)) { return precedence(atermpp::down_cast<trans_or_nil>(x)); }
  return core::detail::max_precedence;
}

// only defined for binary operators
inline bool is_left_associative(const seq&)  { return false; }
inline bool is_left_associative(const alt&)  { return true; }
inline bool is_left_associative(const regular_formula& x)
{
  if (is_seq(x))      { return is_left_associative(atermpp::down_cast<seq>(x)); }
  else if (is_alt(x)) { return is_left_associative(atermpp::down_cast<alt>(x)); }
  return false;
}

inline bool is_right_associative(const seq&)  { return true; }
inline bool is_right_associative(const alt&)  { return true; }
inline bool is_right_associative(const regular_formula& x)
{
  if (is_seq(x))      { return is_right_associative(atermpp::down_cast<seq>(x)); }
  else if (is_alt(x)) { return is_right_associative(atermpp::down_cast<alt>(x)); }
  return false;
}

namespace detail
{

template <typename Derived>
struct printer: public regular_formulas::add_traverser_sort_expressions<action_formulas::detail::printer, Derived>
{
  typedef regular_formulas::add_traverser_sort_expressions<action_formulas::detail::printer, Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::derived;
  using super::print_action_declarations;

  using super::print_expression;
  using super::print_unary_left_operation;
  using super::print_unary_right_operation;
  using super::print_binary_operation;

  void apply(const regular_formulas::seq& x)
  {
    derived().enter(x);
    print_binary_operation(x, " . ");
    derived().leave(x);
  }

  void apply(const regular_formulas::alt& x)
  {
    derived().enter(x);
    print_binary_operation(x, " + ");
    derived().leave(x);
  }

  void apply(const regular_formulas::trans& x)
  {
    derived().enter(x);
    print_unary_right_operation(x, "+");
    derived().leave(x);
  }

  void apply(const regular_formulas::trans_or_nil& x)
  {
    derived().enter(x);
    print_unary_right_operation(x, "*");
    derived().leave(x);
  }

  void apply(const regular_formulas::untyped_regular_formula& x)
  {
    derived().enter(x);
    print_expression(x.left(), false);
    derived().print(" " + std::string(x.name()) + " ");
    print_expression(x.right(), false);
    derived().leave(x);
  }
};

} // namespace detail

/// \brief Prints the object t to a stream.
template <typename T>
void pp(const T& t, std::ostream& out, bool precendence_aware)
{
  core::detail::apply_printer<regular_formulas::detail::printer> printer(out, precendence_aware);
  printer.apply(t);
}

/// \brief Returns a string representation of the object t.
template <typename T>
std::string pp(const T& t, bool precendence_aware = true)
{
  std::ostringstream out;
  regular_formulas::pp(t, out, precendence_aware);
  return out.str();
}

} // namespace regular_formulas

namespace state_formulas {

constexpr inline int precedence(const mu&)                 { return 41; }
constexpr inline int precedence(const nu&)                 { return 41; }
constexpr inline int precedence(const forall&)             { return 42; }
constexpr inline int precedence(const exists&)             { return 42; }
constexpr inline int precedence(const infimum&)            { return 42; }
constexpr inline int precedence(const supremum&)           { return 42; }
constexpr inline int precedence(const sum&)                { return 42; }
constexpr inline int precedence(const imp&)                { return 45; }
constexpr inline int precedence(const or_&)                { return 46; }
constexpr inline int precedence(const and_&)               { return 47; }
constexpr inline int precedence(const plus&)               { return 43; }
constexpr inline int precedence(const const_multiply&)     { return 44; }
constexpr inline int precedence(const const_multiply_alt&) { return 44; }
constexpr inline int precedence(const must&)               { return 48; }
constexpr inline int precedence(const may&)                { return 48; }
constexpr inline int precedence(const not_&)               { return 50; }
constexpr inline int precedence(const minus&)              { return 50; }
inline int precedence(const state_formula& x)
{
  if      (is_mu(x))                 { return precedence(atermpp::down_cast<mu>(x)); }
  else if (is_nu(x))                 { return precedence(atermpp::down_cast<nu>(x)); }
  else if (is_forall(x))             { return precedence(atermpp::down_cast<forall>(x)); }
  else if (is_exists(x))             { return precedence(atermpp::down_cast<exists>(x)); }
  else if (is_infimum(x))            { return precedence(atermpp::down_cast<infimum>(x)); }
  else if (is_supremum(x))           { return precedence(atermpp::down_cast<supremum>(x)); }
  else if (is_sum(x))                { return precedence(atermpp::down_cast<sum>(x)); }
  else if (is_imp(x))                { return precedence(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))                 { return precedence(atermpp::down_cast<or_>(x)); }
  else if (is_and(x))                { return precedence(atermpp::down_cast<and_>(x)); }
  else if (is_plus(x))               { return precedence(atermpp::down_cast<plus>(x)); }
  else if (is_const_multiply(x))     { return precedence(atermpp::down_cast<const_multiply>(x)); }
  else if (is_const_multiply_alt(x)) { return precedence(atermpp::down_cast<const_multiply_alt>(x)); }
  else if (is_must(x))               { return precedence(atermpp::down_cast<must>(x)); }
  else if (is_may(x))                { return precedence(atermpp::down_cast<may>(x)); }
  else if (is_not(x))                { return precedence(atermpp::down_cast<not_>(x)); }
  else if (is_minus(x))                { return precedence(atermpp::down_cast<minus>(x)); }
  return core::detail::max_precedence;
}

// only defined for binary operators
inline bool is_left_associative(const imp&)  { return false; }
inline bool is_left_associative(const or_&)  { return true; }
inline bool is_left_associative(const and_&) { return true; }
inline bool is_left_associative(const plus&) { return true; }
inline bool is_left_associative(const const_multiply&) { return true; }
inline bool is_left_associative(const const_multiply_alt&) { return true; }
inline bool is_left_associative(const state_formula& x)
{
  if (is_imp(x))      { return is_left_associative(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))  { return is_left_associative(atermpp::down_cast<or_>(x)); }
  else if (is_and(x)) { return is_left_associative(atermpp::down_cast<and_>(x)); }
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
inline bool is_right_associative(const state_formula& x)
{
  if (is_imp(x))      { return is_right_associative(atermpp::down_cast<imp>(x)); }
  else if (is_or(x))  { return is_right_associative(atermpp::down_cast<or_>(x)); }
  else if (is_and(x)) { return is_right_associative(atermpp::down_cast<and_>(x)); }
  else if (is_plus(x)) { return is_right_associative(atermpp::down_cast<plus>(x)); }
  else if (is_const_multiply(x)) { return is_right_associative(atermpp::down_cast<const_multiply>(x)); }
  else if (is_const_multiply_alt(x)) { return is_right_associative(atermpp::down_cast<const_multiply_alt>(x)); }
  return false;
}

namespace detail
{

template <typename Derived>
struct printer: public state_formulas::add_traverser_sort_expressions<regular_formulas::detail::printer, Derived>
{
  typedef state_formulas::add_traverser_sort_expressions<regular_formulas::detail::printer, Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::derived;
  using super::print_abstraction;
  using super::print_variables;
  using super::print_list;
  using super::print_action_declarations;

  using super::print_expression;
  using super::print_unary_left_operation;
  using super::print_binary_operation;
  using super::print_unary_operand;

  // Determines whether or not data expressions should be wrapped inside 'val'.
  std::vector<bool> val;

  void disable_val()
  {
    val.push_back(false);
  }

  void enable_val()
  {
    assert(!val.empty());
    val.pop_back();
  }

  void apply(const data::data_expression& x)
  {
    bool print_val = val.empty();
    derived().enter(x);
    if (print_val)
    {
      disable_val();
      derived().print("val(");
    }
    super::apply(x);
    if (print_val)
    {
      derived().print(")");
      enable_val();
    }
    derived().leave(x);
  }

  void apply(const state_formulas::true_& x)
  {
    derived().enter(x);
    derived().print("true");
    derived().leave(x);
  }

  void apply(const state_formulas::false_& x)
  {
    derived().enter(x);
    derived().print("false");
    derived().leave(x);
  }

  void apply(const state_formulas::not_& x)
  {
    derived().enter(x);
    print_unary_left_operation(x, "!");
    derived().leave(x);
  }

  void apply(const state_formulas::minus& x)
  {
    derived().enter(x);
    print_unary_left_operation(x, "-");
    derived().leave(x);
  }

  void apply(const state_formulas::and_& x)
  {
    derived().enter(x);
    print_binary_operation(x, " && ");
    derived().leave(x);
  }

  void apply(const state_formulas::or_& x)
  {
    derived().enter(x);
    print_binary_operation(x, " || ");
    derived().leave(x);
  }

  void apply(const state_formulas::imp& x)
  {
    derived().enter(x);
    print_binary_operation(x, " => ");
    derived().leave(x);
  }

  void apply(const state_formulas::plus& x)
  {
    derived().enter(x);
    print_binary_operation(x, " + ");
    derived().leave(x);
  }

  void apply(const state_formulas::const_multiply& x)
  {
    derived().enter(x);
    derived().apply(x.left());
    derived().print("*"); 
    derived().apply(x.right());
    derived().leave(x);
  }

  void apply(const state_formulas::const_multiply_alt& x)
  {
    derived().enter(x);
    derived().apply(x.left());
    derived().print("*"); 
    derived().apply(x.right());
    derived().leave(x);
  }

  void apply(const state_formulas::forall& x)
  {
    derived().enter(x);
    print_abstraction(x, "forall");
    derived().leave(x);
  }

  void apply(const state_formulas::exists& x)
  {
    derived().enter(x);
    print_abstraction(x, "exists");
    derived().leave(x);
  }

  void apply(const state_formulas::infimum& x)
  {
    derived().enter(x);
    print_abstraction(x, "inf");
    derived().leave(x);
  }

  void apply(const state_formulas::supremum& x)
  {
    derived().enter(x);
    print_abstraction(x, "sup");
    derived().leave(x);
  }

  void apply(const state_formulas::sum& x)
  {
    derived().enter(x);
    print_abstraction(x, "sum");
    derived().leave(x);
  }

  void apply(const state_formulas::must& x)
  {
    derived().enter(x);
    derived().print("[");
    disable_val();
    derived().apply(x.formula());
    enable_val();
    derived().print("]");
    print_unary_operand(x, x.operand());
    derived().leave(x);
  }

  void apply(const state_formulas::may& x)
  {
    derived().enter(x);
    derived().print("<");
    disable_val();
    derived().apply(x.formula());
    enable_val();
    derived().print(">");
    print_unary_operand(x, x.operand());
    derived().leave(x);
  }

  void apply(const state_formulas::yaled& x)
  {
    derived().enter(x);
    derived().print("yaled");
    derived().leave(x);
  }

  void apply(const state_formulas::yaled_timed& x)
  {
    disable_val();
    derived().enter(x);
    derived().print("yaled");
    derived().print(" @ ");
    print_expression(x.time_stamp(), precedence(x.time_stamp()) < core::detail::max_precedence);
    derived().leave(x);
    enable_val();
  }

  void apply(const state_formulas::delay& x)
  {
    derived().enter(x);
    derived().print("delay");
    derived().leave(x);
  }

  void apply(const state_formulas::delay_timed& x)
  {
    disable_val();
    derived().enter(x);
    derived().print("delay");
    derived().print(" @ ");
    print_expression(x.time_stamp(), precedence(x.time_stamp()) < core::detail::max_precedence);
    derived().leave(x);
    enable_val();
  }

  void apply(const state_formulas::variable& x)
  {
    disable_val();
    derived().enter(x);
    derived().apply(x.name());
    print_list(x.arguments(), "(", ")", ", ", false);
    derived().leave(x);
    enable_val();
  }

  void apply(const data::untyped_data_parameter& x)
  {
    disable_val();
    derived().enter(x);
    derived().apply(x.name());
    print_list(x.arguments(), "(", ")", ", ", false);
    derived().leave(x);
    enable_val();
  }

  // TODO: merge this function with the version in data/print.h (?)
  void print_assignments(const data::assignment_list& assignments)
  {
    if (assignments.empty())
    {
      return;
    }
    disable_val();
    derived().print("(");
    for (auto i = assignments.begin(); i != assignments.end(); ++i)
    {
      if (i != assignments.begin())
      {
        derived().print(", ");
      }
      derived().apply(i->lhs());
      derived().print(": ");
      derived().apply(i->lhs().sort());
      derived().print(" = ");
      derived().apply(i->rhs());
    }
    derived().print(")");
    enable_val();
  }

  void apply(const state_formulas::nu& x)
  {
    derived().enter(x);
    derived().print("nu ");
    derived().apply(x.name());
    print_assignments(x.assignments());
    derived().print(". ");
    derived().apply(x.operand());
    derived().leave(x);
  }

  void apply(const state_formulas::mu& x)
  {
    derived().enter(x);
    derived().print("mu ");
    derived().apply(x.name());
    print_assignments(x.assignments());
    derived().print(". ");
    derived().apply(x.operand());
    derived().leave(x);
  }

  void apply(const state_formulas::state_formula_specification& x)
  {
    derived().enter(x);
    derived().apply(x.data());
    print_action_declarations(x.action_labels(), "act  ",";\n\n", ";\n     ");
    derived().print("form ");
    derived().apply(x.formula());
    derived().print(";\n");
    derived().leave(x);
  }
};

} // namespace detail

/// \brief Prints the object t to a stream.
template <typename T>
void pp(const T& t, std::ostream& out, bool precendence_aware)
{
  core::detail::apply_printer<state_formulas::detail::printer> printer(out, precendence_aware);
  printer.apply(t);
}

/// \brief Returns a string representation of the object t.
template <typename T>
std::string pp(const T& t, bool precendence_aware = true)
{
  std::ostringstream out;
  state_formulas::pp(t, out, precendence_aware);
  return out.str();
}

} // namespace state_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULA_PRINT_H

