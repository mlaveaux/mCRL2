// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/print.h
/// \brief Functions for pretty printing ATerms.

#ifndef MCRL2_CORE_PRINT_H
#define MCRL2_CORE_PRINT_H

#include <cctype>
#include "mcrl2/core/print_format.h"
#include "mcrl2/core/traverser.h"

namespace mcrl2
{
namespace core
{

/// \cond INTERNAL_DOCS
namespace detail
{

const int max_precedence = 10000;

template <typename T>
int precedence(const T&)
{
  return max_precedence;
}

template <typename Derived>
struct printer: public core::traverser<Derived>
{
  typedef core::traverser<Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;

  std::ostream* m_out;
  bool m_precedence_aware {true}; /// \brief Indicates whether the printer should take the precedence of operators into account.

  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }

  std::ostream& out()
  {
    return *m_out;
  }

  void print(const std::string& s)
  {
    out() << s;
  }

  template <typename T>
  void print_expression(const T& x, bool needs_parentheses)
  {
    if (!m_precedence_aware || needs_parentheses)
    {
      derived().print("(");
    }
    derived().apply(x);
    if (!m_precedence_aware || needs_parentheses)
    {
      derived().print(")");
    }
  }

  template <typename T, typename U>
  void print_unary_operand(const T& x, const U& operand)
  {
    print_expression(operand, !m_precedence_aware || precedence(operand) < precedence(x));
  }

  template <typename T>
  void print_unary_left_operation(const T& x, const std::string& op)
  {
    derived().print(op);
    print_unary_operand(x, x.operand());
  }

  template <typename T>
  void print_unary_right_operation(const T& x, const std::string& op)
  {
    print_unary_operand(x, x.operand());
    derived().print(op);
  }

  template <typename T>
  void print_binary_operation(const T& x, const std::string& op)
  {
    const auto& x1 = x.left();
    const auto& x2 = x.right();
    auto p = precedence(x);
    auto p1 = precedence(x1);
    auto p2 = precedence(x2);
    print_expression(x1, !m_precedence_aware || (p1 < p) || (p1 == p && !is_left_associative(x)));
    derived().print(op);
    print_expression(x2, !m_precedence_aware || (p2 < p) || (p2 == p && !is_right_associative(x)));
  }

  template <typename Container>
  void print_list(const Container& container,
                  const std::string& opener = "(",
                  const std::string& closer = ")",
                  const std::string& separator = ", ",
                  bool print_empty_container = false
                 )
  {
    if (container.empty() && !print_empty_container)
    {
      return;
    }
    derived().print(opener);
    for (auto i = container.begin(); i != container.end(); ++i)
    {
      if (i != container.begin())
      {
        derived().print(separator);
      }
      derived().apply(*i);
    }
    derived().print(closer);
  }

  template <typename T>
  void apply(const std::list<T>& x)
  {
    derived().enter(x);
    print_list(x, "", "", ", ");
    derived().leave(x);
  }

  template <typename T>
  void apply(const atermpp::term_list<T>& x)
  {
    derived().enter(x);
    print_list(x, "", "", ", ");
    derived().leave(x);
  }

  template <typename T>
  void apply(const std::set<T>& x)
  {
    derived().enter(x);
    print_list(x, "", "", ", ");
    derived().leave(x);
  }

  void apply(const core::identifier_string& x)
  {
    derived().enter(x);
    if (x == core::identifier_string())
    {
      derived().print("@NoValue");
    }
    else
    {
      derived().print(std::string(x));
    }
    derived().leave(x);
  }

  void apply(const atermpp::aterm_list& x)
  {
    derived().enter(x);
    derived().print(utilities::to_string(x));
    derived().leave(x);
  }

  void apply(const atermpp::aterm& x)
  {
    derived().enter(x);
    derived().print(utilities::to_string(x));
    derived().leave(x);
  }

  void apply(const atermpp::aterm_int& x)
  {
    derived().enter(x);
    derived().print(utilities::to_string(x));
    derived().leave(x);
  }
};

template <template <class> class Traverser>
struct apply_printer: public Traverser<apply_printer<Traverser>>
{
  typedef Traverser<apply_printer<Traverser>> super;

  using super::enter;
  using super::leave;
  using super::apply;

  /// \brief precedence_aware Indicates whether the printer should take precendence into account when printing expressions.
  apply_printer(std::ostream& out, bool precedence_aware)
  {
    typedef printer<apply_printer<Traverser> > Super;
    static_cast<Super&>(*this).m_out = &out;
    static_cast<Super&>(*this).m_precedence_aware = precedence_aware;
  }

};

} // namespace detail
/// \endcond

/// \brief Prints the object x to a stream.
struct stream_printer
{
  template <typename T>
  void operator()(const T& x, std::ostream& out, bool precedence_aware)
  {
    core::detail::apply_printer<core::detail::printer> printer(out, precedence_aware);
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

} // namespace core

} // namespace mcrl2

#endif // MCRL2_CORE_PRINT_H
