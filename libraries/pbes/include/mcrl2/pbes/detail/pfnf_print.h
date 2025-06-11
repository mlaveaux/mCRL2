// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/detail/pfnf_print.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_DETAIL_PFNF_PRINT_H
#define MCRL2_PBES_DETAIL_PFNF_PRINT_H

#include "mcrl2/pbes/detail/is_pfnf.h"
#include "mcrl2/pbes/print.h"

namespace mcrl2 {

namespace pbes_system {

namespace detail {

template <typename Derived>
struct pfnf_printer: public pbes_system::detail::printer<Derived>
{
  typedef pbes_system::detail::printer<Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::derived;
  using super::print_abstraction;
  using super::print_list;
  using super::print_variables;

  bool is_abstraction(const pbes_system::pbes_expression& x)
  {
    return is_forall(x) || is_exists(x);
  }

  template <typename Abstraction>
  std::string abstraction_operator(const Abstraction& x) const
  {
    if (is_forall(x))
    {
      return "forall";
    }
    else if (is_exists(x))
    {
      return "exists";
    }
    else
    {
      throw mcrl2::runtime_error("error: unknown abstraction!");
    }
    return "";
  }

  template <typename Abstraction>
  pbes_expression abstraction_body(const Abstraction& x) const
  {
    if (is_forall(x))
    {
      return forall(x).body();
    }
    else if (is_exists(x))
    {
      return exists(x).body();
    }
    else
    {
      throw mcrl2::runtime_error("error: unknown abstraction!");
    }
    return pbes_expression();
  }

  template <typename Abstraction>
  data::variable_list abstraction_variables(const Abstraction& x) const
  {
    if (is_forall(x))
    {
      return forall(x).variables();
    }
    else if (is_exists(x))
    {
      return exists(x).variables();
    }
    else
    {
      throw mcrl2::runtime_error("error: unknown abstraction!");
    }
    return data::variable_list();
  }

  template <typename Abstraction>
  void print_pbes_abstraction(const Abstraction& x)
  {
    std::string op = abstraction_operator(x);
    derived().enter(x);
    derived().print(op + " ");
    print_variables(abstraction_variables(x), true, true, false, "", "", ", ");
    derived().print(". ");
    pbes_expression body = abstraction_body(x);
    if (is_abstraction(body))
    {
      print_pbes_abstraction(body);
    }
    else
    {
      std::vector<pbes_expression> implications = pfnf_implications(body);
      print_list(implications, "\n       (\n         ", "\n       )", "\n      && ", false);
    }
    derived().leave(x);
  }

  void apply(const pbes_system::pbes_expression& x)
  {
    derived().enter(x);
    if (is_abstraction(x))
    {
      print_pbes_abstraction(x);
    }
    else
    {
      super::apply(x);
    }
    derived().leave(x);
  }
};

} // namespace detail

/// \brief Prints the object x to a stream.
struct pfnf_stream_printer
{
  template <typename T>
  void operator()(const T& x, std::ostream& out, bool precedence_aware)
  {
    core::detail::apply_printer<pbes_system::detail::pfnf_printer> printer(out, precedence_aware);
    printer.apply(x);
  }
};

/// \brief Returns a PFNF string representation of the object x.
template <typename T>
std::string pfnf_pp(const T& x, bool precedence_aware = true)
{
  std::ostringstream out;
  pfnf_stream_printer()(x, out, precedence_aware);
  return out.str();
}

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_DETAIL_PFNF_PRINT_H
