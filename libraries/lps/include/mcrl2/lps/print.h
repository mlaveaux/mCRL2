// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/print.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_PRINT_H
#define MCRL2_LPS_PRINT_H

#include "mcrl2/data/parse.h"
#include "mcrl2/lps/traverser.h"
#include "mcrl2/process/print.h"

namespace mcrl2
{

namespace lps
{

namespace detail
{

template <typename Derived>
struct printer: public lps::add_traverser_sort_expressions<process::detail::printer, Derived>
{
  typedef lps::add_traverser_sort_expressions<process::detail::printer, Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::derived;
  using super::print_action_declarations;
  using super::print_assignments;
  using super::print_condition;
  using super::print_expression;
  using super::print_list;
  using super::print_variables;

  bool m_print_summand_numbers;

  printer()
    : m_print_summand_numbers(false)
  {}

  bool& print_summand_numbers()
  {
  	return m_print_summand_numbers;
  }

  template <typename Container>
  void print_numbered_list(const Container& container,
                           const std::string& separator = ", ",
                           const std::string& number_separator = "",
                           std::size_t index = 0,
                           bool print_start_separator = false,
                           bool print_empty_container = false
                          )
  {
    if (container.empty() && !print_empty_container)
    {
      return;
    }
    for (auto i = container.begin(); i != container.end(); ++i)
    {
      derived().print("\n");
      derived().print(number_separator);
      derived().print("%");
      derived().print(utilities::number2string(index++));

      derived().print("\n");
      if (i == container.begin() && !print_start_separator)
      {
        derived().print(number_separator);
      }
      else
      {
        derived().print(separator);
      }
      derived().apply(*i);
    }
  }

  void apply(const lps::deadlock& x)
  {
    derived().enter(x);
    derived().print("delta");
    if (x.has_time())
    {
      derived().print(" @ ");
      print_expression(x.time(), precedence(x.time()) < core::detail::max_precedence);
    }
    derived().leave(x);
  }

  void apply(const lps::multi_action& x)
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
    if (x.has_time())
    {
      derived().print(" @ ");
      print_expression(x.time(), precedence(x.time()) < core::detail::max_precedence);
    }
    derived().leave(x);
  }

  void apply(const lps::deadlock_summand& x)
  {
    derived().enter(x);
    print_variables(x.summation_variables(), true, true, false, "sum ", ".\n         ", ",");
    print_condition(x.condition(), " ->\n         ");
    derived().apply(x.deadlock());
    derived().leave(x);
  }

  void apply(const lps::stochastic_distribution& x)
  {
    derived().enter(x);
    if (x.variables().empty()) // do not print the empty distribution
    {
      return;
    }
    derived().print("dist ");
    print_variables(x.variables(), true, true, false, "", "");
    derived().print("[");
    derived().apply(x.distribution());
    derived().print("]");
    derived().leave(x);
  }

  void print_distribution(const lps::action_summand&)
  { }

  void print_distribution(const lps::stochastic_action_summand& x)
  {
    if (x.distribution().is_defined())
    {
      derived().apply(x.distribution());
      derived().print(" .\n         ");
    }
  }

  template <typename ActionSummand>
  void print_action_summand(const ActionSummand& x)
  {
    derived().enter(x);
    print_variables(x.summation_variables(), true, true, false, "sum ", ".\n         ", ",");
    print_condition(x.condition(), " ->\n         ");
    derived().apply(x.multi_action());
    derived().print(" .\n         ");
    print_distribution(x);
    derived().print("P(");
    print_assignments(x.assignments(), true, "", "", ", ");
    derived().print(")");
    derived().leave(x);
  }

  void apply(const lps::action_summand& x)
  {
    print_action_summand(x);
  }

  void apply(const lps::stochastic_action_summand& x)
  {
    print_action_summand(x);
  }

  void apply(const lps::process_initializer& x)
  {
    derived().enter(x);
    derived().print("init P");
    print_variables(x.expressions(), false);
    derived().print(";");
    derived().leave(x);
  }

  void apply(const lps::stochastic_process_initializer& x)
  {
    derived().enter(x);
    derived().print("init ");
    if (x.distribution().is_defined())
    {
      derived().apply(x.distribution());
      derived().print(" . ");
    }
    derived().print("P");
    print_variables(x.expressions(), false);
    derived().print(";");
    derived().leave(x);
  }

  // this overload is enabled for linear_process and stochastic_linear_process
  template <typename LinearProcess>
  void print_linear_process(const LinearProcess& x)
  {
    derived().enter(x);
    derived().print("proc P");
    print_variables(x.process_parameters(), true, true, false, "(", ")", ", ");

    if (m_print_summand_numbers)
    {
      derived().print(" =");

      std::string separator        = "     + ";
      std::string number_separator = "       ";

      // print action summands
      print_numbered_list(x.action_summands(), separator, number_separator, 1, false);

      // print deadlock summands
      print_numbered_list(x.deadlock_summands(), separator, number_separator, x.action_summands().size() + 1, true);

      // print delta if there are no summands
      if (x.action_summands().empty() && (x.deadlock_summands().empty()))
      {
        deadlock_summand_vector v;
        v.push_back(deadlock_summand(data::variable_list(), data::sort_bool::true_(), lps::deadlock(data::parse_data_expression("0"))));
        print_numbered_list(v, separator, number_separator, 1, true);
      }
    }
    else
    {
      derived().print(" =\n       ");

      // print action summands
      std::string opener = "";
      std::string closer = "";
      std::string separator = "\n     + ";
      print_list(x.action_summands(), opener, closer, separator);

      // print deadlock summands
      if (!x.action_summands().empty())
      {
        opener = separator;
      }
      print_list(x.deadlock_summands(), opener, closer, separator);

      // print delta if there are no summands
      if (x.action_summands().empty() && (x.deadlock_summands().empty()))
      {
        deadlock_summand_vector v;
        v.push_back(deadlock_summand(data::variable_list(), data::sort_bool::true_(), lps::deadlock(data::parse_data_expression("0"))));
        print_list(v, opener, closer, separator);
      }
    }

    derived().print(";\n");
    derived().leave(x);
  }

  void apply(const linear_process& x)
  {
    print_linear_process(x);
  }

  void apply(const stochastic_linear_process& x)
  {
    print_linear_process(x);
  }

  template <typename Specification>
  void print_specification(const Specification& x)
  {
    derived().enter(x);
    derived().apply(x.data());
    print_action_declarations(x.action_labels(), "act  ",";\n\n", ";\n     ");
    print_variables(x.global_variables(), true, true, true, "glob ", ";\n\n", ";\n     ");
    derived().apply(x.process());
    derived().print("\n");
    derived().apply(x.initial_process());
    derived().print("\n");
    derived().leave(x);
  }

  void apply(const specification& x)
  {
    print_specification(x);
  }

  void apply(const stochastic_specification& x)
  {
    print_specification(x);
  }

  /* void apply(const lps::state &x)
  {
    derived().enter(x);
    derived().print("state(");
    bool first = true;
    for (lps::state::const_iterator i = x.begin(); i != x.end(); i++)
    {
      if (!first)
      {
        derived().print(", ");
      }
      first = false;
      print_expression(*i);
    }
    derived().print(")");
    derived().leave(x);
  } */
};

} // namespace detail

/// \brief Prints the object x to a stream.
struct stream_printer
{
  template <typename T>
  void operator()(const T& x, std::ostream& out, bool precedence_aware)
  {
    core::detail::apply_printer<lps::detail::printer> printer(out, precedence_aware);
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

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_PRINT_H
