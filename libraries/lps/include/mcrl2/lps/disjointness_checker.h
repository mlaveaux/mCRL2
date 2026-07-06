// Author(s): Luc Engelen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/disjointness_checker.h
/// \brief Checks disjointness of read/write parameters of summands in a linear process.

// Interface to class Disjointness_Checker
// file: disjointness_checker.h

#ifndef MCRL2_LPS_DISJOINTNESS_CHECKER_H
#define MCRL2_LPS_DISJOINTNESS_CHECKER_H

#include "mcrl2/data/assignment.h"
#include "mcrl2/lps/linear_process.h"

/// \brief Class that can determine if two summands are syntactically disjoint.
/// Two summands are syntactically disjoint if the following conditions hold:
/// - The set of variables used by one summand is disjoint from the set of variables changed by the other summand and
///   vice versa.
/// - The set of variables changed by one summand is disjoint from the set of variables changed by the other summand.
///
/// An instance of the class Disjointness_Checker is created using the constructor
/// Disjointness_Checker::Disjointness_Checker. The parameter a_process_equations is used to pass the summands to be
/// checked for disjointness. The function Disjointness_Checker::disjoint indicates whether the two summands with numbers
/// n_1 and n_2 are syntactically disjoint.

namespace mcrl2::lps::detail
{

class Disjointness_Checker
{
  private:
    /// \brief The number of summands of the LPS passed as argument of the constructor.
    std::size_t f_number_of_summands;

    /// \brief A two dimensional array, indicating which parameters a summand uses, for each of the summands.
    std::vector<std::set<data::variable> > f_used_parameters_per_summand;

    /// \brief A two dimensional array, indicating which parameters a summand changes, for each of the summands.
    std::vector<std::set<data::variable> > f_changed_parameters_per_summand;

    /// \brief Updates the array Disjointness_Checker::f_used_parameters_per_summand, given the expression a_expression.
    void process_data_expression(std::size_t n, const data::data_expression& x);

    /// \brief Updates the array Disjointness_Checker::f_used_parameters_per_summand, given the multiaction a.
    void process_multi_action(std::size_t n, const multi_action& a);

    /// \brief Updates the arrays Disjointness_Checker::f_changed_parameters_per_summand and
    /// \brief Disjointness_Checker::f_used_parameters_per_summand, given the summand s.
    void process_summand(std::size_t n, const action_summand& s);

  public:
    /// \brief Constructor that initializes the sets Disjointness_Checker::f_used_parameters_per_summand and
    /// \brief Disjointness_Checker::f_changed_parameters_per_summand, and the indexed set
    /// \brief Disjointness_Checker::f_parameter_set.
    /// precondition: the argument passed as parameter a_process_equations is a specification of process equations in mCRL2
    /// format
    /// precondition: the arguments passed as parameters n_1 and n_2 correspond to summands in
    /// the proces equations passed as parameter a_process_equations. They lie in the interval from and including 1 upto and
    /// including the highest summand number
    template <typename LinearProcess>
    Disjointness_Checker(const LinearProcess& a_process_equation);

    /// \brief Indicates whether or not the summands with number n_1 and n_2 are disjoint.
    bool disjoint(std::size_t n1, std::size_t n2);
};

inline
void Disjointness_Checker::process_data_expression(std::size_t n, const data::data_expression& x)
{
  // This should probably once be replaced by a visitor.
  if (data::is_machine_number(x))
  {
    // Do nothing.
  }
  else if (data::is_variable(x))
  {
    f_used_parameters_per_summand[n].insert(data::variable(x));
  }
  else if (data::is_where_clause(x))
  {
    const data::where_clause& t = atermpp::down_cast<data::where_clause>(x);
    process_data_expression(n, t.body());
    const data::assignment_list& assignments = t.assignments();
    for (const data::assignment& assignment : assignments)
    {
      process_data_expression(n, assignment.rhs());
    }
  }
  else if (data::is_function_symbol(x))
  {
    // do nothing
  }
  else if (data::is_application(x))
  {
    const data::application& t = atermpp::down_cast<data::application>(x);
    process_data_expression(n, t.head());
    for (const auto& i : t)
    {
      process_data_expression(n, i);
    }
  }
  else if (data::is_abstraction(x))
  {
    const data::abstraction& t = atermpp::down_cast<data::abstraction>(x);
    process_data_expression(n, t.body());
  }
  else
  {
    throw mcrl2::runtime_error("disjointness checker encountered unknown term " + data::pp(x));
  }
}

inline
void Disjointness_Checker::process_multi_action(std::size_t n, const multi_action& a)
{
  if (a.has_time())
  {
    process_data_expression(n, a.time());
  }

  const process::action_list& v_actions = a.actions();
  for (const auto & v_action : v_actions)
  {
    const data::data_expression_list v_expressions=v_action.arguments();

    for (const auto & v_expression : v_expressions)
    {
      process_data_expression(n, v_expression);
    }
  }
}

inline
void Disjointness_Checker::process_summand(std::size_t n, const action_summand& s)
{
  // variables used in condition
  process_data_expression(n, s.condition());

  // variables used in multiaction
  process_multi_action(n, s.multi_action());

  // variables used and changed in assignments
  const data::assignment_list& v_assignments = s.assignments();
  for (const auto & v_assignment : v_assignments)
  {
    // variables changed in the assignment
    f_changed_parameters_per_summand[n].insert(v_assignment.lhs());

    // variables used in assignment
    process_data_expression(n, v_assignment.rhs());
  }
}

template <typename LinearProcess>
Disjointness_Checker::Disjointness_Checker(const LinearProcess& a_process_equation)
{
  const auto& v_summands = a_process_equation.action_summands();
  std::size_t v_summand_number = 1;

  f_number_of_summands = v_summands.size();
  f_used_parameters_per_summand =
    std::vector<std::set<data::variable> >(f_number_of_summands + 1, std::set<data::variable>());
  f_changed_parameters_per_summand =
    std::vector<std::set<data::variable> >(f_number_of_summands + 1, std::set<data::variable>());

  for (const auto & v_summand : v_summands)
  {
    process_summand(v_summand_number, v_summand);
    v_summand_number++;
  }
}

inline
bool Disjointness_Checker::disjoint(std::size_t n1, std::size_t n2)
{
  using utilities::detail::has_empty_intersection;
  assert(n1 <= f_number_of_summands && n2 <= f_number_of_summands);
  bool v_used_1_changed_2 = has_empty_intersection(f_used_parameters_per_summand[n1], f_changed_parameters_per_summand[n2]);
  bool v_used_2_changed_1 = has_empty_intersection(f_used_parameters_per_summand[n2], f_changed_parameters_per_summand[n1]);
  bool v_changed_1_changed_2 = has_empty_intersection(f_changed_parameters_per_summand[n1], f_changed_parameters_per_summand[n2]);
  return v_used_1_changed_2 && v_used_2_changed_1 && v_changed_1_changed_2;
}

} // namespace mcrl2::lps::detail

#endif
