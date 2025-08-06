// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/is_well_typed.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_IS_WELL_TYPED_H
#define MCRL2_LPS_IS_WELL_TYPED_H

#include <iostream>

#include "mcrl2/data/detail/sequence_algorithm.h"
#include "mcrl2/lps/detail/action_utility.h"
#include "mcrl2/lps/stochastic_specification.h"
#include <boost/iterator/transform_iterator.hpp>

namespace mcrl2::lps::detail
{

/// \brief Function object for applying a substitution to LPS data types.
struct lps_well_typed_checker
{
  // The result of the last well typedness check.
  bool result = false;

  // Error message are written to the output stream error.
  mutable std::ostringstream error;

  /// \brief Checks if the sort of t has type real
  bool check_time(const data::data_expression& t, const std::string& type) const
  {
    if (!data::sort_real::is_real(t.sort()))
    {
      error << "is_well_typed(" << type << ") failed: time " << t << " doesn't have sort real." << std::endl;
      return false;
    }
    return true;
  }

  /// \brief Checks if the sort of t has type bool
  bool check_condition(const data::data_expression& t, const std::string& type) const
  {
    if (!data::sort_bool::is_bool(t.sort()))
    {
      error << "is_well_typed(" << type << ") failed: condition " << t << " doesn't have sort bool." << std::endl;
      return false;
    }
    return true;
  }

  /// \brief Checks if the assignments are well typed and have unique left hand sides
  bool check_assignments(const data::assignment_list& l, const std::string& type) const
  {
    if (!is_well_typed_container(l))
    {
      error << "is_well_typed(" << type << ") failed: the assignments " << l << " are not well typed." << std::endl;
      return false;
    }
    auto lhs = [](const data::assignment& a) { return a.lhs(); };
    if (data::detail::sequence_contains_duplicates(
          boost::make_transform_iterator(l.begin(), lhs),
          boost::make_transform_iterator(l.end()  , lhs)
        )
       )
    {
      error << "is_well_typed(" << type << ") failed: data assignments " << l << " don't have unique left hand sides." << std::endl;
      return false;
    }
    return true;
  }

  /// \brief Checks well typedness of the elements of a container
  template <typename Container>
  bool is_well_typed_container(const Container& c) const
  {
    for (auto i = c.begin(); i != c.end(); ++i)
    {
      if (!is_well_typed(*i))
      {
        return false;
      }
    }
    return true;
  }

  /// \brief Checks well typedness of a sort expression.
  /// \param d A sort expression.
  bool is_well_typed(const data::sort_expression& d) const
  {
    (void)d; // Suppress an unused variable warning.
    return true;
  }

  /// \brief Checks well typedness of a variable.
  /// \param d A variable.
  bool is_well_typed(const data::variable& d) const
  {
    (void)d; // Suppress an unused variable warning.
    return true;
  }

  /// \brief Checks well typedness of a data expression
  /// \param d A data expression
  bool is_well_typed(const data::data_expression& d) const
  {
    (void)d; // Suppress an unused variable warning.
    return true;
  }

  /// \brief Traverses an assignment
  /// \param a An assignment
  bool is_well_typed(const data::assignment& a) const
  {
    if (a.lhs().sort() != a.rhs().sort())
    {
      std::clog << "is_well_typed(data_assignment) failed: the left and right hand sides "
                << a.lhs() << " and " << a.rhs() << " have different sorts." << std::endl;
      return false;
    }
    return true;
  }

  /// \brief Traverses an action label.
  /// \param d An action label.
  bool is_well_typed(const process::action_label& d) const
  {
    (void)d; // Suppress an unused variable warning.
    return true;
  }

  /// \brief Traverses an action.
  /// \param a An action.
  bool is_well_typed(const process::action& a) const
  {
    (void)a; // Suppress an unused variable warning.
    return true;
  }

  /// \brief Checks well typedness of a deadlock
  /// \param d A deadlock
  /// \return Returns true if
  /// <ul>
  /// <li>the (optional) time has sort Real</li>
  /// </ul>
  bool is_well_typed(const deadlock& d) const
  {
    if (d.has_time())
    {
      check_time(d.time(), "deadlock");
    }
    return true;
  }

  /// \brief Checks well typedness of a multi-action
  /// \param a A multi-action
  /// \return Returns true if
  /// <ul>
  /// <li>the (optional) time has sort Real</li>
  /// </ul>
  bool is_well_typed(const multi_action& a) const
  {
    if (a.has_time())
    {
      check_time(a.time(), "multi_action");
    }
    return true;
  }

  /// \brief Checks well typedness of a summand
  /// \param s An action summand
  bool is_well_typed(const action_summand& s) const
  {
    if (!data::detail::unique_names(s.summation_variables()))
    {
      error << "is_well_typed(action_summand) failed: summation variables " << core::detail::print_list(s.summation_variables()) << " don't have unique names." << std::endl;
      return false;
    }
    if (!check_condition(s.condition(), "action_summand"))
    {
      return false;
    }
    if (!is_well_typed(s.multi_action()))
    {
      return false;
    }
    if (!check_assignments(s.assignments(), "action_summand"))
    {
      return false;
    }
    return true;
  }

  /// \brief Checks well typedness of a summand
  /// \param s A summand
  bool is_well_typed(const deadlock_summand& s) const
  {
    if (!check_condition(s.condition(), "deadlock_summand"))
    {
      return false;
    }
    if (!is_well_typed(s.deadlock()))
    {
      return false;
    }
    return true;
  }

  /// \brief Checks well typedness of a linear process
  /// \param p A linear_process
  /// \return True if
  /// <ul>
  /// <li>the process parameters have unique names</li>
  /// <li>process parameters and summation variables have different names</li>
  /// <li>the left hand sides of the assignments of summands are contained in the process parameters</li>
  /// <li>the summands are well typed</li>
  /// </ul>
  template <typename ActionSummand>
  bool is_well_typed(const linear_process_base<ActionSummand>& p) const
  {
    // check 2)
    if (!data::detail::unique_names(p.process_parameters()))
    {
      error << "is_well_typed(linear_process) failed: process parameters " << core::detail::print_list(p.process_parameters()) << " don't have unique names." << std::endl;
      return false;
    }

    // check 4)
    std::set<core::identifier_string> names;
    for (auto i = p.process_parameters().begin(); i != p.process_parameters().end(); ++i)
    {
      names.insert(i->name());
    }
    for (auto i = p.action_summands().begin(); i != p.action_summands().end(); ++i)
    {
      if (!data::detail::check_variable_names(i->summation_variables(), names))
      {
        error << "is_well_typed(linear_process) failed: some of the names of the summation variables " << core::detail::print_list(i->summation_variables()) << " also appear as process parameters." << std::endl;
        return false;
      }
    }

    // check 5)
    for (auto i = p.action_summands().begin(); i != p.action_summands().end(); ++i)
    {
      if (!data::detail::check_assignment_variables(i->assignments(), p.process_parameters()))
      {
        error << "is_well_typed(linear_process) failed: some left hand sides of the assignments " << core::detail::print_list(i->assignments()) << " do not appear as process parameters." << std::endl;
        return false;
      }
    }

    // check 6)
    if (!is_well_typed_container(p.action_summands()))
    {
      return false;
    }
    if (!is_well_typed_container(p.deadlock_summands()))
    {
      return false;
    }
    return true;
  }

  /// \brief Checks well typedness of a linear process specification.
  /// \param spec A linear process specification.
  /// \param free_variables Free variables that can be used.
  /// \return True if
  /// <ul>
  /// <li>the sorts occurring in the summation variables are declared in the data specification</li>
  /// <li>the sorts occurring in the process parameters are declared in the data specification </li>
  /// <li>the sorts occurring in the free variables are declared in the data specification     </li>
  /// <li>the sorts occurring in the action labels are declared in the data specification      </li>
  /// <li>the action labels occurring in the process are contained in action_labels()          </li>
  /// <li>the process is well typed                                                            </li>
  /// <li>the data specification is well typed                                                 </li>
  /// <li>the initial process is well typed                                                    </li>
  /// <li>the free variables occurring in the linear process are declared in the global variable specification</li>
  /// <li>the free variables occurring in the initial process are declared in the global variable specification</li>
  /// <li>the global variables have unique names</li>
  /// </ul>
  template <typename LinearProcess, typename InitialProcessExpression>
  bool is_well_typed(const specification_base<LinearProcess, InitialProcessExpression>& spec,
                     const std::set<data::variable>& free_variables) const
  {
    std::set<data::sort_expression> declared_sorts = data::detail::make_set(spec.data().sorts());
    std::set<process::action_label> declared_labels = data::detail::make_set(spec.action_labels());
    auto const& action_summands = spec.process().action_summands();

    // check 1)
    for (auto i = action_summands.begin(); i != action_summands.end(); ++i)
    {
      if (!(data::detail::check_variable_sorts(i->summation_variables(), declared_sorts)))
      {
        error << "is_well_typed(specification) failed: some of the sorts of the summation variables " << core::detail::print_list(i->summation_variables()) << " are not declared in the data specification " << core::detail::print_list(spec.data().sorts()) << std::endl;
        return false;
      }
    }

    // check 2)
    if (!(data::detail::check_variable_sorts(spec.process().process_parameters(), declared_sorts)))
    {
      error << "is_well_typed(specification) failed: some of the sorts of the process parameters " << core::detail::print_list(spec.process().process_parameters()) << " are not declared in the data specification " << core::detail::print_list(spec.data().sorts()) << std::endl;
      return false;
    }

    // check 3)
    if (!(data::detail::check_variable_sorts(spec.global_variables(), declared_sorts)))
    {
      error << "is_well_typed(specification) failed: some of the sorts of the free variables " << core::detail::print_list(spec.global_variables()) << " are not declared in the data specification " << core::detail::print_list(spec.data().sorts()) << std::endl;
      return false;
    }

    // check 4)
    if (!(detail::check_action_label_sorts(spec.action_labels(), declared_sorts)))
    {
      error << "is_well_typed(specification) failed: some of the sorts occurring in the action labels " << core::detail::print_list(spec.action_labels()) << " are not declared in the data specification " << core::detail::print_list(spec.data().sorts()) << std::endl;
      return false;
    }

    // check 5)
    for (const action_summand& s: action_summands)
    {
      if (!(detail::check_action_labels(s.multi_action().actions(), declared_labels)))
      {
        error << "is_well_typed(specification) failed: some of the labels occurring in the actions " << core::detail::print_list(s.multi_action().actions()) << " are not declared in the action specification " << core::detail::print_list(spec.action_labels()) << std::endl;
        return false;
      }
    }
    if (!is_well_typed(spec.process()))
    {
      return false;
    }
    if (!spec.data().is_well_typed())
    {
      return false;
    }
    if (!free_variables.empty())
    {
      error << "is_well_typed(specification) failed: some of the free variables were not declared\n";
      error << "declared global variables: " << core::detail::print_list(spec.global_variables()) << std::endl;
      error << "occurring free variables: " << core::detail::print_list(free_variables) << std::endl;
      return false;
    }

    // check 3)
    if (!data::detail::unique_names(spec.global_variables()))
    {
      error << "is_well_typed(specification) failed: global variables " << core::detail::print_list(spec.global_variables()) << " don't have unique names." << std::endl;
      return false;
    }

    return true;
  }

  bool is_well_typed(const specification& spec) const
  {
    std::set<data::variable> free_variables = lps::find_free_variables(spec);
    return is_well_typed(spec, free_variables);
  }

  bool is_well_typed(const stochastic_specification& spec) const
  {
    std::set<data::variable> free_variables = lps::find_free_variables(spec);
    return is_well_typed(spec, free_variables);
  }

  template <typename Term>
  bool operator()(const Term& t) const
  {
    return is_well_typed(t);
  }
};

/// \brief Checks well typedness of an LPS object.
template <typename T>
bool is_well_typed(const T& x)
{
  lps::detail::lps_well_typed_checker checker;
  return checker(x);
}

/// \brief Checks well typedness of an LPS object, and will print error messages to stderr.
template <typename T>
bool check_well_typedness(const T& x)
{
  lps::detail::lps_well_typed_checker checker;
  bool result = checker(x);
  if (!result)
  {
    mCRL2log(log::error) << checker.error.str();
  }
  return result;
}

} // namespace mcrl2::lps::detail

#endif // MCRL2_LPS_IS_WELL_TYPED_H
