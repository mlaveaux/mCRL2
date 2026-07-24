// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/typecheck.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_TYPECHECK_H
#define MCRL2_MODAL_FORMULA_TYPECHECK_H

#include "mcrl2/lps/typecheck.h"
#include "mcrl2/modal_formula/detail/state_variable_context.h"
#include "mcrl2/modal_formula/is_monotonous.h"
#include "mcrl2/modal_formula/normalize_sorts.h"

namespace mcrl2
{

namespace action_formulas
{

namespace detail
{

struct typecheck_builder: public action_formula_builder<typecheck_builder>
{
  using super = action_formula_builder<typecheck_builder>;
  using super::apply;

  data::data_type_checker& m_data_type_checker;
  data::detail::variable_context m_variable_context;
  const process::detail::action_context& m_action_context;

  typecheck_builder(data::data_type_checker& data_typechecker,
                    const data::detail::variable_context& variable_context,
                    const process::detail::action_context& action_context
                   )
    : m_data_type_checker(data_typechecker),
      m_variable_context(variable_context),
      m_action_context(action_context)
  {}

  process::action typecheck_action(const core::identifier_string& name, const data::data_expression_list& parameters)
  {
    return process::typecheck_action(name, parameters, m_data_type_checker, m_variable_context, m_action_context);
  }

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    result = atermpp::down_cast<T>(m_data_type_checker.typecheck_data_expression(x, data::sort_bool::bool_(), m_variable_context));
  }

  template <class T>
  void apply(T& result, const action_formulas::at& x)
  {
    data::data_expression new_time = m_data_type_checker.typecheck_data_expression(x.time_stamp(), data::sort_real::real_(), m_variable_context);
    make_at(result, [&](action_formula& r){ (*this).apply(r, x.operand()); }, new_time);
  }

  template <class T>
  void apply(T& result, const process::untyped_multi_action& x)
  {
    // If x has size 1, first try to type check it as a data expression.
    if (x.actions().size() == 1)
    {
      const data::untyped_data_parameter& y = x.actions().front();
      try
      {
        result = atermpp::down_cast<T>(data::typecheck_untyped_data_parameter(m_data_type_checker, y.name(), y.arguments(), data::sort_bool::bool_(), m_variable_context));
        return;
      }
      // NOLINTNEXTLINE(bugprone-empty-catch) intentional: fall through to multi-action typecheck below
      catch (mcrl2::runtime_error& )
      {
        // skip
      }
    }
    // Type check it as a multi action
    process::action_list new_arguments;
    for (const data::untyped_data_parameter& a: x.actions())
    {
      new_arguments.push_front(typecheck_action(a.name(), a.arguments()));
    }
    new_arguments = atermpp::reverse(new_arguments);
    result = atermpp::down_cast<T>(action_formulas::multi_action(new_arguments));
  }

  template <class T>
  void apply(T& result, const action_formulas::forall& x)
  {
    try
    {
      auto m_variable_context_copy = m_variable_context;
      m_variable_context.add_context_variables(x.variables(), m_data_type_checker);
      action_formula body;
      (*this).apply(body, x.body());
      m_variable_context = m_variable_context_copy;
      result = atermpp::down_cast<T>(forall(x.variables(), body));
    }
    catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
    {
      throw mcrl2::runtime_error(std::string(e.what()) + "\nwhile typechecking " + action_formulas::pp(x));
    }
  }

  template <class T>
  void apply(T& result,  const action_formulas::exists& x)
  {
    try
    {
      auto m_variable_context_copy = m_variable_context;
      m_variable_context.add_context_variables(x.variables(), m_data_type_checker);
      action_formula body;
      (*this).apply(body, x.body());
      m_variable_context = m_variable_context_copy;
      result = atermpp::down_cast<T>(exists(x.variables(), body));
    }
    catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
    {
      throw mcrl2::runtime_error(std::string(e.what()) + "\nwhile typechecking " + action_formulas::pp(x));
    }
  }
};

inline
typecheck_builder make_typecheck_builder(
                    data::data_type_checker& data_typechecker,
                    const data::detail::variable_context& variables,
                    const process::detail::action_context& actions
                   )
{
  return typecheck_builder(data_typechecker, variables, actions);
}

} // namespace detail

template <typename ActionLabelContainer = std::vector<state_formulas::variable>, typename VariableContainer = std::vector<data::variable> >
action_formula typecheck_action_formula(const action_formula& x,
                                        const data::data_specification& dataspec,
                                        const VariableContainer& variables,
                                        const ActionLabelContainer& actions
                                       )
{
  data::data_type_checker data_typechecker(dataspec);
  data::detail::variable_context variable_context;
  variable_context.add_context_variables(variables, data_typechecker);
  process::detail::action_context action_context;
  action_context.add_context_action_labels(actions, data_typechecker);
  action_formula result;
  detail::make_typecheck_builder(data_typechecker, variable_context, action_context).
             apply(result, action_formulas::normalize_sorts(x, data_typechecker.typechecked_data_specification()));
  return result;
}

inline
action_formula typecheck_action_formula(const action_formula& x, const lps::stochastic_specification& lpsspec)
{
  return typecheck_action_formula(x, lpsspec.data(), lpsspec.global_variables(), lpsspec.action_labels());
}

} // namespace action_formulas

namespace regular_formulas
{

namespace detail
{

struct typecheck_builder: public regular_formula_builder<typecheck_builder>
{
  using super = regular_formula_builder<typecheck_builder>;
  using super::apply;

  data::data_type_checker& m_data_type_checker;
  const data::detail::variable_context& m_variable_context;
  const process::detail::action_context& m_action_context;

  typecheck_builder(data::data_type_checker& data_typechecker,
                    const data::detail::variable_context& variables,
                    const process::detail::action_context& actions
                   )
    : m_data_type_checker(data_typechecker),
      m_variable_context(variables),
      m_action_context(actions)
  {}

  data::data_expression make_fbag_union(const data::data_expression& left, const data::data_expression& right)
  {
    const data::sort_expression& s = atermpp::down_cast<data::application>(left).head().sort();
    const data::container_sort& cs = atermpp::down_cast<data::container_sort>(s);
    return data::sort_fbag::union_(cs.element_sort(), left, right);
  }

  data::data_expression make_bag_union(const data::data_expression& left, const data::data_expression& right)
  {
    const data::sort_expression& s = atermpp::down_cast<data::application>(left).head().sort();
    const data::container_sort& cs = atermpp::down_cast<data::container_sort>(s);
    return data::sort_bag::union_(cs.element_sort(), left, right);
  }

  data::data_expression make_fset_union(const data::data_expression& left, const data::data_expression& right)
  {
    const data::sort_expression& s = atermpp::down_cast<data::application>(left).head().sort();
    const data::container_sort& cs = atermpp::down_cast<data::container_sort>(s);
    return data::sort_fset::union_(cs.element_sort(), left, right);
  }

  data::data_expression make_set_union(const data::data_expression& left, const data::data_expression& right)
  {
    const data::sort_expression& s = atermpp::down_cast<data::application>(left).head().sort();
    const data::container_sort& cs = atermpp::down_cast<data::container_sort>(s);
    return data::sort_set::union_(cs.element_sort(), left, right);
  }

  data::data_expression make_plus(const data::data_expression& left, const data::data_expression& right)
  {
    if (data::sort_real::is_real(left.sort()) || data::sort_real::is_real(right.sort()))
    {
      return data::sort_real::plus(left, right);
    }
    else if (data::sort_int::is_int(left.sort()) || data::sort_int::is_int(right.sort()))
    {
      return data::sort_int::plus(left, right);
    }
    else if (data::sort_nat::is_nat(left.sort()) || data::sort_nat::is_nat(right.sort()))
    {
      return data::sort_nat::plus(left, right);
    }
    else if (data::sort_pos::is_pos(left.sort()) || data::sort_pos::is_pos(right.sort()))
    {
      return data::sort_pos::plus(left, right);
    }
    else if (data::sort_bag::is_union_application(left) || data::sort_bag::is_union_application(right))
    {
      return make_bag_union(left, right);
    }
    else if (data::sort_fbag::is_union_application(left) || data::sort_fbag::is_union_application(right))
    {
      return make_fbag_union(left, right);
    }
    else if (data::sort_set::is_union_application(left) || data::sort_set::is_union_application(right))
    {
      return make_set_union(left, right);
    }
    else if (data::sort_fset::is_union_application(left) || data::sort_fset::is_union_application(right))
    {
      return make_fset_union(left, right);
    }
    throw mcrl2::runtime_error("could not typecheck " + data::pp(left) + " + " + data::pp(right));
  }

  data::data_expression make_element_at(const data::data_expression& left, const data::data_expression& right) const
  {
    const data::sort_expression& s = atermpp::down_cast<data::application>(left).head().sort();
    const data::container_sort& cs = atermpp::down_cast<data::container_sort>(s);
    return data::sort_list::element_at(cs.element_sort(), left, right);
  }

  template <class T>
  void apply(T& result, const regular_formulas::untyped_regular_formula& x)
  {
    regular_formula left;
    (*this).apply(left, x.left());
    regular_formula right;
    (*this).apply(right, x.right());
    if (data::is_data_expression(left) && data::is_data_expression(right))
    {
      if (x.name() == core::identifier_string("."))
      {
        result = atermpp::down_cast<T>(make_element_at(atermpp::down_cast<data::data_expression>(left), atermpp::down_cast<data::data_expression>(right)));
      }
      else
      {
        result = atermpp::down_cast<T>(make_plus(atermpp::down_cast<data::data_expression>(left), atermpp::down_cast<data::data_expression>(right)));
      }
    }
    else
    {
      if (x.name() == core::identifier_string("."))
      {
        result = seq(left, right);
      }
      else
      {
        result = alt(left, right);
      }
    }
  }

  /* Has to be a regular formula as result.... */
  // template <class T>A
  void apply(regular_formula& result, const action_formulas::action_formula& x)
  {
    action_formulas::detail::make_typecheck_builder(m_data_type_checker, m_variable_context, m_action_context).apply(result, x);
  }
};

inline
typecheck_builder make_typecheck_builder(
                    data::data_type_checker& data_typechecker,
                    const data::detail::variable_context& variables,
                    const process::detail::action_context& actions
                   )
{
  return typecheck_builder(data_typechecker, variables, actions);
}

} // namespace detail

template <typename ActionLabelContainer = std::vector<state_formulas::variable>, typename VariableContainer = std::vector<data::variable> >
regular_formula typecheck_regular_formula(const regular_formula& x,
                                          const data::data_specification& dataspec,
                                          const VariableContainer& variables,
                                          const ActionLabelContainer& actions
                                         )
{
  data::data_type_checker data_typechecker(dataspec);
  data::detail::variable_context variable_context;
  variable_context.add_context_variables(variables, data_typechecker);
  process::detail::action_context action_context;
  action_context.add_context_action_labels(actions, data_typechecker);
  regular_formula result;
  detail::make_typecheck_builder(data_typechecker, variable_context, action_context).
           apply(result, regular_formulas::normalize_sorts(x, data_typechecker.typechecked_data_specification()));
  return result;
}

inline
regular_formula typecheck_regular_formula(const regular_formula& x, const lps::stochastic_specification& lpsspec)
{
  return typecheck_regular_formula(x, lpsspec.data(), lpsspec.global_variables(), lpsspec.action_labels());
}

} // namespace regular_formulas

namespace state_formulas
{

namespace detail
{

struct typecheck_builder: public state_formula_builder<typecheck_builder>
{
  using super = state_formula_builder<typecheck_builder>;
  using super::apply;

  data::data_type_checker& m_data_type_checker;
  data::detail::variable_context m_variable_context;
  const process::detail::action_context& m_action_context;
  detail::state_variable_context m_state_variable_context;
  const bool m_formula_is_quantitative;        // If true, allow real values, otherwise restrict to a boolean formula.  

  typecheck_builder(data::data_type_checker& data_typechecker,
                    const data::detail::variable_context& variable_context,
                    const process::detail::action_context& action_context,
                    const detail::state_variable_context& state_variable_context,
                    const bool formula_is_quantitative
                   )
    : m_data_type_checker(data_typechecker),
      m_variable_context(variable_context),
      m_action_context(action_context),
      m_state_variable_context(state_variable_context),
      m_formula_is_quantitative(formula_is_quantitative)
  {}

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    if (m_formula_is_quantitative)
    {
      try
      {
        result = atermpp::down_cast<T>(m_data_type_checker.typecheck_data_expression(x, data::sort_real::real_(), m_variable_context));
      }
      catch (mcrl2::runtime_error& e1)
      {
        try 
        {
          result = atermpp::down_cast<T>(m_data_type_checker.typecheck_data_expression(x, data::sort_bool::bool_(), m_variable_context));
        }
        catch (mcrl2::runtime_error& e2) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
        {
          throw mcrl2::runtime_error("Fail to type data expression as bool or real: " + pp(x) + ".\n" +
                                     e1.what() + "\n" +
                                     e2.what());
        }
      }
    }
    else
    {
      result = atermpp::down_cast<T>(m_data_type_checker.typecheck_data_expression(x, data::sort_bool::bool_(), m_variable_context));
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::forall& x)
  {
    if (m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("Forall is not allowed to capture values in a quantitative modal formula " + state_formulas::pp(x) + ". Use inf for infimum instead. ");
    }
    else 
    {
      try
      {
        auto m_variable_context_copy = m_variable_context;
        m_variable_context.add_context_variables(x.variables(), m_data_type_checker);
        state_formula body;
        (*this).apply(body, x.body());
        m_variable_context = m_variable_context_copy;
        result = forall(x.variables(), body);
      }
      catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
      {
        throw mcrl2::runtime_error(std::string(e.what()) + "\nwhile typechecking " + state_formulas::pp(x));
      }
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::exists& x)
  {
    if (m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("Exists is not allowed to capture values in a quantitative modal formula " + state_formulas::pp(x) + ". Use sup for supremum instead. ");
    }
    else 
    {
      try
      {
        auto m_variable_context_copy = m_variable_context;
        m_variable_context.add_context_variables(x.variables(), m_data_type_checker);
        state_formula body;
        (*this).apply(body, x.body());
        m_variable_context = m_variable_context_copy;
        result = exists(x.variables(), body);
      }
      catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
      {
        throw mcrl2::runtime_error(std::string(e.what()) + "\nwhile typechecking " + state_formulas::pp(x));
      }
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::infimum& x)
  {
    if (!m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("Infimum is not allowed in an ordinary, non quantitative modal formula " + state_formulas::pp(x) + ". Use forall instead. ");
    }
    else 
    {
      try
      {
        auto m_variable_context_copy = m_variable_context;
        m_variable_context.add_context_variables(x.variables(), m_data_type_checker);
        state_formula body;
        (*this).apply(body, x.body());
        m_variable_context = m_variable_context_copy;
        result = infimum(x.variables(), body);
      }
      catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
      {
        throw mcrl2::runtime_error(std::string(e.what()) + "\nwhile typechecking " + state_formulas::pp(x));
      }
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::supremum& x)
  {
    if (!m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("Supremum is not allowed in an ordinary, non quantitative modal formula " + state_formulas::pp(x) + ". Use exists instead. ");
    }
    else 
    {
      try
      {
        auto m_variable_context_copy = m_variable_context;
        m_variable_context.add_context_variables(x.variables(), m_data_type_checker);
        state_formula body;
        (*this).apply(body, x.body());
        m_variable_context = m_variable_context_copy;
        result = supremum(x.variables(), body);
      }
      catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
      {
        throw mcrl2::runtime_error(std::string(e.what()) + "\nwhile typechecking " + state_formulas::pp(x));
      }
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::sum& x)
  {
    if (!m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("Sum is not allowed in an ordinary, non quantitative modal formula " + state_formulas::pp(x) + ". ");
    }
    else
    {
      try
      {
        auto m_variable_context_copy = m_variable_context;
        m_variable_context.add_context_variables(x.variables(), m_data_type_checker);
        state_formula body;
        (*this).apply(body, x.body());
        m_variable_context = m_variable_context_copy;
        result = sum(x.variables(), body);
      }
      catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
      {
        throw mcrl2::runtime_error(std::string(e.what()) + "\nwhile typechecking " + state_formulas::pp(x));
      }
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::may& x)
  {
    regular_formulas::regular_formula formula;
    regular_formulas::detail::make_typecheck_builder(m_data_type_checker, m_variable_context, m_action_context).apply(formula, x.formula());
    state_formula operand;
    (*this).apply(operand, x.operand());
    make_may(result, formula, operand);
  }

  template <class T>
  void apply(T& result, const state_formulas::must& x)
  {
    regular_formulas::regular_formula formula;
    regular_formulas::detail::make_typecheck_builder(m_data_type_checker, m_variable_context, m_action_context).apply(formula, x.formula());
    state_formula operand;
    (*this).apply(operand, x.operand());
    make_must(result, formula, operand);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay_timed& x)
  {
    data::data_expression new_time = m_data_type_checker.typecheck_data_expression(x.time_stamp(), data::sort_real::real_(), m_variable_context);
    make_delay_timed(result, new_time);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled_timed& x)
  {
    data::data_expression new_time = m_data_type_checker.typecheck_data_expression(x.time_stamp(), data::sort_real::real_(), m_variable_context);
    make_yaled_timed(result, new_time);
  }

  state_formula apply_untyped_parameter(const core::identifier_string& name, const data::data_expression_list& arguments)
  {
    data::sort_expression_list expected_sorts = m_state_variable_context.matching_state_variable_sorts(name, arguments);
    data::data_expression_list new_arguments;
    auto q1 = expected_sorts.begin();
    auto q2 = arguments.begin();
    for (; q1 != expected_sorts.end(); ++q1, ++q2)
    {
      new_arguments.push_front(m_data_type_checker.typecheck_data_expression(*q2, *q1, m_variable_context));
    }
    new_arguments = atermpp::reverse(new_arguments);
    return state_formulas::variable(name, new_arguments);
  }

  template <class T>
  void apply(T& result, const state_formulas::variable& x)
  {
    result = apply_untyped_parameter(x.name(), x.arguments());
  }

  template <class T>
  void apply(T& result, const data::untyped_data_parameter& x)
  {
    result = apply_untyped_parameter(x.name(), x.arguments());
  }

  data::variable_list assignment_variables(const data::assignment_list& x) const
  {
    std::vector<data::variable> result;
    for (const data::assignment& a: x)
    {
      result.push_back(a.lhs());
    }
    return data::variable_list(result.begin(), result.end());
  }

  template <typename MuNuFormula>
  state_formula apply_mu_nu(const MuNuFormula& x, bool is_mu)
  {
    data::assignment_list new_assignments = m_data_type_checker.typecheck_assignment_list(x.assignments(), m_variable_context);

    // add the assignment variables to the context
    auto m_variable_context_copy = m_variable_context;
    data::variable_list x_variables = assignment_variables(x.assignments());
    m_variable_context.add_context_variables(x_variables, m_data_type_checker);

    // add x to the state variable context
    auto m_state_variable_context_copy = m_state_variable_context;
    m_state_variable_context.add_state_variable(x.name(), x_variables, m_data_type_checker);

    // typecheck the operand
    state_formula new_operand;
    (*this).apply(new_operand, x.operand());

    // restore the context
    m_variable_context = m_variable_context_copy;
    m_state_variable_context = m_state_variable_context_copy;

    if (is_mu)
    {
      return mu(x.name(), new_assignments, new_operand);
    }
    else
    {
      return nu(x.name(), new_assignments, new_operand);
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::nu& x)
  {
    result = apply_mu_nu(x, false);
  }

  template <class T>
  void apply(T& result, const state_formulas::mu& x)
  {
    result = apply_mu_nu(x, true);
  }

  template <class T>
  void apply(T& result, const state_formulas::not_& x)
  {
    if (m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("No ! allowed in a quantitative modal formula " + state_formulas::pp(x) + ".");
    }
    else
    {
      apply(result, static_cast<not_>(x).operand());
      make_not(result, result);
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::minus& x)
  {
    if (!m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("No unary minus (-) allowed in a boolean modal formula " + state_formulas::pp(x) + ".");
    }
    else
    {
      apply(result, static_cast<minus>(x).operand());
      make_minus(result, result);
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::plus& x)
  {
    if (!m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("No addition (+) allowed in a boolean modal formula " + state_formulas::pp(x) + ".");
    }
    else
    {
      apply(result, static_cast<plus>(x).left());
      state_formula rhs;
      apply(rhs, static_cast<plus>(x).right());
      make_plus(result, result, rhs);
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply& x)
  {
    if (!m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("No multiplication with a constant allowed in a boolean modal formula " + state_formulas::pp(x) + ".");
    }
    else
    {
      data::data_expression constant = m_data_type_checker.typecheck_data_expression(static_cast<const_multiply>(x).left(), data::sort_real::real_(), m_variable_context);
      state_formula rhs;
      apply(rhs, static_cast<const_multiply>(x).right());
      make_const_multiply(result, constant, rhs);
    }
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply_alt& x)
  {
    if (!m_formula_is_quantitative)
    {
      throw mcrl2::runtime_error("No multiplication with a constant allowed in a boolean modal formula " + state_formulas::pp(x) + ".");
    }
    else
    {
      state_formula lhs;
      apply(lhs, static_cast<const_multiply_alt>(x).left());
      data::data_expression constant = m_data_type_checker.typecheck_data_expression(static_cast<const_multiply_alt>(x).right(), data::sort_real::real_(), m_variable_context);
      make_const_multiply_alt(result, lhs, constant);
    }
  }
};

inline
typecheck_builder make_typecheck_builder(
                    data::data_type_checker& data_typechecker,
                    const data::detail::variable_context& variable_context,
                    const process::detail::action_context& action_context,
                    const detail::state_variable_context& state_variable_context,
                    const bool formula_is_quantitative
                   )
{
  return typecheck_builder(data_typechecker, variable_context, action_context, state_variable_context, formula_is_quantitative);
}

} // namespace detail

class state_formula_type_checker
{
  protected:
    data::data_type_checker m_data_type_checker;
    data::detail::variable_context m_variable_context;
    process::detail::action_context m_action_context;
    detail::state_variable_context m_state_variable_context;
    const bool m_formula_is_quantitative;

  public:
    /** \brief     Constructor for a state_formula type checker.
     *  \param[in] dataspec The data specification against which state formulas are checked.
     *  \param[in] action_labels The data labels that can occur in a state formula.
     *  \param[in] variables A container containing variables that can occur in state formulas.
     **/
    template <typename ActionLabelContainer = std::vector<state_formulas::variable>, typename VariableContainer = std::vector<data::variable> >
    state_formula_type_checker(const data::data_specification& dataspec,
                               const bool formula_is_quantitative,
                               const ActionLabelContainer& action_labels = ActionLabelContainer(),
                               const VariableContainer& variables = VariableContainer()
                              )
      : m_data_type_checker(dataspec),
        m_formula_is_quantitative(formula_is_quantitative)
    {
      m_variable_context.add_context_variables(variables, m_data_type_checker);
      m_action_context.add_context_action_labels(action_labels, m_data_type_checker);
    }

    //check correctness of the state formula as follows:
    //1) determine the types of actions according to the definitions
    //   in spec
    //2) determine the types of data expressions according to the
    //   definitions in spec
    //3) check for name conflicts of data variable declarations in
    //   forall, exists, mu and nu quantifiers
    //4) check for monotonicity of fixpoint variables
    state_formula typecheck_state_formula(const state_formula& x)
    {
      mCRL2log(log::verbose) << "type checking state formula..." << std::endl;
      state_formula result;
      detail::make_typecheck_builder(m_data_type_checker, m_variable_context, m_action_context, m_state_variable_context, m_formula_is_quantitative).
               apply(result, state_formulas::normalize_sorts(x, m_data_type_checker.typechecked_data_specification()));
      return result;
    }
};

/** \brief     Type check a state formula.
 *  Throws an exception if something went wrong.
 *  \param[in] x A state formula that has not been type checked.
 *  \param[in] dataspec The data specification against which the formulas is checked.
 *  \param[in] action_labels A declaration of action labels that can be used in the state formulas.
 *  \param[in] variables A container with global data variables.
 *  \post      formula is type checked.
 **/
template <typename VariableContainer, typename ActionLabelContainer>
state_formula typecheck_state_formula(const state_formula& x,
                                      const bool formula_is_quantitative,
                                      const data::data_specification& dataspec = data::data_specification(),
                                      const ActionLabelContainer& action_labels = ActionLabelContainer(),
                                      const VariableContainer& variables = VariableContainer()
                                     )
{
  try
  {
    state_formula_type_checker type_checker(dataspec, formula_is_quantitative, action_labels, variables);
    return type_checker.typecheck_state_formula(x);
  }
  catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
  {
    throw mcrl2::runtime_error(std::string(e.what()) + "\ncould not type check modal formula " + state_formulas::pp(x));
  }
}

/** \brief     Type check a state formula.
 *  Throws an exception if something went wrong.
 *  \param[in] x A state formula that has not been type checked.
 *  \param[in] lpsspec A linear process specifications containing data, action and global variable declarations
 *                     to be used when typechecking the formula.
 *  \post      formula is type checked.
 **/
inline
state_formula typecheck_state_formula(const state_formula& x,
                                      const lps::stochastic_specification& lpsspec,
                                      const bool formula_is_quantitative
                                     )
{
  return typecheck_state_formula(x, formula_is_quantitative, lpsspec.data(), lpsspec.action_labels(), lpsspec.global_variables());
}

/// \brief Typecheck the state formula specification formspec. It is assumed that the formula is self contained,
/// i.e. all actions and sorts must be declared.
inline
void typecheck_state_formula_specification(state_formula_specification& formspec, const bool formula_is_quantitative)
{
  try
  {
    data::data_type_checker checker(formspec.data());
    data::data_specification dataspec = checker.typechecked_data_specification();
    state_formulas::normalize_sorts(formspec, dataspec);
    state_formula_type_checker type_checker(dataspec, formula_is_quantitative, formspec.action_labels(), std::vector<data::variable>());
    formspec.formula() = type_checker.typecheck_state_formula(formspec.formula());
    formspec.data() = checker.typechecked_data_specification();
    formspec.data().translate_user_notation();
  }
  catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
  {
    throw mcrl2::runtime_error(std::string(e.what()) + "\ncould not type check modal formula specification " + state_formulas::pp(formspec));
  }
}

/// \brief Typecheck the state formula specification formspec. It is assumed that the formula is not self contained,
/// i.e. some of the actions and sorts may be declared in lpsspec.
inline
void typecheck_state_formula_specification(state_formula_specification& formspec, const lps::stochastic_specification& lpsspec, const bool formula_is_quantitative)
{
  try
  {
    data::data_type_checker checker(formspec.data());
    data::data_specification dataspec = checker.typechecked_data_specification();
    state_formulas::normalize_sorts(formspec, dataspec);
    state_formula_type_checker type_checker(dataspec, formula_is_quantitative, lpsspec.action_labels() + formspec.action_labels(), lpsspec.global_variables());
    formspec.formula() = type_checker.typecheck_state_formula(formspec.formula());
    formspec.data() = checker.typechecked_data_specification();
    formspec.data().translate_user_notation();
  }
  catch (mcrl2::runtime_error& e) // NOLINT(bugprone-empty-catch) rethrows a wrapped exception with context
  {
    throw mcrl2::runtime_error(std::string(e.what()) + "\ncould not type check modal formula specification " + state_formulas::pp(formspec));
  }
}

} // namespace state_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULA_TYPECHECK_H
