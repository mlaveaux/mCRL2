// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/builder.h
/// \brief Builder class for linear process specifications.

#ifndef MCRL2_LPS_BUILDER_H
#define MCRL2_LPS_BUILDER_H

#include "mcrl2/lps/stochastic_specification.h"
#include "mcrl2/process/builder.h"

namespace mcrl2::lps
{

// Adds sort expression traversal to a builder
//--- start generated add_sort_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_sort_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  void update(lps::deadlock& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::data_expression result_time;
    static_cast<Derived&>(*this).apply(result_time, x.time());
    x.time() = result_time;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_multi_action(result, [&](process::action_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::deadlock_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_summation_variables;
    static_cast<Derived&>(*this).apply(result_summation_variables, x.summation_variables());
    x.summation_variables() = result_summation_variables;
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    static_cast<Derived&>(*this).update(x.deadlock());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::action_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_summation_variables;
    static_cast<Derived&>(*this).apply(result_summation_variables, x.summation_variables());
    x.summation_variables() = result_summation_variables;
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    lps::multi_action result_multi_action;
    static_cast<Derived&>(*this).apply(result_multi_action, x.multi_action());
    x.multi_action() = result_multi_action;
    data::assignment_list result_assignments;
    static_cast<Derived&>(*this).apply(result_assignments, x.assignments());
    x.assignments() = result_assignments;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::process_initializer& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_process_initializer(result, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.expressions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::linear_process& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_process_parameters;
    static_cast<Derived&>(*this).apply(result_process_parameters, x.process_parameters());
    x.process_parameters() = result_process_parameters;
    static_cast<Derived&>(*this).update(x.deadlock_summands());
    static_cast<Derived&>(*this).update(x.action_summands());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::action_label_list result_action_labels;
    static_cast<Derived&>(*this).apply(result_action_labels, x.action_labels());
    x.action_labels() = result_action_labels;
    static_cast<Derived&>(*this).update(x.global_variables());
    static_cast<Derived&>(*this).update(x.process());
    process_initializer result_initial_process;
    static_cast<Derived&>(*this).apply(result_initial_process, x.initial_process());
    x.initial_process() = result_initial_process;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::stochastic_distribution& x)
  {
    static_cast<Derived&>(*this).enter(x);
    result = x; if (x.is_defined()) { lps::make_stochastic_distribution(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); }); }
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_action_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_summation_variables;
    static_cast<Derived&>(*this).apply(result_summation_variables, x.summation_variables());
    x.summation_variables() = result_summation_variables;
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    lps::multi_action result_multi_action;
    static_cast<Derived&>(*this).apply(result_multi_action, x.multi_action());
    x.multi_action() = result_multi_action;
    data::assignment_list result_assignments;
    static_cast<Derived&>(*this).apply(result_assignments, x.assignments());
    x.assignments() = result_assignments;
    stochastic_distribution result_distribution;
    static_cast<Derived&>(*this).apply(result_distribution, x.distribution());
    x.distribution() = result_distribution;
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_linear_process& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_process_parameters;
    static_cast<Derived&>(*this).apply(result_process_parameters, x.process_parameters());
    x.process_parameters() = result_process_parameters;
    static_cast<Derived&>(*this).update(x.deadlock_summands());
    static_cast<Derived&>(*this).update(x.action_summands());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::action_label_list result_action_labels;
    static_cast<Derived&>(*this).apply(result_action_labels, x.action_labels());
    x.action_labels() = result_action_labels;
    static_cast<Derived&>(*this).update(x.global_variables());
    static_cast<Derived&>(*this).update(x.process());
    stochastic_process_initializer result_initial_process;
    static_cast<Derived&>(*this).apply(result_initial_process, x.initial_process());
    x.initial_process() = result_initial_process;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::stochastic_process_initializer& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_stochastic_process_initializer(result, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.expressions()); }, [&](stochastic_distribution& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); });
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct sort_expression_builder: public add_sort_expressions<process::sort_expression_builder, Derived>
{
};
//--- end generated add_sort_expressions code ---//

// Adds data expression traversal to a builder
//--- start generated add_data_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_data_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  void update(lps::deadlock& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::data_expression result_time;
    static_cast<Derived&>(*this).apply(result_time, x.time());
    x.time() = result_time;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_multi_action(result, [&](process::action_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::deadlock_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    static_cast<Derived&>(*this).update(x.deadlock());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::action_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    lps::multi_action result_multi_action;
    static_cast<Derived&>(*this).apply(result_multi_action, x.multi_action());
    x.multi_action() = result_multi_action;
    data::assignment_list result_assignments;
    static_cast<Derived&>(*this).apply(result_assignments, x.assignments());
    x.assignments() = result_assignments;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::process_initializer& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_process_initializer(result, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.expressions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::linear_process& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.deadlock_summands());
    static_cast<Derived&>(*this).update(x.action_summands());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.process());
    process_initializer result_initial_process;
    static_cast<Derived&>(*this).apply(result_initial_process, x.initial_process());
    x.initial_process() = result_initial_process;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::stochastic_distribution& x)
  {
    static_cast<Derived&>(*this).enter(x);
    result = x; if (x.is_defined()) { lps::make_stochastic_distribution(result, x.variables(), [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); }); }
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_action_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    lps::multi_action result_multi_action;
    static_cast<Derived&>(*this).apply(result_multi_action, x.multi_action());
    x.multi_action() = result_multi_action;
    data::assignment_list result_assignments;
    static_cast<Derived&>(*this).apply(result_assignments, x.assignments());
    x.assignments() = result_assignments;
    stochastic_distribution result_distribution;
    static_cast<Derived&>(*this).apply(result_distribution, x.distribution());
    x.distribution() = result_distribution;
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_linear_process& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.deadlock_summands());
    static_cast<Derived&>(*this).update(x.action_summands());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.process());
    stochastic_process_initializer result_initial_process;
    static_cast<Derived&>(*this).apply(result_initial_process, x.initial_process());
    x.initial_process() = result_initial_process;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::stochastic_process_initializer& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_stochastic_process_initializer(result, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.expressions()); }, [&](stochastic_distribution& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); });
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct data_expression_builder: public add_data_expressions<process::data_expression_builder, Derived>
{
};
//--- end generated add_data_expressions code ---//

//--- start generated add_variables code ---//
template <template <class> class Builder, class Derived>
struct add_variables: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  void update(lps::deadlock& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::data_expression result_time;
    static_cast<Derived&>(*this).apply(result_time, x.time());
    x.time() = result_time;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_multi_action(result, [&](process::action_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::deadlock_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_summation_variables;
    static_cast<Derived&>(*this).apply(result_summation_variables, x.summation_variables());
    x.summation_variables() = result_summation_variables;
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    static_cast<Derived&>(*this).update(x.deadlock());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::action_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_summation_variables;
    static_cast<Derived&>(*this).apply(result_summation_variables, x.summation_variables());
    x.summation_variables() = result_summation_variables;
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    lps::multi_action result_multi_action;
    static_cast<Derived&>(*this).apply(result_multi_action, x.multi_action());
    x.multi_action() = result_multi_action;
    data::assignment_list result_assignments;
    static_cast<Derived&>(*this).apply(result_assignments, x.assignments());
    x.assignments() = result_assignments;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::process_initializer& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_process_initializer(result, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.expressions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::linear_process& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_process_parameters;
    static_cast<Derived&>(*this).apply(result_process_parameters, x.process_parameters());
    x.process_parameters() = result_process_parameters;
    static_cast<Derived&>(*this).update(x.deadlock_summands());
    static_cast<Derived&>(*this).update(x.action_summands());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.global_variables());
    static_cast<Derived&>(*this).update(x.process());
    process_initializer result_initial_process;
    static_cast<Derived&>(*this).apply(result_initial_process, x.initial_process());
    x.initial_process() = result_initial_process;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::stochastic_distribution& x)
  {
    static_cast<Derived&>(*this).enter(x);
    result = x; if (x.is_defined()) { lps::make_stochastic_distribution(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); }); }
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_action_summand& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_summation_variables;
    static_cast<Derived&>(*this).apply(result_summation_variables, x.summation_variables());
    x.summation_variables() = result_summation_variables;
    data::data_expression result_condition;
    static_cast<Derived&>(*this).apply(result_condition, x.condition());
    x.condition() = result_condition;
    lps::multi_action result_multi_action;
    static_cast<Derived&>(*this).apply(result_multi_action, x.multi_action());
    x.multi_action() = result_multi_action;
    data::assignment_list result_assignments;
    static_cast<Derived&>(*this).apply(result_assignments, x.assignments());
    x.assignments() = result_assignments;
    stochastic_distribution result_distribution;
    static_cast<Derived&>(*this).apply(result_distribution, x.distribution());
    x.distribution() = result_distribution;
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_linear_process& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::variable_list result_process_parameters;
    static_cast<Derived&>(*this).apply(result_process_parameters, x.process_parameters());
    x.process_parameters() = result_process_parameters;
    static_cast<Derived&>(*this).update(x.deadlock_summands());
    static_cast<Derived&>(*this).update(x.action_summands());
    static_cast<Derived&>(*this).leave(x);
  }

  void update(lps::stochastic_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.global_variables());
    static_cast<Derived&>(*this).update(x.process());
    stochastic_process_initializer result_initial_process;
    static_cast<Derived&>(*this).apply(result_initial_process, x.initial_process());
    x.initial_process() = result_initial_process;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const lps::stochastic_process_initializer& x)
  {
    static_cast<Derived&>(*this).enter(x);
    lps::make_stochastic_process_initializer(result, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.expressions()); }, [&](stochastic_distribution& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); });
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct variable_builder: public add_variables<process::data_expression_builder, Derived>
{
};
//--- end generated add_variables code ---//

} // namespace mcrl2::lps

#endif // MCRL2_LPS_BUILDER_H
