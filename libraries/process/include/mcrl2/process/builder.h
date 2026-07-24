// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/builder.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_BUILDER_H
#define MCRL2_PROCESS_BUILDER_H

#include "mcrl2/data/builder.h"
#include "mcrl2/process/process_specification.h"
#include "mcrl2/process/untyped_multi_action.h"

namespace mcrl2::process
{

/// \brief Base builder class for processes.
template <typename Derived>
struct process_builder_base: public core::builder<Derived>
{
  using super = core::builder<Derived>;
  using super::apply;

  template <class T>
  void apply(T& result, const data::untyped_data_parameter& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }
};


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

  template <class T>
  void apply(T& result, const process::action_label& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_action_label(result, x.name(), [&](data::sort_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.sorts()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(process::process_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::action_label_list result_action_labels;
    static_cast<Derived&>(*this).apply(result_action_labels, x.action_labels());
    x.action_labels() = result_action_labels;
    static_cast<Derived&>(*this).update(x.global_variables());
    static_cast<Derived&>(*this).update(x.equations());
    process_expression result_init;
    static_cast<Derived&>(*this).apply(result_init, x.init());
    x.init() = result_init;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_identifier& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_identifier(result, x.name(), [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_equation(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.formal_parameters()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.expression()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::untyped_multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_untyped_multi_action(result, [&](data::untyped_data_parameter_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_action(result, [&](action_label& result){ static_cast<Derived&>(*this).apply(result, x.label()); }, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_instance(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.actual_parameters()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_instance_assignment(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::delta& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::tau& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sum(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::block& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_block(result, x.block_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::hide& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_hide(result, x.hide_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::rename& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_rename(result, x.rename_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::comm& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_comm(result, x.comm_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::allow& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_allow(result, x.allow_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::sync& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sync(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_at(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_seq(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then_else& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then_else(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.else_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::bounded_init& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_bounded_init(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::left_merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_left_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::choice& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_choice(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::stochastic_operator& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_stochastic_operator(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::untyped_process_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_untyped_process_assignment(result, x.name(), [&](data::untyped_identifier_assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (process::is_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::action>(x));
    }
    else if (process::is_process_instance(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance>(x));
    }
    else if (process::is_process_instance_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance_assignment>(x));
    }
    else if (process::is_delta(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::delta>(x));
    }
    else if (process::is_tau(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::tau>(x));
    }
    else if (process::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sum>(x));
    }
    else if (process::is_block(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::block>(x));
    }
    else if (process::is_hide(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::hide>(x));
    }
    else if (process::is_rename(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::rename>(x));
    }
    else if (process::is_comm(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::comm>(x));
    }
    else if (process::is_allow(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::allow>(x));
    }
    else if (process::is_sync(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sync>(x));
    }
    else if (process::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::at>(x));
    }
    else if (process::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::seq>(x));
    }
    else if (process::is_if_then(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then>(x));
    }
    else if (process::is_if_then_else(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then_else>(x));
    }
    else if (process::is_bounded_init(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::bounded_init>(x));
    }
    else if (process::is_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::merge>(x));
    }
    else if (process::is_left_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::left_merge>(x));
    }
    else if (process::is_choice(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::choice>(x));
    }
    else if (process::is_stochastic_operator(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::stochastic_operator>(x));
    }
    else if (process::is_untyped_process_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_process_assignment>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct sort_expression_builder: public add_sort_expressions<data::sort_expression_builder, Derived>
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

  void update(process::process_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.equations());
    process_expression result_init;
    static_cast<Derived&>(*this).apply(result_init, x.init());
    x.init() = result_init;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_equation(result, x.identifier(), x.formal_parameters(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.expression()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::untyped_multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_untyped_multi_action(result, [&](data::untyped_data_parameter_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_action(result, x.label(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_instance(result, x.identifier(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.actual_parameters()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_instance_assignment(result, x.identifier(), [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::delta& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::tau& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sum(result, x.variables(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::block& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_block(result, x.block_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::hide& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_hide(result, x.hide_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::rename& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_rename(result, x.rename_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::comm& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_comm(result, x.comm_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::allow& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_allow(result, x.allow_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::sync& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sync(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_at(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_seq(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then_else& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then_else(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.else_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::bounded_init& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_bounded_init(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::left_merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_left_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::choice& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_choice(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::stochastic_operator& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_stochastic_operator(result, x.variables(), [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::untyped_process_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_untyped_process_assignment(result, x.name(), [&](data::untyped_identifier_assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (process::is_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::action>(x));
    }
    else if (process::is_process_instance(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance>(x));
    }
    else if (process::is_process_instance_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance_assignment>(x));
    }
    else if (process::is_delta(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::delta>(x));
    }
    else if (process::is_tau(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::tau>(x));
    }
    else if (process::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sum>(x));
    }
    else if (process::is_block(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::block>(x));
    }
    else if (process::is_hide(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::hide>(x));
    }
    else if (process::is_rename(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::rename>(x));
    }
    else if (process::is_comm(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::comm>(x));
    }
    else if (process::is_allow(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::allow>(x));
    }
    else if (process::is_sync(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sync>(x));
    }
    else if (process::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::at>(x));
    }
    else if (process::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::seq>(x));
    }
    else if (process::is_if_then(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then>(x));
    }
    else if (process::is_if_then_else(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then_else>(x));
    }
    else if (process::is_bounded_init(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::bounded_init>(x));
    }
    else if (process::is_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::merge>(x));
    }
    else if (process::is_left_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::left_merge>(x));
    }
    else if (process::is_choice(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::choice>(x));
    }
    else if (process::is_stochastic_operator(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::stochastic_operator>(x));
    }
    else if (process::is_untyped_process_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_process_assignment>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct data_expression_builder: public add_data_expressions<data::data_expression_builder, Derived>
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

  void update(process::process_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.global_variables());
    static_cast<Derived&>(*this).update(x.equations());
    process_expression result_init;
    static_cast<Derived&>(*this).apply(result_init, x.init());
    x.init() = result_init;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_identifier& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_identifier(result, x.name(), [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_equation(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.formal_parameters()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.expression()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::untyped_multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_untyped_multi_action(result, [&](data::untyped_data_parameter_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_action(result, x.label(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_instance(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.actual_parameters()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_instance_assignment(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::delta& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::tau& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sum(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::block& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_block(result, x.block_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::hide& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_hide(result, x.hide_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::rename& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_rename(result, x.rename_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::comm& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_comm(result, x.comm_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::allow& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_allow(result, x.allow_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::sync& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sync(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_at(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_seq(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then_else& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then_else(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.else_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::bounded_init& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_bounded_init(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::left_merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_left_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::choice& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_choice(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::stochastic_operator& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_stochastic_operator(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.distribution()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::untyped_process_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_untyped_process_assignment(result, x.name(), [&](data::untyped_identifier_assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (process::is_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::action>(x));
    }
    else if (process::is_process_instance(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance>(x));
    }
    else if (process::is_process_instance_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance_assignment>(x));
    }
    else if (process::is_delta(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::delta>(x));
    }
    else if (process::is_tau(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::tau>(x));
    }
    else if (process::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sum>(x));
    }
    else if (process::is_block(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::block>(x));
    }
    else if (process::is_hide(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::hide>(x));
    }
    else if (process::is_rename(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::rename>(x));
    }
    else if (process::is_comm(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::comm>(x));
    }
    else if (process::is_allow(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::allow>(x));
    }
    else if (process::is_sync(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sync>(x));
    }
    else if (process::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::at>(x));
    }
    else if (process::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::seq>(x));
    }
    else if (process::is_if_then(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then>(x));
    }
    else if (process::is_if_then_else(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then_else>(x));
    }
    else if (process::is_bounded_init(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::bounded_init>(x));
    }
    else if (process::is_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::merge>(x));
    }
    else if (process::is_left_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::left_merge>(x));
    }
    else if (process::is_choice(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::choice>(x));
    }
    else if (process::is_stochastic_operator(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::stochastic_operator>(x));
    }
    else if (process::is_untyped_process_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_process_assignment>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct variable_builder: public add_variables<data::data_expression_builder, Derived>
{
};
//--- end generated add_variables code ---//

// Adds process expression traversal to a builder
//--- start generated add_process_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_process_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  void update(process::process_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.equations());
    process_expression result_init;
    static_cast<Derived&>(*this).apply(result_init, x.init());
    x.init() = result_init;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_equation(result, x.identifier(), x.formal_parameters(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.expression()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::delta& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::tau& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sum(result, x.variables(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::block& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_block(result, x.block_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::hide& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_hide(result, x.hide_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::rename& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_rename(result, x.rename_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::comm& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_comm(result, x.comm_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::allow& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_allow(result, x.allow_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::sync& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sync(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_at(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, x.time_stamp());
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_seq(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then(result, x.condition(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then_else& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then_else(result, x.condition(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.else_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::bounded_init& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_bounded_init(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::left_merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_left_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::choice& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_choice(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::stochastic_operator& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_stochastic_operator(result, x.variables(), x.distribution(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::untyped_process_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::process_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (process::is_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::action>(x));
    }
    else if (process::is_process_instance(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance>(x));
    }
    else if (process::is_process_instance_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance_assignment>(x));
    }
    else if (process::is_delta(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::delta>(x));
    }
    else if (process::is_tau(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::tau>(x));
    }
    else if (process::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sum>(x));
    }
    else if (process::is_block(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::block>(x));
    }
    else if (process::is_hide(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::hide>(x));
    }
    else if (process::is_rename(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::rename>(x));
    }
    else if (process::is_comm(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::comm>(x));
    }
    else if (process::is_allow(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::allow>(x));
    }
    else if (process::is_sync(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sync>(x));
    }
    else if (process::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::at>(x));
    }
    else if (process::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::seq>(x));
    }
    else if (process::is_if_then(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then>(x));
    }
    else if (process::is_if_then_else(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then_else>(x));
    }
    else if (process::is_bounded_init(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::bounded_init>(x));
    }
    else if (process::is_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::merge>(x));
    }
    else if (process::is_left_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::left_merge>(x));
    }
    else if (process::is_choice(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::choice>(x));
    }
    else if (process::is_stochastic_operator(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::stochastic_operator>(x));
    }
    else if (process::is_untyped_process_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_process_assignment>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct process_expression_builder: public add_process_expressions<process::process_builder_base, Derived>
{
};
//--- end generated add_process_expressions code ---//

//--- start generated add_process_identifiers code ---//
template <template <class> class Builder, class Derived>
struct add_process_identifiers: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  void update(process::process_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.equations());
    process_expression result_init;
    static_cast<Derived&>(*this).apply(result_init, x.init());
    x.init() = result_init;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_identifier& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::process_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_equation(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, x.formal_parameters(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.expression()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_instance(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, x.actual_parameters());
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::process_instance_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_process_instance_assignment(result, [&](process_identifier& result){ static_cast<Derived&>(*this).apply(result, x.identifier()); }, x.assignments());
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::delta& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::tau& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sum(result, x.variables(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::block& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_block(result, x.block_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::hide& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_hide(result, x.hide_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::rename& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_rename(result, x.rename_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::comm& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_comm(result, x.comm_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::allow& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_allow(result, x.allow_set(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::sync& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_sync(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_at(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, x.time_stamp());
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_seq(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then(result, x.condition(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::if_then_else& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_if_then_else(result, x.condition(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.then_case()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.else_case()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::bounded_init& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_bounded_init(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::left_merge& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_left_merge(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::choice& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_choice(result, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::stochastic_operator& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::make_stochastic_operator(result, x.variables(), x.distribution(), [&](process_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const process::untyped_process_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const process::process_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (process::is_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::action>(x));
    }
    else if (process::is_process_instance(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance>(x));
    }
    else if (process::is_process_instance_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::process_instance_assignment>(x));
    }
    else if (process::is_delta(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::delta>(x));
    }
    else if (process::is_tau(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::tau>(x));
    }
    else if (process::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sum>(x));
    }
    else if (process::is_block(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::block>(x));
    }
    else if (process::is_hide(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::hide>(x));
    }
    else if (process::is_rename(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::rename>(x));
    }
    else if (process::is_comm(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::comm>(x));
    }
    else if (process::is_allow(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::allow>(x));
    }
    else if (process::is_sync(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::sync>(x));
    }
    else if (process::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::at>(x));
    }
    else if (process::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::seq>(x));
    }
    else if (process::is_if_then(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then>(x));
    }
    else if (process::is_if_then_else(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::if_then_else>(x));
    }
    else if (process::is_bounded_init(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::bounded_init>(x));
    }
    else if (process::is_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::merge>(x));
    }
    else if (process::is_left_merge(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::left_merge>(x));
    }
    else if (process::is_choice(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::choice>(x));
    }
    else if (process::is_stochastic_operator(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::stochastic_operator>(x));
    }
    else if (process::is_untyped_process_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_process_assignment>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct process_identifier_builder: public add_process_identifiers<process::process_builder_base, Derived>
{
};
//--- end generated add_process_identifiers code ---//

} // namespace mcrl2::process

#endif // MCRL2_PROCESS_BUILDER_H
