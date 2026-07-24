// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/builder.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_BUILDER_H
#define MCRL2_DATA_BUILDER_H

#include "mcrl2/core/builder.h"
#include "mcrl2/data/alias.h"

#ifdef MCRL2_ENABLE_MACHINENUMBERS
#include "mcrl2/data/bag64.h"
#else
#include "mcrl2/data/bag.h"
#endif

#include "mcrl2/data/bag_comprehension.h"
#include "mcrl2/data/exists.h"
#include "mcrl2/data/lambda.h"
#include "mcrl2/data/set_comprehension.h"
#include "mcrl2/data/structured_sort.h"
#include "mcrl2/data/machine_number.h"
#include "mcrl2/data/untyped_data_parameter.h"
#include "mcrl2/data/untyped_possible_sorts.h"
#include "mcrl2/data/untyped_set_or_bag_comprehension.h"
#include "mcrl2/data/untyped_sort_variable.h"
#include "mcrl2/data/where_clause.h"




namespace mcrl2::data
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

  template <class T>
  void apply(T& result, const data::variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_variable(result, x.name(), [&](sort_expression& result){ static_cast<Derived&>(*this).apply(result, x.sort()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::function_symbol& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_function_symbol(result, x.name(), [&](sort_expression& result){ static_cast<Derived&>(*this).apply(result, x.sort()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::application& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_application(result,
       x.head(),
       x.begin(),
       x.end(),
       [&](data_expression& result, const data::data_expression& t){ static_cast<Derived&>(*this).apply(result,t);} );
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::where_clause& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_where_clause(result, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); }, [&](assignment_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.declarations()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::machine_number& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_identifier& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_assignment(result, [&](variable& result){ static_cast<Derived&>(*this).apply(result, x.lhs()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_identifier_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_identifier_assignment(result, x.lhs(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::basic_sort& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::container_sort& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_container_sort(result, x.container_name(), [&](sort_expression& result){ static_cast<Derived&>(*this).apply(result, x.element_sort()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::structured_sort& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_structured_sort(result, [&](structured_sort_constructor_list& result){ static_cast<Derived&>(*this).apply(result, x.constructors()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::function_sort& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_function_sort(result, [&](sort_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.domain()); }, [&](sort_expression& result){ static_cast<Derived&>(*this).apply(result, x.codomain()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_sort& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_possible_sorts& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_possible_sorts(result, [&](sort_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.sorts()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_sort_variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_forall(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_exists(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::lambda& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_lambda(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::set_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_set_comprehension(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::bag_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_bag_comprehension(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_set_or_bag_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_set_or_bag_comprehension(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::structured_sort_constructor_argument& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_structured_sort_constructor_argument(result, x.name(), [&](sort_expression& result){ static_cast<Derived&>(*this).apply(result, x.sort()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::structured_sort_constructor& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_structured_sort_constructor(result, x.name(), [&](structured_sort_constructor_argument_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); }, x.recogniser());
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::alias& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_alias(result, x.name(), [&](sort_expression& result){ static_cast<Derived&>(*this).apply(result, x.reference()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::data_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_data_equation(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.lhs()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_data_parameter& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_data_parameter(result, x.name(), [&](data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_abstraction(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::abstraction>(x));
    }
    else if (data::is_variable(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::variable>(x));
    }
    else if (data::is_function_symbol(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::function_symbol>(x));
    }
    else if (data::is_application(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::application>(x));
    }
    else if (data::is_where_clause(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::where_clause>(x));
    }
    else if (data::is_machine_number(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::machine_number>(x));
    }
    else if (data::is_untyped_identifier(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_identifier>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::assignment_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::assignment>(x));
    }
    else if (data::is_untyped_identifier_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_identifier_assignment>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::sort_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_basic_sort(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::basic_sort>(x));
    }
    else if (data::is_container_sort(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::container_sort>(x));
    }
    else if (data::is_structured_sort(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::structured_sort>(x));
    }
    else if (data::is_function_sort(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::function_sort>(x));
    }
    else if (data::is_untyped_sort(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_sort>(x));
    }
    else if (data::is_untyped_possible_sorts(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_possible_sorts>(x));
    }
    else if (data::is_untyped_sort_variable(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_sort_variable>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::abstraction& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::forall>(x));
    }
    else if (data::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::exists>(x));
    }
    else if (data::is_lambda(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::lambda>(x));
    }
    else if (data::is_set_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::set_comprehension>(x));
    }
    else if (data::is_bag_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::bag_comprehension>(x));
    }
    else if (data::is_untyped_set_or_bag_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_set_or_bag_comprehension>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct sort_expression_builder: public add_sort_expressions<core::builder, Derived>
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

  template <class T>
  void apply(T& result, const data::variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::function_symbol& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::application& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_application(result,
       x.head(),
       x.begin(),
       x.end(),
       [&](data_expression& result, const data::data_expression& t){ static_cast<Derived&>(*this).apply(result,t);} );
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::where_clause& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_where_clause(result, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); }, [&](assignment_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.declarations()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::machine_number& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_identifier& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_assignment(result, x.lhs(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_identifier_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_identifier_assignment(result, x.lhs(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_forall(result, x.variables(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_exists(result, x.variables(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::lambda& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_lambda(result, x.variables(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::set_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_set_comprehension(result, x.variables(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::bag_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_bag_comprehension(result, x.variables(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_set_or_bag_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_set_or_bag_comprehension(result, x.variables(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::data_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_data_equation(result, x.variables(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.lhs()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_data_parameter& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_data_parameter(result, x.name(), [&](data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_abstraction(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::abstraction>(x));
    }
    else if (data::is_variable(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::variable>(x));
    }
    else if (data::is_function_symbol(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::function_symbol>(x));
    }
    else if (data::is_application(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::application>(x));
    }
    else if (data::is_where_clause(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::where_clause>(x));
    }
    else if (data::is_machine_number(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::machine_number>(x));
    }
    else if (data::is_untyped_identifier(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_identifier>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::assignment_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::assignment>(x));
    }
    else if (data::is_untyped_identifier_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_identifier_assignment>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::abstraction& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::forall>(x));
    }
    else if (data::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::exists>(x));
    }
    else if (data::is_lambda(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::lambda>(x));
    }
    else if (data::is_set_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::set_comprehension>(x));
    }
    else if (data::is_bag_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::bag_comprehension>(x));
    }
    else if (data::is_untyped_set_or_bag_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_set_or_bag_comprehension>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct data_expression_builder: public add_data_expressions<core::builder, Derived>
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

  template <class T>
  void apply(T& result, const data::variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::function_symbol& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::application& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_application(result,
       x.head(),
       x.begin(),
       x.end(),
       [&](data_expression& result, const data::data_expression& t){ static_cast<Derived&>(*this).apply(result,t);} );
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::where_clause& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_where_clause(result, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); }, [&](assignment_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.declarations()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::machine_number& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_identifier& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_assignment(result, [&](variable& result){ static_cast<Derived&>(*this).apply(result, x.lhs()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_identifier_assignment& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_identifier_assignment(result, x.lhs(), [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_forall(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_exists(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::lambda& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_lambda(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::set_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_set_comprehension(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::bag_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_bag_comprehension(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_set_or_bag_comprehension& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_set_or_bag_comprehension(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::data_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_data_equation(result, [&](variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.condition()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.lhs()); }, [&](data_expression& result){ static_cast<Derived&>(*this).apply(result, x.rhs()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_data_parameter& x)
  {
    static_cast<Derived&>(*this).enter(x);
    data::make_untyped_data_parameter(result, x.name(), [&](data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_abstraction(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::abstraction>(x));
    }
    else if (data::is_variable(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::variable>(x));
    }
    else if (data::is_function_symbol(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::function_symbol>(x));
    }
    else if (data::is_application(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::application>(x));
    }
    else if (data::is_where_clause(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::where_clause>(x));
    }
    else if (data::is_machine_number(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::machine_number>(x));
    }
    else if (data::is_untyped_identifier(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_identifier>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::assignment_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::assignment>(x));
    }
    else if (data::is_untyped_identifier_assignment(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_identifier_assignment>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const data::abstraction& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::forall>(x));
    }
    else if (data::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::exists>(x));
    }
    else if (data::is_lambda(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::lambda>(x));
    }
    else if (data::is_set_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::set_comprehension>(x));
    }
    else if (data::is_bag_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::bag_comprehension>(x));
    }
    else if (data::is_untyped_set_or_bag_comprehension(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_set_or_bag_comprehension>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct variable_builder: public add_variables<core::builder, Derived>
{
};
//--- end generated add_variables code ---//

} // namespace mcrl2::data



#endif // MCRL2_DATA_BUILDER_H
