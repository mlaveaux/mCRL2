// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/builder.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_BUILDER_H
#define MCRL2_MODAL_FORMULA_BUILDER_H

#include "mcrl2/lps/builder.h"
#include "mcrl2/modal_formula/state_formula_specification.h"

namespace mcrl2
{

namespace action_formulas
{

/// \brief Base class for action_formula_builder.
template <typename Derived>
struct action_formula_builder_base: public core::builder<Derived>
{
  using super = core::builder<Derived>;
  using super::apply;

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const data::untyped_data_parameter& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }
};

//--- start generated action_formulas::add_sort_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_sort_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const action_formulas::true_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::false_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_not(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_and(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_or(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_imp(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_forall(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_exists(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_at(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_multi_action(result, [&](process::action_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::action_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (action_formulas::is_true(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::true_>(x));
    }
    else if (action_formulas::is_false(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::false_>(x));
    }
    else if (action_formulas::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::not_>(x));
    }
    else if (action_formulas::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::and_>(x));
    }
    else if (action_formulas::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::or_>(x));
    }
    else if (action_formulas::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::imp>(x));
    }
    else if (action_formulas::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::forall>(x));
    }
    else if (action_formulas::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::exists>(x));
    }
    else if (action_formulas::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::at>(x));
    }
    else if (action_formulas::is_multi_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::multi_action>(x));
    }
    else if (process::is_untyped_multi_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_multi_action>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct sort_expression_builder: public add_sort_expressions<lps::sort_expression_builder, Derived>
{
};
//--- end generated action_formulas::add_sort_expressions code ---//

//--- start generated action_formulas::add_data_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_data_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const action_formulas::true_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::false_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_not(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_and(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_or(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_imp(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_forall(result, x.variables(), [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_exists(result, x.variables(), [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_at(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_multi_action(result, [&](process::action_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::action_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (action_formulas::is_true(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::true_>(x));
    }
    else if (action_formulas::is_false(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::false_>(x));
    }
    else if (action_formulas::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::not_>(x));
    }
    else if (action_formulas::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::and_>(x));
    }
    else if (action_formulas::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::or_>(x));
    }
    else if (action_formulas::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::imp>(x));
    }
    else if (action_formulas::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::forall>(x));
    }
    else if (action_formulas::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::exists>(x));
    }
    else if (action_formulas::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::at>(x));
    }
    else if (action_formulas::is_multi_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::multi_action>(x));
    }
    else if (process::is_untyped_multi_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_multi_action>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct data_expression_builder: public add_data_expressions<lps::data_expression_builder, Derived>
{
};
//--- end generated action_formulas::add_data_expressions code ---//

//--- start generated action_formulas::add_variables code ---//
template <template <class> class Builder, class Derived>
struct add_variables: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const action_formulas::true_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::false_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_not(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_and(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_or(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_imp(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_forall(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_exists(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_at(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_multi_action(result, [&](process::action_list& result){ static_cast<Derived&>(*this).apply(result, x.actions()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::action_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (action_formulas::is_true(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::true_>(x));
    }
    else if (action_formulas::is_false(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::false_>(x));
    }
    else if (action_formulas::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::not_>(x));
    }
    else if (action_formulas::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::and_>(x));
    }
    else if (action_formulas::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::or_>(x));
    }
    else if (action_formulas::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::imp>(x));
    }
    else if (action_formulas::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::forall>(x));
    }
    else if (action_formulas::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::exists>(x));
    }
    else if (action_formulas::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::at>(x));
    }
    else if (action_formulas::is_multi_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::multi_action>(x));
    }
    else if (process::is_untyped_multi_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_multi_action>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct variable_builder: public add_variables<lps::data_expression_builder, Derived>
{
};
//--- end generated action_formulas::add_variables code ---//

//--- start generated action_formulas::add_action_formula_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_action_formula_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const action_formulas::true_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::false_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_not(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_and(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_or(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_imp(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_forall(result, x.variables(), [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_exists(result, x.variables(), [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::at& x)
  {
    static_cast<Derived&>(*this).enter(x);
    action_formulas::make_at(result, [&](action_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); }, x.time_stamp());
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::multi_action& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const action_formulas::action_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (action_formulas::is_true(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::true_>(x));
    }
    else if (action_formulas::is_false(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::false_>(x));
    }
    else if (action_formulas::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::not_>(x));
    }
    else if (action_formulas::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::and_>(x));
    }
    else if (action_formulas::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::or_>(x));
    }
    else if (action_formulas::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::imp>(x));
    }
    else if (action_formulas::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::forall>(x));
    }
    else if (action_formulas::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::exists>(x));
    }
    else if (action_formulas::is_at(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::at>(x));
    }
    else if (action_formulas::is_multi_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::multi_action>(x));
    }
    else if (process::is_untyped_multi_action(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<process::untyped_multi_action>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct action_formula_builder: public add_action_formula_expressions<action_formulas::action_formula_builder_base, Derived>
{
};
//--- end generated action_formulas::add_action_formula_expressions code ---//

} // namespace action_formulas

namespace regular_formulas
{

/// \brief Builder class for regular_formula_builder. Used as a base class for pbes_expression_builder.
template <typename Derived>
struct regular_formula_builder_base: public core::builder<Derived>
{
  using super = core::builder<Derived>;
  using super::apply;

  template <class T>
  void apply(T& result, const data::data_expression& x)
  // void apply(regular_formula& result, const data::data_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }
 
  template <class T>
  void apply(T& result, const action_formulas::action_formula& x)
  // void apply(regular_formula& result, const action_formulas::action_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }
};

//--- start generated regular_formulas::add_sort_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_sort_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const regular_formulas::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_seq(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::alt& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_alt(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::trans& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_trans(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::trans_or_nil& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_trans_or_nil(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::untyped_regular_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_untyped_regular_formula(result, x.name(), [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::regular_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (action_formulas::is_action_formula(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::action_formula>(x));
    }
    else if (regular_formulas::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::seq>(x));
    }
    else if (regular_formulas::is_alt(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::alt>(x));
    }
    else if (regular_formulas::is_trans(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::trans>(x));
    }
    else if (regular_formulas::is_trans_or_nil(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::trans_or_nil>(x));
    }
    else if (regular_formulas::is_untyped_regular_formula(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::untyped_regular_formula>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct sort_expression_builder: public add_sort_expressions<action_formulas::sort_expression_builder, Derived>
{
};
//--- end generated regular_formulas::add_sort_expressions code ---//

//--- start generated regular_formulas::add_data_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_data_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const regular_formulas::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_seq(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::alt& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_alt(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::trans& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_trans(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::trans_or_nil& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_trans_or_nil(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::untyped_regular_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_untyped_regular_formula(result, x.name(), [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::regular_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (action_formulas::is_action_formula(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::action_formula>(x));
    }
    else if (regular_formulas::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::seq>(x));
    }
    else if (regular_formulas::is_alt(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::alt>(x));
    }
    else if (regular_formulas::is_trans(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::trans>(x));
    }
    else if (regular_formulas::is_trans_or_nil(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::trans_or_nil>(x));
    }
    else if (regular_formulas::is_untyped_regular_formula(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::untyped_regular_formula>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct data_expression_builder: public add_data_expressions<action_formulas::data_expression_builder, Derived>
{
};
//--- end generated regular_formulas::add_data_expressions code ---//

//--- start generated regular_formulas::add_variables code ---//
template <template <class> class Builder, class Derived>
struct add_variables: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const regular_formulas::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_seq(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::alt& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_alt(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::trans& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_trans(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::trans_or_nil& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_trans_or_nil(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::untyped_regular_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_untyped_regular_formula(result, x.name(), [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::regular_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (action_formulas::is_action_formula(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::action_formula>(x));
    }
    else if (regular_formulas::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::seq>(x));
    }
    else if (regular_formulas::is_alt(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::alt>(x));
    }
    else if (regular_formulas::is_trans(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::trans>(x));
    }
    else if (regular_formulas::is_trans_or_nil(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::trans_or_nil>(x));
    }
    else if (regular_formulas::is_untyped_regular_formula(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::untyped_regular_formula>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct variable_builder: public add_variables<action_formulas::data_expression_builder, Derived>
{
};
//--- end generated regular_formulas::add_variables code ---//

//--- start generated regular_formulas::add_regular_formula_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_regular_formula_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const regular_formulas::seq& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_seq(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::alt& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_alt(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::trans& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_trans(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::trans_or_nil& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_trans_or_nil(result, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::untyped_regular_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    regular_formulas::make_untyped_regular_formula(result, x.name(), [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const regular_formulas::regular_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (action_formulas::is_action_formula(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<action_formulas::action_formula>(x));
    }
    else if (regular_formulas::is_seq(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::seq>(x));
    }
    else if (regular_formulas::is_alt(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::alt>(x));
    }
    else if (regular_formulas::is_trans(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::trans>(x));
    }
    else if (regular_formulas::is_trans_or_nil(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::trans_or_nil>(x));
    }
    else if (regular_formulas::is_untyped_regular_formula(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<regular_formulas::untyped_regular_formula>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct regular_formula_builder: public add_regular_formula_expressions<regular_formulas::regular_formula_builder_base, Derived>
{
};
//--- end generated regular_formulas::add_regular_formula_expressions code ---//

} // namespace regular_formulas

namespace state_formulas
{

/// \brief Builder class for pbes_expressions. Used as a base class for pbes_expression_builder.
template <typename Derived>
struct state_formula_builder_base: public core::builder<Derived>
{
  using super = core::builder<Derived>;
  using super::apply;

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  } 

  template <class T>
  void apply(T& result, const data::untyped_data_parameter& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  } 

};

//--- start generated state_formulas::add_sort_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_sort_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const state_formulas::true_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::false_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_not(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::minus& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_minus(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_and(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_or(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_imp(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::plus& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_plus(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_const_multiply(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply_alt& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_const_multiply_alt(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_forall(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_exists(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::infimum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_infimum(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::supremum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_supremum(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_sum(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::must& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_must(result, [&](regular_formulas::regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.formula()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::may& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_may(result, [&](regular_formulas::regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.formula()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled_timed& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_yaled_timed(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay_timed& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_delay_timed(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_variable(result, x.name(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::nu& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_nu(result, x.name(), [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::mu& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_mu(result, x.name(), [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(state_formulas::state_formula_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    process::action_label_list result_action_labels;
    static_cast<Derived&>(*this).apply(result_action_labels, x.action_labels());
    x.action_labels() = result_action_labels;
    state_formula result_formula;
    static_cast<Derived&>(*this).apply(result_formula, x.formula());
    x.formula() = result_formula;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::state_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (state_formulas::is_true(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::true_>(x));
    }
    else if (state_formulas::is_false(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::false_>(x));
    }
    else if (state_formulas::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::not_>(x));
    }
    else if (state_formulas::is_minus(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::minus>(x));
    }
    else if (state_formulas::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::and_>(x));
    }
    else if (state_formulas::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::or_>(x));
    }
    else if (state_formulas::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::imp>(x));
    }
    else if (state_formulas::is_plus(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::plus>(x));
    }
    else if (state_formulas::is_const_multiply(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::const_multiply>(x));
    }
    else if (state_formulas::is_const_multiply_alt(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::const_multiply_alt>(x));
    }
    else if (state_formulas::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::forall>(x));
    }
    else if (state_formulas::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::exists>(x));
    }
    else if (state_formulas::is_infimum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::infimum>(x));
    }
    else if (state_formulas::is_supremum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::supremum>(x));
    }
    else if (state_formulas::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::sum>(x));
    }
    else if (state_formulas::is_must(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::must>(x));
    }
    else if (state_formulas::is_may(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::may>(x));
    }
    else if (state_formulas::is_yaled(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::yaled>(x));
    }
    else if (state_formulas::is_yaled_timed(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::yaled_timed>(x));
    }
    else if (state_formulas::is_delay(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::delay>(x));
    }
    else if (state_formulas::is_delay_timed(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::delay_timed>(x));
    }
    else if (state_formulas::is_variable(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::variable>(x));
    }
    else if (state_formulas::is_nu(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::nu>(x));
    }
    else if (state_formulas::is_mu(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::mu>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct sort_expression_builder: public add_sort_expressions<regular_formulas::sort_expression_builder, Derived>
{
};
//--- end generated state_formulas::add_sort_expressions code ---//

//--- start generated state_formulas::add_data_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_data_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const state_formulas::true_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::false_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_not(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::minus& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_minus(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_and(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_or(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_imp(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::plus& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_plus(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_const_multiply(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply_alt& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_const_multiply_alt(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_forall(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_exists(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::infimum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_infimum(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::supremum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_supremum(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_sum(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::must& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_must(result, [&](regular_formulas::regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.formula()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::may& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_may(result, [&](regular_formulas::regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.formula()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled_timed& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_yaled_timed(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay_timed& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_delay_timed(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_variable(result, x.name(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::nu& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_nu(result, x.name(), [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::mu& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_mu(result, x.name(), [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(state_formulas::state_formula_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formula result_formula;
    static_cast<Derived&>(*this).apply(result_formula, x.formula());
    x.formula() = result_formula;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::state_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (state_formulas::is_true(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::true_>(x));
    }
    else if (state_formulas::is_false(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::false_>(x));
    }
    else if (state_formulas::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::not_>(x));
    }
    else if (state_formulas::is_minus(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::minus>(x));
    }
    else if (state_formulas::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::and_>(x));
    }
    else if (state_formulas::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::or_>(x));
    }
    else if (state_formulas::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::imp>(x));
    }
    else if (state_formulas::is_plus(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::plus>(x));
    }
    else if (state_formulas::is_const_multiply(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::const_multiply>(x));
    }
    else if (state_formulas::is_const_multiply_alt(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::const_multiply_alt>(x));
    }
    else if (state_formulas::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::forall>(x));
    }
    else if (state_formulas::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::exists>(x));
    }
    else if (state_formulas::is_infimum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::infimum>(x));
    }
    else if (state_formulas::is_supremum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::supremum>(x));
    }
    else if (state_formulas::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::sum>(x));
    }
    else if (state_formulas::is_must(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::must>(x));
    }
    else if (state_formulas::is_may(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::may>(x));
    }
    else if (state_formulas::is_yaled(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::yaled>(x));
    }
    else if (state_formulas::is_yaled_timed(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::yaled_timed>(x));
    }
    else if (state_formulas::is_delay(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::delay>(x));
    }
    else if (state_formulas::is_delay_timed(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::delay_timed>(x));
    }
    else if (state_formulas::is_variable(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::variable>(x));
    }
    else if (state_formulas::is_nu(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::nu>(x));
    }
    else if (state_formulas::is_mu(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::mu>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct data_expression_builder: public add_data_expressions<regular_formulas::data_expression_builder, Derived>
{
};
//--- end generated state_formulas::add_data_expressions code ---//

//--- start generated state_formulas::add_variables code ---//
template <template <class> class Builder, class Derived>
struct add_variables: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const state_formulas::true_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::false_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_not(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::minus& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_minus(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_and(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_or(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_imp(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::plus& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_plus(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_const_multiply(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply_alt& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_const_multiply_alt(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_forall(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_exists(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::infimum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_infimum(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::supremum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_supremum(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_sum(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::must& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_must(result, [&](regular_formulas::regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.formula()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::may& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_may(result, [&](regular_formulas::regular_formula& result){ static_cast<Derived&>(*this).apply(result, x.formula()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled_timed& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_yaled_timed(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay_timed& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_delay_timed(result, [&](data::data_expression& result){ static_cast<Derived&>(*this).apply(result, x.time_stamp()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_variable(result, x.name(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.arguments()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::nu& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_nu(result, x.name(), [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::mu& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_mu(result, x.name(), [&](data::assignment_list& result){ static_cast<Derived&>(*this).apply(result, x.assignments()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(state_formulas::state_formula_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formula result_formula;
    static_cast<Derived&>(*this).apply(result_formula, x.formula());
    x.formula() = result_formula;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::state_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (state_formulas::is_true(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::true_>(x));
    }
    else if (state_formulas::is_false(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::false_>(x));
    }
    else if (state_formulas::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::not_>(x));
    }
    else if (state_formulas::is_minus(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::minus>(x));
    }
    else if (state_formulas::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::and_>(x));
    }
    else if (state_formulas::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::or_>(x));
    }
    else if (state_formulas::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::imp>(x));
    }
    else if (state_formulas::is_plus(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::plus>(x));
    }
    else if (state_formulas::is_const_multiply(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::const_multiply>(x));
    }
    else if (state_formulas::is_const_multiply_alt(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::const_multiply_alt>(x));
    }
    else if (state_formulas::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::forall>(x));
    }
    else if (state_formulas::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::exists>(x));
    }
    else if (state_formulas::is_infimum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::infimum>(x));
    }
    else if (state_formulas::is_supremum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::supremum>(x));
    }
    else if (state_formulas::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::sum>(x));
    }
    else if (state_formulas::is_must(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::must>(x));
    }
    else if (state_formulas::is_may(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::may>(x));
    }
    else if (state_formulas::is_yaled(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::yaled>(x));
    }
    else if (state_formulas::is_yaled_timed(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::yaled_timed>(x));
    }
    else if (state_formulas::is_delay(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::delay>(x));
    }
    else if (state_formulas::is_delay_timed(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::delay_timed>(x));
    }
    else if (state_formulas::is_variable(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::variable>(x));
    }
    else if (state_formulas::is_nu(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::nu>(x));
    }
    else if (state_formulas::is_mu(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::mu>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct variable_builder: public add_variables<regular_formulas::data_expression_builder, Derived>
{
};
//--- end generated state_formulas::add_variables code ---//

//--- start generated state_formulas::add_state_formula_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_state_formula_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  template <class T>
  void apply(T& result, const state_formulas::true_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::false_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_not(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::minus& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_minus(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_and(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_or(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_imp(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::plus& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_plus(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_const_multiply(result, x.left(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::const_multiply_alt& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_const_multiply_alt(result, [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, x.right());
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_forall(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_exists(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::infimum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_infimum(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::supremum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_supremum(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::sum& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_sum(result, x.variables(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::must& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_must(result, x.formula(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::may& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_may(result, x.formula(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::yaled_timed& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::delay_timed& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::nu& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_nu(result, x.name(), x.assignments(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::mu& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formulas::make_mu(result, x.name(), x.assignments(), [&](state_formula& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(state_formulas::state_formula_specification& x)
  {
    static_cast<Derived&>(*this).enter(x);
    state_formula result_formula;
    static_cast<Derived&>(*this).apply(result_formula, x.formula());
    x.formula() = result_formula;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const state_formulas::state_formula& x)
  {
    static_cast<Derived&>(*this).enter(x);
    if (data::is_data_expression(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::data_expression>(x));
    }
    else if (data::is_untyped_data_parameter(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<data::untyped_data_parameter>(x));
    }
    else if (state_formulas::is_true(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::true_>(x));
    }
    else if (state_formulas::is_false(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::false_>(x));
    }
    else if (state_formulas::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::not_>(x));
    }
    else if (state_formulas::is_minus(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::minus>(x));
    }
    else if (state_formulas::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::and_>(x));
    }
    else if (state_formulas::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::or_>(x));
    }
    else if (state_formulas::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::imp>(x));
    }
    else if (state_formulas::is_plus(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::plus>(x));
    }
    else if (state_formulas::is_const_multiply(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::const_multiply>(x));
    }
    else if (state_formulas::is_const_multiply_alt(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::const_multiply_alt>(x));
    }
    else if (state_formulas::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::forall>(x));
    }
    else if (state_formulas::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::exists>(x));
    }
    else if (state_formulas::is_infimum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::infimum>(x));
    }
    else if (state_formulas::is_supremum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::supremum>(x));
    }
    else if (state_formulas::is_sum(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::sum>(x));
    }
    else if (state_formulas::is_must(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::must>(x));
    }
    else if (state_formulas::is_may(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::may>(x));
    }
    else if (state_formulas::is_yaled(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::yaled>(x));
    }
    else if (state_formulas::is_yaled_timed(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::yaled_timed>(x));
    }
    else if (state_formulas::is_delay(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::delay>(x));
    }
    else if (state_formulas::is_delay_timed(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::delay_timed>(x));
    }
    else if (state_formulas::is_variable(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::variable>(x));
    }
    else if (state_formulas::is_nu(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::nu>(x));
    }
    else if (state_formulas::is_mu(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<state_formulas::mu>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct state_formula_builder: public add_state_formula_expressions<state_formulas::state_formula_builder_base, Derived>
{
};
//--- end generated state_formulas::add_state_formula_expressions code ---//

} // namespace state_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULA_BUILDER_H
