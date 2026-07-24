// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/builder.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_BUILDER_H
#define MCRL2_PBES_BUILDER_H

#include "mcrl2/data/builder.h"
#include "mcrl2/pbes/pbes.h"

namespace mcrl2::pbes_system
{

/// \brief Builder class
template <typename Derived>
struct pbes_expression_builder_base: public core::builder<Derived>
{
  using super = core::builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;
  	
  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    result = atermpp::down_cast<T>(x);
  }
  template <class T>
  void apply(T& result, const data::untyped_data_parameter& x)
  {
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
  void apply(T& result, const pbes_system::propositional_variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_propositional_variable(result, x.name(), [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.parameters()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(pbes_system::pbes_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    propositional_variable result_variable;
    static_cast<Derived&>(*this).apply(result_variable, x.variable());
    x.variable() = result_variable;
    pbes_expression result_formula;
    static_cast<Derived&>(*this).apply(result_formula, x.formula());
    x.formula() = result_formula;
    static_cast<Derived&>(*this).leave(x);
  }

  void update(pbes_system::pbes& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.global_variables());
    static_cast<Derived&>(*this).update(x.equations());
    propositional_variable_instantiation result_initial_state;
    static_cast<Derived&>(*this).apply(result_initial_state, x.initial_state());
    x.initial_state() = result_initial_state;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::propositional_variable_instantiation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_propositional_variable_instantiation(result, x.name(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.parameters()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_not(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_and(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_or(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_imp(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_forall(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_exists(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::pbes_expression& x)
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
    else if (pbes_system::is_propositional_variable_instantiation(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::propositional_variable_instantiation>(x));
    }
    else if (pbes_system::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::not_>(x));
    }
    else if (pbes_system::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::and_>(x));
    }
    else if (pbes_system::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::or_>(x));
    }
    else if (pbes_system::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::imp>(x));
    }
    else if (pbes_system::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::forall>(x));
    }
    else if (pbes_system::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::exists>(x));
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

  void update(pbes_system::pbes_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_expression result_formula;
    static_cast<Derived&>(*this).apply(result_formula, x.formula());
    x.formula() = result_formula;
    static_cast<Derived&>(*this).leave(x);
  }

  void update(pbes_system::pbes& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.equations());
    propositional_variable_instantiation result_initial_state;
    static_cast<Derived&>(*this).apply(result_initial_state, x.initial_state());
    x.initial_state() = result_initial_state;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::propositional_variable_instantiation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_propositional_variable_instantiation(result, x.name(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.parameters()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_not(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_and(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_or(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_imp(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_forall(result, x.variables(), [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_exists(result, x.variables(), [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::pbes_expression& x)
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
    else if (pbes_system::is_propositional_variable_instantiation(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::propositional_variable_instantiation>(x));
    }
    else if (pbes_system::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::not_>(x));
    }
    else if (pbes_system::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::and_>(x));
    }
    else if (pbes_system::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::or_>(x));
    }
    else if (pbes_system::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::imp>(x));
    }
    else if (pbes_system::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::forall>(x));
    }
    else if (pbes_system::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::exists>(x));
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

  template <class T>
  void apply(T& result, const pbes_system::propositional_variable& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_propositional_variable(result, x.name(), [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.parameters()); });
    static_cast<Derived&>(*this).leave(x);
  }

  void update(pbes_system::pbes_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    propositional_variable result_variable;
    static_cast<Derived&>(*this).apply(result_variable, x.variable());
    x.variable() = result_variable;
    pbes_expression result_formula;
    static_cast<Derived&>(*this).apply(result_formula, x.formula());
    x.formula() = result_formula;
    static_cast<Derived&>(*this).leave(x);
  }

  void update(pbes_system::pbes& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.global_variables());
    static_cast<Derived&>(*this).update(x.equations());
    propositional_variable_instantiation result_initial_state;
    static_cast<Derived&>(*this).apply(result_initial_state, x.initial_state());
    x.initial_state() = result_initial_state;
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::propositional_variable_instantiation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_propositional_variable_instantiation(result, x.name(), [&](data::data_expression_list& result){ static_cast<Derived&>(*this).apply(result, x.parameters()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_not(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_and(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_or(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_imp(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_forall(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_exists(result, [&](data::variable_list& result){ static_cast<Derived&>(*this).apply(result, x.variables()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::pbes_expression& x)
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
    else if (pbes_system::is_propositional_variable_instantiation(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::propositional_variable_instantiation>(x));
    }
    else if (pbes_system::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::not_>(x));
    }
    else if (pbes_system::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::and_>(x));
    }
    else if (pbes_system::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::or_>(x));
    }
    else if (pbes_system::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::imp>(x));
    }
    else if (pbes_system::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::forall>(x));
    }
    else if (pbes_system::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::exists>(x));
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

// Adds pbes expression traversal to a builder
//--- start generated add_pbes_expressions code ---//
template <template <class> class Builder, class Derived>
struct add_pbes_expressions: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  void update(pbes_system::pbes_equation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_expression result_formula;
    static_cast<Derived&>(*this).apply(result_formula, x.formula());
    x.formula() = result_formula;
    static_cast<Derived&>(*this).leave(x);
  }

  void update(pbes_system::pbes& x)
  {
    static_cast<Derived&>(*this).enter(x);
    static_cast<Derived&>(*this).update(x.equations());
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::propositional_variable_instantiation& x)
  {
    static_cast<Derived&>(*this).enter(x);
    // skip
    static_cast<Derived&>(*this).leave(x);
    result = atermpp::down_cast<T>(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::not_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_not(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.operand()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::and_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_and(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::or_& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_or(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::imp& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_imp(result, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.left()); }, [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.right()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::forall& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_forall(result, x.variables(), [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::exists& x)
  {
    static_cast<Derived&>(*this).enter(x);
    pbes_system::make_exists(result, x.variables(), [&](pbes_expression& result){ static_cast<Derived&>(*this).apply(result, x.body()); });
    static_cast<Derived&>(*this).leave(x);
  }

  template <class T>
  void apply(T& result, const pbes_system::pbes_expression& x)
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
    else if (pbes_system::is_propositional_variable_instantiation(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::propositional_variable_instantiation>(x));
    }
    else if (pbes_system::is_not(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::not_>(x));
    }
    else if (pbes_system::is_and(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::and_>(x));
    }
    else if (pbes_system::is_or(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::or_>(x));
    }
    else if (pbes_system::is_imp(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::imp>(x));
    }
    else if (pbes_system::is_forall(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::forall>(x));
    }
    else if (pbes_system::is_exists(x))
    {
      static_cast<Derived&>(*this).apply(result, atermpp::down_cast<pbes_system::exists>(x));
    }
    static_cast<Derived&>(*this).leave(x);
  }

};

/// \\brief Builder class
template <typename Derived>
struct pbes_expression_builder: public add_pbes_expressions<pbes_system::pbes_expression_builder_base, Derived>
{
};
//--- end generated add_pbes_expressions code ---//

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_BUILDER_H
