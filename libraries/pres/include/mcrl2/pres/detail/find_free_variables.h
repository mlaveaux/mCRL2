// Author(s): Jan Friso Groote. Based on pbes/detail/find_free_variables.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/detail/find_free_variables.h
/// \brief Add your file description here.

#ifndef MCRL2_PRES_DETAIL_FIND_FREE_VARIABLES_H
#define MCRL2_PRES_DETAIL_FIND_FREE_VARIABLES_H

#include "mcrl2/pres/traverser.h"

namespace mcrl2::pres_system::detail
{

struct find_free_variables_traverser: public pres_expression_traverser<find_free_variables_traverser>
{
  using super = pres_expression_traverser<find_free_variables_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  data::variable_list bound_variables;
  std::vector<data::variable_list> quantifier_stack;
  std::set<data::variable> result;
  bool search_propositional_variables;

  find_free_variables_traverser(bool search_propositional_variables_ = true)
    : search_propositional_variables(search_propositional_variables_)
  {}

  find_free_variables_traverser(const data::variable_list& bound_variables_, bool search_propositional_variables_ = true)
    : bound_variables(bound_variables_), search_propositional_variables(search_propositional_variables_)
  {}

  /// \brief Returns true if v is an element of bound_variables or quantifier_stack
  /// \param v A data variable
  /// \return True if v is an element of bound_variables or quantifier_stack
  bool is_bound(const data::variable& v) const
  {
    if (std::find(bound_variables.begin(), bound_variables.end(), v) != bound_variables.end())
    {
      return true;
    }
    for (const data::variable_list& vars: quantifier_stack)
    {
      if (std::find(vars.begin(), vars.end(), v) != vars.end())
      {
        return true;
      }
    }
    return false;
  }

  /// \brief Pushes v on the stack of quantifier variables
  /// \param v A sequence of data variables
  void push(const data::variable_list& v)
  {
    quantifier_stack.push_back(v);
  }

  /// \brief Pops the stack of quantifier variables
  void pop()
  {
    quantifier_stack.pop_back();
  }

  void enter(const infimum& x)
  {
    push(x.variables());
  }

  void leave(const infimum&)
  {
    pop();
  }

  void enter(const supremum& x)
  {
    push(x.variables());
  }

  void leave(const supremum&)
  {
    pop();
  }

  void enter(const sum& x)
  {
    push(x.variables());
  }

  void leave(const sum&)
  {
    pop();
  }

  void enter(const propositional_variable& x)
  {
    if (search_propositional_variables)
    {
      for (const data::variable& v: data::find_free_variables(x.parameters()))
      {
        if (!is_bound(v))
        {
          result.insert(v);
        }
      }
    }
  }

  void enter(const data::data_expression& x)
  {
    for (const data::variable& v: data::find_free_variables(x))
    {
      if (!is_bound(v))
      {
        result.insert(v);
      }
    }
  }
};

inline
std::set<data::variable> find_free_variables(const pres_expression& x, const data::variable_list& bound_variables = data::variable_list(), bool search_propositional_variables = true)
{
  find_free_variables_traverser f(bound_variables, search_propositional_variables);
  f.apply(x);
  return f.result;
}

} // namespace mcrl2::pres_system::detail

#endif // MCRL2_PRES_DETAIL_FIND_FREE_VARIABLES_H
