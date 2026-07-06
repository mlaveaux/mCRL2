// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/has_name_clashes.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_HAS_NAME_CLASHES_H
#define MCRL2_MODAL_FORMULA_HAS_NAME_CLASHES_H

#include "mcrl2/modal_formula/traverser.h"

namespace mcrl2::state_formulas {

namespace detail
{

/// \brief Traverser that checks for name clashes in nested mu's/nu's.
class state_variable_name_clash_checker: public state_formulas::state_formula_traverser<state_variable_name_clash_checker>
{
  public:
    using super = state_formulas::state_formula_traverser<state_variable_name_clash_checker>;

    using super::apply;
    using super::enter;
    using super::leave;

    /// \brief The stack of names.
    std::vector<core::identifier_string> m_name_stack;

    /// \brief Pops the stack
    void pop()
    {
      m_name_stack.pop_back();
    }

    /// \brief Pushes name on the stack.
    void push(const core::identifier_string& name)
    {
      using utilities::detail::contains;
      if (contains(m_name_stack, name))
      {
        throw mcrl2::runtime_error("nested propositional variable " + std::string(name) + " clashes");
      }
      m_name_stack.push_back(name);
    }

    void enter(const mu& x)
    {
      push(x.name());
    }

    void leave(const mu&)
    {
      pop();
    }

    void enter(const nu& x)
    {
      push(x.name());
    }

    void leave(const nu&)
    {
      pop();
    }
};

/// \brief Traverser that checks for name clashes in parameters of nested mu's/nu's and forall/exists.
class state_formula_data_variable_name_clash_checker: public state_formulas::state_formula_traverser<state_formula_data_variable_name_clash_checker>
{
  public:
    using super = state_formulas::state_formula_traverser<state_formula_data_variable_name_clash_checker>;

    using super::apply;
    using super::enter;
    using super::leave;

    std::set<core::identifier_string> m_names;

    // throws an exception if name was already present in m_names
    void insert(const core::identifier_string& name, const state_formula& x)
    {
      auto p = m_names.insert(name);
      if (!p.second)
      {
        throw mcrl2::runtime_error("Data parameter " + data::pp(name) + " in subformula " + state_formulas::pp(x) + " clashes with a data parameter in an enclosing formula.");
      }
    }

    void erase(const core::identifier_string& name)
    {
      m_names.erase(name);
    }

    void enter(const mu& x)
    {
      for (const data::assignment& a: x.assignments())
      {
        insert(a.lhs().name(), x);
      }
    }

    void leave(const mu& x)
    {
      for (const data::assignment& a: x.assignments())
      {
        erase(a.lhs().name());
      }
    }

    void enter(const nu& x)
    {
      for (const data::assignment& a: x.assignments())
      {
        insert(a.lhs().name(), x);
      }
    }

    void leave(const nu& x)
    {
      for (const data::assignment& a: x.assignments())
      {
        erase(a.lhs().name());
      }
    }

    void enter(const forall& x)
    {
      for (const data::variable& v: x.variables())
      {
        insert(v.name(), x);
      }
    }

    void leave(const forall& x)
    {
      for (const data::variable& v: x.variables())
      {
        erase(v.name());
      }
    }

    void enter(const exists& x)
    {
      for (const data::variable& v: x.variables())
      {
        insert(v.name(), x);
      }
    }

    void leave(const exists& x)
    {
      for (const data::variable& v: x.variables())
      {
        erase(v.name());
      }
    }
};

} // namespace detail

/// \brief Throws an exception if the formula contains name clashes
inline
void check_state_variable_name_clashes(const state_formula& x)
{
  detail::state_variable_name_clash_checker checker;
  checker.apply(x);
}

/// \brief Returns true if the formula contains name clashes
inline
bool has_state_variable_name_clashes(const state_formula& x)
{
  try
  {
    check_state_variable_name_clashes(x);
  }
  catch (const mcrl2::runtime_error&)
  {
    return true;
  }
  return false;
}

/// \brief Throws an exception if the formula contains name clashes in the parameters of mu/nu/exists/forall
inline
void check_data_variable_name_clashes(const state_formula& x)
{
  detail::state_formula_data_variable_name_clash_checker checker;
  checker.apply(x);
}

/// \brief Returns true if the formula contains parameter name clashes
inline
bool has_data_variable_name_clashes(const state_formula& x)
{
  try
  {
    check_data_variable_name_clashes(x);
  }
  catch (const mcrl2::runtime_error&)
  {
    return true;
  }
  return false;
}

} // namespace mcrl2::state_formulas

#endif // MCRL2_MODAL_FORMULA_HAS_NAME_CLASHES_H
