// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/resolve_name_clashes.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_RESOLVE_NAME_CLASHES_H
#define MCRL2_MODAL_FORMULA_RESOLVE_NAME_CLASHES_H

#include "mcrl2/modal_formula/builder.h"
#include "mcrl2/utilities/detail/container_utility.h"

namespace mcrl2::state_formulas {

namespace detail
{

class mapvector_substitution
{
  using vector_map = std::map<data::variable, std::vector<data::variable>>;

  vector_map m_substitutions;

public:
  static constexpr bool is_identity_substitution=false;

  mapvector_substitution(const vector_map& substitutions)
   : m_substitutions(substitutions)
  {}

  data::data_expression operator()(const data::variable& v) const
  {
    auto i = m_substitutions.find(v);
    if (i == m_substitutions.end())
    {
      return v;
    }
    return i->second.back();
  };
};

template <typename Derived>
class state_variable_name_clash_resolver: public state_formulas::state_formula_builder<Derived>
{
  public:
    using super = state_formulas::state_formula_builder<Derived>;

    using super::enter;
    using super::leave;
    using super::apply;

    using name_map = std::map<core::identifier_string, std::vector<core::identifier_string>>;

    /// \brief The stack of names.
    name_map m_names;

    /// \brief Generator for fresh variable names.
    utilities::number_postfix_generator m_generator;

    /// \brief Pops the name of the stack
    void pop(const core::identifier_string& name)
    {
      m_names[name].pop_back();
    }

    /// \brief Pushes name on the stack.
    void push(const core::identifier_string& name)
    {
      std::vector<core::identifier_string>& names = m_names[name];
      if (names.empty())
      {
        names.push_back(name);
      }
      else
      {
        names.emplace_back(m_generator(std::string(name) + "_"));
      }
    }

    void enter(const mu& x)
    {
      push(x.name());
    }

    void leave(const mu& x)
    {
      pop(x.name());
    }

    void enter(const nu& x)
    {
      push(x.name());
    }

    void leave(const nu& x)
    {
      pop(x.name());
    }

    // Rename variable
    template <class T>
    void apply(T& result, const mu& x)
    {
      enter(x);
      // N.B. If the two lines below are replace by
      //   state_formula result = mu(m_names[x.name()].back(), x.assignments(), (*this)(x.operand()));
      // a memory error occurs with the clang and intel compilers!
      core::identifier_string name = m_names[x.name()].back();
      state_formula f;
      apply(f, x.operand());
      result = mu(name, x.assignments(), f);
      leave(x);
    }

    // Rename variable
    template <class T>
    void apply(T& result, const nu& x)
    {
      enter(x);
      // N.B. If the two lines below are replace by
      //   state_formula result = nu(m_names[x.name()].back(), x.assignments(), (*this)(x.operand()));
      // a memory error occurs with the clang and intel compilers!
      core::identifier_string name = m_names[x.name()].back();
      state_formula f;
      apply(f, x.operand());
      result = nu(name, x.assignments(), f);
      leave(x);
    }

    // Rename variable
    template <class T>
    void apply(T& result, const variable& x)
    {
      result = variable(m_names[x.name()].back(), x.arguments());
    }
};

class state_formula_data_variable_name_clash_resolver: public state_formulas::data_expression_builder<state_formula_data_variable_name_clash_resolver>
{
  public:
    using super = state_formulas::data_expression_builder<state_formula_data_variable_name_clash_resolver>;

    using super::enter;
    using super::leave;
    using super::apply;

    data::set_identifier_generator& generator;
    std::multiset<data::variable> bound_variables;
    std::map<data::variable, std::vector<data::variable>> substitutions;

    explicit state_formula_data_variable_name_clash_resolver(data::set_identifier_generator& generator_)
     : generator(generator_)
    {}

    void insert(const data::variable& v)
    {
      if (utilities::detail::contains(bound_variables, v))
      {
        substitutions[v].emplace_back(generator(v.name()), v.sort());
      }
      bound_variables.insert(v);
    }

    void erase(const data::variable& v)
    {
      // Remove one variable, which must exist. Do not use bound_variables.erase(v) as it removes all variables v. 
      std::multiset<data::variable>::const_iterator var_iter=bound_variables.find(v);
      assert(var_iter!=bound_variables.end()); 
      bound_variables.erase(var_iter); 

      auto i = substitutions.find(v);
      if (i != substitutions.end())
      {
        i->second.pop_back();
        if (i->second.empty())
        {
          substitutions.erase(i);
        }
      }
    }

    data::assignment_list apply_assignments(const data::assignment_list& x)
    {
      mapvector_substitution sigma(substitutions); 
      return data::assignment_list(x.begin(), 
                                   x.end(), 
                                   [&](const data::assignment& a)
                                   {
                                     return data::assignment(atermpp::down_cast<data::variable>(sigma(a.lhs())), data::replace_free_variables(a.rhs(), sigma));
                                   });
    }

    data::variable_list apply_variables(const data::variable_list& x)
    {
      mapvector_substitution sigma(substitutions); 

      return data::variable_list(x.begin(), x.end(), [&](const data::variable& v)
                                   {
                                     return atermpp::down_cast<data::variable>(sigma(v));
                                   }
      );
    }

    template <class T>
    void apply(T& result, const mu& x)
    {
      for (const data::assignment& a: x.assignments())
      {
        insert(a.lhs());
      }
      apply(result, x.operand());
      result = mu(x.name(), apply_assignments(x.assignments()), result);
      for (const data::assignment& a: x.assignments())
      {
        erase(a.lhs());
      }
    }

    template <class T>
    void apply(T& result, const nu& x)
    {
      for (const data::assignment& a: x.assignments())
      {
        insert(a.lhs());
      }
      apply(result, x.operand());
      result = nu(x.name(), apply_assignments(x.assignments()), result);
      for (const data::assignment& a: x.assignments())
      {
        erase(a.lhs());
      }
    }

    template <class T>
    void apply(T& result, const forall& x)
    {
      for (const data::variable& v: x.variables())
      {
        insert(v);
      }
      apply(result,x.body());
      result = forall(apply_variables(x.variables()), result);
      for (const data::variable& v: x.variables())
      {
        erase(v);
      }
    }

    template <class T>
    void apply(T& result, const exists& x)
    {
      for (const data::variable& v: x.variables())
      {
        insert(v);
      }
      apply(result, x.body());
      result = exists(apply_variables(x.variables()), result);
      for (const data::variable& v: x.variables())
      {
        erase(v);
      }
    }

    template <class T>
    void apply(T& result, const action_formulas::forall& x)
    {
      for (const data::variable& v: x.variables())
      {
        insert(v);
      }
      action_formulas::action_formula body;
      apply(body, x.body());
      result = atermpp::down_cast<T>(action_formulas::forall(apply_variables(x.variables()), body));
      for (const data::variable& v: x.variables())
      {
        erase(v);
      }
    }

    template <class T>
    void apply(T& result, const action_formulas::exists& x)
    {
      for (const data::variable& v: x.variables())
      {
        insert(v);
      }
      action_formulas::action_formula body;
      apply(body, x.body());
      result = atermpp::down_cast<T>(action_formulas::exists(apply_variables(x.variables()), body));
      for (const data::variable& v: x.variables())
      {
        erase(v);
      }
    }

    template <class T>
    void apply(T& result, const data::data_expression& x)
    {
      mapvector_substitution sigma(substitutions); 

      result=atermpp::down_cast<T>(data::replace_free_variables(x, sigma));
    }
};

} // namespace detail

/// \brief Resolves name clashes in state variables of formula x
inline
state_formula resolve_state_variable_name_clashes(const state_formula& x)
{
  state_formula result;
  core::make_apply_builder<detail::state_variable_name_clash_resolver>().apply(result, x);
  return result;
}

/// \brief Resolves name clashes in data variables of formula x
inline
state_formula resolve_state_formula_data_variable_name_clashes(const state_formula& x, const std::set<core::identifier_string>& context_ids = std::set<core::identifier_string>())
{
  data::set_identifier_generator generator;
  generator.add_identifiers(state_formulas::find_identifiers(x));
  generator.add_identifiers(context_ids);
  detail::state_formula_data_variable_name_clash_resolver f(generator);
  state_formula result;
  f.apply(result, x);
  return result;
}

} // namespace mcrl2::state_formulas

#endif // MCRL2_MODAL_FORMULA_RESOLVE_NAME_CLASHES_H
