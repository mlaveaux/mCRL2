// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/negate_variables.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_NEGATE_VARIABLES_H
#define MCRL2_MODAL_FORMULA_NEGATE_VARIABLES_H

#include "mcrl2/modal_formula/builder.h"

namespace mcrl2::state_formulas
{

namespace detail
{

/// Visitor that negates propositional variable instantiations with a given name.
template <typename Derived>
struct state_variable_negator: public state_formulas::state_formula_builder<Derived>
{
  using super = state_formulas::state_formula_builder<Derived>;
  using super::apply;

  core::identifier_string m_name;
  bool m_quantitative;

  state_variable_negator(const core::identifier_string& name, bool quantitative)
    : m_name(name),
      m_quantitative(quantitative)
  {}

  /// \brief Visit variable node.
  /// \param x A term.
  /// \return The result of visiting the node.
  template <class T>
  void apply(T& result, const variable& x)
  {
    if (x.name() == m_name)
    {
      if (m_quantitative)
      {
        state_formulas::make_minus(result, x);
      }
      else
      {
        state_formulas::make_not(result, x);
      }
      return;
    }
    result = x;
  }
};

} // namespace detail

inline
/// \brief Negates variable instantiations in a state formula with a given name.
/// \param name The name of the variables that should be negated.
/// \param x The state formula. 
state_formula negate_variables(const core::identifier_string& name, bool quantitative, const state_formula& x)
{
  state_formula result;
  core::make_apply_builder_arg2<detail::state_variable_negator>(name, quantitative).apply(result, x);
  return result;
}

} // namespace mcrl2::state_formulas

#endif // MCRL2_MODAL_FORMULA_NEGATE_VARIABLES_H
