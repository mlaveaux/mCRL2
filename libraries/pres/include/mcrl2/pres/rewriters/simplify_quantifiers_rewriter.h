// Author(s): Jan Friso Groote. Based on pbes/rewriters/simplify_quantifiers_rewriter.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/rewriters/simplify_quantifiers_rewriter.h
/// \brief add your file description here.

#ifndef MCRL2_PRES_REWRITERS_SIMPLIFY_QUANTIFIERS_REWRITER_H
#define MCRL2_PRES_REWRITERS_SIMPLIFY_QUANTIFIERS_REWRITER_H

#include "mcrl2/pres/rewriters/data_rewriter.h"
#include "mcrl2/pres/rewriters/simplify_rewriter.h"

namespace mcrl2::pres_system
{

namespace detail {

template <template <class> class Builder, class Derived>
struct add_simplify_quantifiers: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::apply;

  template <class T>
  void apply(T& result, const infimum& x)
  {
    pres_expression body;
    super::apply(body, x.body());
    const data::variable_list& variables = x.variables();

    if (is_minus(body))
    {
      make_optimized_supremum(result, variables, atermpp::down_cast<minus>(body).operand());
      make_optimized_minus(result, result);
    }
    if (is_and(body))
    {
      const pres_expression& left = atermpp::down_cast<and_>(body).left();
      const pres_expression& right = atermpp::down_cast<and_>(body).right();
      make_optimized_infimum(result, variables, left); 
      pres_expression result_right;
      make_optimized_infimum(result_right, variables, right);
      make_optimized_and(result, result, result_right);
    }
    else if (is_or(body))
    {
      const pres_expression& left = atermpp::down_cast<or_>(body).left();
      const pres_expression& right = atermpp::down_cast<or_>(body).right();
      data::variable_list lv = data::detail::set_intersection(variables, find_free_variables(left));
      data::variable_list rv = data::detail::set_intersection(variables, find_free_variables(right));
      if (lv.empty())
      {
        make_optimized_infimum(result, rv, right);
        make_optimized_or(result, left, result);
      }
      else if (rv.empty())
      {
        make_optimized_infimum(result, lv, left);
        make_optimized_or(result, result, right);
      }
      else
      {
        make_optimized_infimum(result, variables, body);
      }
    }
    else
    {
      make_optimized_infimum(result, variables, body);
    }
  }

  template <class T>
  void apply(T& result, const supremum& x)
  {
    pres_expression body;
    super::apply(body, x.body());
    const data::variable_list& variables = x.variables();

    if (is_minus(body))
    {
      make_optimized_infimum(result, variables, atermpp::down_cast<minus>(body).operand());
      make_optimized_minus(result, result);
    }
    if (is_or(body))
    {
      const pres_expression& left = atermpp::down_cast<or_>(body).left();
      const pres_expression& right = atermpp::down_cast<or_>(body).right();
      make_optimized_supremum(result, variables, left);
      pres_expression result_right;
      make_optimized_supremum(result_right, variables, right);
      make_optimized_or(result, result, result_right);
    }
    else if (is_and(body))
    {
      const pres_expression& left = atermpp::down_cast<and_>(body).left();
      const pres_expression& right = atermpp::down_cast<and_>(body).right();
      data::variable_list lv = data::detail::set_intersection(variables, find_free_variables(left));
      data::variable_list rv = data::detail::set_intersection(variables, find_free_variables(right));
      if (lv.empty())
      {
        make_optimized_supremum(result, rv, right);
        make_optimized_and(result, left, result);
      }
      else if (rv.empty())
      {
        make_optimized_supremum(result, lv, left);
        make_optimized_and(result, right, result);
      }
      else
      {
        make_optimized_supremum(result, variables, body);
      }
    }
    else
    {
      make_optimized_supremum(result, variables, body);
    }
  }

  template <class T>
  void apply(T& result, const sum& x)
  {
    pres_expression body;
    super::apply(body, x.body());
    const data::variable_list& variables = x.variables();
    // TODO: This could be optimised. 
    make_sum(result, variables, body);
  }

};

template <typename Derived>
struct simplify_quantifiers_builder: public add_simplify_quantifiers<pres_system::detail::simplify_builder, Derived>
{ };

template <typename Derived, typename DataRewriter, typename SubstitutionFunction>
struct simplify_quantifiers_data_rewriter_builder: public add_data_rewriter<pres_system::detail::simplify_quantifiers_builder, Derived, DataRewriter, SubstitutionFunction>
{
  using super = add_data_rewriter<pres_system::detail::simplify_quantifiers_builder,
      Derived,
      DataRewriter,
      SubstitutionFunction>;
  using super::enter;
  using super::leave;

  simplify_quantifiers_data_rewriter_builder(const data::data_specification&, const DataRewriter& R, SubstitutionFunction& sigma)
    : super(R, sigma)
  {}
};

} // namespace detail

/// \brief A rewriter that simplifies boolean expressions and quantifiers.
struct simplify_quantifiers_rewriter
{
  using term_type = pres_expression;
  using variable_type = data::variable;

  pres_expression operator()(const pres_expression& x) const
  {
    pres_expression result;
    core::make_apply_builder<detail::simplify_quantifiers_builder>().apply(result, x);
    return result;
  }
};

/// \brief A rewriter that simplifies boolean expressions and quantifiers, and rewrites data expressions.
template <typename DataRewriter>
struct simplify_quantifiers_data_rewriter
{
  using term_type = pres_expression;
  using variable_type = data::variable;

  const data::data_specification& m_data_spec;
  const DataRewriter& m_R;

  simplify_quantifiers_data_rewriter(const data::data_specification& data_spec, const DataRewriter& R)
    : m_data_spec(data_spec),
      m_R(R)
  {}

  pres_expression operator()(const pres_expression& x) const
  {
    data::no_substitution sigma;
    pres_expression result;
    detail::make_apply_rewriter_builder<detail::simplify_quantifiers_data_rewriter_builder>(m_data_spec, m_R, sigma).apply(result, x);
    return result;
  }

  template <typename SubstitutionFunction>
  pres_expression operator()(const pres_expression& x, SubstitutionFunction& sigma) const
  {
    pres_expression result;
    detail::make_apply_rewriter_builder<detail::simplify_quantifiers_data_rewriter_builder>(m_data_spec, m_R, sigma).apply(result, x);
    return result;
  }
};

} // namespace mcrl2::pres_system

#endif // MCRL2_PRES_REWRITERS_SIMPLIFY_QUANTIFIERS_REWRITER_H
