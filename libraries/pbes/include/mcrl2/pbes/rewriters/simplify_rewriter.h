// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/rewriters/simplify_rewriter.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_REWRITERS_SIMPLIFY_REWRITER_H
#define MCRL2_PBES_REWRITERS_SIMPLIFY_REWRITER_H

#include "mcrl2/data/detail/enumerator_identifier_generator.h"
#include "mcrl2/data/replace_capture_avoiding_with_an_identifier_generator.h"
#include "mcrl2/data/optimized_boolean_operators.h"
#include "mcrl2/pbes/pbes_expression.h"
#include "mcrl2/pbes/replace.h"
#include "mcrl2/pbes/rewriters/data_rewriter.h"
#include "mcrl2/pbes/substitutions.h"

namespace mcrl2::pbes_system {

namespace detail {

template <template <class> class Builder, class Derived>
struct add_simplify: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::apply;

  template <class T>
  void apply(T& result, const not_& x)
  {
    assert(&result!=&x);  // Result is used as temporary store and cannot match x.
    apply(result, x.operand());
    if (is_false(result))
    {
      result = true_();
      return;
    }
    if (is_true(result))
    {
      result = false_();
      return;
    }
    if (is_not(result))
    {
      result = atermpp::down_cast<not_>(result).operand();
      return;
    }
    make_not(result, result);
  }

  template <class T>
  void apply(T& result, const and_& x)
  {
    assert(&result!=&x);  // Result is used as temporary store and cannot match x.
    apply(result, x.left());
    if (is_false(result))
    {
      result = false_();
      return;
    }
    if (is_true(result))
    {
      apply(result, x.right());
      return;
    }
    pbes_expression right;
    apply(right, x.right());
    if (is_false(right))
    {
      result = false_();
      return;
    }
    if (is_true(right) || result==right)
    {
      return;
    }
    if (data::is_data_expression(result) && data::is_data_expression(right))
    {
      data::sort_bool::make_and(atermpp::assign_cast<data::data_expression>(result),
                                     atermpp::down_cast<data::data_expression>(result),
                                     atermpp::down_cast<data::data_expression>(right));
      return;
    }
    make_and(result, result, right);
  }

  template <class T>
  void apply(T& result, const or_& x)
  {
    assert(&result!=&x);  // Result is used as temporary store and cannot match x.
    apply(result, x.left());
    if (is_true(result))
    {
      result = true_();
      return;
    }
    if (is_false(result))
    {
      apply(result, x.right());
      return;
    }
    pbes_expression right;
    apply(right, x.right());
    if (is_true(right))
    {
      result = true_();
      return;
    }
    if (is_false(right) || result==right)
    {
      return;
    }
    if (data::is_data_expression(result) && data::is_data_expression(right))
    {
      data::sort_bool::make_or(atermpp::assign_cast<data::data_expression>(result),
                                    atermpp::down_cast<data::data_expression>(result),
                                    atermpp::down_cast<data::data_expression>(right));
      return;
    }
    make_or(result, result, right);
  }

  template <class T>
  void apply(T& result, const imp& x)
  {
    assert(&result!=&x);  // Result is used as temporary store and cannot match x.
    if (is_false(x.left()))  // This test is cheap.
    {
      result = true_();
      return;
    }
    apply(result, x.right());
    if (is_true(result))
    {
      result = true_();
      return;
    }
    if (is_false(result))
    {
      apply(result, x.left());
      if (is_not(result))
      {
        result = atermpp::down_cast<not_>(result).operand();
        return;
      }
      if (is_true(result))
      {
        result = false_();
        return;
      }
      if (is_false(result))
      {
        result = true_();
        return;
      }
      make_not(result, result);
      return;
    }
    pbes_expression left;
    apply(left, x.left());
    if (is_true(left))
    {
      return;
    }
    if (is_false(left) || result==left)
    {
      result = true_();
      return;
    }
    make_imp(result, left, result);
  }

  template <class T>
  void apply(T& result, const forall& x)
  {
    apply(result, x.body());
    pbes_system::make_optimized_forall(result, x.variables(), result, true);
  }

  template <class T>
  void apply(T& result, const exists& x)
  {
    apply(result, x.body());
    pbes_system::make_optimized_exists(result, x.variables(), result, true);
  }
};

template <typename Derived>
struct simplify_builder: public add_simplify<pbes_system::pbes_expression_builder, Derived>
{ };

template <typename Derived, typename DataRewriter, typename MutableSubstitution, typename PbesSubstitution>
struct simplify_data_rewriter_builder : public add_data_rewriter<pbes_system::detail::simplify_builder, Derived, DataRewriter, MutableSubstitution>
{
  using super = add_data_rewriter<pbes_system::detail::simplify_builder, Derived, DataRewriter, MutableSubstitution>;
  using super::enter;
  using super::leave;
  using super::apply;

  using substitution_administration_type = data::detail::substitution_updater_with_an_identifier_generator<MutableSubstitution,
                                                                  data::enumerator_identifier_generator>;
  substitution_administration_type substitution_administration;
  const PbesSubstitution& m_sigma_pbes;

  simplify_data_rewriter_builder(const DataRewriter& R, MutableSubstitution& sigma)
    : super(R, sigma),
      substitution_administration(sigma, const_cast<data::enumerator_identifier_generator&>(R.identifier_generator()))
  {}

  simplify_data_rewriter_builder(const DataRewriter& R, MutableSubstitution& sigma, const PbesSubstitution& sigma_pbes)
    : super(R, sigma),
      substitution_administration(sigma, const_cast<data::enumerator_identifier_generator&>(R.identifier_generator())),
      m_sigma_pbes(sigma_pbes)
  {}

  MutableSubstitution& substitution()
  {
    return substitution_administration.substitution();
  }


  template<class T, typename SubstitutionType>
  void apply_substitution(T& result, const propositional_variable_instantiation& x, const SubstitutionType& sigma)
  {
    result = sigma(x);
  }

  /// Overload to have a trivial function in case the substitution is no_substitution.
  template<class T>
  void apply_substitution(T&, const propositional_variable_instantiation&, const no_substitution&)
  {
    return;
  }

  template<class T>
  void apply(T& result, const propositional_variable_instantiation& x)
  {
    super::apply(result, x);
    apply_substitution(result, atermpp::down_cast<propositional_variable_instantiation>(result), m_sigma_pbes);
  }

  template <class T>
  void apply(T& result, const forall& x)
  {
    const data::variable_list vl = substitution_administration.push(x.variables());
    super::apply(result, x.body());
    pbes_system::make_optimized_forall(result, vl, result, true);
    substitution_administration.pop(vl);
  }

  template <class T>
  void apply(T& result, const exists& x)
  {
    const data::variable_list vl = substitution_administration.push(x.variables());
    super::apply(result, x.body());
    pbes_system::make_optimized_exists(result, vl, result, true);
    substitution_administration.pop(vl);
  }

};

template<template<class, class, class, class> class Builder,
  class DataRewriter,
  class MutableSubstitution,
  class PbesSubstitution = no_substitution>
struct apply_data_rewriter_with_pbes_substitution_builder
  : public Builder<
      apply_data_rewriter_with_pbes_substitution_builder<Builder, DataRewriter, MutableSubstitution, PbesSubstitution>,
      DataRewriter,
      MutableSubstitution,
      PbesSubstitution>
{
  using super = Builder<
    apply_data_rewriter_with_pbes_substitution_builder<Builder, DataRewriter, MutableSubstitution, PbesSubstitution>,
    DataRewriter,
    MutableSubstitution,
    PbesSubstitution>;
  using super::enter;
  using super::leave;

  apply_data_rewriter_with_pbes_substitution_builder(const DataRewriter& R, MutableSubstitution& sigma, const PbesSubstitution& sigma_pbes)
    : super(R, sigma, sigma_pbes)
  {}

  apply_data_rewriter_with_pbes_substitution_builder(const DataRewriter& R,
    MutableSubstitution& sigma)
    : super(R, sigma, no_substitution())
  {}
};

template<template<class, class, class, class> class Builder,
  class DataRewriter,
  class MutableSubstitution,
  class PbesSubstitution>
apply_data_rewriter_with_pbes_substitution_builder<Builder, DataRewriter, MutableSubstitution, PbesSubstitution>
make_apply_data_rewriter_with_pbes_substitution_builder(const DataRewriter& R,
  MutableSubstitution& sigma,
  const PbesSubstitution& sigma_pbes)
{
  return apply_data_rewriter_with_pbes_substitution_builder<Builder, DataRewriter, MutableSubstitution, PbesSubstitution>(R,
    sigma,
    sigma_pbes);
}

template<template<class, class, class, class> class Builder,
  class DataRewriter,
  class MutableSubstitution>
apply_data_rewriter_with_pbes_substitution_builder<Builder, DataRewriter, MutableSubstitution, no_substitution>
make_apply_data_rewriter_with_pbes_substitution_builder(const DataRewriter& R,
  MutableSubstitution& sigma)
{
  return apply_data_rewriter_with_pbes_substitution_builder<Builder,
    DataRewriter,
    MutableSubstitution,
    no_substitution>(R, sigma);
}

} // namespace detail

/// \brief A rewriter that simplifies boolean expressions in a term.
struct simplify_rewriter
{
  using term_type = pbes_expression;
  using variable_type = data::variable;

  pbes_expression operator()(const pbes_expression& x) const
  {
    pbes_expression result;
    core::make_apply_builder<detail::simplify_builder>().apply(result, x);
    return result;
  }

  void operator()(pbes_expression& result, const pbes_expression& x) const
  {
    core::make_apply_builder<detail::simplify_builder>().apply(result, x);
  }

  // sigma maps propositional variable instantiations to PBES expression.
  template<class PbesSubstitution>
  pbes_expression operator()(const pbes_expression& x, const PbesSubstitution& sigma) const
  {
    pbes_expression result;
    pbes_system::detail::make_replace_propositional_variables_builder<detail::simplify_builder>(sigma).apply(result, x);
    return result;
  }

  // sigma maps propositional variable instantiations to PBES expression.
  template <class PbesSubstitution>
  void operator()(pbes_expression& result, const pbes_expression& x, const PbesSubstitution& sigma) const
  {
    assert(&result != &x); // precondition for simplify_builder.apply
    pbes_system::detail::make_replace_propositional_variables_builder<detail::simplify_builder>(sigma)
      .apply(result, x);
  }
};

/// \brief A rewriter that simplifies boolean expressions in a term, and rewrites data expressions using DataRewriter.
template <typename DataRewriter>
struct simplify_data_rewriter
{
  using term_type = pbes_expression;
  using variable_type = data::variable;

  const DataRewriter& R;

  explicit simplify_data_rewriter(const DataRewriter& R_)
    : R(R_)
  {}

  pbes_expression operator()(const pbes_expression& x) const
  {
    data::no_substitution sigma;
    pbes_expression result;
    detail::make_apply_data_rewriter_with_pbes_substitution_builder<
      pbes_system::detail::simplify_data_rewriter_builder>(R, sigma)
      .apply(result, x);
    return result;
  }

  template <typename Substitution>
  pbes_expression operator()(const pbes_expression& x, Substitution& sigma) const
  {
    pbes_expression result;
    detail::make_apply_data_rewriter_with_pbes_substitution_builder<
      pbes_system::detail::simplify_data_rewriter_builder>(R, sigma)
      .apply(result, x);
    return result;
  }

  template <typename Substitution>
  void operator()(pbes_expression& result, const pbes_expression& x, Substitution& sigma) const
  {
    detail::make_apply_data_rewriter_with_pbes_substitution_builder<
      pbes_system::detail::simplify_data_rewriter_builder>(R, sigma)
      .apply(result, x);
  }

  template<typename Substitution, typename PbesSubstitution>
  pbes_expression operator()(const pbes_expression& x,
    Substitution& sigma,
    const PbesSubstitution& sigma_pbes) const
  {
    pbes_expression result;
    detail::make_apply_data_rewriter_with_pbes_substitution_builder<
      pbes_system::detail::simplify_data_rewriter_builder>(R, sigma, sigma_pbes)
      .apply(result, x);
    return result;
  }

  template<typename Substitution, typename PbesSubstitution>
  void operator()(pbes_expression& result, const pbes_expression& x, Substitution& sigma, const PbesSubstitution& sigma_pbes) const
  {
    detail::make_apply_data_rewriter_with_pbes_substitution_builder<pbes_system::detail::simplify_data_rewriter_builder>(R, sigma, sigma_pbes)
      .apply(result, x);
  }
};

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_REWRITERS_SIMPLIFY_REWRITER_H
