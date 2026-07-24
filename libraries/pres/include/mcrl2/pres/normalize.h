// Author(s): Jan Friso Groote. Based on pbes/normalize.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/normalize.h
/// \brief Normalization of pres expressions.

#ifndef MCRL2_PRES_NORMALIZE_H
#define MCRL2_PRES_NORMALIZE_H

#include "mcrl2/data/consistency.h"
#include "mcrl2/pres/builder.h"
#include "mcrl2/pres/traverser.h"



namespace mcrl2::pres_system
{

/// \cond INTERNAL_DOCS
// \brief Visitor for checking if a pres expression is normalized.
struct is_normalized_traverser: public pres_expression_traverser<is_normalized_traverser>
{
  using super = pres_expression_traverser<is_normalized_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  bool result = true;

  /// \brief Visit not node
  void enter(const minus& /* x */)
  {
    result = false;
  }

  /// \brief Visit imp node
  void enter(const imp& /* x */)
  {
    result = false;
  }
};
/// \endcond

/// \cond INTERNAL_DOCS
// \brief Visitor for checking if a pres expression is normalized.
struct normalize_builder: public pres_expression_builder<normalize_builder>
{
  using super = pres_expression_builder<normalize_builder>;
  using super::apply;

  bool negated = false;

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    result =  atermpp::down_cast<T>(negated ? data::not_(x) : x);
  }

  template <class T>
  void apply(T& result, const minus& x)
  {
    negated = !negated;
    super::apply(result, x.operand());
    negated = !negated;
  }

  template <class T>
  void apply(T& result, const and_& x)
  {
    pres_expression left;
    super::apply(left, x.left());
    pres_expression right;
    super::apply(right, x.right());
    if (negated)
    {
      make_or(result, left, right);
    }
    else
    {
      make_and(result, left, right);
    }
  }

  template <class T>
  void apply(T& result, const or_& x)
  {
    pres_expression left;
    super::apply(left, x.left());
    pres_expression right;
    super::apply(right, x.right());
    if (negated)
    {
      make_and(result, left, right);
    }
    else
    {
      make_or(result, left, right);
    }
  }

  template <class T>
  void apply(T& result, const imp& x)
  {
    negated = !negated;
    pres_expression left;
    super::apply(left, x.left());
    negated = !negated;
    pres_expression right;
    super::apply(right, x.right());
    if (negated)
    {
      make_and(result, left, right);
    }
    else
    {
      make_or(result, left, right);
    }
  }

  template <class T>
  void apply(T& result, const infimum& x)
  {
    pres_expression body;
    super::apply(body, x.body());
    result = negated ? optimized_supremum(x.variables(), body) : optimized_infimum(x.variables(), body);
  }

  template <class T>
  void apply(T& result, const supremum& x)
  {
    pres_expression body;
    super::apply(body, x.body());
    result = negated ? optimized_infimum(x.variables(), body) : optimized_supremum(x.variables(), body);
  }

  template <class T>
  void apply(T& result, const propositional_variable_instantiation& x)
  {
    if (negated)
    {
      throw mcrl2::runtime_error(std::string("normalize error: illegal argument ") + pp(x));
    }
    result = x;
  }
};
/// \endcond

/// \brief Checks if a pres expression is normalized
/// \param x A PRES expression
/// \return True if the pres expression is normalized
template <typename T>
bool is_normalized(const T& x)
{
  is_normalized_traverser f;
  f.apply(x);
  return f.result;
}

/// \brief The function normalize brings (embedded) pres expressions into positive normal form,
/// i.e. a formula without any occurrences of ! or =>.
/// \param x an object containing pres expressions
template <typename T>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void normalize(T& x)
{
  normalize_builder f;
  f.update(x);
}

/// \brief The function normalize brings (embedded) pres expressions into positive normal form,
/// i.e. a formula without any occurrences of ! or =>.
/// \param x an object containing pres expressions
template <typename T>
  requires(std::is_base_of_v<atermpp::aterm, T>)
T normalize(const T& x)
{
  T result;
  normalize_builder f;
  f.apply(result, x);
  return result;
}

} // namespace mcrl2::pres_system



#endif // MCRL2_PRES_NORMALIZE_H
