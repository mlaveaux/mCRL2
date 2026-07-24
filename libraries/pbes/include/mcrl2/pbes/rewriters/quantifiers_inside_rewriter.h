// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/rewriters/quantifiers_inside_rewriter.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_REWRITERS_QUANTIFIERS_INSIDE_REWRITER_H
#define MCRL2_PBES_REWRITERS_QUANTIFIERS_INSIDE_REWRITER_H

#include "mcrl2/data/rewriters/quantifiers_inside_rewriter.h"
#include "mcrl2/pbes/builder.h"
#include "mcrl2/pbes/join.h"



namespace mcrl2::pbes_system {

namespace detail {

pbes_expression quantifiers_inside(const pbes_expression& x);
pbes_expression quantifiers_inside_forall(const std::set<data::variable>& variables, const pbes_expression& x);
pbes_expression quantifiers_inside_exists(const std::set<data::variable>& variables, const pbes_expression& x);

template <typename BinaryOperation>
std::tuple<pbes_expression, pbes_expression> compute_Phi_Psi(const std::vector<pbes_expression>& X, const std::set<data::variable>& V, BinaryOperation op, pbes_expression empty_sequence_result)
{
  using utilities::detail::set_difference;
  using utilities::detail::set_intersection;

  std::vector<std::set<data::variable>> vars; // free variables
  for (const pbes_expression& x_j: X)
  {
    vars.push_back(set_intersection(find_free_variables(x_j), V));
  }

  std::vector<std::set<data::variable>>::const_iterator j = std::min_element(vars.begin(), 
                                                                             vars.end(),
                                                                             [&](const std::set<data::variable>& x, const std::set<data::variable>& y)
                                                                             {
                                                                               return x.size() < y.size();
                                                                             }
                                                                            );
  const std::set<data::variable>& Z = *j;

  std::vector<pbes_expression> phi;
  std::vector<pbes_expression> psi;
  for (std::size_t i = 0; i < X.size(); i++)
  {
    if (std::includes(Z.begin(), Z.end(), vars[i].begin(), vars[i].end()))
    {
      phi.push_back(X[i]);
    }
    else
    {
      psi.push_back(X[i]);
    }
  }
  pbes_expression Phi = utilities::detail::join(phi.begin(), phi.end(), op, empty_sequence_result);
  pbes_expression Psi = utilities::detail::join(psi.begin(), psi.end(), op, empty_sequence_result);
  return { Phi, Psi };
}

struct quantifiers_inside_builder: public pbes_expression_builder<quantifiers_inside_builder>
{
  using super = pbes_expression_builder<quantifiers_inside_builder>;
  using super::apply;

  template <class T>
  void apply(T& result, const forall& x)
  {
    const pbes_expression& phi = x.body();
    std::set<data::variable> W = data::detail::make_variable_set(x.variables());
    apply(result, phi);
    result = quantifiers_inside_forall(W, result);
  }

  template <class T>
  void apply(T& result, const exists& x)
  {
    const pbes_expression& phi = x.body();
    std::set<data::variable> W = data::detail::make_variable_set(x.variables());
    apply(result, phi);
    result = quantifiers_inside_exists(W, result);
  }

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    atermpp::assign_cast<data::data_expression>(result) = data::quantifiers_inside_rewrite(x);
  }

};

struct quantifiers_inside_forall_builder: public data_expression_builder<quantifiers_inside_forall_builder>
{
  using super = data_expression_builder<quantifiers_inside_forall_builder>;
  using super::apply;

  const std::set<data::variable>& V;

  explicit quantifiers_inside_forall_builder(const std::set<data::variable>& V_)
      : V(V_)
  {}

  static void make_forall_(pbes_expression& result, const data::variable_list& vars, const pbes_expression& body)
  {
    if (vars.empty())
    {
      result=body;
    }
    else 
    {
      pbes_system::make_forall(result, vars, body);
    }
  }

  // default case
  template <typename T>
  void apply_default(pbes_expression& result, const T& x)
  {
    using utilities::detail::set_intersection;
    std::set<data::variable> W = set_intersection(V, find_free_variables(x));
    make_forall_(result, data::variable_list(W.begin(), W.end()), x);
  }

  template <class T>
  void apply(T& result, const forall& x)
  {
    using utilities::detail::set_union;
    std::set<data::variable> W = data::detail::make_variable_set(x.variables());
    result = quantifiers_inside_forall(set_union(V, W), x.body());
  }

  template <class T>
  void apply(T& result, const not_& x)
  {
    const pbes_expression& phi = x.operand();
    result = not_(quantifiers_inside_exists(V, phi));
  }

  template <class T>
  void apply(T& result, const and_& x)
  {
    const pbes_expression& phi = x.left();
    const pbes_expression& psi = x.right();
    make_and(result, quantifiers_inside_forall(V, phi), quantifiers_inside_forall(V, psi));
  }

  template <class T>
  void apply(T& result, const or_& x)
  {
    using utilities::detail::set_difference;
    using utilities::detail::set_intersection;
    using tr = core::term_traits<pbes_expression>;

    std::vector<pbes_expression> X;
    utilities::detail::split(x, std::back_inserter(X), tr::is_or, tr::left, tr::right);
    const auto [Phi, Psi] = compute_Phi_Psi(X, V, tr::or_, tr::false_());
    if (is_false(Psi))
    {
      make_forall_(result, data::detail::make_variable_list(set_intersection(V, find_free_variables(x))), x);
      return;
    }
    std::set<data::variable> vars_Phi = find_free_variables(Phi);
    std::set<data::variable> vars_Psi = find_free_variables(Psi);
    make_optimized_or(result,
                      quantifiers_inside_forall(set_difference(set_intersection(V, vars_Phi), vars_Psi), Phi),
                      quantifiers_inside_forall(set_difference(set_intersection(V, vars_Psi), vars_Phi), Psi));
    make_forall_(result, 
                 data::detail::make_variable_list(set_intersection(V, set_intersection(vars_Phi, vars_Psi))),
                 result);
  }

  template <class T>
  void apply(T& result, const imp& x)
  {
    result = quantifiers_inside_forall(V, or_(not_(x.left()), x.right()));
  }

  template <class T>
  void apply(T& result, const propositional_variable_instantiation& x)
  {
    apply_default(result, x);
  }

  template <class T>
  void apply(T& result, const exists& x)
  {
    apply_default(result, x);
  }

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    result = atermpp::down_cast<T>(data::detail::quantifiers_inside_forall(V, x));
  }
};

struct quantifiers_inside_exists_builder: public pbes_expression_builder<quantifiers_inside_exists_builder>
{
  using super = pbes_expression_builder<quantifiers_inside_exists_builder>;
  using super::apply;

  const std::set<data::variable>& V;

  explicit quantifiers_inside_exists_builder(const std::set<data::variable>& V_)
      : V(V_)
  {
  }

  static void make_exists_(pbes_expression& result, const data::variable_list& vars, const pbes_expression& body)
  {
    if (vars.empty())
    { 
      result=body;
    }
    else 
    {
      pbes_system::make_exists(result, vars, body);
    }
  }

  // default case
  template <typename T>
  void apply_default(pbes_expression& result, const T& x)
  {
    using utilities::detail::set_intersection;
    std::set<data::variable> W = set_intersection(V, find_free_variables(x));
    make_exists_(result, data::variable_list(W.begin(), W.end()), x);
  }

  template <class T>
  void apply(T& result, const exists& x)
  {
    using utilities::detail::set_union;
    std::set<data::variable> W = data::detail::make_variable_set(x.variables());
    result = quantifiers_inside_exists(set_union(V, W), x.body());
  }

  template <class T>
  void apply(T& result, const not_& x)
  {
    const pbes_expression& phi = x.operand();
    result = not_(quantifiers_inside_forall(V, phi));
  }

  template <class T>
  void apply(T& result, const or_& x)
  {
    const pbes_expression& phi = x.left();
    const pbes_expression& psi = x.right();
    result = or_(quantifiers_inside_exists(V, phi), quantifiers_inside_exists(V, psi));
  }

  template <class T>
  void apply(T& result, const and_& x)
  {
    using utilities::detail::set_difference;
    using utilities::detail::set_intersection;
    using tr = core::term_traits<pbes_expression>;

    std::vector<pbes_expression> X;
    utilities::detail::split(x, std::back_inserter(X), tr::is_and, tr::left, tr::right);
    const auto [Phi, Psi] = compute_Phi_Psi(X, V, tr::and_, tr::true_());
    if (is_true(Psi))
    {
      make_exists_(result, data::detail::make_variable_list(set_intersection(V, find_free_variables(x))), x);
      return;
    }
    std::set<data::variable> vars_Phi = find_free_variables(Phi);
    std::set<data::variable> vars_Psi = find_free_variables(Psi);
    make_optimized_and(result,
                       quantifiers_inside_exists(set_difference(set_intersection(V, vars_Phi), vars_Psi), Phi),
                       quantifiers_inside_exists(set_difference(set_intersection(V, vars_Psi), vars_Phi), Psi));
    make_exists_(result, data::detail::make_variable_list(set_intersection(V, set_intersection(vars_Phi, vars_Psi))),
                          result);
  }

  template <class T>
  void apply(T& result, const imp& x)
  {
    result = quantifiers_inside_exists(V, and_(not_(x.left()), x.right()));
  }

  template <class T>
  void apply(T& result, const propositional_variable_instantiation& x)
  {
    apply_default(result, x);
  }

  template <class T>
  void apply(T& result, const forall& x)
  {
    apply_default(result, x);
  }

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    result = atermpp::down_cast<T>(data::detail::quantifiers_inside_exists(V, x));
  }
};

inline
pbes_expression quantifiers_inside(const pbes_expression& x)
{
  quantifiers_inside_builder f;
  pbes_expression result;
  f.apply(result, x);
  return result;
}

inline
pbes_expression quantifiers_inside_forall(const std::set<data::variable>& variables, const pbes_expression& x)
{
  quantifiers_inside_forall_builder f(variables);
  pbes_expression result;
  f.apply(result, x);
  return result;
}

inline
pbes_expression quantifiers_inside_exists(const std::set<data::variable>& variables, const pbes_expression& x)
{
  quantifiers_inside_exists_builder f(variables);
  pbes_expression result;
  f.apply(result, x);
  return result;
}

} // namespace detail

/// \brief A rewriter that pushes quantifiers inside in a PBES expression.
class quantifiers_inside_rewriter
{
  public:
    /// \brief The term type
    using term_type = pbes_expression;

    /// \brief The variable type
    using variable_type = data::variable;

    /// \brief Rewrites a pbes expression.
    /// \param x A term
    /// \return The rewrite result.
    pbes_expression operator()(const pbes_expression& x) const
    {
      return detail::quantifiers_inside(x);
    }
};

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_REWRITERS_QUANTIFIERS_INSIDE_REWRITER_H
