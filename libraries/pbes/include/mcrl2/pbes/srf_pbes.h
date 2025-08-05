// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/srf_pbes.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_SRF_PBES_H
#define MCRL2_PBES_SRF_PBES_H

#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/core/detail/print_utility.h"
#include "mcrl2/data/data_expression.h"
#include "mcrl2/pbes/detail/pbes_remove_counterexample_info.h"
#include "mcrl2/pbes/find.h"
#include "mcrl2/pbes/join.h"
#include "mcrl2/pbes/pbes_equation.h"
#include "mcrl2/pbes/pbes_expression.h"
#include "mcrl2/pbes/pbes_functions.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/rewriters/pbes2data_rewriter.h"
#include "mcrl2/utilities/logger.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2
{

namespace pbes_system
{

namespace detail
{

inline pbes_expression make_not(const pbes_expression& x)
{
  if (data::is_data_expression(x) && data::sort_bool::is_not_application(x))
  {
    return data::sort_bool::arg(atermpp::down_cast<data::data_expression>(x));
  }
  else if (is_not(x))
  {
    return accessors::arg(x);
  }
  return not_(x);
}

template <typename Iterator>
pbes_expression make_conjunction(Iterator first, Iterator last)
{
  if (first == last)
  {
    return true_();
  }

  Iterator it = first;
  pbes_expression result = *it;
  ++it;

  while (it != last)
  {
    make_and_(result, result, *it);
    ++it;
  }

  return result;
}

template <typename Iterator>
pbes_expression make_disjunction(Iterator first, Iterator last)
{
  if (first == last)
  {
    return false_();
  }

  Iterator it = first;
  pbes_expression result = *it;
  ++it;

  while (it != last)
  {
    make_or_(result, result, *it);
    ++it;
  }

  return result;
}

/// @brief Assuming that x is of the form (x1 && x2 && ... xn), generates
///        expression !x1 || !x2 || ... !xn
/// @param x a pbes_expression that is conjunctive
/// @return disjunctive expression that is equivalent to !x.
inline pbes_expression distribute_not_over_and(const pbes_expression& x)
{
  std::set<pbes_expression> conjuncts = split_and(x);
  std::set<pbes_expression> disjuncts;
  // negate each of the expressions
  for (const pbes_expression& conjunct : conjuncts)
  {
    disjuncts.insert(make_not(conjunct));
  }
  return make_disjunction(disjuncts.begin(), disjuncts.end());
}

/// A summand in a srf equation, if allow_ce is true then counter example information is
/// allowed in the condition. This means that the condition contains propositional variables.
template <bool allow_ce>
class pre_srf_summand
{
  using condition_type = std::conditional_t<allow_ce, pbes_expression, data::data_expression>;
protected:
  data::variable_list m_parameters;
  condition_type m_condition;
  propositional_variable_instantiation m_X;

public:
  pre_srf_summand(data::variable_list parameters,
      const pbes_expression& condition,
      propositional_variable_instantiation X)
      : m_parameters(std::move(parameters)),
        m_X(std::move(X))
  {
    if constexpr (!allow_ce)
    {
      m_condition = detail::pbes2data(condition);
    }
    else
    {
      m_condition = condition;
    }
  }

  const data::variable_list& parameters() const { return m_parameters; }

  data::variable_list& parameters() { return m_parameters; }

  const condition_type& condition() const { return m_condition; }

  condition_type& condition() { return m_condition; }

  const propositional_variable_instantiation& variable() const { return m_X; }

  propositional_variable_instantiation& variable() { return m_X; }

  void add_variables(const data::variable_list& variables) { m_parameters = variables + m_parameters; }

  void add_condition(const pbes_expression& f)
  {
    if constexpr (allow_ce)
    {
      m_condition = pbes_system::and_(f, m_condition);
    }
    else
    {
      m_condition = data::and_(pbes2data(f), m_condition);
    }
  }

  pbes_expression to_pbes(bool conjunctive) const
  {
    if (conjunctive)
    {
      // if counterexample information is allowed, we may have generated conjuncts
      // with guards of the form !X, where X is a variable that encodes counterexample
      // information. As other algorithms expect expressions in positive form, we need
      // to push the negation inside.
      return make_forall_(m_parameters, or_(distribute_not_over_and(m_condition), m_X));
    }
    else
    {
      return make_exists_(m_parameters, and_(m_condition, m_X));
    }
  }
};

template <bool allow_ce>
inline std::ostream& operator<<(std::ostream& out, const pre_srf_summand<allow_ce>& summand)
{
  return out << "variables = " << core::detail::print_list(summand.parameters()) << " f = " << summand.condition()
             << " X = " << summand.variable();
}

template <bool allow_ce>
class pre_srf_equation
{
protected:
  fixpoint_symbol m_sigma;
  propositional_variable m_variable;
  std::vector<pre_srf_summand<allow_ce>> m_summands;
  bool m_conjunctive;

public:
  explicit pre_srf_equation(const fixpoint_symbol& sigma,
      const propositional_variable& variable,
      std::vector<pre_srf_summand<allow_ce>> summands,
      bool conjunctive)
      : m_sigma(sigma),
        m_variable(variable),
        m_summands(std::move(summands)),
        m_conjunctive(conjunctive)
  {}

  const fixpoint_symbol& symbol() const { return m_sigma; }

  fixpoint_symbol& symbol() { return m_sigma; }

  const propositional_variable& variable() const { return m_variable; }

  propositional_variable& variable() { return m_variable; }

  const std::vector<pre_srf_summand<allow_ce>>& summands() const { return m_summands; }

  std::vector<pre_srf_summand<allow_ce>>& summands() { return m_summands; }

  bool& is_conjunctive() { return m_conjunctive; }

  bool is_conjunctive() const { return m_conjunctive; }

  pbes_equation to_pbes() const
  {
    if (is_counter_example_variable(variable()))
    {
      if (is_counter_example_positive(variable().name()))
      {
        return pbes_equation(m_sigma, m_variable, true_());
      }
      else
      {
        return pbes_equation(m_sigma, m_variable, false_());
      }
    }

    std::vector<pbes_expression> v;
    for (const auto& summand : m_summands)
    {
      v.push_back(summand.to_pbes(m_conjunctive));
    }
    pbes_expression rhs = m_conjunctive ? join_and(v.begin(), v.end()) : join_or(v.begin(), v.end());
    return pbes_equation(m_sigma, m_variable, rhs);
  }

  /// Ensures that the equation is total, by adding the summands corresponding to true and false
  void make_total(const pre_srf_summand<allow_ce>& true_summand, const pre_srf_summand<allow_ce>& false_summand)
  {
    if (m_conjunctive)
    {
      m_summands.push_back(true_summand);
    }
    else
    {
      m_summands.push_back(false_summand);
    }
  }
};

template <bool allow_ce>
inline std::ostream& operator<<(std::ostream& out, const pre_srf_equation<allow_ce>& eqn)
{
  out << "srf equation" << std::endl;
  for (const auto& summand : eqn.summands())
  {
    out << summand << std::endl;
  }
  return out;
}

template <bool allow_ce>
std::vector<pre_srf_summand<allow_ce>> srf_or(const pbes_expression& phi,
    std::deque<pbes_equation>& equations,
    const pbes_equation& eqn,
    const data::variable_list& V,
    data::set_identifier_generator& id_generator,
    const core::identifier_string& X_true,
    const core::identifier_string& X_false,
    std::vector<pre_srf_equation<allow_ce>>& result,
    bool merge_simple_expressions);

template <bool allow_ce>
struct srf_or_traverser : public pbes_expression_traverser<srf_or_traverser<allow_ce>>
{
  using super = pbes_expression_traverser<srf_or_traverser>;
  using super::apply;
  using super::enter;
  using super::leave;

  // The remaining PBES equations
  std::deque<pbes_equation>& equations;

  // The current equation
  const pbes_equation& eqn;

  const data::variable_list& V;

  // Used for creating new equations
  data::set_identifier_generator& id_generator;

  // The names of the true and false equations
  const core::identifier_string& X_true;
  const core::identifier_string& X_false;

  // The equations of the resulting srf PBES
  std::vector<pre_srf_equation<allow_ce>>& result;

  // The summands of the generated equation
  std::vector<pre_srf_summand<allow_ce>> summands;

  /// If true do not introduce a new PBES equation for simple expressions, and instead add them to every summand
  /// separately.
  bool m_merge_simple_expressions = false;

  srf_or_traverser(std::deque<pbes_equation>& equations_,
      const pbes_equation& eqn_,
      const data::variable_list& V_,
      data::set_identifier_generator& id_generator_,
      const core::identifier_string& X_true_,
      const core::identifier_string& X_false_,
      std::vector<pre_srf_equation<allow_ce>>& result_,
      bool merge_simple_expressions)
      : equations(equations_),
        eqn(eqn_),
        V(V_),
        id_generator(id_generator_),
        X_true(X_true_),
        X_false(X_false_),
        result(result_),
        m_merge_simple_expressions(merge_simple_expressions)
  {}

  void apply(const and_& x)
  {
    if (m_merge_simple_expressions && is_simple_expression(x.left(), allow_ce))
    {
      std::size_t size = summands.size();
      apply(x.right());
      for (auto i = summands.begin() + size; i != summands.end(); ++i)
      {
        i->add_condition(x.left());
      }
    }
    else if (m_merge_simple_expressions && is_simple_expression(x.right(), allow_ce))
    {
      std::size_t size = summands.size();
      apply(x.left());
      for (auto i = summands.begin() + size; i != summands.end(); ++i)
      {
        i->add_condition(x.right());
      }
    }
    else
    {
      const propositional_variable& X = eqn.variable();
      propositional_variable X1(id_generator(X.name()), V);
      const pbes_expression& f = true_();
      equations.emplace_front(eqn.symbol(), X1, x);
      summands.emplace_back(data::variable_list(),
          f,
          propositional_variable_instantiation(X1.name(), data::make_data_expression_list(V)));
    }
  }


  void apply(const or_& x)
  {
    if (!allow_ce)
    {
      return super::apply(x);
    }

    mCRL2log(log::trace) << "Or expression " << x << "\n";

    // Special case for pre-SRF with simple expressions.
    std::set<pbes_expression> clauses = split_or(x, false);

    std::vector<pbes_expression> simple_clauses;
    for (const auto& clause : clauses)
    {
      if (m_merge_simple_expressions && is_simple_expression(clause, false))
      {
        simple_clauses.emplace_back(clause);
        apply(clause);
      }
    }

    if (simple_clauses.empty())
    {
      // No simple clauses, so we can just apply the or_ operator.
      super::apply(x);
      return;
    }

    pbes_expression condition = make_disjunction(simple_clauses.begin(), simple_clauses.end());
    mCRL2log(log::trace) << "Simple condition " << condition << "\n";

    for (const auto& clause : clauses)
    {
      if (!m_merge_simple_expressions || !is_simple_expression(clause, false))
      {
        mCRL2log(log::trace) << "Clause " << clause << "\n";
        std::size_t size = summands.size();
        apply(clause);
        for (auto i = summands.begin() + size; i != summands.end(); ++i)
        {
          i->add_condition(detail::make_not(condition));
        }
      }
    }
  }


  void apply(const exists& x)
  {
    std::vector<pre_srf_summand<allow_ce>> body_summands = srf_or(x.body(),
        equations,
        eqn,
        V + x.variables(),
        id_generator,
        X_true,
        X_false,
        result,
        m_merge_simple_expressions);
    for (auto& summand : body_summands)
    {
      summand.add_variables(x.variables());
    }
    summands.insert(summands.end(), body_summands.begin(), body_summands.end());
  }

  void apply(const forall& x)
  {
    if (is_simple_expression(x.body(), allow_ce))
    {
      const pbes_expression& f = x.body();
      const propositional_variable_instantiation& X = propositional_variable_instantiation(X_true, {});
      summands.emplace_back(x.variables(), f, X);
    }
    else
    {
      const propositional_variable& X = eqn.variable();
      propositional_variable X1(id_generator(X.name()), V);
      const pbes_expression& f = true_();
      equations.emplace_front(eqn.symbol(), X1, x);
      summands.emplace_back(data::variable_list(),
          f,
          propositional_variable_instantiation(X1.name(), data::make_data_expression_list(V)));
    }
  }

  void apply(const propositional_variable_instantiation& x)
  {
    const pbes_expression& f = true_();
    summands.emplace_back(data::variable_list(), f, x);
  }

  void apply(const pbes_expression& x)
  {
    if (is_simple_expression(x, allow_ce))
    {
      const propositional_variable_instantiation& X = propositional_variable_instantiation(X_true, {});
      const pbes_expression& f = x;
      summands.emplace_back(data::variable_list(), f, X);
    }
    else
    {
      super::apply(x);
    }
  }

  void apply(const not_& x) {
    throw mcrl2::runtime_error("srf_or_traverser::apply(not_) unsupported term " + pbes_system::pp(x));
  }

  void apply(const imp& x) {
    throw mcrl2::runtime_error("srf_or_traverser::apply(imp) unsupported term " + pbes_system::pp(x));
  }
};

template<bool allow_ce>
inline std::vector<pre_srf_summand<allow_ce>> srf_or(const pbes_expression& phi,
    std::deque<pbes_equation>& equations,
    const pbes_equation& eqn,
    const data::variable_list& V,
    data::set_identifier_generator& id_generator,
    const core::identifier_string& X_true,
    const core::identifier_string& X_false,
    std::vector<pre_srf_equation<allow_ce>>& result,
    bool merge_simple_expressions)
{
  srf_or_traverser<allow_ce> f(equations, eqn, V, id_generator, X_true, X_false, result, merge_simple_expressions);
  f.apply(phi);
  return std::move(f.summands);
}

template<bool allow_ce>
std::vector<pre_srf_summand<allow_ce>> srf_and(const pbes_expression& phi,
    std::deque<pbes_equation>& equations,
    const pbes_equation& eqn,
    const data::variable_list& V,
    data::set_identifier_generator& id_generator,
    const core::identifier_string& X_true,
    const core::identifier_string& X_false,
    std::vector<pre_srf_equation<allow_ce>>& result,
    bool merge_simple_expressions);

template<bool allow_ce>
struct srf_and_traverser : public pbes_expression_traverser<srf_and_traverser<allow_ce>>
{
  using super = pbes_expression_traverser<srf_and_traverser>;
  using super::apply;
  using super::enter;
  using super::leave;

  // The remaining PBES equations
  std::deque<pbes_equation>& equations;

  // The current equation
  const pbes_equation& eqn;

  const data::variable_list& V;

  // Used for creating new equations
  data::set_identifier_generator& id_generator;

  // The names of the true and false equations
  const core::identifier_string& X_true;
  const core::identifier_string& X_false;

  // The equations of the resulting srf PBES
  std::vector<pre_srf_equation<allow_ce>>& result;

  // The summands of the generated equation
  std::vector<pre_srf_summand<allow_ce>> summands;

  /// If true do not introduce a new PBES equation for simple expressions, and instead add them to every summand
  /// separately.
  bool m_merge_simple_expressions = false;

  srf_and_traverser(std::deque<pbes_equation>& equations_,
      const pbes_equation& eqn_,
      const data::variable_list& V_,
      data::set_identifier_generator& id_generator_,
      const core::identifier_string& X_true_,
      const core::identifier_string& X_false_,
      std::vector<pre_srf_equation<allow_ce>>& result_,
      bool merge_simple_expressions)
      : equations(equations_),
        eqn(eqn_),
        V(V_),
        id_generator(id_generator_),
        X_true(X_true_),
        X_false(X_false_),
        result(result_),
        m_merge_simple_expressions(merge_simple_expressions)
  {}

  void apply(const or_& x)
  {
    if (m_merge_simple_expressions && is_simple_expression(x.left(), allow_ce))
    {
      std::size_t size = summands.size();
      apply(x.right());
      for (auto i = summands.begin() + size; i != summands.end(); ++i)
      {
        i->add_condition(detail::make_not(x.left()));
      }
    }
    else if (m_merge_simple_expressions && is_simple_expression(x.right(), allow_ce))
    {
      std::size_t size = summands.size();
      apply(x.left());
      for (auto i = summands.begin() + size; i != summands.end(); ++i)
      {
        i->add_condition(detail::make_not(x.right()));
      }
    }
    else
    {
      const propositional_variable& X = eqn.variable();
      propositional_variable X1(id_generator(X.name()), V);
      const pbes_expression& f = true_();
      equations.emplace_front(eqn.symbol(), X1, x);
      summands.emplace_back(data::variable_list(),
          f,
          propositional_variable_instantiation(X1.name(), data::make_data_expression_list(V)));
    }
  }

  void apply(const and_& x)
  {
    if (!allow_ce)
    {
      super::apply(x);
      return;
    }

    mCRL2log(log::trace) << "Expression " << x << "\n";

    // Special case for pre-SRF with simple expressions.
    std::set<pbes_expression> clauses = split_and(x, false);

    std::vector<pbes_expression> simple_clauses;
    // Collect simple expressions to strengthen conjuncts of result of
    // recursive calls. We do not include counterexample information to
    // ensure that the PBES remains in positive form.
    for (const auto& clause : clauses)
    {
      if (m_merge_simple_expressions && is_simple_expression(clause, false))
      {
        simple_clauses.emplace_back(clause);
        apply(clause);
      }
    }

    if (simple_clauses.empty())
    {
      // No simple clauses, so we can just apply the or_ operator.
      super::apply(x);
      return;
    }

    // condition used for strengthening guards of dependencies.
    pbes_expression condition = make_conjunction(simple_clauses.begin(), simple_clauses.end());
    mCRL2log(log::trace) << "Simple condition " << condition << "\n";

    // Recursively apply (pre)SRF transformation.
    for (const auto& clause : clauses)
    {
      if (!m_merge_simple_expressions || !is_simple_expression(clause, false))
      {
        mCRL2log(log::trace) << "Clause " << clause << "\n";
        std::size_t size = summands.size();
        apply(clause);
        for (auto i = summands.begin() + size; i != summands.end(); ++i)
        {
          i->add_condition(condition);
        }
      }
    }
  }

  void apply(const forall& x)
  {
    std::vector<pre_srf_summand<allow_ce>> body_summands = srf_and(x.body(),
        equations,
        eqn,
        V + x.variables(),
        id_generator,
        X_true,
        X_false,
        result,
        m_merge_simple_expressions);
    for (auto& summand : body_summands)
    {
      summand.add_variables(x.variables());
    }
    summands.insert(summands.end(), body_summands.begin(), body_summands.end());
  }

  void apply(const exists& x)
  {
    if (is_simple_expression(x.body(), allow_ce))
    {
      const pbes_expression& f = x.body();
      const propositional_variable_instantiation& X = propositional_variable_instantiation(X_true, {});
      summands.emplace_back(x.variables(), f, X);
    }
    else
    {
      const propositional_variable& X = eqn.variable();
      propositional_variable X1(id_generator(X.name()), V);
      const pbes_expression& f = true_();
      equations.emplace_front(eqn.symbol(), X1, x);
      summands.emplace_back(data::variable_list(),
          f,
          propositional_variable_instantiation(X1.name(), data::make_data_expression_list(V)));
    }
  }

  void apply(const propositional_variable_instantiation& x)
  {
    const pbes_expression& f = true_();
    summands.emplace_back(data::variable_list(), f, x);
  }

  void apply(const pbes_expression& x)
  {
    if (is_simple_expression(x, allow_ce))
    {
      const propositional_variable_instantiation& X = propositional_variable_instantiation(X_false, {});
      const pbes_expression& f = x;
      summands.emplace_back(data::variable_list(), detail::make_not(f), X);
    }
    else
    {
      super::apply(x);
    }
  }

  void apply(const not_& x) {
    throw mcrl2::runtime_error("srf_and_traverser::apply(not_) unsupported term " + pbes_system::pp(x));
  }

  void apply(const imp& x) {
    throw mcrl2::runtime_error("srf_and_traverser::apply(imp) unsupported term " + pbes_system::pp(x));
  }
};

template<bool allow_ce>
inline std::vector<pre_srf_summand<allow_ce>> srf_and(const pbes_expression& phi,
    std::deque<pbes_equation>& equations,
    const pbes_equation& eqn,
    const data::variable_list& V,
    data::set_identifier_generator& id_generator,
    const core::identifier_string& X_true,
    const core::identifier_string& X_false,
    std::vector<pre_srf_equation<allow_ce>>& result,
    bool merge_simple_expressions)
{
  srf_and_traverser f(equations, eqn, V, id_generator, X_true, X_false, result, merge_simple_expressions);
  f.apply(phi);
  return std::move(f.summands);
}

inline bool is_conjunctive(const pbes_expression& phi, bool allow_ce)
{
  if (is_simple_expression(phi, allow_ce))
  {
    return false;
  }
  else if (is_propositional_variable_instantiation(phi))
  {
    return false;
  }
  else if (is_or(phi))
  {
    const auto& phi_ = atermpp::down_cast<or_>(phi);
    return (is_simple_expression(phi_.left(), allow_ce) && is_propositional_variable_instantiation(phi_.right()))
           || (is_simple_expression(phi_.right(), allow_ce) && is_propositional_variable_instantiation(phi_.left()));
  }
  else if (is_and(phi))
  {
    const auto& phi_ = atermpp::down_cast<and_>(phi);
    bool result = !((is_simple_expression(phi_.left(), allow_ce) && is_propositional_variable_instantiation(phi_.right()))
             || (is_simple_expression(phi_.right(), allow_ce) && is_propositional_variable_instantiation(phi_.left())));
    return result;
  }
  else if (is_exists(phi))
  {
    return false;
  }
  else if (is_forall(phi))
  {
    return true;
  }
  throw mcrl2::runtime_error("is_conjunctive: unexpected case " + pbes_system::pp(phi));
}

/// explicit representation of a pbes in SRF format
///
/// If allow_ce is true, then counter example expressions are allowed in the
/// PBES. We refer to this as kind of `pre_srf`. When the counter example
/// information is replaced by `true` or `false` in every condition, we obtain a
/// `srf`.
template <bool allow_ce>
class pre_srf_pbes
{
protected:
  data::data_specification m_dataspec;
  std::vector<pre_srf_equation<allow_ce>> m_equations;
  propositional_variable_instantiation m_initial_state;

public:
  pre_srf_pbes() = default;

  pre_srf_pbes(const data::data_specification& dataspec,
      std::vector<pre_srf_equation<allow_ce>> equations,
      propositional_variable_instantiation initial_state)
      : m_dataspec(dataspec),
        m_equations(std::move(equations)),
        m_initial_state(std::move(initial_state))
  {}

  const std::vector<pre_srf_equation<allow_ce>>& equations() const { return m_equations; }

  std::vector<pre_srf_equation<allow_ce>>& equations() { return m_equations; }

  const propositional_variable_instantiation& initial_state() const { return m_initial_state; }

  propositional_variable_instantiation& initial_state() { return m_initial_state; }

  const data::data_specification& data() const { return m_dataspec; }

  data::data_specification& data() { return m_dataspec; }

  pbes to_pbes() const
  {
    std::vector<pbes_equation> v;
    for (const auto& eqn : equations())
    {
      v.push_back(eqn.to_pbes());
    }
    return pbes(m_dataspec, std::set<data::variable>(), v, m_initial_state);
  }

  // Adds extra clauses to the equations to enforce that the PBES is in total SRF format
  // Precondition: the last two equations must be the equations corresponding to false and true
  void make_total()
  {
    // TODO: Remove this hack.
    std::size_t N = m_equations.size();
    const auto& false_summand = m_equations[N - 2].summands().front();
    const auto& true_summand = m_equations[N - 1].summands().front();
    for (std::size_t i = 0; i < N - 2; i++)
    {
      m_equations[i].make_total(true_summand, false_summand);
    }
  }
};

/// \brief Converts a PBES into standard recursive form
/// \pre The pbes p must be normalized
/// \param allow_ce If true, then counter example expressions are allowed in the simple expressions of the PBES.
/// \param merge_simple_expressions If true, then simple expressions are merged into the summands of the SRF equations.
///
/// \return A pre-SRF PBES, which is a PBES in SRF format
template<bool allow_ce>
inline detail::pre_srf_pbes<allow_ce> pbes2pre_srf(const pbes& p, bool merge_simple_expressions = true)
{
  data::set_identifier_generator id_generator;
  for (const core::identifier_string& id : pbes_system::find_identifiers(p))
  {
    id_generator.add_identifier(id);
  }

  core::identifier_string X_false = id_generator("X_false");
  core::identifier_string X_true = id_generator("X_true");
  pbes_equation eqn_false(fixpoint_symbol::mu(),
      propositional_variable(X_false, {}),
      or_(data::sort_bool::false_(), propositional_variable_instantiation(X_false, {})));
  pbes_equation eqn_true(fixpoint_symbol::nu(),
      propositional_variable(X_true, {}),
      propositional_variable_instantiation(X_true, {}));

  const auto& p_equations = p.equations();
  std::deque<pbes_equation> equations(p_equations.begin(), p_equations.end());
  equations.emplace_back(eqn_false);
  equations.emplace_back(eqn_true);

  std::vector<detail::pre_srf_equation<allow_ce>> srf_equations;
  while (!equations.empty())
  {
    pbes_equation eqn = equations.front();
    equations.pop_front();
    bool is_conjunctive = detail::is_conjunctive(eqn.formula(), allow_ce);
    std::vector<detail::pre_srf_summand<allow_ce>> summands = is_conjunctive ? detail::srf_and(eqn.formula(),
        equations,
        eqn,
        eqn.variable().parameters(),
        id_generator,
        X_true,
        X_false,
        srf_equations,
        merge_simple_expressions)
      : detail::srf_or(eqn.formula(),
        equations,
        eqn,
        eqn.variable().parameters(),
        id_generator,
        X_true,
        X_false,
        srf_equations,
        merge_simple_expressions);
    srf_equations.emplace_back(eqn.symbol(), eqn.variable(), summands, is_conjunctive);
  }

  auto result = detail::pre_srf_pbes<allow_ce>(p.data(),
      std::vector<detail::pre_srf_equation<allow_ce>>(srf_equations.begin(), srf_equations.end()),
      p.initial_state());
  return result;
}

} // namespace detail

using srf_summand = detail::pre_srf_summand<false>;

using srf_equation = detail::pre_srf_equation<false>;

/// This is a PBES in SRF format
using srf_pbes = detail::pre_srf_pbes<false>;

/// This is a PBES in pre-SRF format, where counter example expressions, i.e., propositional variables are allowed in
/// simple expressions.
using srf_pbes_with_ce = detail::pre_srf_pbes<true>;


/// \brief Converts a pre-SRF PBES into standard recursive form. Note that the
/// counter example information of the pre_srf_pbes is removed since otherwise
/// the result is not in SRF.
inline srf_pbes pre_srf2srfpbes(const srf_pbes_with_ce& p)
{
  // Used to remove counter example information
  detail::subsitute_counterexample f(true, true);
  simplify_rewriter simplify;

  std::vector<detail::pre_srf_equation<false>> equations;
  for (const auto& equation : p.equations())
  {
    if (!detail::is_counter_example_equation(equation.to_pbes()))
    {
      std::vector<detail::pre_srf_summand<false>> summands;
      for (const auto& summand : equation.summands())
      {
        pbes_expression result;
        f.apply(result, summand.variable());
        propositional_variable_instantiation variable = atermpp::down_cast<propositional_variable_instantiation>(simplify(result));

        f.apply(result, summand.condition());

        summands.emplace_back(summand.parameters(), simplify(result), variable);
      }

      equations.emplace_back(equation.symbol(), equation.variable(), summands, equation.is_conjunctive());
    }
  }

  return srf_pbes(p.data(), equations, p.initial_state());
}


/// \brief Converts a PBES into pre standard recursive form
/// \pre The pbes p must be normalized
inline srf_pbes_with_ce pbes2pre_srf(const pbes& p, bool merge_simple_expressions = true)
{
  return detail::pbes2pre_srf<true>(p, merge_simple_expressions);
}

/// \brief Converts a PBES into standard recursive form
/// \pre The pbes p must be normalized
inline srf_pbes pbes2srf(const pbes& p, bool merge_simple_expressions = true)
{
  return detail::pbes2pre_srf<false>(p, merge_simple_expressions);
}

} // namespace pbes_system

inline bool is_srf(const pbes_system::pbes& pbes, bool merge_simple_expressions = true)
{
  return pbes_system::pbes2srf(pbes, merge_simple_expressions).equations().size() == pbes.equations().size() + 2;
}

} // namespace mcrl2

#endif // MCRL2_PBES_SRF_PBES_H
