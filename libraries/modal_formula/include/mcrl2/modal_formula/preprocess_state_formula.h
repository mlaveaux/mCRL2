// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/preprocess_state_formula.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_PREPROCESS_STATE_FORMULA_H
#define MCRL2_MODAL_FORMULA_PREPROCESS_STATE_FORMULA_H

#include "mcrl2/data/detail/find.h"
#include "mcrl2/data/xyz_identifier_generator.h"
#include "mcrl2/modal_formula/has_name_clashes.h"
#include "mcrl2/modal_formula/is_monotonous.h"
#include "mcrl2/modal_formula/normalize.h"
#include "mcrl2/modal_formula/resolve_name_clashes.h"
#include "mcrl2/modal_formula/state_formula_rename.h"

namespace mcrl2::state_formulas
{

namespace detail
{

struct count_modal_operator_nesting_traverser: public state_formula_traverser<count_modal_operator_nesting_traverser>
{
  using super = state_formula_traverser<count_modal_operator_nesting_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  std::size_t result = 0;
  std::vector<std::size_t> nesting_depth;

  count_modal_operator_nesting_traverser()
  {
    enter_scope();
  }

  void enter_scope()
  {
    nesting_depth.push_back(0);
  }

  void leave_scope()
  {
    nesting_depth.pop_back();
  }

  void increase_nesting_depth()
  {
    nesting_depth.back()++;
    if (nesting_depth.back() > result)
    {
      result = nesting_depth.back();
    }
  }

  void decrease_nesting_depth()
  {
    nesting_depth.back()--;
  }

  void enter(const must&)
  {
    increase_nesting_depth();
  }

  void leave(const must&)
  {
    decrease_nesting_depth();
  }

  void enter(const may&)
  {
    increase_nesting_depth();
  }

  void leave(const may&)
  {
    decrease_nesting_depth();
  }

  void enter(const mu&)
  {
    enter_scope();
  }

  void leave(const mu&)
  {
    leave_scope();
  }

  void enter(const nu&)
  {
    enter_scope();
  }

  void leave(const nu&)
  {
    leave_scope();
  }

  void enter(const forall&)
  {
    enter_scope();
  }

  void leave(const forall&)
  {
    leave_scope();
  }

  void enter(const exists&)
  {
    enter_scope();
  }

  void leave(const exists&)
  {
    leave_scope();
  }
};

inline
std::size_t count_modal_operator_nesting(const state_formula& x)
{
  count_modal_operator_nesting_traverser f;
  f.apply(x);
  return f.result;
}

/// \brief Traverser that detects unscoped modal operators.
struct has_unscoped_modal_formula_traverser: public state_formula_traverser<has_unscoped_modal_formula_traverser>
{
  using super = state_formula_traverser<has_unscoped_modal_formula_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  bool result = false;

  /// \brief Stack of the fixpoint operators (mu/nu) that enclose the current node.
  /// When it is empty, any modal operator that is visited is unscoped.
  std::vector<state_formula> fixpoints;

  void push(const state_formula& x)
  {
    fixpoints.push_back(x);
  }

  void pop()
  {
    fixpoints.pop_back();
  }

  void enter(const must&)
  {
    if (fixpoints.empty())
    {
      result = true;
    }
  }

  void enter(const may&)
  {
    if (fixpoints.empty())
    {
      result = true;
    }
  }

  void enter(const mu& x)
  {
    push(x);
  }

  void leave(const mu&)
  {
    pop();
  }

  void enter(const nu& x)
  {
    push(x);
  }

  void leave(const nu&)
  {
    pop();
  }
};

/// \brief Checks whether a state formula contains an unscoped modal operator.
/// A modal operator (must '[a]' or may '<a>') is <em>unscoped</em> when, on the path from the root of
/// the formula to that operator, it is not preceded by a fixpoint operator (mu or nu). In other words,
/// the modal operator is not contained in the scope of any fixpoint. For example '[a]true' and
/// 'true => \<a\>true' contain an unscoped modal operator, whereas in 'mu X. [a]\<b\>true' both modal
/// operators are scoped by the surrounding 'mu X'.
inline
bool has_unscoped_modal_formulas(const state_formula& x)
{
  has_unscoped_modal_formula_traverser f;
  f.apply(x);
  return f.result;
}

/// Visitor that transforms state formulas. This can be useful if the state formula contains nested modal operators.
template <typename IdentifierGenerator>
struct state_formula_preprocess_nested_modal_operators_builder: public state_formulas::state_formula_builder<state_formula_preprocess_nested_modal_operators_builder<IdentifierGenerator> >
{
  using super = state_formulas::state_formula_builder<
      state_formula_preprocess_nested_modal_operators_builder<IdentifierGenerator>>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  /// \brief An identifier generator
  IdentifierGenerator& generator;
  std::vector<state_formula> fixpoints;

  /// \brief Constructor
  /// \param generator A generator for fresh identifiers
  state_formula_preprocess_nested_modal_operators_builder(IdentifierGenerator& generator)
    : generator(generator)
  {}

  void push(const state_formula& x)
  {
    fixpoints.push_back(x);
  }

  void pop()
  {
    fixpoints.pop_back();
  }

  bool is_mu()
  {
    if (fixpoints.empty())
    {
      return true;
    }
    return state_formulas::is_mu(fixpoints.back());
  }

  void enter(const mu& x)
  {
    push(x);
  }

  void leave(const mu&)
  {
    pop();
  }

  void enter(const nu& x)
  {
    push(x);
  }

  void leave(const nu&)
  {
    pop();
  }

  template <class T>
  void apply(T& result, const must& x)
  {
    state_formula operand;
    apply(operand, x.operand());
    if (!has_unscoped_modal_formulas(operand))
    {
      make_must(result, x.formula(), operand);
      return;
    }
    core::identifier_string X = generator("X");
    if (is_mu())
    {
      make_must(result, x.formula(), mu(X, {}, operand));
      return;
    }
    else
    {
      make_must(result, x.formula(), nu(X, {}, operand));
      return;
    }
  }

  template <class T>
  void apply(T& result, const may& x)
  {
    state_formula operand;
    apply(operand, x.operand());
    if (!has_unscoped_modal_formulas(operand))
    {
      make_may(result, x.formula(), operand);
      return;
    }
    core::identifier_string X = generator("X");
    if (is_mu())
    {
      make_may(result, x.formula(), mu(X, {}, operand));
      return;
    }
    else
    {
      make_may(result, x.formula(), nu(X, {}, operand));
      return;
    }
  }
};

/// \brief Utility function for creating a state_formula_preprocess_nested_modal_operators_builder.
/// \param generator A generator for fresh identifiers
/// \return a state_formula_preprocess_nested_modal_operators_builder
template <typename IdentifierGenerator>
state_formula_preprocess_nested_modal_operators_builder<IdentifierGenerator> make_state_formula_preprocess_nested_modal_operators_builder(IdentifierGenerator& generator)
{
  return state_formula_preprocess_nested_modal_operators_builder<IdentifierGenerator>(generator);
}

} // namespace detail

/// \brief Preprocesses a state formula that contains (nested) modal operators
/// \param x A modal formula
inline
state_formula preprocess_nested_modal_operators(const state_formula& x)
{
  data::set_identifier_generator generator;
  std::set<core::identifier_string> ids = state_formulas::find_identifiers(x);
  generator.add_identifiers(ids);
  state_formula result;
  detail::make_state_formula_preprocess_nested_modal_operators_builder(generator).apply(result, x);
  return result;
}

/// \brief Renames data variables and predicate variables in the formula \p f, and
///        wraps the formula inside a 'nu' if needed. This is needed as a preprocessing
///        step for the algorithm.
/// \param formula A modal formula.
/// \param context_ids A set of identifier strings.
/// \param preprocess_modal_operators A boolean indicating that dummy fixed point symbols can be
///                                   inserted which makes subsequent handling easier.
/// \param warn_for_modal_operator_nesting A boolean enabling warnings for modal operator nesting.
/// \return The preprocessed formula.
inline
state_formulas::state_formula preprocess_state_formula(const state_formulas::state_formula& formula,
                                                       const std::set<core::identifier_string>& context_ids,
                                                       bool preprocess_modal_operators,
                                                       bool quantitative = false,
                                                       bool warn_for_modal_operator_nesting = true
                                                      )
{
  state_formulas::state_formula f = formula;

  if (!state_formulas::is_monotonous(f))
  {
    throw mcrl2::runtime_error("The formula " + state_formulas::pp(f) + " is not monotonous!");
  }

  if (!preprocess_modal_operators && warn_for_modal_operator_nesting && state_formulas::detail::count_modal_operator_nesting(formula) >= 3)
  {
    mCRL2log(log::info) <<
      "Warning: detected nested modal operators. This may result in a long execution time.\n"
      "Use the option -m (for lps2pbes/lps2pres), -p (for lts2pbes/lts2pres) or insert dummy fix \n"
      "point operators in between manually to speed up the transformation." << std::endl;
  }

  mCRL2log(log::debug) << "Formula before preprocessing: " << f << ".\n";

  // rename data variables in f, to prevent name clashes with data variables in the context
  std::set<core::identifier_string> ids = state_formulas::find_identifiers(f);
  ids.insert(context_ids.begin(), context_ids.end());
  f = state_formulas::rename_variables(f, ids);

  // rename predicate variables in f, to prevent name clashes
  data::xyz_identifier_generator xyz_generator;
  xyz_generator.add_identifiers(state_formulas::find_identifiers(f));
  f = rename_predicate_variables(f, xyz_generator);

  mCRL2log(log::debug) << "Formula after renaming variables: " << f << ".\n";

  // add dummy fixpoints between nested modal operators
  if (preprocess_modal_operators)
  {
    state_formula f0 = f;
    f = state_formulas::preprocess_nested_modal_operators(f);
    if (f0 != f)
    {
      mCRL2log(log::debug) << "Formula after inserting dummy fix points between modal operators:  " << f << ".\n";
    }
  }

  // remove occurrences of ! and =>
  if (!state_formulas::is_normalized(f))
  {
    f = state_formulas::normalize(f, quantitative); // true indicates that the formula is quantitative.
    mCRL2log(log::debug) << "Formula after normalization:  " << f << ".\n";
    assert(state_formulas::is_normalized(f));
  }

  // wrap the formula inside a 'nu' if needed
  if (!state_formulas::is_mu(f) && !state_formulas::is_nu(f))
  {
    data::set_identifier_generator generator;
    generator.add_identifiers(state_formulas::find_identifiers(f));
    core::identifier_string X = generator("X");
    f = state_formulas::nu(X, data::assignment_list(), f);
    mCRL2log(log::debug) << "Formula after wrapping the formula inside a 'nu':  " << f << ".\n";
  }

  // resolve name clashes like mu X(n: Nat). forall n: Nat
  if (state_formulas::has_data_variable_name_clashes(f))
  {
    f = resolve_state_formula_data_variable_name_clashes(f, context_ids);
    mCRL2log(log::debug) << "formula after removing data variable name clashes:  " << f << std::endl;
  }

  mCRL2log(log::debug) << "formula after preprocessing:  " << f << std::endl;

  return f;
}

} // namespace mcrl2::state_formulas

#endif // MCRL2_MODAL_FORMULA_PREPROCESS_STATE_FORMULA_H
