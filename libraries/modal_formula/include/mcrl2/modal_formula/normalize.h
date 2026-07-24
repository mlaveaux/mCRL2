// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/normalize.h
/// \brief Add your file description here.

#ifndef MCRL2_MODAL_FORMULA_NORMALIZE_H
#define MCRL2_MODAL_FORMULA_NORMALIZE_H

#include "mcrl2/modal_formula/negate_variables.h"
#include "mcrl2/modal_formula/traverser.h"

namespace mcrl2::state_formulas
{

/// \cond INTERNAL_DOCS
// \brief Visitor for checking if a state formula is normalized.
struct is_normalized_traverser: public state_formula_traverser<is_normalized_traverser>
{
  using super = state_formula_traverser<is_normalized_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  bool result = true;
  
  /// \brief Visit not node
  void enter(const not_& /* x */)
  {
    result = false;
  }

  /// \brief Visit minus node
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

// \brief Visitor for normalizing a state formula.
struct normalize_builder: public state_formula_builder<normalize_builder>
{
  using super = state_formula_builder<normalize_builder>;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

  bool m_quantitative;
  bool m_negated;

  // negated indicates that the formulas is under a negation.
  // quantitative indicates that the formula yields a real value, 
  // in which case forall, exists and not are replaced by supremum, infimum and minus. 
  normalize_builder(bool quantitative, bool negated)
    : m_quantitative(quantitative),
      m_negated(negated)
  {}

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    if (m_quantitative)
    {
      if (m_negated)
      {
        if (x.sort()==data::sort_bool::bool_()) 
        {
          result = atermpp::down_cast<T>(data::sort_bool::not_(x));
        }
        else
        {
          result = atermpp::down_cast<T>(data::sort_real::negate(x));
        }
      }
      else
      {
        result = atermpp::down_cast<T>(x);
      }
    }
    else // x is an ordinary modal formula.
    {
      assert(x.sort()==data::sort_bool::bool_());
      result = atermpp::down_cast<T>(m_negated ? data::sort_bool::not_(x) : x);
    }
  }

  template <class T>
  void apply(T& result, const true_& /*x*/)
  {
    if (m_negated)
    {
      result = false_();
    }
    else
    {
      result = true_();
    }
  }

  template <class T>
  void apply(T& result, const false_& /*x*/)
  {
    if (m_negated)
    {
      result = true_();
    }
    else
    {
      result = false_();
    }
  }

  template <class T>
  void apply(T& result, const not_& x)
  {
    assert(!m_quantitative);
    m_negated=!m_negated;
    apply(result, x.operand());
    m_negated=!m_negated;
  }

  template <class T>
  void apply(T& result, const minus& x)
  {
    assert(m_quantitative);
    m_negated=!m_negated;
    apply(result, x.operand());
    m_negated=!m_negated;
  }

  template <class T>
  void apply(T& result, const and_& x)
  {
    state_formula left;
    state_formula right;
    apply(left, x.left());
    apply(right, x.right()); 
    if (m_negated)
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
    state_formula left;
    state_formula right;
    apply(left, x.left());
    apply(right, x.right()); 
    if (m_negated)
    {
      make_and(result, left, right);
    }
    else
    {
      make_or(result, left, right);
    }
  }

  template <class T>
  void apply(T& result, const plus& x)
  {
    state_formula left;
    state_formula right;
    apply(left, x.left());
    apply(right, x.right()); 
    if (m_negated)
    {
      // The plus operator can be inverted if for instance the operator eq_mininf(x) is present, being equal to -infinity
      // if x is minus infinity, and being plus infinity otherwise. 
      // The dual plus then is: e1 dual+ e2 = eq_mininf(e1) && eq_mininf(e2) && (e1+e2).
      throw mcrl2::runtime_error("Cannot negate the plus operator in quantitative modal formulas: " + pp(x) + ".");
    }
    else
    {
      make_plus(result, left, right);
    }
  }

  template <class T>
  void apply(T& result, const imp& x)
  {
    state_formula y = m_quantitative ?
                        or_(minus(x.left()), x.right()) :
                        or_(not_(x.left()), x.right());
    apply(result, y);
  }

  template <class T>
  void apply(T& result, const forall& x)
  {
    state_formula body;
    apply(body, x.body());
    if (m_negated)
    {
      state_formulas::make_exists(result, x.variables(), body); 
    }
    else
    {
      state_formulas::make_forall(result, x.variables(), body);
    }
  }

  template <class T>
  void apply(T& result, const exists& x)
  {
    state_formula body;
    apply(body, x.body());
    if (m_negated)
    {
      state_formulas::make_forall(result, x.variables(), body); 
    }
    else
    {
      state_formulas::make_exists(result, x.variables(), body);
    }
  }

  template <class T>
  void apply(T& result, const supremum& x)
  {
    state_formula body;
    apply(body, x.body());
    if (m_negated)
    {
      make_infimum(result, x.variables(), body);
    }
    else
    {
      make_supremum(result, x.variables(), body); 
    }
  }

  template <class T>
  void apply(T& result, const infimum& x)
  {
    state_formula body;
    apply(body, x.body());
    if (m_negated)
    {
      make_supremum(result, x.variables(), body); 
    }
    else
    {
      make_infimum(result, x.variables(), body);
    }
  }

  template <class T>
  void apply(T& result, const sum& x)
  {
    state_formula body;
    apply(body, x.body());
    make_sum(result, x.variables(), body);
  }

  template <class T>
  void apply(T& result, const variable& x)
  {
    if (m_negated)
    {
      throw mcrl2::runtime_error(std::string("normalize error: illegal argument ") + pp(x));
    }
    result = x;
  }

  template <class T>
  void apply(T& result, const must& x)
  {
    state_formula operand;
    apply(operand, x.operand());
    if (m_negated)
    {
      make_may(result, x.formula(), operand);
    }
    else
    {
      make_must(result, x.formula(), operand);
    }
  }

  template <class T>
  void apply(T& result, const may& x)
  {
    state_formula operand;
    apply(operand, x.operand());
    if (m_negated)
    {
      make_must(result, x.formula(), operand);
    }
    else
    {
      make_may(result, x.formula(), operand);
    }
  }

  template <class T>
  void apply(T& result, const mu& x)
  {
    state_formula operand;
    if (m_negated)
    {
      apply(operand, negate_variables(x.name(), m_quantitative, x.operand()));
      make_nu(result, x.name(), x.assignments(), operand);
    }
    else
    {
      apply(operand, x.operand());
      make_mu(result, x.name(), x.assignments(), operand);
    }
  }

  template <class T>
  void apply(T& result, const nu& x)
  {
    state_formula operand;
    if (m_negated)
    {
      apply(operand, negate_variables(x.name(), m_quantitative, x.operand()));
      make_mu(result, x.name(), x.assignments(), operand);
    }
    else
    {
      apply(operand, x.operand());
      make_nu(result, x.name(), x.assignments(), operand);
    }
  }

  template <class T>
  void apply(T& result, const delay& x)
  {
    if (m_negated)
    {
      result=yaled();
    }
    else
    {
      result=x;
    }
  }

  template <class T>
  void apply(T& result, const delay_timed& x)
  {
    if (m_negated)
    {
      make_yaled_timed(result, x.time_stamp());
    }
    else
    {
      result=x;
    }
  }

  template <class T>
  void apply(T& result, const yaled& x)
  {
    if (m_negated)
    {
      result=delay();
    }
    else
    {
      result=x;
    }
  }

  template <class T>
  void apply(T& result,  const yaled_timed& x)
  {
    if (m_negated)
    {
      make_delay_timed(result, x.time_stamp());
    }
    else
    {
      result=x;
    }
  }
};
/// \endcond

/// \brief Checks if a state formula is normalized.
/// \param x A PBES expression.
/// \return True if the state formula is normalized.
template <typename T>
bool is_normalized(const T& x)
{
  is_normalized_traverser f;
  f.apply(x);
  return f.result;
}

/// \brief The function normalize brings (embedded) state formulas into positive normal form,
/// i.e. a formula without any occurrences of ! or =>.
/// \param x an object containing state formulas.
/// \param quantitative Indication whether the formula is a quantitative boolean formula.
/// \param negated Indication whether the formula must be interpreted as being negated.
template <typename T>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void normalize(T& x,
    bool quantitative = false,
    bool negated = false)
{
  normalize_builder f(quantitative, negated);
  f.update(x);
}

/// \brief The function normalize brings (embedded) state formulas into positive normal form,
/// i.e. a formula without any occurrences of ! or =>.
/// \param x an object containing state formulas
/// \param quantitative Indication whether the formula is a quantitative boolean formula.
/// \param negated Indication whether the formula must be interpreted as being negated.
template <typename T>
  requires(std::is_base_of_v<atermpp::aterm, T>)
T normalize(const T& x,
    bool quantitative = false,
    bool negated = false)
{
  T result;
  normalize_builder f(quantitative, negated);
  f.apply(result, x);
  return result;
}

} // namespace mcrl2::state_formulas

#endif // MCRL2_MODAL_FORMULA_NORMALIZE_H
