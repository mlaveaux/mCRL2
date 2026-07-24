// Author(s): Jan Friso Groote. Based on pbes/rewriters/simplify_data_rewriter by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/rewriters/simplify_rewriter.h
/// \brief add your file description here.

#ifndef MCRL2_PRES_REWRITERS_SIMPLIFY_REWRITER_H
#define MCRL2_PRES_REWRITERS_SIMPLIFY_REWRITER_H

#include "mcrl2/data/standard.h"
#include "mcrl2/data/detail/enumerator_identifier_generator.h"
#include "mcrl2/data/replace_capture_avoiding_with_an_identifier_generator.h"
#include "mcrl2/pres/rewriters/data_rewriter.h"

namespace mcrl2::pres_system {

namespace detail {

template <template <class> class Builder, class Derived>
struct add_simplify: public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::apply;

protected:
  template <class T>
  void const_multiply_helper(T& result, const data::data_expression& f, const pres_expression& x)
  {
    if (data::sort_real::is_zero(f))
    {
      result=atermpp::down_cast<T>(f);
      return;
    }
    if (data::sort_real::is_one(f))
    {
      apply(result, x);
      return;
    }
    apply(result, x);
    assert(!is_const_multiply_alt(result));
    if (data::sort_real::is_zero(result))
    {
      return;
    }
    if (data::sort_real::is_one(result))
    {
      apply(result, f);
      return;
    }
    if (data::is_data_expression(result))
    { 
      const data::data_expression& d=atermpp::down_cast<data::data_expression>(result);
      if (d.sort()==data::sort_real::real_())
      {
        assert(!data::sort_real::is_one(d));
        result=atermpp::down_cast<T>(data::sort_real::times(f, d));
        return;
      }
      assert(d.sort()==data::sort_bool::bool_());
      return;
    }

    if (is_const_multiply(result))
    {
      const const_multiply& cm=atermpp::down_cast<const_multiply>(result);
      apply(result, const_multiply(data::sort_real::times(f, cm.left()), cm.right()));
      return;
    }
    if ((is_true(result) || is_false(result) || is_eqinf(result) || is_eqninf(result)) &&
        data::sort_real::is_larger_zero(f))
    {
      return;
    }

    make_const_multiply(result, f, result);
    return;  
  }

public:
  add_simplify()=default;

  template <class T>
  void apply(T& result, const minus& x)
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
    if (is_minus(result))
    {
      result = atermpp::down_cast<minus>(result).operand();
      return;
    }
    make_minus(result, result);
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
    pres_expression right;
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
    if (data::is_data_expression(result) && 
        data::is_data_expression(right))
    {
      const data::data_expression& resulta = atermpp::down_cast<data::data_expression>(result);
      const data::data_expression& righta = atermpp::down_cast<data::data_expression>(right);
      if (resulta.sort() == data::sort_real::real_() &&
          righta.sort() == data::sort_real::real_())
      {
        right = atermpp::down_cast<pres_expression>(data::sort_real::minimum(resulta, righta));
        apply(result, right);
        return;
      }
    }
    make_and(result,result, right);
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
    pres_expression right;
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
      const data::data_expression& resulta = atermpp::down_cast<data::data_expression>(result);
      const data::data_expression& righta = atermpp::down_cast<data::data_expression>(right);
      if (resulta.sort() == data::sort_real::real_() &&
          righta.sort() == data::sort_real::real_())
      {
        right = atermpp::down_cast<pres_expression>(data::sort_real::maximum(resulta, righta));
        apply(result, right);
        return;
      }
    }
    make_or(result,result, right);
  }

  template <class T>
  void apply(T& result, const plus& x)
  {
    assert(&result!=&x);  // Result is used as temporary store and cannot match x. 
    apply(result, x.left());
    if (is_true(result))
    {
      result = true_();
      return;
    }
    if (data::sort_real::is_zero(result))
    {
      apply(result, x.right());
      return;
    }
    pres_expression right;
    apply(right, x.right());
    if (is_true(right))
    {
      result = true_();
      return;
    }
    if (is_false(right))
    {
      result=right;
      return;
    }
    if (is_false(result) || data::sort_real::is_zero(right))
    {
      return;
    }
    if (data::is_data_expression(result) && data::is_data_expression(right))
    {
      const data::data_expression& resulta = atermpp::down_cast<data::data_expression>(result);
      const data::data_expression& righta = atermpp::down_cast<data::data_expression>(right);
      if (resulta.sort() == data::sort_real::real_() &&
          righta.sort() == data::sort_real::real_())

      {
        right = atermpp::down_cast<pres_expression>(data::sort_real::plus(resulta, righta));
        apply(result, right);
        return;
      }
    }
    make_plus(result, result, right);
    
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
      if (is_minus(result))
      { 
        result = atermpp::down_cast<minus>(result).operand();
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
      make_minus(result, result);
      return;
    }
    pres_expression left;
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
  void apply(T& result, const infimum& x)
  {
    apply(result, x.body());
    make_optimized_infimum(result, x.variables(), result, true);
  }

  template <class T>
  void apply(T& result, const supremum& x)
  {
    apply(result, x.body());
    make_optimized_supremum(result, x.variables(), result, true);
  }

  template <class T>
  void apply(T& result, const sum& x)
  {
    pres_expression body;
    apply(body, x.body());
    make_sum(result, x.variables(), result);
  }

  template <class T>
  void apply(T& result, const const_multiply& x)
  {
    const_multiply_helper(result, x.left(), x.right());
  }

  template <class T>
  void apply(T& result, const const_multiply_alt& x)
  {
    const_multiply_helper(result, x.right(), x.left());
  }

  template <class T>
  void apply(T& result, const condeq& x)
  {
    apply(result, x.arg1());  
    if (is_true(result))
    {
      apply(result, x.arg3());
      return;
    }
    if (is_false(result)) 
    {
      apply(result, and_(x.arg2(), x.arg3()));
      return;
    }

    pres_expression result1;
    pres_expression result2;
    apply(result1, x.arg2());
    apply(result2, x.arg3());
    if (result1==result2)
    {
      result=result1;
      return;
    }
    make_condeq(result, result, result1, result2);
  }

  template <class T>
  void apply(T& result, const condsm& x)
  {

    apply(result, x.arg1());   
    if (is_true(result))
    {
      apply(result, or_(x.arg2(), x.arg3()));
      return;
    }
    if (is_false(result))
    {
      apply(result, x.arg2());
      return;
    }

    pres_expression result1;
    pres_expression result2;
    apply(result1, x.arg2());
    apply(result2, x.arg3());
    if (result1==result2)
    {
      result=result1;
      return;
    }
    make_condsm(result, result, result1, result2);
  }

  template <class T>
  void apply(T& result, const eqinf& x)
  {
    apply(result, x.operand());
    if (is_true(result))
    {
      result = true_();
      return;
    }
    if (is_false(result))
    {
      result = true_();
      return;
    }
    if (is_and(result))
    {  
       const and_& resulta = atermpp::down_cast<and_>(result);
       pres_expression t1;
       pres_expression t2;
       make_eqinf(t1, resulta.left());
       apply(t2, t1);
       make_eqinf(t1, resulta.right());
       apply(result, t1);
       make_optimized_and(t1, t2, result);
       result = t1;
       return;
    }
    if (is_or(result))
    {  
       const or_& resulta = atermpp::down_cast<or_>(result);
       pres_expression t1;
       pres_expression t2;
       make_eqinf(t1, resulta.left());
       apply(t2, t1);
       make_eqinf(t1, resulta.right());
       apply(result, t1);
       make_optimized_or(t1, t2, result);
       result = t1;
       return;
    }
    if (is_plus(result))
    {  
       const plus& resulta = atermpp::down_cast<plus>(result);
       pres_expression t1;
       pres_expression t2;
       make_eqinf(t1, resulta.left());
       apply(t2, t1);
       make_eqinf(t1, resulta.right());
       apply(result, t1);
       make_optimized_plus(t1, t2, result);
       result = t1;
       return;
    }
    if (is_const_multiply(result))
    {  
       const const_multiply& resulta = atermpp::down_cast<const_multiply>(result);
       pres_expression t1;
       make_eqinf(t1, resulta.right());
       apply(result, t1);
       return;
    }
    if (is_const_multiply_alt(result))
    {  
       const const_multiply_alt& resulta = atermpp::down_cast<const_multiply_alt>(result);
       pres_expression t1;
       make_eqinf(t1, resulta.left());
       apply(result, t1);
       return;
    }
    if (is_eqinf(result))
    {
      return;
    }
    if (is_eqninf(result))
    {
      return;
    }
    if (data::is_application(result) && atermpp::down_cast<data::application>(result).sort()==data::sort_real::real_())
    {
      result = false_();
      return;
    }
    make_eqinf(result, result);
  }

  template <class T>
  void apply(T& result, const eqninf& x)
  {
    apply(result, x.operand());
    if (is_true(result))
    {
      result = true_();
      return;
    }
    if (is_false(result))
    {
      result = true_();
      return;
    }
    if (is_and(result))
    {
       const and_& resulta = atermpp::down_cast<and_>(result);
       pres_expression t1;
       pres_expression t2;
       make_eqninf(t1, resulta.left());
       apply(t2, t1);
       make_eqninf(t1, resulta.right());
       apply(result, t1);
       make_optimized_and(t1, t2, result);
       result = t1;
       return;
    }
    if (is_or(result))
    {
       const or_& resulta = atermpp::down_cast<or_>(result);
       pres_expression t1;
       pres_expression t2;
       make_eqninf(t1, resulta.left());
       apply(t2, t1);
       make_eqninf(t1, resulta.right());
       apply(result, t1);
       make_optimized_or(t1, t2, result);
       result = t1;
       return;
    }
    if (is_plus(result))
    {
       const plus& resulta = atermpp::down_cast<plus>(result);
       pres_expression t;
       pres_expression t1;
       pres_expression t2;
       pres_expression t3;
       pres_expression t4;
       make_eqninf(t, resulta.left());
       apply(t1, t);
       make_eqinf(t, resulta.right());
       apply(t2, t);
       make_eqinf(t, resulta.left());
       apply(t3, t);
       make_eqninf(t, resulta.right());
       apply(t4, t);
       make_optimized_or(t, t1, t2);
       make_optimized_or(t1, t3, t4);
       make_optimized_and(result, t, t1);
       return;
    }
    if (is_const_multiply(result))
    {
       const const_multiply& resulta = atermpp::down_cast<const_multiply>(result);
       pres_expression t1;
       make_eqninf(t1, resulta.right());
       apply(result, t1);
       return;
    }
    if (is_const_multiply_alt(result))
    {
       const const_multiply_alt& resulta = atermpp::down_cast<const_multiply_alt>(result);
       pres_expression t1;
       make_eqninf(t1, resulta.left());
       apply(result, t1);
       return;
    }

    if (is_eqinf(result))
    {
      return;
    }
    if (is_eqninf(result))
    {
      return;
    }
    if (data::is_application(result) && atermpp::down_cast<data::application>(result).sort()==data::sort_real::real_())
    {
      result = true_();
      return;
    }
    make_eqninf(result, result);
  }

};

template <typename Derived>
struct simplify_builder: public add_simplify<pres_system::pres_expression_builder, Derived>
{ };

template <typename Derived, typename DataRewriter, typename Substitution>
struct simplify_data_rewriter_builder : public add_data_rewriter < pres_system::detail::simplify_builder, Derived, DataRewriter, Substitution >
{
  using super = add_data_rewriter<pres_system::detail::simplify_builder, Derived, DataRewriter, Substitution >;

  using super::apply;
  const data::data_specification& m_data_spec;

  using substitution_administration_type = data::detail::substitution_updater_with_an_identifier_generator<Substitution,
                                                                  data::enumerator_identifier_generator>;
  substitution_administration_type substitution_administration;


  
protected:
  template <class T>
  void const_multiply_helper(T& result, const data::data_expression& f, const pres_expression& x)
  {
    data::data_expression data_result;
    apply(data_result, f);
    if (data::sort_real::is_zero(data_result))
    {
      result=atermpp::down_cast<T>(data_result);
      return;
    }
    apply(result, x);
    assert(!is_const_multiply_alt(result));
    assert(!data::sort_real::is_zero(data_result));
    if (data::sort_real::is_one(data_result))
    {
      return;
    }
    if (is_const_multiply(result))
    {
      data::data_expression data_result1;
      const const_multiply& result_cm=atermpp::down_cast<const_multiply>(result);
      apply(data_result1, data::sort_real::times(data_result, result_cm.left()));
      assert(!data::sort_real::is_zero(data_result1));
      if (!data::sort_real::is_one(data_result1))
      {
        make_const_multiply(result, data_result1, result_cm.right());
      }
      else 
      {
        result = result_cm.right();
      }
      return;
    }
    if (data::is_data_expression(result))
    {
      const data::data_expression& d=atermpp::down_cast<data::data_expression>(result);
      if (d.sort()==data::sort_real::real_())
      {
        if (data::sort_real::is_one(d))
        {
          result=atermpp::down_cast<T>(data_result);
          return;
        }
        if (data::sort_real::is_zero(d))
        {
          return;
        }
        data::data_expression data_result1;
        apply(result, data::sort_real::times(data_result, d));

        return;
      }
      else if (data::sort_real::is_larger_zero(data_result))
      {
        assert(d.sort()==data::sort_bool::bool_());
        return;
      }
    }
    if (is_true(result) || is_false(result) || is_eqinf(result) || is_eqninf(result))
    {
      if (data::sort_real::is_larger_zero(data_result))
      {
        assert(data::sort_real::value<double>(data_result)>0);
        return;
      }
      make_const_multiply(result, data_result, result);
      return;
    }

    make_const_multiply(result, data_result, result);
    return;  
  }

public:
  simplify_data_rewriter_builder(const data::data_specification& data_spec, const DataRewriter& R, Substitution& sigma)
    : super(R, sigma),
      m_data_spec(data_spec),
      substitution_administration(sigma, const_cast<data::enumerator_identifier_generator&>(R.identifier_generator()))
  {}

  template <class T>
  void apply(T& result, const infimum& x)
  { 
    const data::variable_list vl = substitution_administration.push(x.variables());
    apply(result, x.body());
    make_optimized_infimum(result, vl, result, true);
    substitution_administration.pop(vl);
  }

  template <class T>
  void apply(T& result, const supremum& x)
  {
    const data::variable_list vl = substitution_administration.push(x.variables());
    apply(result, x.body());
    make_optimized_supremum(result, vl, result, true);
    substitution_administration.pop(vl);
  }   

  template <class T>
  void apply(T& result, const sum& x)
  {
    const data::variable_list vl = substitution_administration.push(x.variables());
    pres_expression body;
    apply(body, x.body());
    make_optimized_sum(result, vl, body, m_data_spec, super::R);
    substitution_administration.pop(vl);
  }

  template <class T>
  void apply(T& result, const const_multiply& x)
  {
    const_multiply_helper(result, x.left(), x.right());
  }

  template <class T>
  void apply(T& result, const const_multiply_alt& x)
  {
    const_multiply_helper(result, x.right(), x.left());
  }

  template <class T>
  void apply(T& result, const condeq& x)
  {
    apply(result, x.arg1());  
    if (is_true(result))
    {
      apply(result, x.arg3());
      return;
    }
    if (is_false(result)) 
    {
      apply(result, and_(x.arg2(), x.arg3()));
      return;
    }
    if (data::is_data_expression(result))
    {
      const data::data_expression& d=atermpp::down_cast<data::data_expression>(result);
      if (d.sort()==data::sort_real::real_())
      {
        pres_expression aux;
        apply(aux, static_cast<pres_expression>(data::less(data::sort_real::real_zero(),d)));   // NOTE: arg1() is sometimes rewritten twice.
     
        if (is_true(aux))
        {
          apply(result, x.arg3());
          return;
        }
        if (is_false(aux))
        {
          apply(result, and_(x.arg2(), x.arg3()));
          return;
        }
      }
    } 

    pres_expression result1;
    pres_expression result2;
    apply(result1, x.arg2());
    apply(result2, x.arg3());
    if (result1==result2)
    {
      result=result1;
      return;
    }
    make_condeq(result, result, result1, result2);
  }

  template <class T>
  void apply(T& result, const condsm& x)
  {

    apply(result, x.arg1());   
    if (is_true(result))
    {
      apply(result, or_(x.arg2(), x.arg3()));
      return;
    }
    if (is_false(result))
    {
      apply(result, x.arg2());
      return;
    }
    if (data::is_data_expression(result))
    { 
      const data::data_expression& d=atermpp::down_cast<data::data_expression>(result);
      if (d.sort()==data::sort_real::real_())
      { 
        pres_expression aux;
        apply(aux, static_cast<pres_expression>(less(d,data::sort_real::real_zero())));   // NOTE: arg1() is sometimes rewritten twice.
        
        if (is_true(aux))
        { 
          apply(result, x.arg2());
          return;
        }
        if (is_false(aux))
        { 
          apply(result, or_(x.arg2(), x.arg3()));
          return;
        }
      }
    } 

    pres_expression result1;
    pres_expression result2;
    apply(result1, x.arg2());
    apply(result2, x.arg3());
    if (result1==result2)
    {
      result=result1;
      return;
    }
    make_condsm(result, result, result1, result2);
  }

};

} // namespace detail

/// \brief A rewriter that simplifies boolean expressions in a term.
struct simplify_rewriter
{
  using term_type = pres_expression;
  using variable_type = data::variable;

  pres_expression operator()(const pres_expression& x) const
  {
    pres_expression result;
    core::make_apply_builder<detail::simplify_builder>().apply(result, x);
    return result;
  }

  void operator()(pres_expression& result, const pres_expression& x) const
  {
    core::make_apply_builder<detail::simplify_builder>().apply(result, x);
  }
};

/// \brief A rewriter that simplifies boolean expressions in a term, and rewrites data expressions.
template <typename DataRewriter>
struct simplify_data_rewriter
{
  using term_type = pres_expression;
  using variable_type = data::variable;

  const DataRewriter& R;
  const data::data_specification& m_dataspec;

  explicit simplify_data_rewriter(const data::data_specification& dataspec, const DataRewriter& R_)
    : R(R_),
      m_dataspec(dataspec)
  {}

  pres_expression operator()(const pres_expression& x) const
  {
    data::no_substitution sigma;
    pres_expression result;
    detail::make_apply_rewriter_builder<pres_system::detail::simplify_data_rewriter_builder>(m_dataspec, R, sigma).apply(result,x);
    return result;
  }

  template <typename Substitution>
  pres_expression operator()(const pres_expression& x, Substitution& sigma) const
  {
    pres_expression result;
    detail::make_apply_rewriter_builder<pres_system::detail::simplify_data_rewriter_builder>(m_dataspec, R, sigma).apply(result, x);
    return result;
  }

  template <typename Substitution>
  void operator()(pres_expression& result, const pres_expression& x, Substitution& sigma) const
  {
    detail::make_apply_rewriter_builder<pres_system::detail::simplify_data_rewriter_builder>(R, sigma).apply(result, x);
  }
};

} // namespace mcrl2::pres_system



#endif // MCRL2_PRES_REWRITERS_SIMPLIFY_REWRITER_H
