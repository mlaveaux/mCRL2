// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/expression_traits.h
/// \brief Contains term traits for data_expression.

#ifndef MCRL2_DATA_EXPRESSION_TRAITS_H
#define MCRL2_DATA_EXPRESSION_TRAITS_H

#include "mcrl2/core/term_traits.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/exists.h"
#include "mcrl2/data/forall.h"
#include "mcrl2/data/detail/data_sequence_algorithm.h"

namespace mcrl2
{

namespace core
{

/// \brief Contains type information for data expressions.
template <>
struct term_traits<data::data_expression>
{
  /// \brief The term type
  using term_type = data::data_expression;

  /// \brief The variable type
  using variable_type = data::variable;

  /// \brief The variable sequence type
  using variable_sequence_type = data::variable_list;

  /// \brief The value true
  /// \return The value true
  static inline
  const term_type& true_()
  {
    return data::sort_bool::true_();
  }

  /// \brief The value false
  /// \return The value false
  static inline
  const term_type& false_()
  {
    return data::sort_bool::false_();
  }

  /// \brief Operator not
  /// \param p A term
  /// \return Operator not applied to p
  static inline
  term_type not_(const term_type& p)
  {
    return data::sort_bool::not_(p);
  }

  /// \brief Operator not
  /// \param result Operator not applied to p
  /// \param p A term
  static inline
  void make_not(term_type& result, const term_type& p)
  {
    data::sort_bool::make_not(result, p);
  }

  /// \brief Operator and
  /// \param p A term
  /// \param q A term
  /// \return Operator and applied to p and q
  static inline
  term_type and_(const term_type& p, const term_type& q)
  {
    return data::sort_bool::and_(p, q);
  }

  /// \brief Operator and
  /// \param result Operator and applied to p and q
  /// \param p A term
  /// \param q A term
  static inline
  void make_and(term_type& result, const term_type& p, const term_type& q)
  {
    data::sort_bool::make_and(result, p, q);
  }

  /// \brief Operator or
  /// \param p A term
  /// \param q A term
  /// \return Operator or applied to p and q
  static inline
  term_type or_(const term_type& p, const term_type& q)
  {
    return data::sort_bool::or_(p, q);
  }

  /// \brief Operator or
  /// \param result Operator or applied to p and q
  /// \param p A term
  /// \param q A term
  static inline
  void make_or(term_type& result, const term_type& p, const term_type& q)
  {
    data::sort_bool::make_or(result, p, q);
  }

  /// \brief Operator imp
  /// \param p A term
  /// \param q A term
  /// \return Operator or applied to p and q
  static inline
  term_type imp(const term_type& p, const term_type& q)
  {
    return data::sort_bool::implies(p, q);
  }

  /// \brief Operator imp
  /// \param result Operator and applied to p and q
  /// \param p A term
  /// \param q A term
  static inline
  void make_imp(term_type& result, const term_type& p, const term_type& q)
  {
    data::sort_bool::make_implies(result, p, q);
  }

  /// \brief Operator forall
  /// \param d A sequence of variables
  /// \param p A term
  /// \return Operator forall applied to d and p
  static inline
  term_type forall(const variable_sequence_type& d, const term_type& p)
  {
    return data::forall(d, p);
  }

  /// \brief Construct a forall.
  /// \param result Place where the forall is put. 
  /// \param d A sequence of variables
  /// \param p A term
  static inline
  void make_forall(term_type& result, const variable_sequence_type& d, const term_type& p)
  {
    data::make_forall(result, d, p);
  }

  /// \brief Operator exists
  /// \param d A sequence of variables
  /// \param p A term
  /// \return Operator exists applied to d and p
  static inline
  term_type exists(const variable_sequence_type& d, const term_type& p)
  {
    return data::exists(d, p);
  }

  /// \brief Construct an exists.
  /// \param result Place where the forall is put. 
  /// \param d A sequence of variables
  /// \param p A term
  static inline
  void make_exists(term_type& result, const variable_sequence_type& d, const term_type& p)
  {
    data::make_exists(result, d, p);
  }

  /// \brief Test for value true
  /// \param t A term
  /// \return True if the term has the value true
  static inline
  bool is_true(const term_type& t)
  {
    return t == data::sort_bool::true_();
  }

  /// \brief Test for value false
  /// \param t A term
  /// \return True if the term has the value false
  static inline
  bool is_false(const term_type& t)
  {
    return t == data::sort_bool::false_();
  }

  /// \brief Test for operator not
  /// \param t A term
  /// \return True if the term is of type not
  static inline
  bool is_not(const term_type& t)
  {
    return data::sort_bool::is_not_application(t);
  }

  /// \brief Test for operator and
  /// \param t A term
  /// \return True if the term is of type and
  static inline
  bool is_and(const term_type& t)
  {
    return data::sort_bool::is_and_application(t);
  }

  /// \brief Test for operator or
  /// \param t A term
  /// \return True if the term is of type or
  static inline
  bool is_or(const term_type& t)
  {
    return data::sort_bool::is_or_application(t);
  }

  /// \brief Test for implication
  /// \param t A term
  /// \return True if the term is an implication
  static inline
  bool is_imp(const term_type& t)
  {
    return data::sort_bool::is_implies_application(t);
  }

  /// \brief Test for universal quantification
  /// \param t A term
  /// \return True if the term is an universal quantification
  static inline
  bool is_forall(const term_type& t)
  {
    return data::is_forall(t);
  }

  /// \brief Test for existential quantification
  /// \param t A term
  /// \return True if the term is an existential quantification
  static inline
  bool is_exists(const term_type& t)
  {
    return data::is_exists(t);
  }

  /// \brief Test for lambda abstraction
  /// \param t A term
  /// \return True if the term is a lambda expression
  static inline
  bool is_lambda(const term_type& t)
  {
    return data::is_lambda(t);
  }

  /// \brief Conversion from variable to term
  /// \param v A variable
  /// \return The converted variable
  static inline
  const term_type& variable2term(const variable_type& v)
  {
    return v;
  }

  /// \brief Test if a term is a variable
  /// \param t A term
  /// \return True if the term is a variable
  static inline
  bool is_variable(const term_type& t)
  {
    return data::is_variable(t);
  }

  /// \brief Get the n-th argument of a data expression, provided it is an application.
  /// \param t A term which is an application.
  /// \param n The index of the argument. The first index has number 0.
  /// \return the n-th argument of t.
  /// \details This function is linear in n.
  static inline
  const term_type& argument(const term_type& t, const std::size_t n)
  {
    assert(data::is_application(t));
    const data::application& a=atermpp::down_cast<data::application>(t);
    assert(a.size()>n);
    data::application::const_iterator i=a.begin();
    for(std::size_t j=0; j<n; ++j, ++i)
    {
      assert(i!=a.end());
    }
    assert(i!=a.end());
    return *i;
  }

  static inline
  const term_type& left(const term_type& t)
  {
    assert(data::is_application(t));
    const data::application& a=atermpp::down_cast<data::application>(t);
    assert(a.size() == 2);
    return *(a.begin());
  }

  static inline
  const term_type& right(const term_type& t)
  {
    assert(data::is_application(t));
    const data::application& a=atermpp::down_cast<data::application>(t);
    assert(a.size() == 2);
    return *(++(a.begin()));
  }

  static inline
  const term_type& not_arg(const term_type& t)
  {
    assert(is_not(t));
    assert(data::is_application(t));
    const data::application& a=atermpp::down_cast<data::application>(t);
    assert(a.size() == 1);
    return *(a.begin());
  }

  /// \brief Pretty print function
  /// \param t A term
  /// \return A pretty print representation of the term
  static inline
  std::string pp(const term_type& t)
  {
    return data::pp(t);
  }

  // \brief A function that calculates the intersection of two variable lists.
  // \param x A list.
  // \param y Another list.
  // \return The intersection of the two lists.
  static inline
  data::variable_list set_intersection(const data::variable_list& x, const std::set<data::variable>& y)
  {
    return data::detail::set_intersection(x, y);
  }

  // \brief Find the free variables in a term.
  // \param x The input term.
  // \return The set of free variables in a term.
  static inline
  std::set<data::variable> find_free_variables(const data::data_expression& x)
  {
    return data::find_free_variables(x);
  }
};

} // namespace core

} // namespace mcrl2

#endif // MCRL2_DATA_EXPRESSION_TRAITS_H
