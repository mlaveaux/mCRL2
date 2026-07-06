// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file linear_inequalities_utilities.h
/// \brief Contains utility functions for linear inequalities.


#ifndef MCRL2_LPSREALELM_DETAIL_LINEAR_INEQUALITY_UTILITIES_H
#define MCRL2_LPSREALELM_DETAIL_LINEAR_INEQUALITY_UTILITIES_H


#include "mcrl2/data/linear_inequalities.h"

namespace mcrl2::data::detail
{

inline data_expression negate_inequality(const data_expression& e)
{
  if (is_equal_to_application(e))
  {
    return not_equal_to(data::binary_left(atermpp::down_cast<application>(e)),data::binary_right(atermpp::down_cast<application>(e)));
  }
  if (is_not_equal_to_application(e))
  {
    return equal_to(data::binary_left(atermpp::down_cast<application>(e)),data::binary_right(atermpp::down_cast<application>(e)));
  }
  else if (is_less_application(e))
  {
    return greater_equal(data::binary_left(atermpp::down_cast<application>(e)),data::binary_right(atermpp::down_cast<application>(e)));
  }
  else if (is_less_equal_application(e))
  {
    return data::greater(data::binary_left(atermpp::down_cast<application>(e)),data::binary_right(atermpp::down_cast<application>(e)));
  }
  else if (is_greater_application(e))
  {
    return less_equal(data::binary_left(atermpp::down_cast<application>(e)),data::binary_right(atermpp::down_cast<application>(e)));
  }
  else if (is_greater_equal_application(e))
  {
    return data::less(data::binary_left(atermpp::down_cast<application>(e)),data::binary_right(atermpp::down_cast<application>(e)));
  }
  else
  {
    throw mcrl2::runtime_error("Expression " + data::pp(e) + " is expected to be an inequality over sort Real");
  }
}


/// \brief Determine whether a data expression is an inequality
/// \param e A data expression
/// \return true iff e is a data application of ==, <, <=, > or >= to
///      two arguments.

inline bool is_inequality(const data_expression& e)
{
  return is_equal_to_application(e) || is_less_application(e) ||
         is_less_equal_application(e) || is_greater_application(e) ||
         is_greater_equal_application(e);
}

// Functions below should have been defined in the data library.
inline const data_expression& condition_part(const data_expression& e)
{
  assert(is_if_application(e));
  const data::application& a = atermpp::down_cast<const application>(e);
  data::application::const_iterator i = a.begin();
  return *i;
}

inline const data_expression& then_part(const data_expression& e)
{
  assert(is_if_application(e));
  const data::application& a = atermpp::down_cast<const application>(e);
  data::application::const_iterator i = a.begin();
  return *(++i);
}

inline const data_expression& else_part(const data_expression& e)
{
  assert(is_if_application(e));
  const data::application& a = atermpp::down_cast<const application>(e);
  data::application::const_iterator i = a.begin();
  return *(++(++i));
}



/// \brief Splits a condition in expressions ranging over reals and the others
/// \details Conceptually, the condition is first transformed to disjunctive
///          normal form. For each disjunct, there will be an entry in both
///          resulting std::vectors, where the real conditions are in "real_conditions",
///          and the others in non_real_conditions. If there are conjuncts with
///          both real and non-real variables an exception is thrown. If negate
///          is true the result will be negated.
/// \param e A data expression of sort bool.
/// \param real_condition Those parts of e with only variables over sort Real.
/// \param non_real_condition Those parts of e with only variables not of sort Real.
/// \param negate A boolean variable that indicates whether the result must be negated.
/// \return True when e contains real variables
/// \pre The parameter e must be of sort Bool.

static bool split_condition_aux(
  const data_expression& e,
  std::vector < data_expression_list >& real_conditions,
  std::vector < data_expression_list >& non_real_conditions,
  const bool negate=false)
{
  assert(real_conditions.empty());
  assert(non_real_conditions.empty());

  // In these three cases, we rewrite the expression and call this function recursively
  // with the rewritten expression
  if (sort_bool::is_implies_application(e))
  {
    data_expression rewritten = sort_bool::or_(sort_bool::not_(sort_bool::left(e)), sort_bool::right(e));
    return split_condition_aux(rewritten, real_conditions, non_real_conditions, negate);
  }
  else if (is_if_application(e))
  {
    return split_condition_aux(sort_bool::or_(sort_bool::and_(condition_part(e),then_part(e)),
                                   sort_bool::and_(sort_bool::not_(condition_part(e)),else_part(e))),
                    real_conditions,non_real_conditions,negate);
  }
  else if (sort_bool::is_not_application(e))
  {
    return split_condition_aux(sort_bool::arg(e),real_conditions,non_real_conditions,!negate);
  }

  if(sort_bool::is_and_application(e) || sort_bool::is_or_application(e))
  {
    // Recursive case
    std::vector<data_expression_list> real_conditions_aux1;
    std::vector<data_expression_list> non_real_conditions_aux1;
    bool left_is_real = split_condition_aux(sort_bool::left(e),real_conditions_aux1,non_real_conditions_aux1,negate);
    std::vector<data_expression_list> real_conditions_aux2;
    std::vector<data_expression_list> non_real_conditions_aux2;
    bool right_is_real = split_condition_aux(sort_bool::right(e),real_conditions_aux2,non_real_conditions_aux2,negate);
    if(!left_is_real && !right_is_real)
    {
      // There are no real variables on either side so we can
      // just store the expression e in non_real_conditions
      real_conditions.emplace_back();
      non_real_conditions.push_back(data_expression_list({ negate ? data_expression(sort_bool::not_(e)) : e }));
    }
    else if ((!negate && sort_bool::is_and_application(e))  || (negate && sort_bool::is_or_application(e)))
    {
      // Combine the recursive results (whis are disjunctiosn of conjunctions)
      // of the left and right sides of e (which is a conjunction or negated disjunction)
      // by using the distributivity of && and || to obtain a result
      // which is again a disjunction of conjuctions.
      for (std::vector < data_expression_list >::const_iterator
           i1r=real_conditions_aux1.begin(), i1n=non_real_conditions_aux1.begin() ;
           i1r!=real_conditions_aux1.end(); ++i1r, ++i1n)
      {
        for (std::vector < data_expression_list >::const_iterator
             i2r=real_conditions_aux2.begin(), i2n=non_real_conditions_aux2.begin() ;
             i2r!=real_conditions_aux2.end(); ++i2r, ++i2n)
        {
          real_conditions.push_back(*i1r + *i2r);
          non_real_conditions.push_back(*i1n + *i2n);
        }
      }
    }
    else
    {
      assert((!negate && sort_bool::is_or_application(e))  || (negate && sort_bool::is_and_application(e)));

      // Combine the recursive results of the left and right sides
      // of the disjunction (or negated conjunction) by concatenating
      // them.
      real_conditions.insert(real_conditions.end(), real_conditions_aux1.begin(), real_conditions_aux1.end());
      real_conditions.insert(real_conditions.end(), real_conditions_aux2.begin(), real_conditions_aux2.end());
      non_real_conditions.insert(non_real_conditions.end(), non_real_conditions_aux1.begin(), non_real_conditions_aux1.end());
      non_real_conditions.insert(non_real_conditions.end(), non_real_conditions_aux2.begin(), non_real_conditions_aux2.end());
    }
    return left_is_real || right_is_real;
  }
  else if (is_inequality(e) && (data::binary_left(atermpp::down_cast<application>(e)).sort() == sort_real::real_() || data::binary_right(atermpp::down_cast<application>(e)).sort() == sort_real::real_()))
  {
    // Base case 1: an inequality over real numbers
    std::set < variable > vars=data::find_all_variables(e);
    for (const variable& v: vars)
    {
      if (v.sort() != sort_real::real_())
      {
        throw  mcrl2::runtime_error("Expression " + data::pp(e) + " contains variable " +
                                    data::pp(v) + " not of sort Real.");
      }
    }
    real_conditions.push_back(data_expression_list({ negate ? negate_inequality(e) : e }));
    non_real_conditions.emplace_back();
    return true;
  }
  else
  {
    // Base case 2: an expression not containing real numbers
    // e is assumed to be a non_real expression.
    std::set < variable > vars=data::find_all_variables(e);
    for (const variable& v: vars)
    {
      if (v.sort() == sort_real::real_())
      {
        throw  mcrl2::runtime_error("Expression " + data::pp(e) + " contains variable " +
                                    data::pp(v) + " of sort Real.");
      }
    }
    non_real_conditions.push_back(data_expression_list({ negate ? data_expression(sort_bool::not_(e)) : e }));
    real_conditions.emplace_back();
    return false;
  }
}

/// \brief This function first splits the given condition e into real conditions and
///        non real conditions.
/// \details This function first uses split_condition_aux to split the condition e. Then
//          it merges equal real conditions by merging the non-real conditions. No further
//          calculations take place with the non-real conditions, but if the non-real conditions
//          lead to unnecessary copying, this may lead to a huge overhead in removing the
//          real conditions.
inline void split_condition(
  const data_expression& e,
  std::vector < data_expression_list >& real_conditions,
  std::vector < data_expression >& non_real_conditions)
{
  std::vector < data_expression_list > aux_real_conditions;
  std::vector < data_expression_list > aux_non_real_conditions;

  split_condition_aux(e,aux_real_conditions, aux_non_real_conditions);
  assert(aux_non_real_conditions.size()==aux_real_conditions.size() && aux_non_real_conditions.size()>0);

  // For every list of real expressions, gather the corresponding non real expressions
  std::map< data_expression_list, data_expression > non_real_expression_map;
  for(std::vector < data_expression_list >::const_iterator i=aux_real_conditions.begin(), j=aux_non_real_conditions.begin();
              i!=aux_real_conditions.end(); ++i, ++j)
  {
    // Find the entry for *i, inserting false if it does not exist yet
    std::map< data_expression_list, data_expression >::iterator insert_result =
        non_real_expression_map.insert(std::make_pair(*i, sort_bool::false_())).first;
    insert_result->second = lazy::or_(insert_result->second, lazy::join_and(j->begin(), j->end()));
  }
  // Convert the map to a pair of vectors
  for(const std::pair<const data_expression_list, data_expression >& expr_pair: non_real_expression_map)
  {
    real_conditions.push_back(expr_pair.first);
    non_real_conditions.push_back(expr_pair.second);
  }
  assert(non_real_conditions.size()==real_conditions.size() && non_real_conditions.size()>0);
}

} // end namespace detail

#endif // MCRL2_LPSREALELM_DETAIL_LINEAR_INEQUALITY_UTILITIES_H
