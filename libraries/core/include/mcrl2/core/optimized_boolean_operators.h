// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/optimized_boolean_operators.h
/// \brief add your file description here.

#ifndef MCRL2_CORE_OPTIMIZED_BOOLEAN_OPERATORS_H
#define MCRL2_CORE_OPTIMIZED_BOOLEAN_OPERATORS_H

#include "mcrl2/core/term_traits.h"

namespace mcrl2::core
{

namespace detail
{
/// \return The value <tt>!arg</tt>
template <typename TermTraits>
inline
void make_optimized_not(typename TermTraits::term_type& result, 
                        const typename TermTraits::term_type& arg, 
                        TermTraits)
{
  using tr = TermTraits;

  if (tr::is_true(arg))
  {
    result=tr::false_();
  }
  else if (tr::is_false(arg))
  {
    result=tr::true_();
  }
  else if (tr::is_not(arg))
  {
    result=tr::not_arg(arg);
  }
  else
  {
    tr::make_not(result, arg);
  }
}

/// \brief Make a conjunction and optimize it if possible.
/// \param left A term
/// \param right A term
/// \return The value <tt>left && right</tt>
template <typename TermTraits>
inline
void make_optimized_and(typename TermTraits::term_type& result, 
                        const typename TermTraits::term_type& left, 
                        const typename TermTraits::term_type& right, 
                        TermTraits)
{
  using tr = TermTraits;

  if (tr::is_true(left))
  {
    result=right;
  }
  else if (tr::is_false(left))
  {
    result=tr::false_();
  }
  else if (tr::is_true(right))
  {
    result=left;
  }
  else if (tr::is_false(right))
  {
    result=tr::false_();
  }
  else if (left == right)
  {
    result=left;
  }
  else
  {
    tr::make_and(result, left, right);
  }
}

/// \brief Make a disjunction
/// \param left A term
/// \param right A term
/// \return The value <tt>left || right</tt>
template <typename TermTraits>
inline
void make_optimized_or(typename TermTraits::term_type& result, 
                       const typename TermTraits::term_type& left, 
                       const typename TermTraits::term_type& right, 
                       TermTraits)
{
  using tr = TermTraits;

  if (tr::is_true(left))
  {
    result=tr::true_();
  }
  else if (tr::is_false(left))
  {
    result=right;
  }
  else if (tr::is_true(right))
  {
    result=tr::true_();
  }
  else if (tr::is_false(right))
  {
    result=left;
  }
  else if (left == right)
  {
    result=left;
  }
  else
  {
    tr::make_or(result, left, right);
  }
}

/// \brief Make an implication
/// \param left A term
/// \param right A term
/// \return The value <tt>left => right</tt>
template <typename TermTraits>
inline
void make_optimized_imp(typename TermTraits::term_type& result, 
                        const typename TermTraits::term_type& left, 
                        const typename TermTraits::term_type& right, 
                        TermTraits t)
{
  using tr = TermTraits;

  if (tr::is_true(left))
  {
    result=right;
  }
  else if (tr::is_false(left))
  {
    result=tr::true_();
  }
  else if (tr::is_true(right))
  {
    result=tr::true_();
  }
  else if (tr::is_false(right))
  {
    make_optimized_not(result,left, t);
  }
  else if (left == right)
  {
    result=tr::true_();
  }
  else
  {
    tr::make_imp(result, left, right);
  }
}

/// \brief Make a universal quantification
/// \param v A sequence of variables
/// \param arg A term
////// \param remove_variables If true, remove bound variables that do not occur in \a arg.
/// \details If the variable sequence v is empty, the result is arg. 
/// \return The universal quantification <tt>forall v.arg</tt>
template <typename TermTraits>
inline
void make_optimized_forall(typename TermTraits::term_type& result, 
                           const typename TermTraits::variable_sequence_type& v, 
                           const typename TermTraits::term_type& arg, 
                           bool remove_variables, 
                           TermTraits)
{
  using tr = TermTraits;

  if (v.empty())
  {
    result = arg;
  }
  else if (tr::is_true(arg))
  {
    result = tr::true_();
  }
  else if (tr::is_false(arg))
  {
    result = tr::false_();
  }
  else
  {
    if (remove_variables)
    {
      typename tr::variable_sequence_type variables = tr::set_intersection(v, tr::find_free_variables(arg));
      if (variables.empty())
      {
        result = arg;
      }
      else
      {
        tr::make_forall(result, variables, arg);
      }
    }
    else
    {
      tr::make_forall(result, v, arg);
    }
  }
}

/// \brief Make an existential quantification
/// \param v A sequence of variables
/// \param arg A term
/// \param remove_variables If true, remove bound variables that do not occur in \a arg.
/// \details If the sequence of variable v is empty, the result is arg. 
/// \return The existential quantification <tt>exists v.arg</tt>
template <typename TermTraits>
inline
void make_optimized_exists(typename TermTraits::term_type& result,
                           const typename TermTraits::variable_sequence_type& v, 
                           const typename TermTraits::term_type& arg, 
                           bool remove_variables, 
                           TermTraits)
{
  using tr = TermTraits;

  if (v.empty())
  {
    result = arg;
  }
  else if (tr::is_true(arg))
  {
    result = tr::true_();
  }
  else if (tr::is_false(arg))
  {
    result = tr::false_();
  }
  else
  {
    if (remove_variables)
    {
       typename tr::variable_sequence_type variables = tr::set_intersection(v, find_free_variables(arg));
      if (variables.empty())
      {
        result = arg;
      }
      else
      {
        result = tr::exists(variables, arg);
      }
    }
    else
    {
      result = tr::exists(v, arg);
    }
  }
}

} // namespace detail

/// \brief Make a negation
/// \param result The placeholder where the negation of the argument is put. This is more efficient than returning as a function result.
/// \param arg A term
template <typename Term>
inline
void make_optimized_not(Term& result, const Term& arg)
{
  detail::make_optimized_not(result, arg, core::term_traits<Term>());
}

/// \brief Make a negation
/// \param arg A term
/// \return The application of not to the argument.
template <typename Term>
inline
Term optimized_not(const Term& arg)
{
  Term result;
  detail::make_optimized_not(result, arg, core::term_traits<Term>());
  return result;
}

/// \brief Make a conjunction, and optimize if possible.
/// \param result Variable to contain the resulting optimized and.
/// \param p A term
/// \param q A term
template <typename Term>
inline
void make_optimized_and(Term& result, const Term& p, const Term& q)
{
  return detail::make_optimized_and(result, p, q, core::term_traits<Term>());
}

/// \brief Make a conjunction, and optimize if possible.
/// \param p A term
/// \param q A term
/// \return The application of and to the arguments.
template <typename Term>
inline
Term optimized_and(const Term& p, const Term& q)
{
  Term result;
  detail::make_optimized_and(result, p, q, core::term_traits<Term>());
  return result;
} 

/// \brief Make a disjunction
/// \param result Variable to contain the resulting optimized or.
/// \param p A term
/// \param q A term
template <typename Term>
inline
void make_optimized_or(Term& result, const Term& p, const Term& q)
{
  detail::make_optimized_or(result, p, q, core::term_traits<Term>());
}

/// \brief Make a disjunction, and optimize if possible.
/// \param p A term
/// \param q A term
/// \return The application of the disjunction to the arguments.
template <typename Term>
inline
Term optimized_or(const Term& p, const Term& q)
{
  Term result;
  detail::make_optimized_or(result, p, q, core::term_traits<Term>());
  return result;
} 

/// \brief Make an implication
/// \param result Variable to contain the resulting optimized implication.
/// \param p A term
/// \param q A term
template <typename Term>
inline
void make_optimized_imp(Term& result, const Term& p, const Term& q)
{
  detail::make_optimized_imp(result, p, q, core::term_traits<Term>());
}

/// \brief Make a conjunction, and optimize if possible.
/// \param p A term
/// \param q A term
/// \return The application of an optimized implication to the arguments.
template <typename Term>
inline
Term optimized_imp(const Term& p, const Term& q)
{
  Term result;
  detail::make_optimized_imp(result, p, q, core::term_traits<Term>());
  return result;
} 

/// \brief Make an optimized universal quantification, applying optimisations when possible.
/// \param result Place where the quantified expressions is placed. If the variable list l is empty, result becomes p.
/// \param l A sequence of variables.
/// \param p A term.
/// \param remove_variables If true, unused quantifier variables are removed.
template <typename Term, typename VariableSequence>
inline
void make_optimized_forall(Term& result, const VariableSequence& l, const Term& p, bool remove_variables = false)
{
  detail::make_optimized_forall(result, l, p, remove_variables, core::term_traits<Term>());
}

/// \brief Make an optimized universal quantification, applying optimisations.
/// \param l A sequence of variables
/// \param p A term
/// \param remove_variables If true, unused quantifier variables are removed
/// \return The application of universal quantification to the arguments. If the l is empty, p is returned. 
template <typename Term, typename VariableSequence>
inline
Term optimized_forall(const VariableSequence& l, const Term& p, bool remove_variables = false)
{
  Term result;
  detail::make_optimized_forall(result, l, p, remove_variables, core::term_traits<Term>());
  return result;
}

/// \brief Make an optimized existential quantification, optimizing the result.
/// \param result Place where the quantified expressions is placed. If l is empty, result becomes p. 
/// \param l A sequence of variables.
/// \param p A term.
/// \param remove_variables If true, unused quantifier variables are removed
template <typename Term, typename VariableSequence>
inline
void make_optimized_exists(Term& result, const VariableSequence& l, const Term& p, bool remove_variables = false)
{
  detail::make_optimized_exists(result, l, p, remove_variables, core::term_traits<Term>());
}

/// \brief Make an optimized existential quantification.
/// \param l A sequence of variables
/// \param p A term
/// \param remove_variables If true, unused quantifier variables are removed
/// \return The application of existential quantification to the arguments. If l is empty, p is returned. 
template <typename Term, typename VariableSequence>
inline
Term optimized_exists(const VariableSequence& l, const Term& p, bool remove_variables = false)
{
  Term result;
  detail::make_optimized_exists(result, l, p, remove_variables, core::term_traits<Term>());
  return result;
}

} // namespace mcrl2::core

#endif // MCRL2_CORE_OPTIMIZED_BOOLEAN_OPERATORS_H
