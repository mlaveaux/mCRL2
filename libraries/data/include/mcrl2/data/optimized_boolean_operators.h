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

#ifndef MCRL2_DATA_OPTIMIZED_BOOLEAN_OPERATORS_H
#define MCRL2_DATA_OPTIMIZED_BOOLEAN_OPERATORS_H

#include "mcrl2/core/term_traits.h"
#include "mcrl2/core/optimized_boolean_operators.h"
#include "mcrl2/data/expression_traits.h"

namespace mcrl2::data
{

/// \brief Make a negation
/// \param result Variable to which the resulting term is assigned.
/// \param arg A term
inline
void make_optimized_not(data_expression& result, const data_expression& arg)
{
  mcrl2::core::make_optimized_not(result, arg);
}

/// \brief Make a negation
/// \param arg A term
/// \return The application of not to the argument.
inline
data_expression optimized_not(const data_expression& arg)
{
  return mcrl2::core::optimized_not(arg);
}

/// \brief Make a conjunction, and optimize if possible.
/// \param result Contains the optimized and.
/// \param p A term
/// \param q A term
inline
void make_optimized_and(data_expression& result, const data_expression& p, const data_expression& q)
{
  mcrl2::core::make_optimized_and(result, p, q);
}

/// \brief Make a conjunction, and optimize if possible.
/// \param p A term
/// \param q A term
/// \return The application of and to the arguments.
inline
data_expression optimized_and(const data_expression& p, const data_expression& q)
{
  return mcrl2::core::optimized_and(p, q);
}

/// \brief Make a disjunction
/// \param result Will contains the optimized or.
/// \param p A term
/// \param q A term
inline
void make_optimized_or(data_expression& result, const data_expression& p, const data_expression& q)
{
  mcrl2::core::make_optimized_or(result, p, q);
}

/// \brief Make a disjunction
/// \param p A term
/// \param q A term
/// \return The application of or to the arguments.
inline
data_expression optimized_or(const data_expression& p, const data_expression& q)
{
  return mcrl2::core::optimized_or(p, q);
}

/// \brief Make an implication and apply optimisations.
/// \param result Will contains the optimized implication.
/// \param p A term
/// \param q A term
/// \return The application of implication to the arguments.
inline
void make_optimized_imp(data_expression& result, const data_expression& p, const data_expression& q)
{
  mcrl2::core::make_optimized_imp(result, p, q);
}

/// \brief Make an implication and apply optimisations.
/// \param p A term
/// \param q A term
/// \return The application of implication to the arguments.
inline
data_expression optimized_imp(const data_expression& p, const data_expression& q)
{
  return mcrl2::core::optimized_imp(p, q);
}

/// \brief Make a universal quantification and apply optimisations.
/// \param result Will contains the optimized universal implication
/// \param l A sequence of variables
/// \param p A term
/// \param remove_variables If true, unused quantifier variables are removed
inline
void make_optimized_forall(data_expression& result, const variable_list& l, const data_expression& p, bool remove_variables = false)
{
  mcrl2::core::make_optimized_forall(result, l, p, remove_variables);
}

/// \brief Make a universal quantification and apply optimisations.
/// \param l A sequence of variables
/// \param p A term
/// \param remove_variables If true, unused quantifier variables are removed
/// \return The application of universal quantification to the arguments.
inline
data_expression optimized_forall(const variable_list& l, const data_expression& p, bool remove_variables = false)
{
  return mcrl2::core::optimized_forall(l, p, remove_variables);
}

/// \brief Make an existential quantification and apply optimisations.
/// \param result Will contains the optimized existential implication.
/// \param l A sequence of variables
/// \param p A term
/// \param remove_variables If true, unused quantifier variables are removed
/// \return The application of existential quantification to the arguments.
inline
void make_optimized_exists(data_expression& result, const variable_list& l, const data_expression& p, bool remove_variables = false)
{
  mcrl2::core::make_optimized_exists(result, l, p, remove_variables);
}

/// \brief Make an existential quantification and apply optimisations.
/// \param l A sequence of variables
/// \param p A term
/// \param remove_variables If true, unused quantifier variables are removed
/// \return The application of existential quantification to the arguments.
inline
data_expression optimized_exists(const variable_list& l, const data_expression& p, bool remove_variables = false)
{
  return mcrl2::core::optimized_exists(l, p, remove_variables);
}

} // namespace mcrl2::data

#endif // MCRL2_DATA_OPTIMIZED_BOOLEAN_OPERATORS_H
