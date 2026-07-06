// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/find.h
/// \brief Search functions of the data library.

#ifndef MCRL2_LPS_FIND_H
#define MCRL2_LPS_FIND_H

#include "mcrl2/lps/add_binding.h"
#include "mcrl2/lps/traverser.h"
#include "mcrl2/process/find.h"

namespace mcrl2::lps
{

//--- start generated lps find code ---//
/// \brief Writes all variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
template <typename T, typename OutputIterator>
void find_all_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_all_variables_traverser<lps::variable_traverser>(o).apply(x);
}

/// \brief Returns all variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All variables that occur in the object x
template <typename T>
std::set<data::variable> find_all_variables(const T& x)
{
  std::set<data::variable> result;
  lps::find_all_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are added.
template <typename T, typename OutputIterator>
void find_free_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_free_variables_traverser<lps::data_expression_traverser, lps::add_data_variable_traverser_binding>(o).apply(x);
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
/// \param[in] bound a container of variables
template <typename T, typename OutputIterator, typename VariableContainer>
void find_free_variables_with_bound(const T& x, OutputIterator o, const VariableContainer& bound)
{
  data::detail::make_find_free_variables_traverser<lps::data_expression_traverser, lps::add_data_variable_traverser_binding>(o, bound).apply(x);
}

/// \brief Returns all free variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All free variables that occur in the object x
template <typename T>
std::set<data::variable> find_free_variables(const T& x)
{
  std::set<data::variable> result;
  lps::find_free_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Returns all free variables that occur in an object.
/// \param[in] x an object containing variables
/// \param[in] bound a container of variables
/// \return All free variables that occur in the object x
template <typename T, typename VariableContainer>
std::set<data::variable> find_free_variables_with_bound(const T& x, VariableContainer const& bound)
{
  std::set<data::variable> result;
  lps::find_free_variables_with_bound(x, std::inserter(result, result.end()), bound);
  return result;
}

/// \brief Writes all identifiers that occur in an object to an output iterator.
/// \param[in] x an object containing identifiers
/// \param[in,out] o an output iterator to which all identifiers occurring in x are written.
template <typename T, typename OutputIterator>
void find_identifiers(const T& x, OutputIterator o)
{
  data::detail::make_find_identifiers_traverser<lps::identifier_string_traverser>(o).apply(x);
}

/// \brief Returns all identifiers that occur in an object.
/// \param[in] x an object containing identifiers
/// \return All identifiers that occur in the object x
template <typename T>
std::set<core::identifier_string> find_identifiers(const T& x)
{
  std::set<core::identifier_string> result;
  lps::find_identifiers(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all sort expressions that occur in an object to an output iterator.
/// \param[in] x an object containing sort expressions
/// \param[in,out] o an output iterator to which all sort expressions occurring in x are written.
template <typename T, typename OutputIterator>
void find_sort_expressions(const T& x, OutputIterator o)
{
  data::detail::make_find_sort_expressions_traverser<lps::sort_expression_traverser>(o).apply(x);
}

/// \brief Returns all sort expressions that occur in an object.
/// \param[in] x an object containing sort expressions
/// \return All sort expressions that occur in the object x
template <typename T>
std::set<data::sort_expression> find_sort_expressions(const T& x)
{
  std::set<data::sort_expression> result;
  lps::find_sort_expressions(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all function symbols that occur in an object to an output iterator.
/// \param[in] x an object containing function symbols
/// \param[in,out] o an output iterator to which all function symbols occurring in x are written.
template <typename T, typename OutputIterator>
void find_function_symbols(const T& x, OutputIterator o)
{
  data::detail::make_find_function_symbols_traverser<lps::data_expression_traverser>(o).apply(x);
}

/// \brief Returns all function symbols that occur in an object.
/// \param[in] x an object containing function symbols
/// \return All function symbols that occur in the object x
template <typename T>
std::set<data::function_symbol> find_function_symbols(const T& x)
{
  std::set<data::function_symbol> result;
  lps::find_function_symbols(x, std::inserter(result, result.end()));
  return result;
}
//--- end generated lps find code ---//

/// \brief Returns true if the term has a given free variable as subterm.
/// \param[in] x an expression
/// \param[in] v a variable
/// \return True if v occurs free in x.
template <typename T>
bool search_free_variable(const T& x, const data::variable& v)
{
  data::detail::search_free_variable_traverser<lps::data_expression_traverser, lps::add_data_variable_traverser_binding> f(v);
  f.apply(x);
  return f.found;
}

/// \brief Returns all action labels that occur in an object
/// \param[in] x an object containing action labels
/// \param[in,out] o an output iterator to which all action labels occurring in x are written.
/// \return All action labels that occur in the term x
template <typename T, typename OutputIterator>
void find_action_labels(const T& x, OutputIterator o)
{
  process::detail::make_find_action_labels_traverser<lps::action_label_traverser>(o).apply(x);
}

/// \brief Returns all action labels that occur in an object
/// \param[in] x an object containing action labels
/// \return All action labels that occur in the object x
template <typename T>
std::set<process::action_label> find_action_labels(const T& x)
{
  std::set<process::action_label> result;
  lps::find_action_labels(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Returns all variables inside a multi-action.
/// \param x A multi-action from which variables are.
/// \param o The output to which the variables are added. 
template <class OutputIterator>
void find_all_variables(const multi_action& x, OutputIterator o)
{
  for(const process::action& a: x.actions())
  {
    process::find_all_variables(a,o);
  }
}

/// \brief Returns all variables inside a multi-action.
/// \param x A multi-action from which variables are.
std::set<data::variable> find_all_variables(const multi_action& x);

} // namespace mcrl2::lps

#endif // MCRL2_LPS_FIND_H
