// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/find.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_FIND_H
#define MCRL2_MODAL_FORMULA_FIND_H

#include "mcrl2/process/find.h"
#include "mcrl2/modal_formula/add_binding.h"
#include "mcrl2/modal_formula/traverser.h"

namespace mcrl2
{

namespace action_formulas
{

//--- start generated action_formulas find code ---//
/// \brief Writes all variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
template <typename T, typename OutputIterator>
void find_all_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_all_variables_traverser<action_formulas::variable_traverser>(o).apply(x);
}

/// \brief Returns all variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All variables that occur in the object x
template <typename T>
std::set<data::variable> find_all_variables(const T& x)
{
  std::set<data::variable> result;
  action_formulas::find_all_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are added.
template <typename T, typename OutputIterator>
void find_free_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_free_variables_traverser<action_formulas::data_expression_traverser, action_formulas::add_data_variable_traverser_binding>(o).apply(x);
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
/// \param[in] bound a container of variables
template <typename T, typename OutputIterator, typename VariableContainer>
void find_free_variables_with_bound(const T& x, OutputIterator o, const VariableContainer& bound)
{
  data::detail::make_find_free_variables_traverser<action_formulas::data_expression_traverser, action_formulas::add_data_variable_traverser_binding>(o, bound).apply(x);
}

/// \brief Returns all free variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All free variables that occur in the object x
template <typename T>
std::set<data::variable> find_free_variables(const T& x)
{
  std::set<data::variable> result;
  action_formulas::find_free_variables(x, std::inserter(result, result.end()));
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
  action_formulas::find_free_variables_with_bound(x, std::inserter(result, result.end()), bound);
  return result;
}

/// \brief Writes all identifiers that occur in an object to an output iterator.
/// \param[in] x an object containing identifiers
/// \param[in,out] o an output iterator to which all identifiers occurring in x are written.
template <typename T, typename OutputIterator>
void find_identifiers(const T& x, OutputIterator o)
{
  data::detail::make_find_identifiers_traverser<action_formulas::identifier_string_traverser>(o).apply(x);
}

/// \brief Returns all identifiers that occur in an object.
/// \param[in] x an object containing identifiers
/// \return All identifiers that occur in the object x
template <typename T>
std::set<core::identifier_string> find_identifiers(const T& x)
{
  std::set<core::identifier_string> result;
  action_formulas::find_identifiers(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all sort expressions that occur in an object to an output iterator.
/// \param[in] x an object containing sort expressions
/// \param[in,out] o an output iterator to which all sort expressions occurring in x are written.
template <typename T, typename OutputIterator>
void find_sort_expressions(const T& x, OutputIterator o)
{
  data::detail::make_find_sort_expressions_traverser<action_formulas::sort_expression_traverser>(o).apply(x);
}

/// \brief Returns all sort expressions that occur in an object.
/// \param[in] x an object containing sort expressions
/// \return All sort expressions that occur in the object x
template <typename T>
std::set<data::sort_expression> find_sort_expressions(const T& x)
{
  std::set<data::sort_expression> result;
  action_formulas::find_sort_expressions(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all function symbols that occur in an object to an output iterator.
/// \param[in] x an object containing function symbols
/// \param[in,out] o an output iterator to which all function symbols occurring in x are written.
template <typename T, typename OutputIterator>
void find_function_symbols(const T& x, OutputIterator o)
{
  data::detail::make_find_function_symbols_traverser<action_formulas::data_expression_traverser>(o).apply(x);
}

/// \brief Returns all function symbols that occur in an object.
/// \param[in] x an object containing function symbols
/// \return All function symbols that occur in the object x
template <typename T>
std::set<data::function_symbol> find_function_symbols(const T& x)
{
  std::set<data::function_symbol> result;
  action_formulas::find_function_symbols(x, std::inserter(result, result.end()));
  return result;
}
//--- end generated action_formulas find code ---//

} // namespace action_formulas

namespace regular_formulas
{

//--- start generated regular_formulas find code ---//
/// \brief Writes all variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
template <typename T, typename OutputIterator>
void find_all_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_all_variables_traverser<regular_formulas::variable_traverser>(o).apply(x);
}

/// \brief Returns all variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All variables that occur in the object x
template <typename T>
std::set<data::variable> find_all_variables(const T& x)
{
  std::set<data::variable> result;
  regular_formulas::find_all_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are added.
template <typename T, typename OutputIterator>
void find_free_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_free_variables_traverser<regular_formulas::data_expression_traverser, regular_formulas::add_data_variable_traverser_binding>(o).apply(x);
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
/// \param[in] bound a container of variables
template <typename T, typename OutputIterator, typename VariableContainer>
void find_free_variables_with_bound(const T& x, OutputIterator o, const VariableContainer& bound)
{
  data::detail::make_find_free_variables_traverser<regular_formulas::data_expression_traverser, regular_formulas::add_data_variable_traverser_binding>(o, bound).apply(x);
}

/// \brief Returns all free variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All free variables that occur in the object x
template <typename T>
std::set<data::variable> find_free_variables(const T& x)
{
  std::set<data::variable> result;
  regular_formulas::find_free_variables(x, std::inserter(result, result.end()));
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
  regular_formulas::find_free_variables_with_bound(x, std::inserter(result, result.end()), bound);
  return result;
}

/// \brief Writes all identifiers that occur in an object to an output iterator.
/// \param[in] x an object containing identifiers
/// \param[in,out] o an output iterator to which all identifiers occurring in x are written.
template <typename T, typename OutputIterator>
void find_identifiers(const T& x, OutputIterator o)
{
  data::detail::make_find_identifiers_traverser<regular_formulas::identifier_string_traverser>(o).apply(x);
}

/// \brief Returns all identifiers that occur in an object.
/// \param[in] x an object containing identifiers
/// \return All identifiers that occur in the object x
template <typename T>
std::set<core::identifier_string> find_identifiers(const T& x)
{
  std::set<core::identifier_string> result;
  regular_formulas::find_identifiers(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all sort expressions that occur in an object to an output iterator.
/// \param[in] x an object containing sort expressions
/// \param[in,out] o an output iterator to which all sort expressions occurring in x are written.
template <typename T, typename OutputIterator>
void find_sort_expressions(const T& x, OutputIterator o)
{
  data::detail::make_find_sort_expressions_traverser<regular_formulas::sort_expression_traverser>(o).apply(x);
}

/// \brief Returns all sort expressions that occur in an object.
/// \param[in] x an object containing sort expressions
/// \return All sort expressions that occur in the object x
template <typename T>
std::set<data::sort_expression> find_sort_expressions(const T& x)
{
  std::set<data::sort_expression> result;
  regular_formulas::find_sort_expressions(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all function symbols that occur in an object to an output iterator.
/// \param[in] x an object containing function symbols
/// \param[in,out] o an output iterator to which all function symbols occurring in x are written.
template <typename T, typename OutputIterator>
void find_function_symbols(const T& x, OutputIterator o)
{
  data::detail::make_find_function_symbols_traverser<regular_formulas::data_expression_traverser>(o).apply(x);
}

/// \brief Returns all function symbols that occur in an object.
/// \param[in] x an object containing function symbols
/// \return All function symbols that occur in the object x
template <typename T>
std::set<data::function_symbol> find_function_symbols(const T& x)
{
  std::set<data::function_symbol> result;
  regular_formulas::find_function_symbols(x, std::inserter(result, result.end()));
  return result;
}
//--- end generated regular_formulas find code ---//

} // namespace regular_formulas

namespace state_formulas
{

//--- start generated state_formulas find code ---//
/// \brief Writes all variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
template <typename T, typename OutputIterator>
void find_all_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_all_variables_traverser<state_formulas::variable_traverser>(o).apply(x);
}

/// \brief Returns all variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All variables that occur in the object x
template <typename T>
std::set<data::variable> find_all_variables(const T& x)
{
  std::set<data::variable> result;
  state_formulas::find_all_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are added.
template <typename T, typename OutputIterator>
void find_free_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_free_variables_traverser<state_formulas::data_expression_traverser, state_formulas::add_data_variable_traverser_binding>(o).apply(x);
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
/// \param[in] bound a container of variables
template <typename T, typename OutputIterator, typename VariableContainer>
void find_free_variables_with_bound(const T& x, OutputIterator o, const VariableContainer& bound)
{
  data::detail::make_find_free_variables_traverser<state_formulas::data_expression_traverser, state_formulas::add_data_variable_traverser_binding>(o, bound).apply(x);
}

/// \brief Returns all free variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All free variables that occur in the object x
template <typename T>
std::set<data::variable> find_free_variables(const T& x)
{
  std::set<data::variable> result;
  state_formulas::find_free_variables(x, std::inserter(result, result.end()));
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
  state_formulas::find_free_variables_with_bound(x, std::inserter(result, result.end()), bound);
  return result;
}

/// \brief Writes all identifiers that occur in an object to an output iterator.
/// \param[in] x an object containing identifiers
/// \param[in,out] o an output iterator to which all identifiers occurring in x are written.
template <typename T, typename OutputIterator>
void find_identifiers(const T& x, OutputIterator o)
{
  data::detail::make_find_identifiers_traverser<state_formulas::identifier_string_traverser>(o).apply(x);
}

/// \brief Returns all identifiers that occur in an object.
/// \param[in] x an object containing identifiers
/// \return All identifiers that occur in the object x
template <typename T>
std::set<core::identifier_string> find_identifiers(const T& x)
{
  std::set<core::identifier_string> result;
  state_formulas::find_identifiers(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all sort expressions that occur in an object to an output iterator.
/// \param[in] x an object containing sort expressions
/// \param[in,out] o an output iterator to which all sort expressions occurring in x are written.
template <typename T, typename OutputIterator>
void find_sort_expressions(const T& x, OutputIterator o)
{
  data::detail::make_find_sort_expressions_traverser<state_formulas::sort_expression_traverser>(o).apply(x);
}

/// \brief Returns all sort expressions that occur in an object.
/// \param[in] x an object containing sort expressions
/// \return All sort expressions that occur in the object x
template <typename T>
std::set<data::sort_expression> find_sort_expressions(const T& x)
{
  std::set<data::sort_expression> result;
  state_formulas::find_sort_expressions(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all function symbols that occur in an object to an output iterator.
/// \param[in] x an object containing function symbols
/// \param[in,out] o an output iterator to which all function symbols occurring in x are written.
template <typename T, typename OutputIterator>
void find_function_symbols(const T& x, OutputIterator o)
{
  data::detail::make_find_function_symbols_traverser<state_formulas::data_expression_traverser>(o).apply(x);
}

/// \brief Returns all function symbols that occur in an object.
/// \param[in] x an object containing function symbols
/// \return All function symbols that occur in the object x
template <typename T>
std::set<data::function_symbol> find_function_symbols(const T& x)
{
  std::set<data::function_symbol> result;
  state_formulas::find_function_symbols(x, std::inserter(result, result.end()));
  return result;
}
//--- end generated state_formulas find code ---//

namespace detail {

// collects state variable names in a set
struct state_variable_name_traverser: public state_formulas::state_formula_traverser<state_variable_name_traverser>
{
  using super = state_formulas::state_formula_traverser<state_variable_name_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  std::set<core::identifier_string> names;

  void apply(const state_formulas::variable& x)
  {
    names.insert(x.name());
  }
};

template <template <class> class Traverser, class OutputIterator>
struct find_state_variables_traverser: public Traverser<find_state_variables_traverser<Traverser, OutputIterator> >
{
  using super = Traverser<find_state_variables_traverser<Traverser, OutputIterator>>;
  using super::enter;
  using super::leave;
  using super::apply;

  OutputIterator out;

  find_state_variables_traverser(OutputIterator out_)
    : out(out_)
  {}

  void apply(const variable& v)
  {
    *out = v;
  }
};

template <template <class> class Traverser, class OutputIterator>
find_state_variables_traverser<Traverser, OutputIterator>
make_find_state_variables_traverser(OutputIterator out)
{
  return find_state_variables_traverser<Traverser, OutputIterator>(out);
}

template <template <class> class Traverser, template <template <class> class, class> class Binder, class OutputIterator>
struct find_free_state_variables_traverser: public Binder<Traverser, find_free_state_variables_traverser<Traverser, Binder, OutputIterator> >
{
  using super = Binder<Traverser, find_free_state_variables_traverser<Traverser, Binder, OutputIterator>>;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::is_bound;
  using super::bound_variables;
  using super::increase_bind_count;

  OutputIterator out;

  find_free_state_variables_traverser(OutputIterator out_)
    : out(out_)
  {}

/*
  template <typename VariableContainer>
  find_free_state_variables_traverser(OutputIterator out_, const VariableContainer& v)
    : out(out_)
  {
    increase_bind_count(v);
  }
*/

  void apply(const variable& v)
  {
    if (!is_bound(v.name()))
    {
      *out = v;
    }
  }
};

template <template <class> class Traverser, template <template <class> class, class> class Binder, class OutputIterator>
find_free_state_variables_traverser<Traverser, Binder, OutputIterator>
make_find_free_state_variables_traverser(OutputIterator out)
{
  return find_free_state_variables_traverser<Traverser, Binder, OutputIterator>(out);
}

template <template <class> class Traverser, template <template <class> class, class> class Binder, class OutputIterator, class VariableContainer>
find_free_state_variables_traverser<Traverser, Binder, OutputIterator>
make_find_free_state_variables_traverser(OutputIterator out, const VariableContainer& v)
{
  return find_free_state_variables_traverser<Traverser, Binder, OutputIterator>(out, v);
}

} // namespace detail

/// \brief Returns the names of the state variables that occur in x.
/// \param[in] x A state formula
inline
std::set<core::identifier_string> find_state_variable_names(const state_formula& x)
{
  detail::state_variable_name_traverser f;
  f.apply(x);
  return f.names;
}

/// \brief Writes all state variables that occur in an object to an output iterator.
/// \param[in] x an object containing state variables
/// \param[in,out] o an output iterator to which all state variables occurring in x are written.
template <typename T, typename OutputIterator>
void find_state_variables(const T& x, OutputIterator o)
{
  state_formulas::detail::make_find_state_variables_traverser<state_formulas::state_variable_traverser>(o).apply(x);
}

/// \brief Returns all state variables that occur in an object
/// \param[in] x an object containing variables
/// \return All state variables that occur in the object x
template <typename T>
std::set<state_formulas::variable> find_state_variables(const T& x)
{
  std::set<state_formulas::variable> result;
  state_formulas::find_state_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all free state variables that occur in an object to an output iterator.
/// \param[in] x an object containing state variables
/// \param[in,out] o an output iterator to which all state variables occurring in x are added.
template <typename T, typename OutputIterator>
void find_free_state_variables(const T& x, OutputIterator o)
{
  state_formulas::detail::make_find_free_state_variables_traverser<state_formulas::state_variable_traverser, state_formulas::add_state_variable_binding>(o).apply(x);
}

/// \brief Returns all free state variables that occur in an object
/// \param[in] x an object containing variables
/// \return All state variables that occur in the object x
template <typename T>
std::set<state_formulas::variable> find_free_state_variables(const T& x)
{
  std::set<state_formulas::variable> result;
  state_formulas::find_free_state_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all action labels that occur in an object to an output iterator.
/// \param[in] x an object containing action labels
/// \param[in,out] o an output iterator to which all action labels occurring in x are written.
template <typename T, typename OutputIterator>
void find_action_labels(const T& x, OutputIterator o)
{
  process::detail::make_find_action_labels_traverser<state_formulas::action_label_traverser>(o).apply(x);
}

/// \brief Returns all action labels that occur in an object
/// \param[in] x an object containing action labels
/// \return All action labels that occur in the object x
template <typename T>
std::set<process::action_label> find_action_labels(const T& x)
{
  std::set<process::action_label> result;
  state_formulas::find_action_labels(x, std::inserter(result, result.end()));
  return result;
}

} // namespace state_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULA_FIND_H
