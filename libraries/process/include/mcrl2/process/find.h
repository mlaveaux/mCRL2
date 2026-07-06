// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/find.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_FIND_H
#define MCRL2_PROCESS_FIND_H

#include "mcrl2/data/find.h"
#include "mcrl2/process/add_binding.h"
#include "mcrl2/process/traverser.h"

namespace mcrl2::process
{

namespace detail
{

/// \cond INTERNAL_DOCS
template <template <class> class Traverser, class OutputIterator>
struct find_action_labels_traverser: public Traverser<find_action_labels_traverser<Traverser, OutputIterator> >
{
  using super = Traverser<find_action_labels_traverser<Traverser, OutputIterator>>;
  using super::enter;
  using super::leave;
  using super::apply;

  OutputIterator out;

  find_action_labels_traverser(OutputIterator out_)
    : out(out_)
  {}

  void apply(const data::data_expression& /* x */)
  {
    // Intentionally empty;
  }

  void apply(const data::untyped_data_parameter& /* x */)
  {
    // Intentionally empty;
  }

  void apply(const process::untyped_multi_action& /* x */)
  {
    // Intentionally empty;
  }

  void apply(const process::action_label& x)
  {
    *out = x;
  }

#if BOOST_MSVC
#include "mcrl2/core/detail/traverser_msvc.inc.h"
#endif
};

template <template <class> class Traverser, class OutputIterator>
find_action_labels_traverser<Traverser, OutputIterator> make_find_action_labels_traverser(OutputIterator out)
{
  return find_action_labels_traverser<Traverser, OutputIterator>(out);
}

struct find_action_names_traverser: public process::action_label_traverser<find_action_names_traverser>
{
  using super = process::action_label_traverser<find_action_names_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  std::set<core::identifier_string> result;

  void apply(const process::action_label& x)
  {
    result.insert(x.name());
  }

  void apply(const process::block& x)
  {
    super::apply(x);
    const core::identifier_string_list& B = x.block_set();
    result.insert(B.begin(), B.end());
  }

  void apply(const process::hide& x)
  {
    super::apply(x);
    const core::identifier_string_list& I = x.hide_set();
    result.insert(I.begin(), I.end());
  }

  void apply(const process::rename& x)
  {
    super::apply(x);
    for (const rename_expression& r: x.rename_set())
    {
      result.insert(r.source());
      result.insert(r.target());
    }
  }

  void apply(const process::comm& x)
  {
    super::apply(x);
    for (const communication_expression& c: x.comm_set())
    {
      core::identifier_string_list names = c.action_name().names();
      result.insert(names.begin(), names.end());
      result.insert(c.name());
    }
  }

  void apply(const process::allow& x)
  {
    super::apply(x);
    for (const action_name_multiset& i: x.allow_set())
    {
      const core::identifier_string_list& names = i.names();
      result.insert(names.begin(), names.end());
    }
  }
};
/// \endcond

} // namespace detail

//--- start generated process find code ---//
/// \brief Writes all variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
template <typename T, typename OutputIterator>
void find_all_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_all_variables_traverser<process::variable_traverser>(o).apply(x);
}

/// \brief Returns all variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All variables that occur in the object x
template <typename T>
std::set<data::variable> find_all_variables(const T& x)
{
  std::set<data::variable> result;
  process::find_all_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are added.
template <typename T, typename OutputIterator>
void find_free_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_free_variables_traverser<process::data_expression_traverser, process::add_data_variable_traverser_binding>(o).apply(x);
}

/// \brief Writes all free variables that occur in an object to an output iterator.
/// \param[in] x an object containing variables
/// \param[in,out] o an output iterator to which all variables occurring in x are written.
/// \param[in] bound a container of variables
template <typename T, typename OutputIterator, typename VariableContainer>
void find_free_variables_with_bound(const T& x, OutputIterator o, const VariableContainer& bound)
{
  data::detail::make_find_free_variables_traverser<process::data_expression_traverser, process::add_data_variable_traverser_binding>(o, bound).apply(x);
}

/// \brief Returns all free variables that occur in an object.
/// \param[in] x an object containing variables
/// \return All free variables that occur in the object x
template <typename T>
std::set<data::variable> find_free_variables(const T& x)
{
  std::set<data::variable> result;
  process::find_free_variables(x, std::inserter(result, result.end()));
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
  process::find_free_variables_with_bound(x, std::inserter(result, result.end()), bound);
  return result;
}

/// \brief Writes all identifiers that occur in an object to an output iterator.
/// \param[in] x an object containing identifiers
/// \param[in,out] o an output iterator to which all identifiers occurring in x are written.
template <typename T, typename OutputIterator>
void find_identifiers(const T& x, OutputIterator o)
{
  data::detail::make_find_identifiers_traverser<process::identifier_string_traverser>(o).apply(x);
}

/// \brief Returns all identifiers that occur in an object.
/// \param[in] x an object containing identifiers
/// \return All identifiers that occur in the object x
template <typename T>
std::set<core::identifier_string> find_identifiers(const T& x)
{
  std::set<core::identifier_string> result;
  process::find_identifiers(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all sort expressions that occur in an object to an output iterator.
/// \param[in] x an object containing sort expressions
/// \param[in,out] o an output iterator to which all sort expressions occurring in x are written.
template <typename T, typename OutputIterator>
void find_sort_expressions(const T& x, OutputIterator o)
{
  data::detail::make_find_sort_expressions_traverser<process::sort_expression_traverser>(o).apply(x);
}

/// \brief Returns all sort expressions that occur in an object.
/// \param[in] x an object containing sort expressions
/// \return All sort expressions that occur in the object x
template <typename T>
std::set<data::sort_expression> find_sort_expressions(const T& x)
{
  std::set<data::sort_expression> result;
  process::find_sort_expressions(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Writes all function symbols that occur in an object to an output iterator.
/// \param[in] x an object containing function symbols
/// \param[in,out] o an output iterator to which all function symbols occurring in x are written.
template <typename T, typename OutputIterator>
void find_function_symbols(const T& x, OutputIterator o)
{
  data::detail::make_find_function_symbols_traverser<process::data_expression_traverser>(o).apply(x);
}

/// \brief Returns all function symbols that occur in an object.
/// \param[in] x an object containing function symbols
/// \return All function symbols that occur in the object x
template <typename T>
std::set<data::function_symbol> find_function_symbols(const T& x)
{
  std::set<data::function_symbol> result;
  process::find_function_symbols(x, std::inserter(result, result.end()));
  return result;
}
//--- end generated process find code ---//

/// \brief Writes all action labels that occur in an object to an output iterator.
/// \param[in] x an object containing action labels
/// \param[in,out] o an output iterator to which all action labels occurring in x are written.
template <typename T, typename OutputIterator>
void find_action_labels(const T& x, OutputIterator o)
{
  process::detail::make_find_action_labels_traverser<process::action_label_traverser>(o)(x);
}

/// \brief Returns all action labels that occur in an object
/// \param[in] x an object containing action labels
/// \return All action labels that occur in the object x
template <typename T>
std::set<process::action_label> find_action_labels(const T& x)
{
  std::set<process::action_label> result;
  process::find_action_labels(x, std::inserter(result, result.end()));
  return result;
}

/// \brief Returns all action names that occur in an object
/// \param[in] x an object containing action names
/// \return All action names that occur in the object x
template <typename T>
std::set<core::identifier_string> find_action_names(const T& x)
{
  detail::find_action_names_traverser f;
  f.apply(x);
  return f.result;
}

/// \brief Finds an equation that corresponds to a process identifier
/// \param[in] equations a sequence of process equations
/// \param[in] id The identifier of the equation that is searched for.
/// \return The equation with the given process identifier. Throws an exception if no such equation was found.
inline
const process_equation& find_equation(const std::vector<process_equation>& equations, const process_identifier& id)
{
  for (const process_equation& equation: equations)
  {
    if (equation.identifier() == id)
    {
      return equation;
    }
  }
  throw mcrl2::runtime_error("unknown process identifier " + process::pp(id));
}

} // namespace mcrl2::process

#endif // MCRL2_PROCESS_FIND_H
