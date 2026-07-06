// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/replace.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_REPLACE_H
#define MCRL2_LPS_REPLACE_H

#include "mcrl2/data/substitutions/mutable_map_substitution.h"
#include "mcrl2/lps/replace_capture_avoiding.h"
#include "mcrl2/process/replace.h"

namespace mcrl2::lps
{

namespace detail {

} // namespace detail

//--- start generated lps replace code ---//
template <typename T, typename Substitution>
void replace_sort_expressions(T& x,
                              const Substitution& sigma,
                              bool innermost
                             )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_replace_sort_expressions_builder<lps::sort_expression_builder>(sigma, innermost).update(x);
}

template <typename T, typename Substitution>
T replace_sort_expressions(const T& x,
                           const Substitution& sigma,
                           bool innermost
                          )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::make_replace_sort_expressions_builder<lps::sort_expression_builder>(sigma, innermost).apply(result, x);
  return result;
}

template <typename T, data::IsSubstitution Substitution>
void replace_data_expressions(T& x,
                              const Substitution& sigma,
                              bool innermost
                             )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_replace_data_expressions_builder<lps::data_expression_builder>(sigma, innermost).update(x);
}

template <typename T, data::IsSubstitution Substitution>
T replace_data_expressions(const T& x,
                           const Substitution& sigma,
                           bool innermost
                          )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::make_replace_data_expressions_builder<lps::data_expression_builder>(sigma, innermost).apply(result, x);
  return result;
}


template <typename T, data::IsSubstitution Substitution>
void replace_variables(T& x,
                       const Substitution& sigma
                      )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  core::make_update_apply_builder<lps::data_expression_builder>(sigma).update(x);
}

template <typename T, data::IsSubstitution Substitution>
T replace_variables(const T& x,
                    const Substitution& sigma
                   )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  core::make_update_apply_builder<lps::data_expression_builder>(sigma).apply(result, x);
  return result;
}

/* Replace all variables, including those in binders and the left hand side of assignments */
template <typename T, data::IsSubstitution Substitution>
void replace_all_variables(T& x,
                           const Substitution& sigma
                          )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  core::make_update_apply_builder<lps::sort_expression_builder>(sigma).update(x);
}

template <typename T, data::IsSubstitution Substitution>
T replace_all_variables(const T& x,
                        const Substitution& sigma
                       )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  core::make_update_apply_builder<lps::sort_expression_builder>(sigma).apply(result, x);
  return result;
}

/// \\brief Applies the substitution sigma to x.
/// \\pre { The substitution sigma must have the property that FV(sigma(x)) is included in {x} for all variables x. }
template <typename T, data::IsSubstitution Substitution>
void replace_free_variables(T& x,
                            const Substitution& sigma
                           )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  assert(data::is_simple_substitution(sigma));
  data::detail::make_replace_free_variables_builder<lps::data_expression_builder, lps::add_data_variable_builder_binding>(sigma).update(x);
}

/// \\brief Applies the substitution sigma to x.
/// \\pre { The substitution sigma must have the property that FV(sigma(x)) is included in {x} for all variables x. }
template <typename T, data::IsSubstitution Substitution>
T replace_free_variables(const T& x,
                         const Substitution& sigma
                        )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  assert(data::is_simple_substitution(sigma));
  T result;
  data::detail::make_replace_free_variables_builder<lps::data_expression_builder, lps::add_data_variable_builder_binding>(sigma).apply(result, x);
  return result;
}

/// \\brief Applies the substitution sigma to x, where the elements of bound_variables are treated as bound variables.
/// \\pre { The substitution sigma must have the property that FV(sigma(x)) is included in {x} for all variables x. }
template <typename T, data::IsSubstitution Substitution, typename VariableContainer>
void replace_free_variables(T& x,
                            const Substitution& sigma,
                            const VariableContainer& bound_variables
                           )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  assert(data::is_simple_substitution(sigma));
  data::detail::make_replace_free_variables_builder<lps::data_expression_builder, lps::add_data_variable_builder_binding>(sigma).update(x, bound_variables);
}

/// \\brief Applies the substitution sigma to x, where the elements of bound_variables are treated as bound variables.
/// \\pre { The substitution sigma must have the property that FV(sigma(x)) is included in {x} for all variables x. }
template <typename T, data::IsSubstitution Substitution, typename VariableContainer>
T replace_free_variables(const T& x,
                         const Substitution& sigma,
                         const VariableContainer& bound_variables
                        )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  assert(data::is_simple_substitution(sigma));
  T result;
  data::detail::make_replace_free_variables_builder<lps::data_expression_builder, lps::add_data_variable_builder_binding>(sigma).apply(result, x, bound_variables);
  return result;
}
//--- end generated lps replace code ---//

namespace detail {

/// \cond INTERNAL_DOCS
template <template <class> class Builder, template <template <class> class, class> class Binder, data::IsSubstitution Substitution>
struct replace_process_parameter_builder: public Binder<Builder, replace_process_parameter_builder<Builder, Binder, Substitution> >
{
  using super = Binder<Builder, replace_process_parameter_builder<Builder, Binder, Substitution>>;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;
  using super::is_bound;
  using super::bound_variables;
  using super::increase_bind_count;

  Substitution sigma;
  std::size_t count; // the bind count

  explicit replace_process_parameter_builder(Substitution sigma_)
    : sigma(sigma_), count(1)
  {}

  template <typename VariableContainer>
  replace_process_parameter_builder(Substitution sigma_, const VariableContainer& bound_variables)
    : sigma(sigma_), count(1)
  {
    increase_bind_count(bound_variables);
  }

  void apply(data::variable& result, const data::variable& x)
  {
    if (bound_variables().count(x) == count)
    {
      result = atermpp::down_cast<data::variable>(sigma(x));
      return;
    }
    result = x;
  }

  template <class T>
  void apply(T& result, const data::variable& x)
  {
    if (bound_variables().count(x) == count)
    {
      result = sigma(x);
      return;
    }
    result = static_cast<data::data_expression>(x);
  }

  template <class T>
  void apply(T& result, const data::assignment& x)
  {
    data::variable lhs;
    apply(lhs, x.lhs());
    data::data_expression rhs;
    apply(rhs, x.rhs());
    data::make_assignment(result, lhs, rhs);
  }

  void update(lps::deadlock_summand& x)
  {
    count = 1;
    super::update(x);
  }

  void update(lps::action_summand& x)
  {
    count = 1;
    super::update(x);
    data::assignment_list assignments;
    super::apply(assignments, x.assignments());
    x.assignments() = assignments;
  }

  void update(lps::linear_process& x)
  {
    super::update(x);
    count = 0;
    data::variable_list process_parameters;
    super::apply(process_parameters, x.process_parameters());
    x.process_parameters() = process_parameters;
  }
};

template <template <class> class Builder, template <template <class> class, class> class Binder, data::IsSubstitution Substitution>
replace_process_parameter_builder<Builder, Binder, Substitution>
make_replace_process_parameters_builder(Substitution sigma)
{
  return replace_process_parameter_builder<Builder, Binder, Substitution>(sigma);
}
/// \endcond

} // namespace detail

/// \brief Applies a substitution to the process parameters of the specification spec.
template <data::IsSubstitution Substitution>
void replace_process_parameters(specification& spec, Substitution sigma)
{
  lps::detail::make_replace_process_parameters_builder<lps::data_expression_builder, lps::add_data_variable_builder_binding>(sigma).update(spec);
}

/// \brief Applies a substitution to the process parameters of the specification spec.
inline
void replace_summand_variables(specification& spec, data::mutable_map_substitution<>& sigma)
{
  data::set_identifier_generator id_generator;
  for (const data::variable& v: data::substitution_variables(sigma))
  {
    id_generator.add_identifier(v.name());
  }

  for (action_summand& s: spec.process().action_summands())
  {
    s.summation_variables() = data::replace_variables(s.summation_variables(), sigma);
    s.condition() = data::replace_variables_capture_avoiding(s.condition(), sigma, id_generator);
    s.multi_action() = lps::replace_variables_capture_avoiding(s.multi_action(), sigma, id_generator);
    s.assignments() = data::replace_variables_capture_avoiding(s.assignments(), sigma, id_generator);
  }
  for (deadlock_summand& s: spec.process().deadlock_summands())
  {
    s.summation_variables() = data::replace_variables(s.summation_variables(), sigma);
    s.condition() = data::replace_variables_capture_avoiding(s.condition(), sigma, id_generator);
    lps::replace_variables_capture_avoiding(s.deadlock(), sigma, id_generator);
  }
}

} // namespace mcrl2::lps

#endif // MCRL2_LPS_REPLACE_H
