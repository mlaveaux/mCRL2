// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/replace_capture_avoiding_with_an_identifier_generator.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_REPLACE_CAPTURE_AVOIDING_WITH_AN_IDENTIFIER_GENERATOR_H
#define MCRL2_DATA_REPLACE_CAPTURE_AVOIDING_WITH_AN_IDENTIFIER_GENERATOR_H

#include "mcrl2/data/add_binding.h"
#include "mcrl2/data/builder.h"

#include <ranges>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm.cxx"
#else
  import atermpp;
#endif

namespace mcrl2::data
{

namespace detail {

template <typename Substitution, typename IdentifierGenerator>
class substitution_updater_with_an_identifier_generator
{
  protected:
    Substitution& m_sigma;
    IdentifierGenerator& m_id_generator;
    std::vector<data::assignment> m_undo; // why not a stack datatype?

  public:
    substitution_updater_with_an_identifier_generator(Substitution& sigma, IdentifierGenerator& id_generator)
      : m_sigma(sigma), m_id_generator(id_generator)
    {}

    Substitution& substitution()
    {
      return m_sigma;
    }

    data::variable bind(const data::variable& v)
    {
      m_undo.push_back(data::assignment(v, m_sigma(v))); // save the old assignment of m_sigma.
      if (m_sigma.variable_occurs_in_a_rhs(v))                    // v notin FV(m_sigma).
      {
        m_sigma[v] = v;                                  // Clear sigma[v].
        return v;
      }
      else
      {
        data::variable w(m_id_generator(v.name()), v.sort());
        m_sigma[v] = w;
        return w;
      }
    }

    data::variable push(const data::variable& v)
    {
      return bind(v);
    }

    void pop(const data::variable& )
    {
      const data::assignment& a = m_undo.back();
      m_sigma[a.lhs()] = a.rhs();
      m_undo.pop_back();
    }

    template <typename VariableContainer>
    VariableContainer push(const VariableContainer& container)
    {
      std::vector<data::variable> result;
      for (const variable& v: container)
      {
        result.push_back(bind(v));
      }
      return VariableContainer(result.begin(), result.end());
    }

    template <typename VariableContainer>
    void pop(const VariableContainer& container)
    {
      for (const variable& v: container)
      {
        pop(v);
      }
    }
};

template <typename Substitution, typename IdentifierGenerator>
data::variable update_substitution(Substitution& sigma, const data::variable& v, const std::multiset<data::variable>& V, IdentifierGenerator& id_generator)
{
  using utilities::detail::contains;
  if (!contains(V, v) && sigma(v) == v)
  {
    return v;
  }
  else
  {
    id_generator.add_identifier(v.name());
    data::variable w(id_generator(v.name()), v.sort());

    while (sigma(w) != w || contains(V, w))
    {
      w = data::variable(id_generator(v.name()), v.sort());
    }
    sigma[v] = w;
    return w;
  }
}

template <typename Substitution, typename IdentifierGenerator, typename VariableContainer>
VariableContainer update_substitution(Substitution& sigma, const VariableContainer& v, const std::multiset<data::variable>& V, IdentifierGenerator& id_generator)
{
  std::vector<data::variable> result;
  for (typename VariableContainer::const_iterator i = v.begin(); i != v.end(); ++i)
  {
    result.push_back(update_substitution(sigma, *i, V, id_generator));
  }
  return VariableContainer(result.begin(), result.end());
}

// Helper code for replace_capture_avoiding_variables_with_an_identifier_generator.
template<template<class> class Builder,
  template<template<class> class, class, class, class> class Binder, class Substitution, class IdentifierGenerator>
struct replace_capture_avoiding_variables_builder_with_an_identifier_generator
  : public Binder<Builder, replace_capture_avoiding_variables_builder_with_an_identifier_generator<Builder, Binder, Substitution, IdentifierGenerator>, Substitution, IdentifierGenerator>
{
  using super = Binder<Builder,
      replace_capture_avoiding_variables_builder_with_an_identifier_generator<Builder,
          Binder,
          Substitution,
          IdentifierGenerator>,
      Substitution,
      IdentifierGenerator>;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  replace_capture_avoiding_variables_builder_with_an_identifier_generator(Substitution& sigma,
                                                                          IdentifierGenerator& id_generator)
    : super(sigma, id_generator)
  {
  }
};

template<template<class> class Builder,
  template<template<class> class, class, class, class> class Binder, class Substitution, class IdentifierGenerator>
replace_capture_avoiding_variables_builder_with_an_identifier_generator<Builder, Binder, Substitution, IdentifierGenerator>
apply_replace_capture_avoiding_variables_builder_with_an_identifier_generator(Substitution& sigma,
                                                                              IdentifierGenerator& id_generator)
{
  return replace_capture_avoiding_variables_builder_with_an_identifier_generator<Builder, Binder, Substitution, IdentifierGenerator>(
    sigma, id_generator);
}

// In the class below, the IdentifierGenerator is expected to generate fresh identifiers, not
// occurring anywhere using the operator (). The enumerator_identifier_generator can do this.
// The substitution has as property that it provides a method "variables_in_rhs" yielding the
// variables occurring in the right hand side of assignments. The mutable_indexed_substitution has this
// method.
template<template<class> class Builder, class Derived, class Substitution, class IdentifierGenerator>
struct add_capture_avoiding_replacement_with_an_identifier_generator : public Builder<Derived>
{
  using super = Builder<Derived>;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

protected:
  substitution_updater_with_an_identifier_generator<Substitution, IdentifierGenerator> update_sigma;

public:
  add_capture_avoiding_replacement_with_an_identifier_generator(Substitution& sigma_,
      IdentifierGenerator& id_generator_)
      : update_sigma(sigma_, id_generator_)
  {}

  template <class T>
  void apply(T& result, const variable& v)
  {
    result = update_sigma.substitution()(v);
    }

    template <class T>
    void apply(T& result, const data::where_clause& x)
    {      
      const auto& assignments = x.declarations() | std::views::transform([](const data::assignment_expression& t) { return atermpp::down_cast<data::assignment>(t); });
      std::vector<data::variable> tmp;
      for (const data::assignment& a: assignments)
      {
        tmp.push_back(a.lhs());
      }
      std::vector<data::variable> v = update_sigma.push(tmp);

      // The updated substitution should be applied to the body.
      data::data_expression new_body;
      apply(new_body, x.body());
      update_sigma.pop(v);

      // The original substitution should be applied to the right hand sides of the assignments.
      std::vector<data::assignment> a;
      auto j = v.begin();
      for (auto i = assignments.begin(); i != assignments.end(); ++i, ++j)
      {
        data::data_expression rhs;
        apply(rhs, (*i).rhs());
        a.emplace_back(*j, rhs);
      }
      data::make_where_clause(result, new_body, assignment_list(a.begin(), a.end()));
    }

    template <class T>
    void apply(T& result, const data::forall& x)
    {
      const data::variable_list v = update_sigma.push(x.variables());
      data::data_expression body;
      apply(body, x.body());
      data::make_forall(result, v, body);
      update_sigma.pop(v);
    }

    template <class T>
    void apply(T& result, const data::exists& x)
    {
      const data::variable_list v = update_sigma.push(x.variables());
      data::data_expression body;
      apply(body, x.body());
      data::make_exists(result, v, body);
      update_sigma.pop(v);
    }

    template <class T>
    void apply(T& result, const data::lambda& x)
    {
      const data::variable_list v = update_sigma.push(x.variables());
      data::data_expression body;
      apply(body, x.body());
      data::make_lambda(result, v, body);
      update_sigma.pop(v);
    }

    template <class T>
    void apply(T& /* result */, data_equation& /* x */)
    {
      throw mcrl2::runtime_error("not implemented yet");
    }
};

} // namespace detail

//--- start generated data replace_capture_avoiding_with_identifier_generator code ---//
/// \\brief Applies sigma as a capture avoiding substitution to x using an identifier generator.
/// \\details This substitution function is much faster than replace_variables_capture_avoiding, but
///          it requires an identifier generator that generates strings for fresh variables. These
///          strings must be unique in the sense that they have not been used for other variables.
/// \\param x The object to which the subsitution is applied.
/// \\param sigma A mutable substitution of which it can efficiently be checked whether a variable occurs in its
///              right hand side. The class maintain_variables_in_rhs is useful for this purpose.
/// \\param id_generator A generator that generates unique strings, not yet used as variable names.

template <typename T, typename Substitution, typename IdentifierGenerator>
void replace_variables_capture_avoiding_with_an_identifier_generator(T& x,
                       Substitution& sigma,
                       IdentifierGenerator& id_generator
                      )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::apply_replace_capture_avoiding_variables_builder_with_an_identifier_generator<data::data_expression_builder, data::detail::add_capture_avoiding_replacement_with_an_identifier_generator>(sigma, id_generator).update(x);
}

/// \\brief Applies sigma as a capture avoiding substitution to x using an identifier generator..
/// \\details This substitution function is much faster than replace_variables_capture_avoiding, but
///          it requires an identifier generator that generates strings for fresh variables. These
///          strings must be unique in the sense that they have not been used for other variables.
/// \\param x The object to which the substiution is applied.
/// \\param sigma A mutable substitution of which it can efficiently be checked whether a variable occurs in its
///              right hand side. The class maintain_variables_in_rhs is useful for this purpose.
/// \\param id_generator A generator that generates unique strings, not yet used as variable names.
/// \\return The result is the term x to which sigma has been applied.
template <typename T, typename Substitution, typename IdentifierGenerator>
T replace_variables_capture_avoiding_with_an_identifier_generator(const T& x,
                    Substitution& sigma,
                    IdentifierGenerator& id_generator
                   )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::apply_replace_capture_avoiding_variables_builder_with_an_identifier_generator<data::data_expression_builder, data::detail::add_capture_avoiding_replacement_with_an_identifier_generator>(sigma, id_generator).apply(result, x);
  return result;
}
//--- end generated data replace_capture_avoiding_with_identifier_generator code ---//

} // namespace mcrl2::data

#endif // MCRL2_DATA_REPLACE_CAPTURE_AVOIDING_WITH_AN_IDENTIFIER_GENERATOR_H
