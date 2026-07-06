// Author(s): Jan Friso Groote. Based on pbes/rewrite.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pres/include/mcrl2/pres/rewrite.h
/// \brief add your file description here.

#ifndef MCRL2_PRES_REWRITE_H
#define MCRL2_PRES_REWRITE_H

#include "mcrl2/data/rewrite.h"
#include "mcrl2/pres/builder.h"



namespace mcrl2::pres_system {

namespace detail
{

/// \cond INTERNAL_DOCS
template <template <class> class Builder, class Rewriter>
struct rewrite_pres_expressions_builder: public Builder<rewrite_pres_expressions_builder<Builder, Rewriter> >
{
  using super = Builder<rewrite_pres_expressions_builder<Builder, Rewriter>>;
  using super::apply;
  using super::update;

  const Rewriter& R;

  rewrite_pres_expressions_builder(const Rewriter& R_)
    : R(R_)
  {}

  template <class T>
  void apply(T& result, const pres_expression& x)
  {
    result = R(x);
  }

  void update(pres_system::pres& x)
  {
    super::update(x);
    // Handle the initial state. It is skipped by the builder because the type is not pres_expression
    pres_expression initial_state;
    apply(initial_state, static_cast<pres_expression>(x.initial_state()));
    x.initial_state() = atermpp::down_cast<propositional_variable_instantiation>(initial_state);
  }
};

template <template <class> class Builder, class Rewriter>
rewrite_pres_expressions_builder<Builder, Rewriter>
make_rewrite_pres_expressions_builder(const Rewriter& R)
{
  return rewrite_pres_expressions_builder<Builder, Rewriter>(R);
}

template <template <class> class Builder, class Rewriter, class Substitution>
struct rewrite_pres_expressions_with_substitution_builder: public Builder<rewrite_pres_expressions_with_substitution_builder<Builder, Rewriter, Substitution> >
{
  using super = Builder<rewrite_pres_expressions_with_substitution_builder<Builder, Rewriter, Substitution>>;
  using super::apply;

  const Rewriter& R;
  Substitution& sigma;

  rewrite_pres_expressions_with_substitution_builder(const Rewriter& R_, Substitution& sigma_)
    : R(R_),
      sigma(sigma_)
  {}

  template <class T>
  void apply(T& result, const pres_expression& x)
  {
    result = R(x, sigma);
  }
};

template <template <class> class Builder, class Rewriter, class Substitution>
rewrite_pres_expressions_with_substitution_builder<Builder, Rewriter, Substitution>
make_rewrite_pres_expressions_with_substitution_builder(const Rewriter& R, Substitution sigma)
{
  return rewrite_pres_expressions_with_substitution_builder<Builder, Rewriter, Substitution>(R, sigma);
}
/// \endcond

} // namespace detail

//--- start generated pres_system rewrite code ---//
/// \\brief Rewrites all embedded expressions in an object x
/// \\param x an object containing expressions
/// \\param R a rewriter
template <typename T, typename Rewriter>
void rewrite(T& x,
             Rewriter R
            )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_rewrite_data_expressions_builder<pres_system::data_expression_builder>(R).update(x);
}

/// \\brief Rewrites all embedded expressions in an object x
/// \\param x an object containing expressions
/// \\param R a rewriter
/// \\return the rewrite result
template <typename T, typename Rewriter>
T rewrite(const T& x,
          Rewriter R
         )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::make_rewrite_data_expressions_builder<pres_system::data_expression_builder>(R).apply(result, x);
  return result;
}

/// \\brief Rewrites all embedded expressions in an object x, and applies a substitution to variables on the fly
/// \\param x an object containing expressions
/// \\param R a rewriter
/// \\param sigma a substitution
template <typename T, typename Rewriter, data::IsSubstitution Substitution>
void rewrite(T& x,
             Rewriter R,
             const Substitution& sigma
            )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_rewrite_data_expressions_with_substitution_builder<pres_system::data_expression_builder>(R, sigma).update(x);
}

/// \\brief Rewrites all embedded expressions in an object x, and applies a substitution to variables on the fly
/// \\param x an object containing expressions
/// \\param R a rewriter
/// \\param sigma a substitution
/// \\return the rewrite result
template <typename T, typename Rewriter, data::IsSubstitution Substitution>
T rewrite(const T& x,
          Rewriter R,
          const Substitution& sigma
         )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result; 
  data::detail::make_rewrite_data_expressions_with_substitution_builder<pres_system::data_expression_builder>(R, sigma).apply(result, x);
  return result;
}
//--- end generated pres_system rewrite code ---//

/// \brief Rewrites all embedded pres expressions in an object x
/// \param x an object containing expressions
/// \param R a pres rewriter
template <typename T, typename Rewriter>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void pres_rewrite(T& x, const Rewriter& R)
{
  pres_system::detail::make_rewrite_pres_expressions_builder<pres_system::pres_expression_builder>(R).update(x);
}

/// \brief Rewrites all embedded pres expressions in an object x
/// \param x an object containing expressions
/// \param R a pres rewriter
/// \return the rewrite result
template <typename T, typename Rewriter>
  requires(std::is_base_of_v<atermpp::aterm, T>)
T pres_rewrite(const T& x, const Rewriter& R)
{
  T result;
  pres_system::detail::make_rewrite_pres_expressions_builder<pres_system::pres_expression_builder>(R).apply(result, x);
  return result;
}

/// \brief Rewrites all embedded pres expressions in an object x, and applies a substitution to variables on the fly
/// \param x an object containing expressions
/// \param R a pres rewriter
/// \param sigma a substitution
template <typename T, typename Rewriter, typename Substitution>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void pres_rewrite(T& x,
    const Rewriter& R,
    Substitution sigma)
{
  pres_system::detail::make_rewrite_pres_expressions_with_substitution_builder<pres_system::pres_expression_builder>(R, sigma).update(x);
}

/// \brief Rewrites all embedded pres expressions in an object x, and applies a substitution to variables on the fly
/// \param x an object containing expressions
/// \param R a pres rewriter
/// \param sigma a substitution
/// \return the rewrite result
template <typename T, typename Rewriter, typename Substitution>
  requires(std::is_base_of_v<atermpp::aterm, T>)
T pres_rewrite(const T& x,
    const Rewriter& R,
    Substitution sigma)
{
  T result;
  pres_system::detail::make_rewrite_pres_expressions_with_substitution_builder<pres_system::pres_expression_builder>(R, sigma).apply(result, x);
  return result;
}

} // namespace mcrl2::pres_system



#endif // MCRL2_PRES_REWRITE_H
