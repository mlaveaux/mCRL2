// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbes/include/mcrl2/pbes/rewrite.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_REWRITE_H
#define MCRL2_PBES_REWRITE_H

#include "mcrl2/data/rewrite.h"
#include "mcrl2/pbes/builder.h"



namespace mcrl2::pbes_system {

namespace detail
{

/// \cond INTERNAL_DOCS
template <template <class> class Builder, class Rewriter>
struct rewrite_pbes_expressions_builder: public Builder<rewrite_pbes_expressions_builder<Builder, Rewriter> >
{
  using super = Builder<rewrite_pbes_expressions_builder<Builder, Rewriter>>;
  using super::apply;
  using super::update;

  const Rewriter& R;

  rewrite_pbes_expressions_builder(const Rewriter& R_)
    : R(R_)
  {}

  template <class T>
  void apply(T& result, const pbes_expression& x)
  {
    result = R(x);
  }

  void update(pbes_system::pbes& x)
  {
    super::update(x);
    // Handle the initial state. It is skipped by the builder because the type is not pbes_expression
    pbes_expression initial_state;
    apply(initial_state, static_cast<pbes_expression>(x.initial_state()));
    x.initial_state() = atermpp::down_cast<propositional_variable_instantiation>(initial_state);
  }
};

template <template <class> class Builder, class Rewriter>
rewrite_pbes_expressions_builder<Builder, Rewriter>
make_rewrite_pbes_expressions_builder(const Rewriter& R)
{
  return rewrite_pbes_expressions_builder<Builder, Rewriter>(R);
}

template <template <class> class Builder, class Rewriter, class Substitution>
struct rewrite_pbes_expressions_with_substitution_builder: public Builder<rewrite_pbes_expressions_with_substitution_builder<Builder, Rewriter, Substitution> >
{
  using super = Builder<rewrite_pbes_expressions_with_substitution_builder<Builder, Rewriter, Substitution>>;
  using super::apply;

  const Rewriter& R;
  Substitution& sigma;

  rewrite_pbes_expressions_with_substitution_builder(const Rewriter& R_, Substitution& sigma_)
    : R(R_),
      sigma(sigma_)
  {}

  template <class T>
  void apply(T& result, const pbes_expression& x)
  {
    result = R(x, sigma);
  }
};

template <template <class> class Builder, class Rewriter, class Substitution>
rewrite_pbes_expressions_with_substitution_builder<Builder, Rewriter, Substitution>
make_rewrite_pbes_expressions_with_substitution_builder(const Rewriter& R, Substitution& sigma)
{
  return rewrite_pbes_expressions_with_substitution_builder<Builder, Rewriter, Substitution>(R, sigma);
}
/// \endcond

} // namespace detail

//--- start generated pbes_system rewrite code ---//
/// \\brief Rewrites all embedded expressions in an object x
/// \\param x an object containing expressions
/// \\param R a rewriter
template <typename T, typename Rewriter>
void rewrite(T& x,
             Rewriter R
            )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_rewrite_data_expressions_builder<pbes_system::data_expression_builder>(R).update(x);
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
  data::detail::make_rewrite_data_expressions_builder<pbes_system::data_expression_builder>(R).apply(result, x);
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
  data::detail::make_rewrite_data_expressions_with_substitution_builder<pbes_system::data_expression_builder>(R, sigma).update(x);
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
  data::detail::make_rewrite_data_expressions_with_substitution_builder<pbes_system::data_expression_builder>(R, sigma).apply(result, x);
  return result;
}
//--- end generated pbes_system rewrite code ---//

/// \brief Rewrites all embedded pbes expressions in an object x
/// \param x an object containing expressions
/// \param R a pbes rewriter
template <typename T, typename Rewriter>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void pbes_rewrite(T& x, const Rewriter& R)
{
  pbes_system::detail::make_rewrite_pbes_expressions_builder<pbes_system::pbes_expression_builder>(R).update(x);
}

/// \brief Rewrites all embedded pbes expressions in an object x
/// \param x an object containing expressions
/// \param R a pbes rewriter
/// \return the rewrite result
template <typename T, typename Rewriter>
  requires(std::is_base_of_v<atermpp::aterm, T>)
T pbes_rewrite(const T& x, const Rewriter& R)
{
  T result;
  pbes_system::detail::make_rewrite_pbes_expressions_builder<pbes_system::pbes_expression_builder>(R).apply(result, x);
  return result;
}

/// \brief Rewrites all embedded pbes expressions in an object x, and applies a substitution to variables on the fly
/// \param x an object containing expressions
/// \param R a pbes rewriter
/// \param sigma a substitution
template <typename T, typename Rewriter, typename Substitution>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void pbes_rewrite(T& x,
    const Rewriter& R,
    Substitution sigma)
{
  pbes_system::detail::make_rewrite_pbes_expressions_with_substitution_builder<pbes_system::pbes_expression_builder>(R, sigma).update(x);
}

/// \brief Rewrites all embedded pbes expressions in an object x, and applies a substitution to variables on the fly
/// \param x an object containing expressions
/// \param R a pbes rewriter
/// \param sigma a substitution
/// \return the rewrite result
template <typename T, typename Rewriter, typename Substitution>
  requires(std::is_base_of_v<atermpp::aterm, T>)
T pbes_rewrite(const T& x,
    const Rewriter& R,
    Substitution sigma)
{
  T result;
  pbes_system::detail::make_rewrite_pbes_expressions_with_substitution_builder<pbes_system::pbes_expression_builder>(R, sigma).apply(result, x);
  return result;
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_REWRITE_H
