// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file data/include/mcrl2/data/rewrite.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_REWRITE_H
#define MCRL2_DATA_REWRITE_H

#include "mcrl2/data/concepts.h"
#include "mcrl2/data/builder.h"

namespace mcrl2::data {

namespace detail
{

template <template <class> class Builder, class Rewriter>
struct rewrite_data_expressions_builder: public Builder<rewrite_data_expressions_builder<Builder, Rewriter> >
{
  using super = Builder<rewrite_data_expressions_builder<Builder, Rewriter>>;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  Rewriter R;

  rewrite_data_expressions_builder(Rewriter R_)
    : R(R_)
  {}

  template <class T>
  void apply(T& result, const data_expression& x)
  {
    result = R(x);
  }
};

template <template <class> class Builder, class Rewriter>
rewrite_data_expressions_builder<Builder, Rewriter>
make_rewrite_data_expressions_builder(Rewriter R)
{
  return rewrite_data_expressions_builder<Builder, Rewriter>(R);
}

template <template <class> class Builder, class Rewriter, IsSubstitution Substitution>
struct rewrite_data_expressions_with_substitution_builder: public Builder<rewrite_data_expressions_with_substitution_builder<Builder, Rewriter, Substitution> >
{
  using super = Builder<rewrite_data_expressions_with_substitution_builder<Builder, Rewriter, Substitution>>;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  Rewriter R;
  Substitution sigma;

  rewrite_data_expressions_with_substitution_builder(Rewriter R_, Substitution sigma_)
    : R(R_),
      sigma(sigma_)
  {}

  template <class T>
  void apply(T& result, const data_expression& x)
  {
    result = R(x, sigma);
  }
};

template <template <class> class Builder, class Rewriter, IsSubstitution Substitution>
rewrite_data_expressions_with_substitution_builder<Builder, Rewriter, Substitution>
make_rewrite_data_expressions_with_substitution_builder(Rewriter R, Substitution sigma)
{
  return rewrite_data_expressions_with_substitution_builder<Builder, Rewriter, Substitution>(R, sigma);
}

} // namespace detail

//--- start generated data rewrite code ---//
/// \\brief Rewrites all embedded expressions in an object x
/// \\param x an object containing expressions
/// \\param R a rewriter
template <typename T, typename Rewriter>
void rewrite(T& x,
             Rewriter R
            )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_rewrite_data_expressions_builder<data::data_expression_builder>(R).update(x);
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
  data::detail::make_rewrite_data_expressions_builder<data::data_expression_builder>(R).apply(result, x);
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
  data::detail::make_rewrite_data_expressions_with_substitution_builder<data::data_expression_builder>(R, sigma).update(x);
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
  data::detail::make_rewrite_data_expressions_with_substitution_builder<data::data_expression_builder>(R, sigma).apply(result, x);
  return result;
}
//--- end generated data rewrite code ---//

} // namespace mcrl2::data



#endif // MCRL2_DATA_REWRITE_H
