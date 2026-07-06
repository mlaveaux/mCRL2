// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/one_point_rule_rewrite.h
/// \brief Applies the one point rule rewriter to the data expressions in an LPS.

#ifndef MCRL2_LPS_ONE_POINT_RULE_REWRITE_H
#define MCRL2_LPS_ONE_POINT_RULE_REWRITE_H

#include "mcrl2/data/rewriters/one_point_rule_rewriter.h"
#include "mcrl2/lps/builder.h"

namespace mcrl2::lps {

namespace detail {



struct one_point_rule_rewrite_builder: public lps::data_expression_builder<one_point_rule_rewrite_builder>
{
  using super = lps::data_expression_builder<one_point_rule_rewrite_builder>;
  using super::apply;
  using super::update;

  data::one_point_rule_rewriter r;

  template <class T>
  void apply(T& result, const data::data_expression& x)
  {
    result = r(x);
  }
};

} // namespace detail

/// \brief Applies the one point rule rewriter to all embedded data expressions in an object x
/// \param x an object containing data expressions
template <typename T>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void one_point_rule_rewrite(T& x)
{
  detail::one_point_rule_rewrite_builder f{};
  f.update(x);
}

/// \brief Applies the one point rule rewriter to all embedded data expressions in an object x
/// \param x an object containing data expressions
/// \return the rewrite result
template <typename T>
  requires(std::is_base_of_v<atermpp::aterm, T>)
T one_point_rule_rewrite(const T& x)
{
  T result;
  detail::one_point_rule_rewrite_builder f{};
  f.apply(result, x);
  return result;
}

} // namespace mcrl2::lps



#endif // MCRL2_LPS_ONE_POINT_RULE_REWRITE_H
