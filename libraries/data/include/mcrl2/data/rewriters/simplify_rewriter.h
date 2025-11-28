// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/rewriters/simplify_rewriter.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_REWRITERS_SIMPLIFY_REWRITER_H
#define MCRL2_DATA_REWRITERS_SIMPLIFY_REWRITER_H

#include "mcrl2/data/builder.h"
#include "mcrl2/data/expression_traits.h"
#include "mcrl2/data/optimized_boolean_operators.h"

namespace mcrl2::data
{

namespace detail {

template <typename Derived>
class simplify_rewrite_builder: public data_expression_builder<Derived>
{
  public:
    using super = data_expression_builder<Derived>;

    using super::apply;

    Derived& derived()
    {
      return static_cast<Derived&>(*this);
    }

    bool is_not(const data_expression& x) const
    {
      return data::sort_bool::is_not_application(x);
    }

    bool is_and(const data_expression& x) const
    {
      return data::sort_bool::is_and_application(x);
    }

    bool is_or(const data_expression& x) const
    {
      return data::sort_bool::is_or_application(x);
    }

    bool is_imp(const data_expression& x) const
    {
      return data::sort_bool::is_implies_application(x);
    }

    bool is_forall(const data_expression& x) const
    {
      return data::is_forall(x);
    }

    bool is_exists(const data_expression& x) const
    {
      return data::is_exists(x);
    }

    template <class T>
    void apply(T& result, const application& x)
    {
      derived().enter(x);
      if (is_not(x)) // x = !y
      {
        data_expression y;
        derived().apply(y, *x.begin());
        data::optimized_not(result, y);
      }
      else if (is_and(x)) // x = y && z
      {
        data_expression y;
        derived().apply(y, binary_left(x));
        data_expression z;
        derived().apply(z, binary_right(x));
        data::optimized_and(result, y, z);
      }
      else if (is_or(x)) // x = y || z
      {
        data_expression y;
        derived().apply(y, binary_left(x));
        data_expression z;
        derived().apply(z, binary_right(x));
        data::optimized_or(result, y, z);
      }
      else if (is_imp(x)) // x = y => z
      {
        data_expression y;
        derived().apply(y, binary_left(x));
        data_expression z;
        derived().apply(z, binary_right(x));
        data::optimized_imp(result, y, z);
      }
      else
      {
        super::apply(result, x);
      }
      derived().leave(x);
    }

    template <class T>
    void apply(T& result, const forall& x) // x = forall d. y
    {
      const variable_list& d = forall(x).variables();
      data_expression y;
      derived().apply(y, forall(x).body());
      data::optimized_forall(result, d, y, true);
    }

    template <class T>
    void apply(T& result, const exists& x) // x = exists d. y
    {
      const variable_list& d = exists(x).variables();
      data_expression y;
      derived().apply(y, exists(x).body());
      data::optimized_exists(result, d, y, true);
    }
};

} // namespace detail

struct simplify_rewriter
{
  using argument_type = data_expression;
  using result_type = data_expression;

  void operator()(data_expression& result, const data_expression& x) const
  {
    core::make_apply_builder<detail::simplify_rewrite_builder>().apply(result, x);
  }

  data_expression operator()(const data_expression& x) const
  {
    data_expression result;
    core::make_apply_builder<detail::simplify_rewrite_builder>().apply(result, x);
    return result;
  }
};

template <typename T>
void simplify(T& x, std::enable_if_t<!std::is_base_of_v<atermpp::aterm, T>>* = 0)
{
  core::make_update_apply_builder<data::data_expression_builder>(simplify_rewriter()).update(x);
}

template <typename T>
T simplify(const T& x, std::enable_if_t<std::is_base_of_v<atermpp::aterm, T>>* = nullptr)
{
  T result;
  core::make_update_apply_builder<data::data_expression_builder>(simplify_rewriter()).apply(result, x);
  return result;
}

} // namespace mcrl2::data

#endif // MCRL2_DATA_REWRITERS_SIMPLIFY_REWRITER_H

