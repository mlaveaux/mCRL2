// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/transformations.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_TRANSFORMATIONS_H
#define MCRL2_PBES_TRANSFORMATIONS_H

#include "mcrl2/pbes/builder.h"

namespace mcrl2::pbes_system {

namespace detail {

struct order_quantified_variables_builder: public pbes_expression_builder<order_quantified_variables_builder>
{
  using super = pbes_expression_builder<order_quantified_variables_builder>;
  using super::apply;

  const data::data_specification& dataspec;

  order_quantified_variables_builder(const data::data_specification& dataspec_)
    : dataspec(dataspec_)
  {}

  template <class T>
  void apply(T& result, const forall& x)
  {
    pbes_expression body;
    apply(body, x.body());
    pbes_system::make_optimized_forall(result, data::order_variables_to_optimise_enumeration(x.variables(), dataspec), body);
  }

  template <class T>
  void apply(T& result, const exists& x)
  {
    pbes_expression body;
    apply(body, x.body());
    pbes_system::make_optimized_exists(result, data::order_variables_to_optimise_enumeration(x.variables(), dataspec), body);
  }
};

} // namespace detail

inline
pbes_expression order_quantified_variables(const pbes_expression& x, const data::data_specification& dataspec)
{
  detail::order_quantified_variables_builder f(dataspec);
  pbes_expression result;
  f.apply(result, x);
  return result;
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_TRANSFORMATIONS_H
