// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/detail/check_well_formed_bes.h
/// \brief This file contains a function that checks whether a pbes expression is a bes in the sense
//         that is does not contain quantifiers or boolean expressions not equal to true and false.

#ifndef MCRL2_PBES_DETAIL_CHECK_WELL_FORMED_BES_H
#define MCRL2_PBES_DETAIL_CHECK_WELL_FORMED_BES_H

#include "mcrl2/pbes/traverser.h"



namespace mcrl2::pbes_system {

/// \cond INTERNAL_DOCS
/// \brief Visitor for checking if a pbes expression is in BES form.
struct is_well_formed_bes_traverser: public pbes_expression_traverser<is_well_formed_bes_traverser>
{
  using super = pbes_expression_traverser<is_well_formed_bes_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  is_well_formed_bes_traverser() = default;

  void enter(const forall& x)
  {
    throw mcrl2::runtime_error("Bes contains a universal quantifier.\n" + pp(x));
  }

  void enter(const exists& x)
  {
    throw mcrl2::runtime_error("Bes contains an existential quantifier.\n" + pp(x));
  }

  void enter(const data::data_expression& x)
  {
    if (x!=data::sort_bool::true_() && x!=data::sort_bool::false_())
    {
      throw mcrl2::runtime_error("Bes contains an expression that is not equal to true or false.\n" + pp(x));
    }
  } 
};
/// \endcond

/// \brief This function checks whether x is a well formed bes expression, without quantifiers and
//         where all data expressions are equal to either true or false.
/// \param x a pbes expression to be checked.
inline void check_whether_argument_is_a_well_formed_bes(const pbes_expression& x)
{
  is_well_formed_bes_traverser f;
  f.apply(x);
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_DETAIL_CHECK_WELL_FORMED_BES_H
