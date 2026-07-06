// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/substitutions/variable_substitution.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_SUBSTITUTIONS_VARIABLE_SUBSTITUTION_H
#define MCRL2_DATA_SUBSTITUTIONS_VARIABLE_SUBSTITUTION_H

#include "mcrl2/data/is_simple_substitution.h"
#include "mcrl2/data/find.h"
#include "mcrl2/data/undefined.h"

namespace mcrl2::data
{

/// \brief Substitution that maps a single variable to a data expression
struct variable_substitution
{
  static constexpr bool is_identity_substitution=false;

  using variable_type = variable;
  using expression_type = data_expression;
  using argument_type = variable;
  using result_type = data_expression;

  const variable& lhs;
  const data_expression& rhs;

  variable_substitution(const variable& lhs_, const data_expression& rhs_)
  : lhs(lhs_),
    rhs(rhs_)
  {}

  const data_expression& operator()(const variable& x) const
  {
    if (x == lhs)
    {
      return rhs;
    }
    return x;
  }
};

inline
std::set<data::variable> substitution_variables(const variable_substitution& sigma)
{
  std::set<data::variable> result;
  data::find_free_variables(sigma.rhs, std::inserter(result, result.end()));
  return result;
}

} // namespace mcrl2::data

#endif // MCRL2_DATA_SUBSTITUTIONS_VARIABLE_SUBSTITUTION_H
