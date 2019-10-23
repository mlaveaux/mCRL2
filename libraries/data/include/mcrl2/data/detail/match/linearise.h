// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_CONSISTENCY_H
#define MCRL2_DATA_CONSISTENCY_H

#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/detail/match/position.h"
#include "mcrl2/data/variable.h"
#include "mcrl2/atermpp/aterm_io.h"

#include <set>
#include <unordered_map>
#include <vector>

namespace mcrl2::data::detail
{

/// \brief The consistency class is a set of variables, must ensure unique entries for optimal performance..
using consistency_class = std::vector<position>;

/// \brief The consistency partition is a set of consistency classes, the underlying vectors should be sets.
using consistency_partition = std::vector<consistency_class>;

/// \brief Given an equation renames multiple occurrences of the same variable to position variables and returns the
///        equivalence class that must be checked for consistency.
inline
std::pair<data_equation, consistency_partition> make_linear(const data_equation& equation)
{
  // Find the variable of the left-hand side.
  std::set<position> var_positions = fringe(equation.lhs());

  // Find the equivalence classes based on these positions.
  std::unordered_map<variable, consistency_class> mapping;

  // Keep track of a variable at position to its position variable.
  mutable_map_substitution<std::map<variable, variable>> sigma;

  // Replace each variable by its position variable (making it a position annotated term).
  data_expression lhs = equation.lhs();
  for (const position& pos : var_positions)
  {
    // This should always be a valid variable by definition of fringe.
    variable var = static_cast<variable>(at_position(equation.lhs(), pos).value());
    assert(is_variable(var));

    // Keep track of the mapping to consistency classes, and the substitution of var to position_variable(pos).
    mapping[var].push_back(pos);
    sigma[var] = position_variable(pos);
    lhs = assign_at_position(lhs, pos, position_variable(pos));
  }

  // Construct the partition and collect the variables.
  variable_list variables;
  consistency_partition partition;
  for (const auto& [var, consistency] : mapping)
  {
    variables.push_front(var);

    // Ignore trivial classes.
    if (consistency.size() >= 2)
    {
      partition.push_back(consistency);
    }
  }

  // Replace the variables in the right-hand side and the condition.
  data_expression rhs = data::replace_variables(equation.rhs(), sigma);
  data_expression condition = data::replace_variables(equation.rhs(), sigma);

  // Create the new data equation.
  data_equation linear_equation(variables, condition, lhs, rhs);

  return std::make_pair(linear_equation, partition);
}

} // namespace mcrl2::data::detail

#endif // MCRL2_DATA_CONSISTENCY_H
