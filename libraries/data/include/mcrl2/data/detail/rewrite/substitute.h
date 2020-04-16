// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_REWRITE_SUBSTITUTE_H
#define MCRL2_DATA_DETAIL_REWRITE_SUBSTITUTE_H

#include "mcrl2/data/abstraction.h"
#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/variable.h"
#include "mcrl2/data/where_clause.h"
#include "mcrl2/utilities/stack_array.h"

#include <assert.h>

namespace mcrl2
{
namespace data
{
namespace detail
{

constexpr bool NaiveSubstitution = true;

/// \brief Renames (by adapting the substitution sigma) the bound variables of the given abstraction.
/// \returns The list of new variables.
template<typename Substitution, typename Generator>
inline data::variable_list rename_bound_variables(const data::abstraction& abstract, Substitution& sigma, Generator& generator)
{
  data::variable_list new_variables;
  if constexpr (NaiveSubstitution)
  {
    // For every variable introduce a fresh variable.
    for (auto& var : abstract.variables())
    {
      variable fresh_variable(generator(), var.sort());
      new_variables.push_front(fresh_variable);
      sigma[var] = fresh_variable;
    }
  }
  else
  {
    for (auto& var : abstract.variables())
    {
      // If the given variable occurs in a right-hand side or the substitution is non-trivial (i.e. sigma(x) != x) then
      // this abstraction binds an existing variable and must be renamed.
      if (sigma.variable_occurs_in_a_rhs(var) || sigma(var) != var)
      {
        variable fresh_variable(generator(), var.sort());
        new_variables.push_front(fresh_variable);
        sigma[var] = fresh_variable;
      }
      else
      {
        new_variables.push_front(var);
      }
    }

  }

  return new_variables;
}

/// \brief Renames (by adapting the substitution sigma) the bound variables of the given where clause.
/// \returns The list of new variables, and updates the substitution to change these variables.
template<typename Substitution, typename Generator>
inline data::variable_list rename_bound_variables(const where_clause& clause, Substitution& sigma, Generator& generator)
{
  data::variable_list new_variables;
  if constexpr (NaiveSubstitution)
  {
    // For every variable introduce a fresh variable.
    for (auto& assignment : clause.assignments())
    {
      const variable& var = assignment.lhs();

      variable fresh_variable(generator(), var.sort());
      new_variables.push_front(fresh_variable);
      sigma[var] = fresh_variable;
    }
  }
  else
  {
    for (auto& assignment : clause.assignments())
    {
      const variable& var = assignment.lhs();

      // If the given variable occurs in a right-hand side or the substitution is non-trivial (i.e. sigma(x) != x) then
      // this abstraction binds an existing variable and must be renamed.
      if (sigma.variable_occurs_in_a_rhs(var) || sigma(var) != var)
      {
        variable fresh_variable(generator(), var.sort());
        new_variables.push_front(fresh_variable);
        sigma[var] = fresh_variable;
      }
      else
      {
        new_variables.push_front(var);
      }
    }

  }

  return new_variables;
}

/// \brief An identity transformer.
inline data_expression identity(const data_expression& expression) { return expression; }

/// \brief A capture-avoiding substitution of sigma applied to the given term.
template<typename Substitution, typename Generator, typename Transformer>
inline data_expression capture_avoiding_substitution(const data_expression& term, Substitution& sigma, Generator& generator, Transformer f)
{
  // C(x, sigma, V) = sigma(x), where x is a variable
  if (is_variable(term))
  {
    const auto& var = static_cast<const variable&>(term);
    return f(sigma(var));
  }
  // C(f, sigma, V) = f, where f is a function symbol.
  else if (is_function_symbol(term))
  {
    return term;
  }
  // C(lambda x . t, sigma, V) = lambda y . C(t, sigma[x := y], V), where x are variables and y are fresh variables.
  else if (is_abstraction(term))
  {
    const auto& abstract = static_cast<const data::abstraction&>(term);
    data::variable_list new_variables = rename_bound_variables(abstract, sigma, generator);
    return abstraction(abstract.binding_operator(), new_variables, capture_avoiding_substitution(abstract.body(), sigma, generator, f));
  }
  // C(t(t_1, ..., t_n, sigma, V) = C(t, sigma, V) ( C(t_1, sigma, V), ..., C(t_n, sigma, V) )
  else
  {
    assert(is_application(term));
    const auto& appl = static_cast<const application&>(term);

    // Substitution of all arguments.
    MCRL2_DECLARE_STACK_ARRAY(arguments, data_expression, appl.size());
    for (std::size_t index = 0; index < appl.size(); ++index)
    {
      arguments[index] = capture_avoiding_substitution(appl[index], sigma, generator, f);
    }

    // Construct the application, also subsituting the head.
    return application(capture_avoiding_substitution(appl.head(), sigma, generator, f), arguments.begin(), arguments.end());
  }
}

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_REWRITE_SUBSTITUTE_H
