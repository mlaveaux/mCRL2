// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_REWRITE_UTILITY_H
#define MCRL2_DATA_DETAIL_REWRITE_UTILITY_H

#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/variable.h"

namespace mcrl2
{
namespace data
{
namespace detail
{

/// \brief A capture-avoiding substitution of sigma applied to the given term.
template<typename Substitution>
static data_expression capture_avoiding_substitution(const data_expression& term, const Substitution& sigma)
{
  // C(x, sigma, V) = sigma(x), where x is a variable
  if (is_variable(term))
  {
    const auto& var = static_cast<const variable&>(term);
    return sigma.at(var);
  }
  // C(f, sigma, V) = f, where f is a function symbol.
  else if (is_function_symbol(term))
  {
    return term;
  }
  // C(lambda x . t, sigma, V) = lambda y . C(t, sigma[x := y], V), where x and y are variables.
  else if (is_abstraction(term))
  {
    const auto& abstraction = static_cast<const class abstraction&>(term);
    assert(false);
  }
  // C(t(t_1, ..., t_n, sigma, V) = C(t, sigma, V) ( C(t_1, sigma, V), ..., C(t_n, sigma, V) )
  else
  {
    assert(is_application(term));
    const auto& appl = static_cast<const application&>(term);

    // Substitution of all arguments.
    std::vector<data_expression> arguments(appl.size());
    for (std::size_t index = 0; index < appl.size(); ++index)
    {
      arguments[index] = capture_avoiding_substitution(appl[index], sigma);
    }

    // Construct the application, also subsituting the head.
    return application(capture_avoiding_substitution(appl.head(), sigma), arguments.begin(), arguments.end());
  }
}

/// \brief Matches a single left-hand side with the given term and creates the substitution.
template<typename Substitution>
static bool match_lhs(const data_expression& term,  const data_expression& lhs, Substitution& sigma)
{
  if (is_function_symbol(lhs))
  {
    return term == lhs;
  }
  else if (is_variable(lhs))
  {
    const auto& var = static_cast<const variable&>(lhs);

    if (sigma.count(var))
    {
      // If the variable was already assigned they must match.
      return sigma.at(var) == term;
    }
    else
    {
      // Else substitute the given term for the current variable.
      sigma[var] = term;
      return true;
    }
  }
  else
  {
    // The term and lhs are applications.
    const application& lhs_appl  = static_cast<const application&>(lhs);
    const application& term_appl = static_cast<const application&>(term);

    // Both must have the same arity, the head symbol must match and their arguments must match.
    if (lhs_appl.size() != term_appl.size())
    {
      return false;
    }

    if (!match_lhs(term_appl.head(), lhs_appl.head(), sigma))
    {
      return false;
    }

    for (std::size_t i = 0; i < term_appl.size(); i++)
    {
      if (!match_lhs(term_appl[i], lhs_appl[i], sigma))
      {
        return false;
      }
    }

    return true;
  }
}

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_REWRITE_UTILITY_H
