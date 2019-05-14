// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_NAIVE_MATCHER_H
#define MCRL2_DATA_DETAIL_NAIVE_MATCHER_H

#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/detail/match/construction_stack.h"
#include "mcrl2/data/substitutions/mutable_indexed_substitution.h"

#include <vector>

namespace mcrl2
{
namespace data
{
namespace detail
{

using data_equation_extended = std::tuple<data_equation, ConstructionStack, ConstructionStack>;

class NaiveMatcher
{
public:
  /// \brief Initialize a naive matcher with a number of equations.
  NaiveMatcher(const data_equation_vector& equations);

  /// \brief The match function defined in the document. However, instead of returning a set of right-hand sides it makes a (arbitrary)
  ///        choice of which right-hand side to return and applies the matching substitution to it. The given term must be in normal form.
  /// \returns A number of equations that match the given term, under the updated (passed through) substitution sigma.
  std::vector<std::reference_wrapper<const data_equation_extended>> match(const data_expression& term, mutable_indexed_substitution<>& matching_sigma);

private:
  /// A mapping from head symbols to rewrite rules and their corresponding construction stacks. A unique index is used for each head symbol to achieve
  /// this mapping without an unordered_map for performance reasons.
  std::vector<std::vector<data_equation_extended>> m_rewrite_system;

  /// The original list of equations to use for matching.
  std::vector<data_equation_extended> m_equations;
};

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_NAIVE_MATCHER_H
