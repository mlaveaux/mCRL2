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

#include "mcrl2/data/detail/match/matcher.h"
#include "mcrl2/data/detail/match/naive_consistency.h"

#include <vector>

namespace mcrl2
{
namespace data
{
namespace detail
{

template<typename Substitution>
class NaiveMatcher final : public Matcher<Substitution>
{
public:
  /// \brief Initialize a naive matcher with a number of equations.
  NaiveMatcher(const data_equation_vector& equations);
  virtual ~NaiveMatcher() {}

  // Matcher interface

  void match(const data_expression& term) override;

  matching_result<Substitution> next() override;

private:

  using variable_partition = std::vector<variable>;

  /// \brief A mapping from head symbols to rewrite rules and their corresponding construction stacks. A unique index is used for each head symbol to achieve
  ///        this mapping without an unordered_map for performance reasons.
  std::vector<std::vector<linear_data_equation>> m_rewrite_system;

  /// \brief The original list of equations to use for matching.
  std::vector<linear_data_equation> m_equations;

  /// \brief The matching substitution computed in next.
  Substitution m_matching_sigma;

  // Information about the matched term.

  std::size_t m_current_index = 0;
  std::size_t m_head_index;
  data_expression m_term;
};

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_NAIVE_MATCHER_H
