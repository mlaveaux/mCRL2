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

/// \brief Implements the naive matching algorithm with several (simple) extensions.
/// \details Can enable head symbol indexing as an optimization.
template<typename Substitution>
class NaiveMatcher
{
public:

  /// \brief A (simplistic) iterator over the matching results.
  class const_iterator
  {
  public:
    const_iterator(NaiveMatcher& matcher,
      const data_expression& term,
      std::size_t head_index,
      std::size_t current_index)
      : m_matcher(matcher),
        m_term(term),
        m_head_index(head_index),
        m_current_index(current_index)
    {}

    void operator++()
    {
      ++m_current_index;
    }

    matching_result<Substitution> operator*()
    {
      return m_matcher.next(m_term, m_head_index, m_current_index);
    }

  private:
    NaiveMatcher& m_matcher;
    const data_expression& m_term;

    std::size_t m_head_index = 0;
    std::size_t m_current_index;
  };

  /// \brief Initialize a naive matcher with a number of equations.
  NaiveMatcher(const data_equation_vector& equations);
  virtual ~NaiveMatcher() {}

  /// \brief Returns an iterator to the set of matching results.
  const_iterator match(const data_expression&);

private:
  /// \brief A function that is used to obtain the next matching result.
  matching_result<Substitution> next(const data_expression& term, std::size_t head_index, std::size_t index);

  using variable_partition = std::vector<variable>;

  /// \brief A mapping from head symbols to rewrite rules and their corresponding construction stacks. A unique index is used for each head symbol to achieve
  ///        this mapping without an unordered_map for performance reasons.
  std::vector<std::vector<linear_data_equation>> m_rewrite_system;

  /// \brief The original list of equations to use for matching.
  std::vector<linear_data_equation> m_equations;

  /// \brief The matching substitution computed in next.
  Substitution m_matching_sigma;
};

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_NAIVE_MATCHER_H
