// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_ADAPTIVE_MATCHER_H
#define MCRL2_DATA_DETAIL_ADAPTIVE_MATCHER_H

#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/detail/match/matcher.h"
#include "mcrl2/data/detail/match/automaton.h"
#include "mcrl2/data/detail/match/position.h"
#include "mcrl2/data/detail/match/automaton_consistency.h"
#include "mcrl2/utilities/indexed_set.h"

#include <vector>

namespace mcrl2
{
namespace data
{
namespace detail
{

/// \brief Represents the unnamed variable omega.
class not_equal : public function_symbol
{
public:
  not_equal()
    : function_symbol(core::identifier_string("@@not_equal"), data::sort_expression())
  {}
};

/// \returns True iff the given data expression is of type omega.
inline bool is_not_equal(const data_expression& expression)
{
  return expression == not_equal();
}

/// \brief An adaptive non-linear pattern matching automata
template<typename Substitution>
class AdaptiveMatcher final
{
public:
  /// \brief Initialize the automaton matcher with a number of equations.
  AdaptiveMatcher(const data_equation_vector& equations);
  virtual ~AdaptiveMatcher() {}

  /// \brief A (simplistic) iterator over the matching results.
  class const_iterator
  {
  public:
    const_iterator(Substitution& sigma)
      : m_matching_sigma(sigma)
    {}

    /// \brief Construct an iterator over the match_set, computing the consistent patterns from the subterm indexing
    const_iterator(std::vector<indexed_linear_data_equation>* match_set,
      Substitution& sigma)
      : m_match_set(match_set),
        m_matching_sigma(sigma)
    {}

    void operator++()
    {
      // Find the next consistent match.
      ++m_current_index;
    }

    const extended_data_equation& operator*()
    {
      return (*m_match_set)[m_current_index];
    }

    bool operator!=(std::nullptr_t) const
    {
      return (m_match_set != nullptr && m_current_index < m_match_set->size());
    }

  private:
    std::vector<indexed_linear_data_equation>* m_match_set = nullptr;
    Substitution& m_matching_sigma;

    std::vector<bool> m_consistent;
    std::size_t m_current_index = 0;
  };

  // Matcher interface.

  /// \details We assume that the given substitution is the identity.
  const_iterator match(const data_expression& term, Substitution& matching_sigma);

private:

  /// \brief Each state in the pattern match automaton is labelled with a set of variables.
  struct apma_state
  {
    bool is_consistency_state() const { return compare.first != compare.second; }

    bool is_matching_state() const { return position != std::numeric_limits<std::size_t>::max(); }

    // Matching states:
    std::size_t position = std::numeric_limits<std::size_t>::max(); ///< L, the index of the position to be inspected.
    std::vector<std::size_t> argument_positions; ///< These are the positions where arguments must be stored in the subterm table.

    // Consistency states:
    std::pair<std::size_t, std::size_t> compare = std::make_pair(0,0); ///< L, the pair of positions that must be compared.

    // Final states:
    std::vector<indexed_linear_data_equation> match_set; ///< L, the equations that matched.
    std::vector<std::pair<data::variable, std::size_t>> assignments; ///< P, the variables that must be assigned and their subterm position.
  };

  // The underlying automaton.
  using Automaton = IndexedAutomaton<apma_state>;

  /// \returns A subautomaton for the given state and prefix.
  Automaton construct_apma(data_expression pref,
    const std::vector<indexed_linear_data_equation>& L,
    const std::set<std::pair<std::size_t, std::size_t> >& E,
    const std::set<std::pair<std::size_t, std::size_t> >& N
    );

  /// \returns A subset of the given positions, filtering out unecessary positions.
  std::set<position> restrict(const std::set<position>& positions, const std::vector<indexed_linear_data_equation>& L);

  /// \returns A single position selected from a set of positions.
  position select(const std::set<position>& positions, const std::vector<indexed_linear_data_equation>& L);

  /// \returns A unique and consistent index for the given position.
  std::size_t position_index(const position& positions);

  // Information about the underlying automata.

  Automaton m_automaton;

  mcrl2::utilities::indexed_set<position> m_positions;

  std::size_t m_dontcare_index = 0; ///< The index of the dontcare data_expression.
  std::size_t m_equal_index = 0;
  std::size_t m_not_equal_index = 0;

  // Store (temporary) information about the match that can be invalidated.

  std::vector<atermpp::unprotected_aterm> m_subterms; ///< A mapping from indices to subterms.
};

} // namespace detail
} // namespace data
} // namespace mcrl2

// Explicit instantiations.


#endif // MCRL2_DATA_DETAIL_ADAPTIVE_MATCHER_H
