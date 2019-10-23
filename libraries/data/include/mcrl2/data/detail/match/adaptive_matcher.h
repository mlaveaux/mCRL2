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

/// \brief An adaptive pattern matching automata based on the construction of "Adaptive Pattern Matching" by
///        R. C. Sekar et al. This has been extended with a postprocessing to obtain the matching subsitution.
template<typename Substitution>
class AdaptiveMatcher final : public Matcher<Substitution>
{
public:
  /// \brief Initialize the automaton matcher with a number of equations.
  AdaptiveMatcher(const data_equation_vector& equations);
  virtual ~AdaptiveMatcher() {}

  // Matcher interface.

  void match(const data_expression& term) override;

  matching_result<Substitution> next() override;

private:

  /// \brief Each state in the pattern match automaton is labelled with a set of variables.
  struct apma_state
  {
    // Other states:
    std::size_t position = std::numeric_limits<std::size_t>::max(); ///< L, the index of the position to be inspected.
    std::vector<std::size_t> argument_positions; ///< These are the positions where arguments must be stored in the subterm table.

    // Final states:
    std::vector<indexed_linear_data_equation> match_set; ///< L, the equations that matched.
    std::vector<std::pair<data::variable, std::size_t>> variables; ///< P, the variables that still must be assigned and their position.
  };

  // The underlying automaton.
  using Automaton = IndexedAutomaton<apma_state>;

  /// \brief Adds states and transitions to the APMA for the given state and prefix.
  Automaton construct_apma(const Automaton& automaton,
    std::size_t s,
    data_expression pref,
    const std::vector<indexed_linear_data_equation>& L);

  /// \returns A subset of the given positions, filtering out unecessary positions.
  std::set<position> restrict(const std::set<position>& positions, const std::vector<indexed_linear_data_equation>& L);

  /// \returns A single position selected from a set of positions.
  position select(const std::set<position>& positions);

  // Information about the underlying automata.

  Automaton m_automaton;

  mcrl2::utilities::indexed_set<position> m_positions;

  std::size_t m_not_equal_index = 0; ///< The index of the not_equal data_expression.

  // Metrics for the automaton.

  std::size_t m_nof_ambiguous_matches = 0; ///< The number of final states with multiple matches.
  std::size_t m_nof_final_states = 0; ///< The number of final states.

  // Store information about the match.

  Substitution m_matching_sigma;

  std::vector<atermpp::unprotected_aterm> m_subterms; ///< A mapping from indices to subterms.

  std::vector<indexed_linear_data_equation>* m_match_set = nullptr;

  std::size_t m_match_index; ///< The index of the equation that should be returned by the call to next.
};

} // namespace detail
} // namespace data
} // namespace mcrl2

// Explicit instantiations.


#endif // MCRL2_DATA_DETAIL_ADAPTIVE_MATCHER_H
