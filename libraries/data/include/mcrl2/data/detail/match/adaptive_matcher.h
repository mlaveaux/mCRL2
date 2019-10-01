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
#include "mcrl2/data/detail/match/consistency.h"
#include "mcrl2/utilities/indexed_set.h"

#include <vector>

namespace mcrl2
{
namespace data
{
namespace detail
{

using linear_data_equation = std::pair<data_equation_extended, partition>;

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

  const data_equation_extended* next(Substitution& matching_sigma) override;

private:

  /// \brief Represents the unnamed variable omega.
  class not_equal : public function_symbol
  {
  public:
    not_equal()
      : function_symbol(core::identifier_string("@@not_equal"), data::sort_expression())
    {}
  };

  /// \returns True iff the given data expression is of type omega.
  static inline bool is_not_equal(const data_expression& expression)
  {
    return expression == not_equal();
  }

  /// \brief Adds states and transitions to the APMA for the given state and prefix.
  void construct_apma(std::size_t s, data_expression pref);

  /// \brief A local copy of the equations to keep the references to equations stable.
  std::deque<linear_data_equation> m_linear_equations;

  // The underlying automaton.

  /// \brief Each state in the pattern match automaton is labelled with a set of variables.
  class apma_state
  {
  public:
    variable variable;
    std::size_t position; ///< The index of the position to be inspected.
    std::vector<std::reference_wrapper<const linear_data_equation>> match_set;
    std::vector<std::size_t> argument_positions;
  };

  // Information about the underlying automata.

  IndexedAutomaton<apma_state> m_automaton;

  mcrl2::utilities::indexed_set<position> m_positions;

  std::size_t m_not_equal_index;

  // Store information about the match.

  Substitution m_matching_sigma;

  std::vector<data_expression> m_subterms; ///< A mapping from indices to subterms.

  std::vector<std::reference_wrapper<const linear_data_equation>>* m_match_set = nullptr;

  std::size_t m_match_index; ///< The index of the equation that should be returned by the call to next.
};

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_ADAPTIVE_MATCHER_H
