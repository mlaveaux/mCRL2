// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_AUTOMATON_MATCHER_H
#define MCRL2_DATA_DETAIL_AUTOMATON_MATCHER_H

#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/detail/match/matcher.h"
#include "mcrl2/data/substitutions/mutable_indexed_substitution.h"
#include "mcrl2/utilities/unordered_map.h"

#include <vector>

namespace mcrl2
{
namespace data
{
namespace detail
{

using Pattern = std::vector<data_expression>;
using PatternSet = std::set<Pattern>;
using equivalence_classes = std::vector<std::vector<variable>>;
using linear_data_equation = std::pair<data_equation_extended, equivalence_classes>;

template<typename Substitution>
class AutomatonMatcher final : public Matcher<Substitution>
{
public:
  /// \brief Initialize the automaton matcher with a number of equations.
  AutomatonMatcher(const data_equation_vector& equations);
  virtual ~AutomatonMatcher() {}

  // Matcher interface.

  void match(const data_expression& term) override;

  const data_equation_extended* next(Substitution& matching_sigma) override;

private:

  /// \brief Each state in the pattern match automaton is labelled with a set of variables.
  class pma_state
  {
  public:
    std::set<variable> variables;
    std::vector<std::reference_wrapper<const linear_data_equation>> match_set;
  };

  /// \brief Each transition is either a function symbol, or labelled with the unnamed variable omega.
  using pma_transition = data_expression;

  /// \brief Represents the unnamed variable omega.
  class omega : public data_expression
  {
  public:
    static atermpp::function_symbol& g_function_symbol()
    {
      static atermpp::function_symbol constant("@@omega@@", 0);
      return constant;
    }

    omega()
      : data_expression(atermpp::aterm_appl(g_function_symbol()))
    {}
  };

  /// \returns True iff the given data expression is of type omega.
  static inline bool is_omega(const data_expression& expression)
  {
    return expression.function() == omega::g_function_symbol();
  }

  /// \brief Represents the epsilon, indicating the end of the string. It is annotated with the associated equation that was matched.
  class end_of_string : public data_expression
  {
  public:
    static atermpp::function_symbol& g_function_symbol()
    {
      static atermpp::function_symbol constant("@@annotated_end_of_string@@", 1);
      return constant;
    }

    end_of_string(const data_equation& equation)
      : data_expression(atermpp::aterm_appl(g_function_symbol(), equation))
    {}

    const data_equation& equation() const
    {
      return static_cast<const data_equation&>(static_cast<const atermpp::aterm_appl&>(*this).operator[](0));
    }
  };

  /// \returns True iff the given data expression is of type end_of_string.
  static inline bool is_end_of_string(const data_expression& expression)
  {
    return expression.function() == end_of_string::g_function_symbol();
  }

  /// \returns A fresh state.
  pma_state* add_fresh_state();

  /// \brief Computes the necessary transitions for the automaton that is stored in this class.
  void construct_rec(const PatternSet& patterns, pma_state* from);

  /// \returns The set of function symbols such that P / f != empty for f in Sigma \ V.
  std::set<data_expression> get_head_symbols(const PatternSet& P) const;

  /// \returns The set { omega^arity(f) p | p in P}
  PatternSet prepend_omegas(const PatternSet& L, const application& f) const;

  /// \brief Print the pattern set as a readable format to the info stream.
  void print_pattern_set(const PatternSet& set) const;

  /// \brief A local copy of the equations to keep the references to equations stable.
  mcrl2::utilities::unordered_map<data_equation, linear_data_equation> m_mapping;

  /// \returns True iff the given equivalence classes are consistent w.r.t. the substitution
  bool is_consistent(const equivalence_classes& classes, const Substitution& sigma);

  // The underlying automaton.

  pma_state* m_root_state;

  std::deque<pma_state> m_states;

  mcrl2::utilities::unordered_map<std::pair<pma_state*, pma_transition>, pma_state*> m_transitions;

  data_expression m_omega = omega();

  // Store information about the match.

  Substitution m_matching_sigma;

  std::vector<std::reference_wrapper<const linear_data_equation>>* m_match_set = nullptr;

  std::size_t m_match_index; ///< The index of the equation that should be returned by the call to next.
};

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_DETAIL_AUTOMATON_MATCHER_H
