// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_INNERMOST_H
#define MCRL2_DATA_DETAIL_INNERMOST_H

#include "mcrl2/data/detail/rewrite.h"
#include "mcrl2/data/detail/match/automaton_matcher.h"
#include "mcrl2/data/detail/match/naive_matcher.h"
#include "mcrl2/utilities/cache_metric.h"
#include "mcrl2/utilities/fixed_size_cache.h"
#include "mcrl2/utilities/unordered_map.h"
#include "mcrl2/utilities/unordered_set.h"

namespace mcrl2
{
namespace data
{
namespace detail
{

/// \brief A left-most innermost rewriter defined in the "Term Rewriting in mCRL2" document.
class InnermostRewriter final : public Rewriter
{
public:
  using substitution_type = Rewriter::substitution_type;

  InnermostRewriter(const data_specification& data_spec, const used_data_equation_selector& selector);

  // Begin of the Rewriter interface.
  data_expression rewrite(const data_expression &term, substitution_type &sigma) override;

  rewrite_strategy getStrategy() override { return innermost; }
  // End of the Rewriter interface.

private:
  /// \brief Applies the substitution sigma to the term, optionally using the construction stack (depending on EnableConstructionStacks).
  template<typename Substitution>
  data_expression apply_substitution(const data_expression& term, Substitution& sigma, const ConstructionStack& stack);

  /// \brief The rewrite function defined in the document. Takes a term t and a substitution sigma and returns the normal
  ///        form of sigma applied to t with respect to the term rewrite system passed in the constructor.
  data_expression rewrite_impl(const data_expression& term, const substitution_type& sigma);

  /// \returns The normal form of a term of the shape lambda x . u.
  data_expression rewrite_abstraction(const abstraction& abstraction, const substitution_type& sigma);

  /// \returns The normal form of a term of the shape h(u_1, ..., u_n).
  data_expression rewrite_application(const application& appl, const substitution_type& sigma);

  /// \returns The result after applying a single rewrite rule (if applicable) to the given expression.
  data_expression rewrite_single(const data_expression& expression);

  /// \brief Prints the various collected performance metrics after a rewrite() call.
  void print_rewrite_metrics();

  /// \returns True if and only if the given term has been marked as being in normal form.
  bool is_normal_form(const data_expression& term) const;

  /// \brief Marks the given term as being in normal form.
  void mark_normal_form(const data_expression& term);

  mcrl2::utilities::unordered_set_large<data_expression> m_normal_forms; ///< Keeps track of terms that are in normal form.

  mcrl2::utilities::fifo_cache<data_expression, data_expression> m_rewrite_cache; ///< Cache the normal forms of certain data expressions.

  mcrl2::data::detail::AutomatonMatcher m_matcher;

  // These are shared data structures to prevent unnecessary reallocations.

  substitution_type m_identity; ///< The identity substitution.

  mutable mutable_indexed_substitution<variable, data_expression> m_local_sigma; ///< A local substitution to prevent reallocations.

  mutable std::vector<data_expression> m_argument_stack; ///< A reused argument stack

  // These members are only used to keep track of underlying metrics.

  std::unordered_map<data_equation, std::size_t, std::hash<atermpp::aterm_appl>> m_application_count; ///< For every rewrite rule, counts the number of times that it was applied.

  mcrl2::utilities::cache_metric m_rewrite_cache_metric;
};

}
}
}

#endif // MCRL2_DATA_DETAIL_INNERMOST_H
