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
#include "mcrl2/data/detail/rewrite/utility.h"
#include "mcrl2/utilities/unordered_set.h"
#include "mcrl2/utilities/unordered_map.h"
#include "mcrl2/utilities/fixed_size_cache.h"

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
  /// \brief The rewrite function defined in the document. Takes a term t and a substitution sigma and returns the normal
  ///        form of sigma applied to t with respect to the term rewrite system passed in the constructor.
  data_expression rewrite_impl(const data_expression& term, const substitution_type& sigma);

  /// \returns The normal form of a function symbol.
  data_expression rewrite_function_symbol(const function_symbol& symbol);

  /// \returns The normal form of a term of the shape lambda x . u.
  data_expression rewrite_abstraction(const abstraction& abstraction, const substitution_type& sigma);

  /// \returns The normal form of a term of the shape h(u_1, ..., u_n).
  data_expression rewrite_application(const application& appl, const substitution_type& sigma);

  /// \returns True if and only if the given term has been marked as being in normal form.
  bool is_normal_form(const data_expression& term) const;

  /// \brief Marks the given term as being in normal form.
  void mark_normal_form(const data_expression& term);

  /// \brief Applies the substitution sigma to the term, optionally using the construction stack (depending on EnableConstructionStacks).
  template<typename Substitution>
  data_expression apply_substitution(const data_expression& term, const Substitution& sigma, const ConstructionStack& stack) const;

  /// \brief The match function defined in the document. However, instead of returning a set of right-hand sides it makes a (arbitrary)
  ///        choice of which right-hand side to return and applies the matching substitution to it. The given term must be in normal form.
  /// \returns A boolean indicated that matching has succeeded (could be replaced by optional).
  bool match(const data_expression& term, data_expression& rhs);

  /// A mapping from head symbols to rewrite rules and their corresponding construction stacks. A unique index is used for each head symbol to achieve
  /// this mapping without an unordered_map for performance reasons.
  std::vector<std::vector<std::tuple<data_equation, ConstructionStack, ConstructionStack>>> m_rewrite_system;

  mcrl2::utilities::unordered_set_large<data_expression> m_normal_forms; ///< Keeps track of terms that are in normal form.

  mcrl2::utilities::fixed_size_cache<data_expression, data_expression> m_rewrite_cache; ///< Cache the normal forms of certain data expressions.

  // These are shared data structures to prevent unnecessary reallocations.

  mutable mcrl2::utilities::unordered_map_large<variable, data_expression> m_local_sigma; ///< A local substitution to prevent reallocations.

  mutable std::vector<data_expression> m_argument_stack; ///< A reused argument stack

  substitution_type m_identity; ///< The identity substitution.
};

}
}
}

#endif // MCRL2_DATA_DETAIL_INNERMOST_H
