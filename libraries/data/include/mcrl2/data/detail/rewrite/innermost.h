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

  rewrite_strategy getStrategy() override { return innermost; };
  // End of the Rewriter interface.

private:
  /// \brief The rewrite function defined in the document. Takes a term t, a substitution sigma and a set of free
  ///        variables V.
  /// \returns The rewritten term according to the term rewrite system passed in the constructor.
  data_expression rewrite_impl(const data_expression& term, substitution_type& sigma);

  /// \brief Rewrites a term of the form lambda x . u to normal form.
  data_expression rewrite_abstraction(const abstraction& abstraction, substitution_type& sigma);

  /// \brief Rewrites a term of the form h(u_1, ..., u_n) to normal form.
  data_expression rewrite_application(const application& appl, substitution_type& sigma);

  /// \brief The match function defined in the document. However, instead of returning a set of right-hand sides it makes a (arbitrary)
  ///        choice of which right-hand side to return and applied the matching substitution to it. The given term must be in normal form.
  /// \returns A boolean indicated that matching has succeeded (could be replaced by optional).
  bool match(const data_expression& term, data_expression& rhs);

  mcrl2::utilities::unordered_map<data_expression, std::vector<data_equation>> m_rewrite_system; ///< A mapping from function symbols to rewrite rules.
};

}
}
}

#endif // MCRL2_DATA_DETAIL_INNERMOST_H
