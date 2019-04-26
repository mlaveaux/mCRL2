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

  data_expression rewrite(const data_expression &term, substitution_type &sigma) override;

  rewrite_strategy getStrategy() override { return innermost; };

private:
  /// \brief The rewrite function defined in the document. Takes a term, a substitution and a set of free
  ///        variables and returns the rewritten term according to the term rewrite system passed in the constructor.
  data_expression rewrite_impl(const data_expression &term, substitution_type &sigma, std::set<data::variable> free_variables);
};

}
}
}

#endif // MCRL2_DATA_DETAIL_INNERMOST_H
