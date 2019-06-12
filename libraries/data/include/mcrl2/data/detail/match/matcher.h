// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_MATCHER_H
#define MCRL2_DATA_DETAIL_MATCHER_H

#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/detail/match/construction_stack.h"
#include "mcrl2/data/substitutions/mutable_indexed_substitution.h"

namespace mcrl2
{
namespace data
{
namespace detail
{

using data_equation_extended = std::tuple<data_equation, ConstructionStack, ConstructionStack>;

class Matcher : public mcrl2::utilities::noncopyable
{
public:
  /// \brief The match function defined in the document. However, instead of returning a set of right-hand sides it makes a (arbitrary)
  ///        choice of which right-hand side to return and applies the matching substitution to it. The given term must be in normal form.
  /// \returns A number of equations that match the given term, under the updated (passed through) substitution sigma.
  virtual std::vector<std::reference_wrapper<const data_equation_extended>> match(const data_expression& term, mutable_indexed_substitution<>& matching_sigma) = 0;
};

}
}
}

#endif // MCRL2_DATA_DETAIL_MATCHER_H
