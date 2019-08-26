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

/// \brief The interface for matching algorithms.
template<typename Substitution>
class Matcher : public mcrl2::utilities::noncopyable
{
public:

  /// \brief Start matching the given term, use next() to obtain the results per index.
  virtual void match(const data_expression& term) = 0;

  /// \returns The matching equation and adapts matching_sigma accordingly.
  virtual const data_equation_extended* next(Substitution& matching_sigma) = 0;
};

}
}
}

#endif // MCRL2_DATA_DETAIL_MATCHER_H
