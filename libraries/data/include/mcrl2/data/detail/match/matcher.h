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
#include "mcrl2/data/detail/match/consistency.h"
#include "mcrl2/data/substitutions/mutable_indexed_substitution.h"

namespace mcrl2
{
namespace data
{
namespace detail
{

/// \brief Extend a data equation with additional constructs used after matching.
class extended_data_equation
{
public:
  extended_data_equation(data_equation equation)
    : m_condition(equation.condition()),
      m_righthandside(equation.rhs()),
      m_equation(equation)
  {}

  const data_equation& equation() const { return m_equation; }

  const ConstructionStack& condition_stack() const { return m_condition; }
  const ConstructionStack& rhs_stack() const { return m_righthandside; }

private:
  ConstructionStack m_condition;
  ConstructionStack m_righthandside;
  data_equation m_equation;
};


/// \brief The combination of a linear data equation and a consistency partition.
struct linear_data_equation : public extended_data_equation
{
public:
  linear_data_equation(data_equation equation, consistency_partition partition)
    : extended_data_equation(equation),
      m_partition(partition)
  {}

  const consistency_partition& partition() const { return m_partition; }
private:
  consistency_partition m_partition;
};

/// \brief The interface for matching algorithms.
template<typename Substitution>
class Matcher : public mcrl2::utilities::noncopyable
{
public:

  /// \brief Start matching the given term, use next() to obtain the results per index.
  virtual void match(const data_expression& term) = 0;

  /// \returns The matching equation and adapts matching_sigma accordingly.
  virtual const extended_data_equation* next(Substitution& matching_sigma) = 0;
};

}
}
}


#endif // MCRL2_DATA_DETAIL_MATCHER_H
