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
  /// \brief Default constructor.
  extended_data_equation() {}

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

/// \brief The result of the matching procedure, nullptr indicates no match.
/// \details The substitution remains valid to the next call to match.
template<typename Substitution>
struct matching_result
{
  const extended_data_equation* equation;
  Substitution& matching_sigma;
};

}
}
}


#endif // MCRL2_DATA_DETAIL_MATCHER_H
