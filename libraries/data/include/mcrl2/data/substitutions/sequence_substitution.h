// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_SUBSTITUTIONS_SEQUENCE_SUBSTITUTION_H
#define MCRL2_DATA_SUBSTITUTIONS_SEQUENCE_SUBSTITUTION_H

#include "mcrl2/data/assignment.h"

#include <vector>

namespace mcrl2 {

namespace data {

/// \brief Substitution that maps data variables to data expressions. The substitution is stored as an
/// assignment_list.
class sequence_substitution
{
public:
  using variable_type = variable;
  using expression_type = data_expression;
  using argument_type = variable;

  expression_type& operator[](const variable& var)
  {
    m_assignments.emplace_back(var, expression_type());
    return m_assignments.back().second;
  }

  const expression_type& operator()(const variable& v) const
  {
    for (const auto&[var, expression] : m_assignments)
    {
      if (var == v)
      {
        return expression;
      }
    }
    return v;
  }

  void clear() { m_assignments.clear(); }

private:
  std::vector<std::pair<variable_type, expression_type>> m_assignments;
};

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_SUBSTITUTIONS_SEQUENCE_SUBSTITUTION_H
