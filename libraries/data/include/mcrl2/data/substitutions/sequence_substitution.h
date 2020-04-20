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

#include <queue>

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

  sequence_substitution(std::vector<std::pair<variable_type, expression_type>>& assignments)
    : m_assignments(assignments)
  {
    offset = m_assignments.size();
  }

  ~sequence_substitution()
  {
    clear();
  }

  expression_type& operator[](const variable& var)
  {
    m_assignments.emplace_back(var, var);
    ++size;
    return m_assignments.back().second;
  }

  const expression_type& operator()(const variable& v) const
  {
    for (std::size_t i = offset; i < offset + size; ++i)
    {
      const auto& [var, expression] = m_assignments[i];

      if (var == v)
      {
        return expression;
      }
    }
    return v;
  }

  void clear() 
  {
    assert(offset + size <= m_assignments.size());
    if (size != 0)
    {
      m_assignments.erase(m_assignments.begin() + offset);
    }
  }

private:
  std::size_t offset = 0;
  std::size_t size = 0;

  std::vector<std::pair<variable_type, expression_type>>& m_assignments;
};

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_SUBSTITUTIONS_SEQUENCE_SUBSTITUTION_H
