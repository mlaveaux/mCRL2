// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_NAIVE_CONSISTENCY_H
#define MCRL2_DATA_DETAIL_NAIVE_CONSISTENCY_H

#include "mcrl2/data/detail/match/matcher.h"
#include "mcrl2/data/detail/match/linearise.h"

namespace mcrl2::data::detail
{

/// \brief Uses variables (in sigma) to check consistency.
using variable_partition = std::vector<std::vector<variable>>;

/// \brief The combination of a linear data equation and a consistency partition.
class linear_data_equation : public extended_data_equation
{
public:
  linear_data_equation(data_equation original_equation, data_equation linear_equation, const consistency_partition& partition)
    : extended_data_equation(linear_equation),
      m_original_equation(original_equation)
  {
    // Convert the positions to their variables.
    for (const consistency_class& positions : partition)
    {
      std::vector<variable> variables;
      for (const position& pos : positions)
      {
        variables.push_back(position_variable(pos));
      }

      m_partition.push_back(variables);
    }
  }

  const variable_partition& partition() const { return m_partition; }

  const data_equation& original_equation() const { return m_original_equation; }
private:
  data_equation m_original_equation;
  variable_partition m_partition;
};

/// \brief Check whether the given substitution sigma is consistent w.r.t. the given equivalence classes.
template<typename Substitution>
inline
bool is_consistent(const linear_data_equation& equation, const Substitution& sigma)
{
   // We also need to check consistency of the matched rule.
   for (const auto& consistency_class : equation.partition())
   {
     const auto& assigned = sigma(consistency_class.front());
     for (const variable& variable : consistency_class)
     {
       if (sigma(variable) != assigned)
       {
         return false;
       }
     }
   }

   return true;
}

} // namespace mcrl2::data::detail


#endif // MCRL2_DATA_DETAIL_NAIVE_CONSISTENCY_H
