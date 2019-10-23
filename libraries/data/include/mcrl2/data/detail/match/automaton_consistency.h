// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_AUTOMATON_CONSISTENCY_H
#define MCRL2_DATA_DETAIL_AUTOMATON_CONSISTENCY_H

#include "mcrl2/data/detail/match/matcher.h"
#include "mcrl2/data/detail/match/linearise.h"
#include "mcrl2/utilities/indexed_set.h"

namespace mcrl2::data::detail
{

/// \brief Uses variables (in sigma) to check consistency.
using indexed_partition = std::vector<std::vector<std::size_t>>;

/// \brief The combination of a linear data equation and a consistency partition.
class indexed_linear_data_equation : public extended_data_equation
{
public:
  /// \brief Default constructor.
  indexed_linear_data_equation() {}

  /// \brief Converts the positions of the partition to indexed position partitions using the given positions set.
  indexed_linear_data_equation(data_equation equation, const consistency_partition& partition, mcrl2::utilities::indexed_set<position>& indices)
    : extended_data_equation(equation)
  {
    // Convert the positions to their variables.
    for (const consistency_class& positions : partition)
    {
      std::vector<std::size_t> variables;
      for (const position& pos : positions)
      {
        variables.push_back(indices.insert(pos).first);
      }

      m_partition.push_back(variables);
    }
  }

  const indexed_partition& partition() const { return m_partition; }
private:
  indexed_partition m_partition;
};

/// \brief Check whether the given substitution sigma is consistent w.r.t. the given equivalence classes.
inline
bool is_consistent(const indexed_linear_data_equation& equation, const std::vector<atermpp::unprotected_aterm>& index)
{
   // We also need to check consistency of the matched rule.
   for (const auto& consistency_class : equation.partition())
   {
     const atermpp::unprotected_aterm& assigned = index[consistency_class.front()];
     for (const std::size_t& pos : consistency_class)
     {
       if (index[pos] != assigned)
       {
         return false;
       }
     }
   }

   return true;
}

} // namespace mcrl2::data::detail

#endif // MCRL2_DATA_DETAIL_AUTOMATON_CONSISTENCY_H
