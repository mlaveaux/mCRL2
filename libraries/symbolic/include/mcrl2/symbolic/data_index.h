// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_SYMBOLIC_DATA_INDEX_H
#define MCRL2_SYMBOLIC_DATA_INDEX_H

#include "mcrl2/core/detail/print_utility.h"
#include "mcrl2/data/data_expression.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/indexed_set.cxx"
#else
  import utilities;
#endif

namespace mcrl2::symbolic
{

/// \brief A bidirectional mapping between data expressions of a given sort and numbers
class data_expression_index : public mcrl2::utilities::indexed_set<data::data_expression>
{
    friend std::ostream& operator<<(std::ostream&, const data_expression_index&);

  protected:
    data::sort_expression m_sort;

  public:
    data_expression_index(const data::sort_expression& sort)
      : m_sort(sort)
    {}

    std::pair<size_type, bool> insert(const key_type& key)
    {
      assert(key.sort() == m_sort);
      return mcrl2::utilities::indexed_set<data::data_expression>::insert(key);
    }
    
    const data::sort_expression& sort() const
    {
      return m_sort;
    }
};

inline
std::ostream& operator<<(std::ostream& out, const data_expression_index& x)
{
  out << "data index with sort = " << x.sort() << " values = " << core::detail::print_list(x);
  return out;
}

} // namespace mcrl2::symbolic

#endif // MCRL2_SYMBOLIC_DATA_INDEX_H
