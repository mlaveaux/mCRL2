// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/detail/split_finite_variables.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_DETAIL_SPLIT_FINITE_VARIABLES_H
#define MCRL2_DATA_DETAIL_SPLIT_FINITE_VARIABLES_H

#include "mcrl2/data/data_specification.h"

namespace mcrl2::data::detail
{

/// \brief Splits the variables in a subset with finite sort, with infinite sort and those that are not used according to a filter.
/// \param variables A sequence of data variables.
/// \param data A data specification.
/// \param finite_variables A sequence of data variables.
/// \param infinite_variables A sequence of data variables.
/// \param unused A sequence of data variables that do not pass the filter.
/// \param filter A predicate that indicates whether a variable is used.
template<typename FILTER>
inline
void split_finite_variables(data::variable_list variables, 
                            const data::data_specification& data, 
                            data::variable_list& finite_variables, 
                            data::variable_list& infinite_variables,
                            data::variable_list& unused,
                            const FILTER filter)
{
  for (data::variable_list::const_reverse_iterator i = variables.rbegin(); i!=variables.rend(); i++)
  {
    if (!filter(*i))
    {
      unused.push_front(*i);
    }
    else if (data.is_certainly_finite(i->sort()))
    {
      finite_variables.push_front(*i);
    }
    else
    {
      infinite_variables.push_front(*i);
    }
  }
}

/// \brief Splits a list of variables in those that are enumerable, non enumerable and not used, according to a filter.
/// \param variables A sequence of data variables.
/// \param data A data specification.
/// \param rewr A data rewriter to determine which elements are enumerable.
/// \param enumerable_variables A sequence of data variables.
/// \param non_enumerable_variables A sequence of data variables.
/// \param unused A sequence of data variables that do not pass the filter.
/// \param filter A predicate that indicates whether a variable is used.
template <typename Rewriter, typename FILTER>
inline
void split_enumerable_variables(data::variable_list variables, 
                                const data::data_specification& data, 
                                const Rewriter& rewr, 
                                data::variable_list& enumerable_variables, 
                                data::variable_list& non_enumerable_variables,
                                data::variable_list& unused,
                                const FILTER filter)
{
  for (data::variable_list::const_reverse_iterator i = variables.rbegin(); i!=variables.rend(); i++)
  {
    if (!filter(*i))
    {
      unused.push_front(*i);
    }
    else if (is_enumerable(data, rewr, i->sort()))
    {
      enumerable_variables.push_front(*i);
    }
    else
    {
      non_enumerable_variables.push_front(*i);
    }
  }
}

} // namespace mcrl2::data::detail

#endif // MCRL2_DATA_DETAIL_SPLIT_FINITE_VARIABLES_H
