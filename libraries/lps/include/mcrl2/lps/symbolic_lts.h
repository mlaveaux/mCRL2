// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_LPS_SYMBOLIC_LTS_H
#define MCRL2_LPS_SYMBOLIC_LTS_H

#ifdef MCRL2_ENABLE_SYLVAN

#include "mcrl2/lps/lps_summand_group.h"
#include "mcrl2/lps/multi_action.h"
#include "mcrl2/symbolic/data_index.h"

#include <sylvan_ldd.hpp>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/indexed_set.cxx"
#else
  import utilities;
#endif

namespace mcrl2::lps
{

/// \brief Represents a symbolic LTS encoded by a disjunctive transition relation and a set of states.
class symbolic_lts
{
  using ldd = sylvan::ldds::ldd;

public:
  data::data_specification data_spec;
  data::variable_list process_parameters;

  sylvan::ldds::ldd states;
  sylvan::ldds::ldd initial_state;

  std::vector<symbolic::data_expression_index> data_index;
  utilities::indexed_set<lps::multi_action> action_index;

  std::vector<lps_summand_group> summand_groups;
  
private:
};

} // namespace mcrl2::lps

#endif // MCRL2_ENABLE_SYLVAN

#endif // MCRL2_LPS_SYMBOLIC_LTS_H