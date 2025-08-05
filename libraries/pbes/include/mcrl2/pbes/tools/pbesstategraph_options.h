// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/tools/pbesstategraph_options.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_TOOLS_PBESSTATEGRAPH_OPTIONS_H
#define MCRL2_PBES_TOOLS_PBESSTATEGRAPH_OPTIONS_H

#include "mcrl2/data/rewrite_strategy.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/execution_timer.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system {

struct pbesstategraph_options
{
  data::rewrite_strategy rewrite_strategy = data::jitty;
  bool simplify = true;
  bool use_global_variant = false;
  bool print_influence_graph = false;
  bool cache_marking_updates = false;
  int marking_algorithm = 0;
  bool use_alternative_lcfp_criterion = false;
  bool use_alternative_gcfp_relation = false;
  bool use_alternative_gcfp_consistency = false;
  utilities::execution_timer* timer = nullptr; // if it is non-zero, it will be used to display timing information

  bool timing_enabled() const
  {
    return timer != nullptr;
  }
};

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_TOOLS_PBESSTATEGRAPH_OPTIONS_H
