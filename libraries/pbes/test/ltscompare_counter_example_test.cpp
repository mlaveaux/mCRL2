// Author(s): Maurice Laveaux, Jan Martens
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include "mcrl2/data/data_specification.h"
#include "mcrl2/lps/stochastic_specification.h"
#include "mcrl2/process/action_label.h"
#include "mcrl2/process/process_expression.h"
#define BOOST_TEST_MODULE ltscompare_counter_example_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/data/parse.h"
#include "mcrl2/lts/detail/lts_convert.h"
#include "mcrl2/lts/lts_algorithm.h"
#include "mcrl2/lts/lts_equivalence.h"
#include "mcrl2/lts/lts_io.h"
#include "mcrl2/lts/lts_lts.h"
#include "mcrl2/pbes/lts2pbes.h"
#include "mcrl2/pbes/pbesinst_structure_graph.h"
#include "mcrl2/pbes/solve_structure_graph.h"

using namespace mcrl2::lps;
using namespace mcrl2::lts;
using namespace mcrl2::pbes_system;
using namespace mcrl2::state_formulas;

static lts_aut_t parse_aut(const std::string& s)
{
  std::stringstream is(s);
  lts_aut_t l;
  l.load(is);
  return l;
}

inline
void compare(const char* s1, const char* s2, const char* counterexample_name)
{    
  mcrl2::lts::lts_aut_t t1 = parse_aut(s1);
  mcrl2::lts::lts_aut_t t2 = parse_aut(s2);

  compare(t1, t2, lts_equivalence::lts_eq_branching_bisim, true, counterexample_name);

  // Read the counter example file and verify that it is distinguishing
  std::ifstream counter_example(counterexample_name);

  mcrl2::process::action_label_list labels;
  labels.push_front(mcrl2::process::action_label("a", {}));
  labels.push_front(mcrl2::process::action_label("b", {}));
  labels.push_front(mcrl2::process::action_label("tau", {}));

  // Convert aut to lts format for lts2pbes
  lts_lts_t t1_lts;
  lts_lts_t t2_lts;

  lts_convert(t1, t1_lts, {}, labels, {}, true);
  lts_convert(t2, t2_lts, {}, labels, {}, true);

  // Parse a formula and make up some action declaration since that is necessary.  
  stochastic_specification spec({}, t1_lts.action_label_declarations(), {}, {}, {});
  mcrl2::state_formulas::state_formula formula = mcrl2::state_formulas::algorithms::parse_state_formula(counter_example, spec, false);

  lts2pbes_algorithm algorithm1(t1_lts);
  auto t1_pbes = algorithm1.run(formula);

  lts2pbes_algorithm algorithm2(t2_lts);
  auto t2_pbes = algorithm2.run(formula);

  pbessolve_options options;
  structure_graph G1;
  pbesinst_structure_graph_algorithm algorithm_inst1(options, t1_pbes, G1);
  algorithm_inst1.run();

  structure_graph G2;
  pbesinst_structure_graph_algorithm algorithm_inst2(options, t2_pbes, G2);
  algorithm_inst2.run();

  solve_structure_graph_algorithm algorithm_solve(true, false);   
  BOOST_CHECK_NE(algorithm_solve.solve(G1), algorithm_solve.solve(G2));  
}

BOOST_AUTO_TEST_CASE(ltscompare_br_hard_test)
{
    auto br_hard1 = R"""(des (0,9,7)
        (0,"tau",1)
        (6, "tau", 1)
        (1, "tau", 2)
        (2, "tau", 3)
        (1, "tau", 4)
        (2, "tau", 5)
        (0, "tau", 3)
        (4, a, 4)
        (5, b, 5)
        )""";

    
    auto br_hard2 = R"""(des (6,9,7)
        (0,"tau",1)
        (6, "tau", 1)
        (1, "tau", 2)
        (2, "tau", 3)
        (1, "tau", 4)
        (2, "tau", 5)
        (0, "tau", 3)
        (4, a, 4)
        (5, b, 5)
        )""";

    compare(br_hard1, br_hard2, "br_hard_counterexample.mcf");
}

BOOST_AUTO_TEST_CASE(ltscompare_tau_dist_test)
{
    auto tau_dist1 = R"""(des (0,3,2)
        (0,"tau",1)
        (0,"a",1)
        (0,"b",1))""";

    
    auto tau_dist2 = R"""(des (0,4,3)
        (0,"a",2)
        (0,"b",2)
        (0, "tau",1)
        (1,"b",2)
        )""";

    compare(tau_dist1, tau_dist2, "tau_dist_counterexample.mcf");
}