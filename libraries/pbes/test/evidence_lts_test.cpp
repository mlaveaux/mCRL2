// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file evidence_lts_test.cpp
/// \brief Tests for solve_structure_graph_with_evidence_lts: the evidence LTS must contain the
/// same information that solve_structure_graph_with_counter_example writes as an evidence LPS.

#define BOOST_TEST_MODULE evidence_lts_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/lps/linearise.h"
#include "mcrl2/lts/lts_algorithm.h"
#include "mcrl2/lts/lts_equivalence.h"
#include "mcrl2/pbes/evidence_lts.h"
#include "mcrl2/pbes/lps2pbes.h"
#include "mcrl2/pbes/pbesinst_structure_graph.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system;

static lts::lts_aut_t parse_aut(const std::string& s)
{
  std::stringstream is(s);
  lts::lts_aut_t l;
  l.load(is);
  return l;
}

// Runs lps2pbes -c on spec_text/formula_text, instantiates and solves the resulting PBES, and
// writes the evidence directly as an .aut file at filename.
static bool run(const std::string& spec_text, const std::string& formula_text, const std::string& filename)
{
  lps::stochastic_specification stochastic_lpsspec = lps::linearise(spec_text);
  lps::specification lpsspec = lps::remove_stochastic_operators(stochastic_lpsspec);

  bool formula_is_quantitative = false;
  state_formulas::state_formula formula
      = state_formulas::algorithms::parse_state_formula(formula_text, stochastic_lpsspec, formula_is_quantitative);

  bool timed = false;
  bool structured = false;
  bool unoptimized = false;
  bool preprocess_modal_operators = false;
  bool generate_counter_example = true;
  pbes p = lps2pbes(stochastic_lpsspec, formula, timed, structured, unoptimized, preprocess_modal_operators, generate_counter_example);

  pbessolve_options options;
  structure_graph G;
  pbesinst_structure_graph_algorithm algorithm(options, p, G);
  algorithm.run();

  return solve_structure_graph_with_evidence_lts(G, lpsspec, p, algorithm.equation_index(), options.rewrite_strategy, lts::lts_aut, filename);
}

// P has a single state (no process parameters), so every evidence LTS for it collapses onto a
// single self-loop: this is a small but real end-to-end exercise of the summand-index lookup,
// the source/action/target slicing, and the BFS-based state numbering.
static const std::string P = "act a, b;\nproc P = a.P + b.P;\ninit P;\n";

BOOST_AUTO_TEST_CASE(test_witness)
{
  bool result = run(P, "<a>true", "evidence_lts_test_witness.aut");
  BOOST_CHECK(result);

  lts::lts_aut_t evidence;
  evidence.load("evidence_lts_test_witness.aut");
  lts::lts_aut_t expected = parse_aut(R"(des (0,1,1)
(0,"a",0)
)");
  BOOST_CHECK(lts::destructive_compare(evidence, expected, lts::lts_eq_bisim));
}

BOOST_AUTO_TEST_CASE(test_counterexample)
{
  // mu X.<a>X has no base case, so it is false everywhere; the counterexample is the same
  // self-loop as the witness above, this time reached via the Zneg side of the encoding.
  bool result = run(P, "mu X.<a>X", "evidence_lts_test_counterexample.aut");
  BOOST_CHECK(!result);

  lts::lts_aut_t evidence;
  evidence.load("evidence_lts_test_counterexample.aut");
  lts::lts_aut_t expected = parse_aut(R"(des (0,1,1)
(0,"a",0)
)");
  BOOST_CHECK(lts::destructive_compare(evidence, expected, lts::lts_eq_bisim));
}

// A process with an actual data parameter, so the evidence LTS has more than one state and the
// source/target slicing in for_each_evidence_transition is exercised on a non-trivial example.
static const std::string Counter
    = "act inc, reset;\n"
      "proc P(n: Nat) = inc.P(n + 1) + (n > 0) -> reset.P(0);\n"
      "init P(0);\n";

BOOST_AUTO_TEST_CASE(test_multi_state_witness)
{
  bool result = run(Counter, "<inc><reset>true", "evidence_lts_test_multi_state.aut");
  BOOST_CHECK(result);

  lts::lts_aut_t evidence;
  evidence.load("evidence_lts_test_multi_state.aut");
  lts::lts_aut_t expected = parse_aut(R"(des (0,2,2)
(0,"inc",1)
(1,"reset",0)
)");
  BOOST_CHECK(lts::destructive_compare(evidence, expected, lts::lts_eq_bisim));
}
