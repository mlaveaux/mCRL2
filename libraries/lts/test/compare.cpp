// Author(s): Muck van Weerdenburg, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/test/included/unit_test_framework.hpp>

#include "mcrl2/lts/lts_algorithm.h"
#include "mcrl2/lts/lts_aut.h"

using namespace mcrl2::lts;

static lts_aut_t parse_aut(const std::string& s)
{
  std::stringstream is(s);
  lts_aut_t l;
  l.load(is);
  return l;
}

// a.(b+c)
const std::string l1 =
  "des(0,3,4)\n"
  "(0,\"a\",1)\n"
  "(1,\"b\",2)\n"
  "(1,\"c\",3)\n";

// a.b+a.c
const std::string l2 =
  "des(0,4,5)\n"
  "(0,\"a\",1)\n"
  "(0,\"a\",2)\n"
  "(1,\"b\",3)\n"
  "(2,\"c\",4)\n";

// a.b+a.c+ a.(b+c)
const std::string l2a =
  "des(0,7,8)\n"
  "(0,\"a\",1)\n"
  "(0,\"a\",2)\n"
  "(1,\"b\",3)\n"
  "(2,\"c\",4)\n"
  "(0,\"a\",5)\n"
  "(5,\"b\",6)\n"
  "(5,\"c\",7)\n";

// a.tau.(b+c)
const std::string l3 =
  "des(0,4,5)\n"
  "(0,\"a\",1)\n"
  "(1,\"tau\",2)\n"
  "(2,\"b\",3)\n"
  "(2,\"c\",4)\n";

// a.(b+b)
const std::string l4 =
  "des(0,3,4)\n"
  "(0,\"a\",1)\n"
  "(1,\"b\",2)\n"
  "(1,\"b\",3)\n";

// a
const std::string a =
  "des (0,1,2)\n"
  "(0,\"a\",1)\n";

// b
const std::string b =
  "des (0,1,2)\n"
  "(0,\"b\",1)\n";

static inline
bool preorder_compare(const std::string& s1, const std::string& s2, lts_preorder pre)
{
  lts_aut_t t1=parse_aut(s1);
  lts_aut_t t2=parse_aut(s2);
  return compare(t1, t2, pre, false, es_breadth);
}

static inline
bool compare(const std::string& s1, const std::string& s2, lts_equivalence eq, bool counterexample=false)
{
  lts_aut_t t1=parse_aut(s1);
  lts_aut_t t2=parse_aut(s2);
  return compare(t1, t2, eq, counterexample);
}

BOOST_AUTO_TEST_CASE(test_reflexive)
{
  BOOST_CHECK(preorder_compare(l2a,l2a,lts_pre_sim));
  BOOST_CHECK(preorder_compare(l2a,l2a,lts_pre_trace));
  BOOST_CHECK(compare(l2a,l2a,lts_eq_sim));
  BOOST_CHECK(compare(l2a,l2a,lts_eq_trace));
  BOOST_CHECK(preorder_compare(l2a,l2a,lts_pre_ready_sim));
  BOOST_CHECK(compare(l2a,l2a,lts_eq_ready_sim));    
}


BOOST_AUTO_TEST_CASE(test_sim_1_2)
{
  BOOST_CHECK(!preorder_compare(l1,l2,lts_pre_sim));  
  BOOST_CHECK(preorder_compare(l2,l1,lts_pre_sim));
  BOOST_CHECK(!compare(l2,l1,lts_eq_sim));
  BOOST_CHECK(!compare(l1,l2,lts_eq_sim));    
}

BOOST_AUTO_TEST_CASE(test_sim_2_2a)
{
  BOOST_CHECK(preorder_compare(l2,l2a,lts_pre_sim));
  BOOST_CHECK(!preorder_compare(l2a,l2,lts_pre_sim));
  BOOST_CHECK(!compare(l2a,l2,lts_eq_sim));
  BOOST_CHECK(!compare(l2,l2a,lts_eq_sim));    
}

BOOST_AUTO_TEST_CASE(test_sim_1_3)
{
  BOOST_CHECK(!preorder_compare(l1,l3,lts_pre_sim));
  BOOST_CHECK(!preorder_compare(l3,l1,lts_pre_sim));
  BOOST_CHECK(!compare(l3,l1,lts_eq_sim));
  BOOST_CHECK(!compare(l1,l3,lts_eq_sim));    
}


BOOST_AUTO_TEST_CASE(test_sim_1_4)
{
  BOOST_CHECK(!preorder_compare(l1,l4,lts_pre_sim));
  BOOST_CHECK(preorder_compare(l4,l1,lts_pre_sim));
  BOOST_CHECK(!compare(l1,l4,lts_eq_sim));
  BOOST_CHECK(!compare(l4,l1,lts_eq_sim));
}


BOOST_AUTO_TEST_CASE(test_ready_sim_1_2)
{
  BOOST_CHECK(!preorder_compare(l1,l2,lts_pre_ready_sim));
  BOOST_CHECK(!preorder_compare(l2,l1,lts_pre_ready_sim));
  BOOST_CHECK(!compare(l2,l1,lts_eq_ready_sim));
  BOOST_CHECK(!compare(l1,l2,lts_eq_ready_sim));  
}

BOOST_AUTO_TEST_CASE(test_ready_sim_2_2a)
{
  BOOST_CHECK(preorder_compare(l2,l2a,lts_pre_ready_sim));      
  BOOST_CHECK(!preorder_compare(l2a,l2,lts_pre_ready_sim));
  BOOST_CHECK(!compare(l2,l2a,lts_eq_ready_sim));      
  BOOST_CHECK(!compare(l2a,l2,lts_eq_ready_sim));  
}

BOOST_AUTO_TEST_CASE(test_ready_sim_1_3)
{
  BOOST_CHECK(!preorder_compare(l1,l3,lts_pre_ready_sim));
  BOOST_CHECK(!preorder_compare(l3,l1,lts_pre_ready_sim));
  BOOST_CHECK(!compare(l3,l1,lts_eq_ready_sim));  
}

BOOST_AUTO_TEST_CASE(test_ready_sim_1_4)
{
  BOOST_CHECK(!preorder_compare(l1,l4,lts_pre_ready_sim));
  BOOST_CHECK(!preorder_compare(l4,l1,lts_pre_ready_sim));
  BOOST_CHECK(!compare(l1,l4,lts_eq_ready_sim));  
}





BOOST_AUTO_TEST_CASE(test_symmetric_trace_1_2)
{
  BOOST_CHECK(compare(l1,l2,lts_eq_trace));
  BOOST_CHECK(compare(l2,l1,lts_eq_trace));
  BOOST_CHECK(!compare(l2,l1,lts_eq_bisim));
  BOOST_CHECK(!compare(l2,l1,lts_eq_bisim_gv));
  // BOOST_CHECK(!compare(l2,l1,lts_eq_bisim_tb)); Only enable if this algorithm is finished and correct.
  BOOST_CHECK(preorder_compare(l1,l2,lts_pre_trace));
  BOOST_CHECK(preorder_compare(l2,l1,lts_pre_trace));
  BOOST_CHECK(preorder_compare(l1,l2,lts_pre_trace_anti_chain));
  BOOST_CHECK(preorder_compare(l2,l1,lts_pre_trace_anti_chain));
  BOOST_CHECK(preorder_compare(l1,l2,lts_pre_failures_refinement));
  BOOST_CHECK(!preorder_compare(l2,l1,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(l1,l2,lts_pre_failures_divergence_refinement));
  BOOST_CHECK(!preorder_compare(l2,l1,lts_pre_failures_divergence_refinement));
}

BOOST_AUTO_TEST_CASE(test_symmetric_trace_1_2a)
{
  BOOST_CHECK(compare(l2,l2a,lts_eq_trace));
  BOOST_CHECK(compare(l2a,l2,lts_eq_trace));
  BOOST_CHECK(!compare(l2a,l2,lts_eq_bisim));
  BOOST_CHECK(!compare(l2a,l2,lts_eq_bisim_gv));
  // BOOST_CHECK(!compare(l2a,l2,lts_eq_bisim_tb)); Only enable if this algorithm is finished and correct. 
  BOOST_CHECK(preorder_compare(l2,l2a,lts_pre_trace));
  BOOST_CHECK(preorder_compare(l2a,l2,lts_pre_trace));
  BOOST_CHECK(preorder_compare(l2,l2a,lts_pre_trace_anti_chain));
  BOOST_CHECK(preorder_compare(l2a,l2,lts_pre_trace_anti_chain));
  BOOST_CHECK(preorder_compare(l2,l2a,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(l2a,l2,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(l2,l2a,lts_pre_failures_divergence_refinement));
  BOOST_CHECK(preorder_compare(l2a,l2,lts_pre_failures_divergence_refinement));
}

BOOST_AUTO_TEST_CASE(test_symmetric_trace_1_3)
{
  BOOST_CHECK(!compare(l1,l3,lts_eq_trace));
  BOOST_CHECK(!compare(l3,l1,lts_eq_trace));
  BOOST_CHECK(compare(l1,l3,lts_eq_weak_trace));
  BOOST_CHECK(compare(l3,l1,lts_eq_weak_trace));
  BOOST_CHECK(!preorder_compare(l1,l3,lts_pre_trace_anti_chain));
  BOOST_CHECK(!preorder_compare(l3,l1,lts_pre_trace_anti_chain));
  BOOST_CHECK(preorder_compare(l1,l3,lts_pre_weak_trace_anti_chain));
  BOOST_CHECK(preorder_compare(l3,l1,lts_pre_weak_trace_anti_chain));
  BOOST_CHECK(!preorder_compare(l1,l3,lts_pre_failures_refinement));
  BOOST_CHECK(!preorder_compare(l3,l1,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(l1,l3,lts_pre_weak_failures_refinement));
  BOOST_CHECK(preorder_compare(l3,l1,lts_pre_weak_failures_refinement));
  BOOST_CHECK(preorder_compare(l1,l3,lts_pre_failures_divergence_refinement));
  BOOST_CHECK(preorder_compare(l3,l1,lts_pre_failures_divergence_refinement));
}

BOOST_AUTO_TEST_CASE(test_symmetric_trace_1_4)
{
  BOOST_CHECK(!compare(l1,l4,lts_eq_trace));
  BOOST_CHECK(!compare(l4,l1,lts_eq_trace));
}

BOOST_AUTO_TEST_CASE(test_symmetric_weak_trace_2_3)
{
  BOOST_CHECK(compare(l2,l3,lts_eq_weak_trace));
  BOOST_CHECK(compare(l3,l2,lts_eq_weak_trace));
}

BOOST_AUTO_TEST_CASE(test_symmetric_weak_trace_3_4)
{
  BOOST_CHECK(!compare(l4,l3,lts_eq_weak_trace));
  BOOST_CHECK(!compare(l3,l4,lts_eq_weak_trace));
}

// Regression test for bug #1082
BOOST_AUTO_TEST_CASE(test_bisim_a_b)
{
  BOOST_CHECK(!compare(a,b,lts_eq_bisim, true));
}

// a.(b.c1+d1) + a.(b.c2+d2)
  const std::string failures_law_left_hand_side =
   "des(0,8,6)\n"
   "(0,\"a\",1)\n"
   "(1,\"b\",2)\n"
   "(2,\"c1\",3)\n"
   "(1,\"d1\",3)\n"
   "(0,\"a\",4)\n"
   "(4,\"b\",5)\n"
   "(5,\"c2\",3)\n"
   "(4,\"d2\",3)\n";

// a.(b.c1+b.c2+d1) + a.(b.c2+b.c1+d2)
  const std::string failures_law_right_hand_side =
   "des(0,10,6)\n"
   "(0,\"a\",1)\n"
   "(1,\"b\",2)\n"
   "(2,\"c1\",3)\n"
   "(1,\"d1\",3)\n"
   "(1,\"b\",5)\n"
   "(0,\"a\",4)\n"
   "(4,\"b\",5)\n"
   "(4,\"b\",2)\n"
   "(5,\"c2\",3)\n"
   "(4,\"d2\",3)\n";

// The typical law for strong failures equivalence is 
// a.(b.x+u)+a.(b.y+v) = a.(b.x+b.y+u)+a.(b.x+b.y+v).
// The validity of this law is tested here, where the
// variables have been replaced by single actions.
BOOST_AUTO_TEST_CASE(laws_for_failures_equivalence)
{
  BOOST_CHECK(preorder_compare(failures_law_left_hand_side,failures_law_right_hand_side,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(failures_law_right_hand_side,failures_law_left_hand_side,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(failures_law_left_hand_side,failures_law_right_hand_side,lts_pre_trace));
  BOOST_CHECK(preorder_compare(failures_law_right_hand_side,failures_law_left_hand_side,lts_pre_trace));
  BOOST_CHECK(preorder_compare(failures_law_left_hand_side,failures_law_right_hand_side,lts_pre_trace_anti_chain));
  BOOST_CHECK(preorder_compare(failures_law_right_hand_side,failures_law_left_hand_side,lts_pre_trace_anti_chain));
  BOOST_CHECK(preorder_compare(failures_law_left_hand_side,failures_law_right_hand_side,lts_pre_failures_divergence_refinement));
  BOOST_CHECK(preorder_compare(failures_law_right_hand_side,failures_law_left_hand_side,lts_pre_failures_divergence_refinement));
}

// a.b + a.(b+c)
  const std::string ababc =
   "des(0,5,6)\n"
   "(0,\"a\",1)\n"
   "(1,\"b\",2)\n"
   "(0,\"a\",3)\n"
   "(3,\"b\",4)\n"
   "(3,\"c\",5)\n";

// a.(b+c) is typically failure included in a.b+a(b+c), but not vice versa.
BOOST_AUTO_TEST_CASE(failures_inclusion_test)
{
  BOOST_CHECK(preorder_compare(l1,ababc,lts_pre_failures_refinement));
  BOOST_CHECK(!preorder_compare(ababc,l1,lts_pre_failures_refinement));
}

// a.(tau.b + tau.c)
  const std::string a_taub_tauc =
   "des(0,5,6)\n"
   "(0,\"a\",1)\n"
   "(1,\"tau\",2)\n"
   "(2,\"b\",3)\n"
   "(1,\"tau\",4)\n"
   "(4,\"c\",5)\n";

// a.(tau.b+tau.c) and a.b + a.c are weak failure equivalent.
BOOST_AUTO_TEST_CASE(weak_failures_equivalence_test)
{
  BOOST_CHECK(!preorder_compare(l2,a_taub_tauc,lts_pre_failures_refinement));
  BOOST_CHECK(!preorder_compare(a_taub_tauc,l2,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(l2,a_taub_tauc,lts_pre_weak_trace));
  BOOST_CHECK(preorder_compare(a_taub_tauc,l2,lts_pre_weak_trace));
  BOOST_CHECK(preorder_compare(l2,a_taub_tauc,lts_pre_weak_trace_anti_chain));
  BOOST_CHECK(preorder_compare(a_taub_tauc,l2,lts_pre_weak_trace_anti_chain));
  BOOST_CHECK(preorder_compare(l2,a_taub_tauc,lts_pre_weak_failures_refinement));
  BOOST_CHECK(preorder_compare(a_taub_tauc,l2,lts_pre_weak_failures_refinement));
}

// a.b + a.(b+c)
  const std::string abc_div =
   "des(0,4,4)\n"
   "(0,\"a\",1)\n"
   "(1,\"tau\",1)\n"
   "(1,\"b\",2)\n"
   "(1,\"c\",3)\n";

// a.tau*(b+c) is typically failures equivalent to a.(b+c).
// But a.tau*.(b+c) is not failures-divergence included in a.(b+c). The reverse however holds.
BOOST_AUTO_TEST_CASE(failures_divergence_inclusion_test)
{
  BOOST_CHECK(!preorder_compare(l1,abc_div,lts_pre_failures_refinement));
  BOOST_CHECK(!preorder_compare(abc_div,l1,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(abc_div,l1,lts_pre_weak_trace));
  BOOST_CHECK(preorder_compare(l1,abc_div,lts_pre_weak_trace));
  BOOST_CHECK(preorder_compare(abc_div,l1,lts_pre_weak_trace_anti_chain));
  BOOST_CHECK(preorder_compare(l1,abc_div,lts_pre_weak_trace_anti_chain));
  BOOST_CHECK(!preorder_compare(l1,abc_div,lts_pre_weak_failures_refinement));
  BOOST_CHECK(preorder_compare(abc_div,l1,lts_pre_weak_failures_refinement));
  BOOST_CHECK(!preorder_compare(l1,abc_div,lts_pre_failures_divergence_refinement));
  BOOST_CHECK(!preorder_compare(abc_div,l1,lts_pre_failures_divergence_refinement));
}

//  Example by Verum showing an error in the weak failures inclusion check. Did not work with 
//  the -c flag (which does not preprocess using branching bisimulation) but did work after
//  applying branching bisimulation.

  const std::string lts_impl =
   "des (1,8,6)\n"
   "(0,f,1)\n"
   "(1,b,2)\n"
   "(1,a,2)\n"
   "(2,tau,0)\n"
   "(2,tau,5)\n"
   "(3,tau,5)\n"
   "(4,tau,0)\n"
   "(5,t,1)\n";

  const std::string lts_spec =
   "des (1,5,3)\n"
   "(0,f,1)\n"
   "(1,b,2)\n"
   "(1,b,0)\n"
   "(1,a,2)\n"
   "(2,t,1)\n";
  
BOOST_AUTO_TEST_CASE(verum_test)
{
  BOOST_CHECK(!preorder_compare(lts_impl,lts_spec,lts_pre_failures_refinement));
  BOOST_CHECK(!preorder_compare(lts_impl,lts_spec,lts_pre_failures_divergence_refinement));
  BOOST_CHECK(!preorder_compare(lts_impl,lts_spec,lts_pre_weak_failures_refinement));
  BOOST_CHECK(!preorder_compare(lts_impl,lts_spec,lts_pre_weak_trace));
  BOOST_CHECK(!preorder_compare(lts_impl,lts_spec,lts_pre_trace));
  BOOST_CHECK(!preorder_compare(lts_impl,lts_spec,lts_pre_weak_trace_anti_chain));
  BOOST_CHECK(!preorder_compare(lts_impl,lts_spec,lts_pre_trace_anti_chain));

  BOOST_CHECK(!preorder_compare(lts_spec,lts_impl,lts_pre_failures_refinement));
  BOOST_CHECK(preorder_compare(lts_spec,lts_impl,lts_pre_failures_divergence_refinement));
  BOOST_CHECK(preorder_compare(lts_spec,lts_impl,lts_pre_weak_failures_refinement));
  BOOST_CHECK(preorder_compare(lts_spec,lts_impl,lts_pre_weak_trace));
  BOOST_CHECK(!preorder_compare(lts_spec,lts_impl,lts_pre_trace));
  BOOST_CHECK(preorder_compare(lts_spec,lts_impl,lts_pre_weak_trace_anti_chain));
  BOOST_CHECK(!preorder_compare(lts_spec,lts_impl,lts_pre_trace_anti_chain));
}




boost::unit_test::test_suite* init_unit_test_suite(int argc, char* argv[])
{
  return nullptr;
}
