// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file eqelm_test.cpp
/// \brief Add your file description here.

#define BOOST_TEST_MODULE eqelm_test
#include <boost/test/included/unit_test.hpp>
#include "mcrl2/pbes/detail/pbessolve_algorithm.h"
#include "mcrl2/pbes/detail/pbes_property_map.h"

#include "mcrl2/pbes/extended_pbes.h"
#include "mcrl2/pbes/eqelm.h"
#include "mcrl2/pbes/rewriter.h"
#include "mcrl2/pbes/txt2pbes.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system;

// Example provided by Frank Stappers
std::string t1 =
  "pbes nu X(m, n: Nat) =                  \n"
  "        val(n == m) && X(m + 1, n + 1); \n"
  "                                        \n"
  " init X(0, 0);                          \n"
  ;
std::string x1 = "binding_variables = X(m: Nat)";
std::string x2 = "binding_variables = X(n: Nat)";

// Example provided by Tim Willemse.
// The parameters n and m are not equivalent, and thus should not
// be removed.
std::string t2 =
  "pbes nu X(n,m: Nat) =          \n"
  "       forall p: Nat. X(m, p); \n"
  "                               \n"
  "init X(0, 0);                  \n"
  ;
std::string x3 = "binding_variables = X(n,m: Nat)";

/// The test succeeds if any element in expected_results is satisfied.
void test_pbes(const std::string& pbes_spec,
               const std::vector<std::string>& expected_results,
               bool /* compute_conditions */,
               bool /* remove_equations = true */,
               const std::string& msg = "")
{
  using my_pbes_rewriter = simplify_data_rewriter<data::rewriter>;

  pbes p = txt2pbes(pbes_spec);
  pbes q = p;

  // data rewriter
  data::rewriter datar(q.data());

  // pbes rewriter
  my_pbes_rewriter pbesr(datar);

  // constelm algorithm
  pbes_eqelm_algorithm<pbes_expression, data::rewriter, my_pbes_rewriter> algorithm(datar, pbesr);

  // run the algorithm
  algorithm.run(q);
  BOOST_CHECK(q.is_well_typed());
  if (!q.is_well_typed())
  {
    std::cerr << pbes_system::pp(q) << std::endl;
  }

  pbes_system::detail::pbes_property_map info1(q);
  std::vector<std::string> failed;

  for (const std::string& expected_result : expected_results)
  {
      pbes_system::detail::pbes_property_map info2(expected_result);
      std::string diff = info1.compare(info2);
      if (!diff.empty())
      {
        std::stringstream result;
        result << "--- expected result" << std::endl;
        result << expected_result << std::endl;
        result << "--- found result" << std::endl;
        result << info1.to_string() << std::endl;
        result << "--- differences" << std::endl;
        result << diff << std::endl;
        failed.push_back(result.str());
      }
  }

  if (failed.size() == expected_results.size())
  {
      std::cerr << "\n------ FAILED TEST ------ " << msg << std::endl;
      std::cerr << "None of the expected results match: \n\n";

      for (const std::string& result : failed)
      {
          std::cerr << result;
      }
  }

  BOOST_CHECK(failed.size() < expected_results.size());
}

void test_eqelm(const std::string& pbes_spec, const bool expected_outcome)
{
  using my_pbes_rewriter = simplify_data_rewriter<data::rewriter>;

  pbes p = txt2pbes(pbes_spec);
  pbes q = p;

  data::rewriter datar(q.data());
  my_pbes_rewriter pbesr(datar);
  pbes_eqelm_algorithm<pbes_expression, data::rewriter, my_pbes_rewriter> algorithm(datar, pbesr);
  algorithm.run(q);


  bool solution_p = pbes_system::detail::pbessolve(extended_pbes {
    .source_pbes = p,
  });
  bool solution_q = pbes_system::detail::pbessolve(extended_pbes {
    .source_pbes = q,
  });
  BOOST_CHECK(solution_p == expected_outcome);
  BOOST_CHECK(solution_q == expected_outcome);
}

std::string random1 =
  "pbes mu X0(b: Bool) =                            \n"
  "       (X1(0, 0) || false) && true;              \n"
  "     mu X1(n,m: Nat) =                           \n"
  "       (X1(m + 1, 0) || false) || val(!(n < 2)); \n"
  "                                                 \n"
  "init X0(true);                                   \n"
  ;

BOOST_AUTO_TEST_CASE(eqelm_test1)
{
  bool compute_conditions = false;
  bool remove_equations = true;
  test_pbes(t1, {x1, x2}, compute_conditions, remove_equations, "test 1");
  test_pbes(t2, {x3}, compute_conditions, remove_equations, "test 2");
  test_eqelm(random1, false);
}
