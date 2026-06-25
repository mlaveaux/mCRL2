// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file srf_pbes_test.cpp
/// \brief Tests for the transformation of a PBES to standard recursive form.

#define BOOST_TEST_MODULE srf_pbes_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/standard_utility.h"
#include "mcrl2/pbes/normalize.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/txt2pbes.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system;

// A simple exists(W, f) inside a conjunctive context is caught by the generic
// apply(pbes_expression) handler in srf_and_traverser, which produces summand
// ([n:Nat], not(exists k. k<n), X_false) — the forall variable becomes the
// summand parameter and the existential is negated into the condition.
// See also: https://github.com/mCRL2org/mCRL2/issues/1913
BOOST_AUTO_TEST_CASE(test_srf_and_exists_simple_body)
{
  const std::string pbes_text =
      "pbes mu X(b: Bool) =\n"
      "  forall n: Nat. (exists k: Nat. val(k < n)) && X(b);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& x_false_eqn = equations[equations.size() - 2];
  const auto& x_true_eqn  = equations[equations.size() - 1];
  BOOST_REQUIRE(x_false_eqn.symbol() == fixpoint_symbol::mu());
  BOOST_REQUIRE(x_true_eqn.symbol()  == fixpoint_symbol::nu());

  const auto& x_false_name = x_false_eqn.variable().name();
  const auto& x_true_name  = x_true_eqn.variable().name();

  const auto& main_eqn = equations[0];
  BOOST_CHECK(main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 2u);

  bool found_x_false_summand = false;
  for (const auto& summand : main_eqn.summands())
  {
    if (summand.variable().name() == x_false_name)
    {
      found_x_false_summand = true;
      BOOST_CHECK_EQUAL(summand.parameters().size(), 1u);
      BOOST_CHECK(data::sort_bool::is_not_application(summand.condition()));
      BOOST_CHECK(summand.variable().name() != x_true_name);
    }
  }
  BOOST_CHECK(found_x_false_summand);
}

// A simple forall(W, f) inside a disjunctive context is caught by the generic
// apply(pbes_expression) handler in srf_or_traverser, which produces summand
// ([], data::forall([n], n<5), X_true) — the universal is preserved as a
// data-level condition with no summand parameters.
// See also: https://github.com/mCRL2org/mCRL2/issues/1913
BOOST_AUTO_TEST_CASE(test_srf_or_forall_simple_body)
{
  const std::string pbes_text =
      "pbes nu X(b: Bool) =\n"
      "  (X(b) || (forall n: Nat. val(n < 5))) || X(b);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& x_true_eqn = equations[equations.size() - 1];
  BOOST_REQUIRE(x_true_eqn.symbol() == fixpoint_symbol::nu());
  const auto& x_true_name = x_true_eqn.variable().name();

  const auto& main_eqn = equations[0];
  BOOST_CHECK(!main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 3u);

  bool found_forall_summand = false;
  for (const auto& summand : main_eqn.summands())
  {
    if (data::is_forall(summand.condition()))
    {
      found_forall_summand = true;
      BOOST_CHECK(summand.parameters().empty());
      BOOST_CHECK_EQUAL(summand.variable().name(), x_true_name);
    }
  }
  BOOST_CHECK(found_forall_summand);
}
