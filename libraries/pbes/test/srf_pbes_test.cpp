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

// A bare PVI in a disjunctive equation produces summand ([], true, PVI).
BOOST_AUTO_TEST_CASE(test_srf_or_pvi)
{
  const std::string pbes_text =
      "pbes nu X(b: Bool) =\n"
      "  X(b);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(!main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 1u);

  const auto& summand = main_eqn.summands()[0];
  BOOST_CHECK(summand.parameters().empty());
  BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
  BOOST_CHECK_EQUAL(summand.variable().name(), main_eqn.variable().name());
}

// A disjunction of two PVIs produces two summands, each with condition=true.
BOOST_AUTO_TEST_CASE(test_srf_or_two_pvis)
{
  const std::string pbes_text =
      "pbes nu X(b: Bool) =\n"
      "  X(true) || X(false);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(!main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 2u);

  for (const auto& summand : main_eqn.summands())
  {
    BOOST_CHECK(summand.parameters().empty());
    BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
    BOOST_CHECK_EQUAL(summand.variable().name(), main_eqn.variable().name());
  }
}

// An exists with a non-simple (PVI) body in a disjunctive equation recurses into the body
// and prepends the bound variable as the summand parameter.
BOOST_AUTO_TEST_CASE(test_srf_or_exists_pvi_body)
{
  const std::string pbes_text =
      "pbes nu X(b: Bool) =\n"
      "  exists n: Nat. X(n == 0);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(!main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 1u);

  const auto& summand = main_eqn.summands()[0];
  BOOST_REQUIRE_EQUAL(summand.parameters().size(), 1u);
  BOOST_CHECK_EQUAL(data::pp(summand.parameters().front().sort()), "Nat");
  BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
  BOOST_CHECK_EQUAL(summand.variable().name(), main_eqn.variable().name());
}

// A forall with a non-simple (PVI) body in a conjunctive equation recurses into the body
// and prepends the bound variable as the summand parameter.
BOOST_AUTO_TEST_CASE(test_srf_and_forall_pvi_body)
{
  const std::string pbes_text =
      "pbes mu X(b: Bool) =\n"
      "  forall n: Nat. X(n == 0);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 1u);

  const auto& summand = main_eqn.summands()[0];
  BOOST_REQUIRE_EQUAL(summand.parameters().size(), 1u);
  BOOST_CHECK_EQUAL(data::pp(summand.parameters().front().sort()), "Nat");
  BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
  BOOST_CHECK_EQUAL(summand.variable().name(), main_eqn.variable().name());
}

// A conjunction of two PVIs is classified as conjunctive; both PVIs become separate
// summands ([], true, PVI) via the srf_and PVI handler.
BOOST_AUTO_TEST_CASE(test_srf_and_two_pvis)
{
  const std::string pbes_text =
      "pbes mu X(b: Bool) =\n"
      "  X(true) && X(false);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 2u);

  for (const auto& summand : main_eqn.summands())
  {
    BOOST_CHECK(summand.parameters().empty());
    BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
    BOOST_CHECK_EQUAL(summand.variable().name(), main_eqn.variable().name());
  }
}

// A forall with a non-simple body in a disjunctive equation cannot be inlined; a fresh
// equation is introduced and the main equation gets a summand pointing to it.
BOOST_AUTO_TEST_CASE(test_srf_or_forall_non_simple_body)
{
  const std::string pbes_text =
      "pbes nu X(b: Bool) =\n"
      "  X(b) || forall n: Nat. X(n == 0);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  // Fresh equation for the forall body adds one extra equation.
  BOOST_REQUIRE_EQUAL(equations.size(), 4u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(!main_eqn.is_conjunctive());
  // One summand for X(b), one for the fresh equation.
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 2u);

  // One summand targets the fresh equation (condition=true, different name from X).
  bool found_fresh = false;
  for (const auto& summand : main_eqn.summands())
  {
    if (summand.variable().name() != main_eqn.variable().name())
    {
      found_fresh = true;
      BOOST_CHECK(summand.parameters().empty());
      BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
    }
  }
  BOOST_CHECK(found_fresh);
}

// An exists with a non-simple body in a conjunctive equation cannot be inlined; a fresh
// equation is introduced and the main equation gets a summand pointing to it.
BOOST_AUTO_TEST_CASE(test_srf_and_exists_non_simple_body)
{
  const std::string pbes_text =
      "pbes mu X(b: Bool) =\n"
      "  X(b) && exists n: Nat. X(n == 0);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  // Fresh equation for the exists body adds one extra equation.
  BOOST_REQUIRE_EQUAL(equations.size(), 4u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(main_eqn.is_conjunctive());
  // One summand for X(b), one for the fresh equation.
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 2u);

  bool found_fresh = false;
  for (const auto& summand : main_eqn.summands())
  {
    if (summand.variable().name() != main_eqn.variable().name())
    {
      found_fresh = true;
      BOOST_CHECK(summand.parameters().empty());
      BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
    }
  }
  BOOST_CHECK(found_fresh);
}

// val(b) && X(b) is classified as disjunctive (simple && PVI does not make an and_ conjunctive).
// srf_or merges the simple left side into the summand condition rather than introducing a
// fresh equation, producing summand ([], b && true, X(b)).
BOOST_AUTO_TEST_CASE(test_srf_or_and_simple_left)
{
  const std::string pbes_text =
      "pbes nu X(b: Bool) =\n"
      "  val(b) && X(b);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(!main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 1u);

  const auto& summand = main_eqn.summands()[0];
  BOOST_CHECK(summand.parameters().empty());
  // Condition is val(b) merged into true: data::and_(b, true).
  BOOST_CHECK(data::sort_bool::is_and_application(summand.condition()));
  BOOST_CHECK_EQUAL(summand.variable().name(), main_eqn.variable().name());
}

// val(b) || X(b) is classified as conjunctive (simple || PVI pattern).
// srf_and merges not(val(b)) into the summand condition rather than introducing a
// fresh equation, producing summand ([], !b && true, X(b)).
BOOST_AUTO_TEST_CASE(test_srf_and_or_simple_left)
{
  const std::string pbes_text =
      "pbes mu X(b: Bool) =\n"
      "  val(b) || X(b);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 1u);

  const auto& summand = main_eqn.summands()[0];
  BOOST_CHECK(summand.parameters().empty());
  // Condition is not(val(b)) merged into true: data::and_(!b, true).
  BOOST_CHECK(data::sort_bool::is_and_application(summand.condition()));
  BOOST_CHECK_EQUAL(summand.variable().name(), main_eqn.variable().name());
}

// X(b) || (X(true) && X(false)): the and_ sub-expression in a disjunctive context has
// neither side simple, so a fresh equation is introduced for it.
BOOST_AUTO_TEST_CASE(test_srf_or_and_both_non_simple)
{
  const std::string pbes_text =
      "pbes nu X(b: Bool) =\n"
      "  X(b) || (X(true) && X(false));\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  // Original X + fresh equation for the and_ + X_false + X_true.
  BOOST_REQUIRE_EQUAL(equations.size(), 4u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(!main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 2u);

  // One summand targets the fresh equation (not X, X_false, or X_true).
  const auto& x_false_name = equations[equations.size() - 2].variable().name();
  const auto& x_true_name  = equations[equations.size() - 1].variable().name();
  bool found_fresh = false;
  for (const auto& summand : main_eqn.summands())
  {
    const auto& name = summand.variable().name();
    if (name != main_eqn.variable().name() && name != x_false_name && name != x_true_name)
    {
      found_fresh = true;
      BOOST_CHECK(summand.parameters().empty());
      BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
    }
  }
  BOOST_CHECK(found_fresh);
}

// X(b) && (X(true) || X(false)): the or_ sub-expression in a conjunctive context has
// neither side simple, so a fresh equation is introduced for it.
BOOST_AUTO_TEST_CASE(test_srf_and_or_both_non_simple)
{
  const std::string pbes_text =
      "pbes mu X(b: Bool) =\n"
      "  X(b) && (X(true) || X(false));\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes s = pbes2srf(p);

  const auto& equations = s.equations();
  // Original X + fresh equation for the or_ + X_false + X_true.
  BOOST_REQUIRE_EQUAL(equations.size(), 4u);

  const auto& main_eqn = equations[0];
  BOOST_CHECK(main_eqn.is_conjunctive());
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 2u);

  const auto& x_false_name = equations[equations.size() - 2].variable().name();
  const auto& x_true_name  = equations[equations.size() - 1].variable().name();
  bool found_fresh = false;
  for (const auto& summand : main_eqn.summands())
  {
    const auto& name = summand.variable().name();
    if (name != main_eqn.variable().name() && name != x_false_name && name != x_true_name)
    {
      found_fresh = true;
      BOOST_CHECK(summand.parameters().empty());
      BOOST_CHECK_EQUAL(summand.condition(), data::sort_bool::true_());
    }
  }
  BOOST_CHECK(found_fresh);
}

// In pre-SRF mode (allow_ce=true), or_ in a disjunctive context with a simple clause:
// the simple clause gets its own summand (→ X_true), while non-simple clauses each
// receive not(simple_condition) merged into their condition (pbes-level and_).
BOOST_AUTO_TEST_CASE(test_pre_srf_or_simple_clause)
{
  const std::string pbes_text =
      "pbes nu X(b: Bool) =\n"
      "  val(b) || X(true) || X(false);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes_with_ce s = pbes2pre_srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& x_true_eqn = equations[equations.size() - 1];
  BOOST_REQUIRE(x_true_eqn.symbol() == fixpoint_symbol::nu());
  const auto& x_true_name = x_true_eqn.variable().name();

  const auto& main_eqn = equations[0];
  BOOST_CHECK(!main_eqn.is_conjunctive());
  // Three summands: val(b) → X_true, X(true) and X(false) with not(val(b)) merged in.
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 3u);

  bool found_simple_summand = false;
  bool found_merged_summand = false;
  for (const auto& summand : main_eqn.summands())
  {
    if (summand.variable().name() == x_true_name)
    {
      found_simple_summand = true;
      // The simple clause val(b) becomes the condition directly.
      BOOST_CHECK(data::is_data_expression(summand.condition()));
    }
    else
    {
      found_merged_summand = true;
      // Non-simple clauses get not(val(b)) merged in: pbes-level and_.
      BOOST_CHECK(pbes_system::is_and(summand.condition()));
    }
  }
  BOOST_CHECK(found_simple_summand);
  BOOST_CHECK(found_merged_summand);
}

// In pre-SRF mode (allow_ce=true), and_ in a conjunctive context with a simple clause:
// the simple clause gets its own summand (→ X_false via not), while non-simple clauses
// receive the simple condition merged in (pbes-level and_) to strengthen their guard.
BOOST_AUTO_TEST_CASE(test_pre_srf_and_simple_clause)
{
  const std::string pbes_text =
      "pbes mu X(b: Bool) =\n"
      "  val(b) && X(true) && X(false);\n"
      "init X(true);\n";

  pbes p = txt2pbes(pbes_text);
  normalize(p);
  srf_pbes_with_ce s = pbes2pre_srf(p);

  const auto& equations = s.equations();
  BOOST_REQUIRE_EQUAL(equations.size(), 3u);

  const auto& x_false_eqn = equations[equations.size() - 2];
  BOOST_REQUIRE(x_false_eqn.symbol() == fixpoint_symbol::mu());
  const auto& x_false_name = x_false_eqn.variable().name();

  const auto& main_eqn = equations[0];
  BOOST_CHECK(main_eqn.is_conjunctive());
  // Three summands: val(b) → X_false via not, X(true) and X(false) with val(b) merged in.
  BOOST_REQUIRE_EQUAL(main_eqn.summands().size(), 3u);

  bool found_simple_summand = false;
  bool found_merged_summand = false;
  for (const auto& summand : main_eqn.summands())
  {
    if (summand.variable().name() == x_false_name)
    {
      found_simple_summand = true;
      // The simple clause val(b) becomes not(val(b)) targeting X_false.
      BOOST_CHECK(pbes_system::is_not(summand.condition()));
    }
    else
    {
      found_merged_summand = true;
      // Non-simple clauses get val(b) merged in as a guard: pbes-level and_.
      BOOST_CHECK(pbes_system::is_and(summand.condition()));
    }
  }
  BOOST_CHECK(found_simple_summand);
  BOOST_CHECK(found_merged_summand);
}

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
