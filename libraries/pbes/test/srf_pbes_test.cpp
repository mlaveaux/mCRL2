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
/// In particular check that simplify_srf_pbes correctly removes summands
/// that have become trivial, and that pre_srf2srfpbes sufficiently simplifies
/// the result.

#define BOOST_TEST_MODULE srf_pbes_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/standard_utility.h"
#include "mcrl2/pbes/pbes_gauss_elimination.h"
#include "mcrl2/pbes/normalize.h"
#include "mcrl2/pbes/resolve_name_clashes.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/txt2pbes.h"
#include "mcrl2/pbes/unify_parameters.h"

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
      // Non-simple clauses get not(val(b)) merged into their guard.
      BOOST_CHECK(pbes_system::is_not(summand.condition()));
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
      // Non-simple clauses get val(b) merged in as a guard.
      BOOST_CHECK(data::is_data_expression(summand.condition()));
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

// Converts p to pre-SRF (with counter example information allowed) and then to SRF,
// which is the transformation under test.
srf_pbes to_srf_via_ce(const pbes& p)
{
  auto pre_srf = detail::pbes2pre_srf<true>(p, true);
  return pre_srf2srfpbes(pre_srf);
}

const srf_equation& find_equation(const srf_pbes& p, const std::string& name)
{
  for (const auto& eqn: p.equations())
  {
    if (std::string(eqn.variable().name()) == name)
    {
      return eqn;
    }
  }
  throw mcrl2::runtime_error("equation " + name + " not found");
}

template <bool allow_ce>
const detail::pre_srf_equation<allow_ce>& find_equation(const detail::pre_srf_pbes<allow_ce>& p, const std::string& name)
{
  for (const auto& eqn: p.equations())
  {
    if (std::string(eqn.variable().name()) == name)
    {
      return eqn;
    }
  }
  throw mcrl2::runtime_error("equation " + name + " not found");
}

// The name of the "always true"/"always false" recursion variables introduced by
// pbes2pre_srf are the last two equations, in that order (false, then true).
core::identifier_string X_false_name(const srf_pbes& p) { return p.equations()[p.equations().size() - 2].variable().name(); }
core::identifier_string X_true_name(const srf_pbes& p) { return p.equations()[p.equations().size() - 1].variable().name(); }

// A propositional variable instantiation with the given name and no parameters, to use
// as a summand target below.
propositional_variable_instantiation var(const std::string& name) { return propositional_variable_instantiation(name); }

// Builds a minimal, well-formed srf_pbes out of the given equations, for directly unit
// testing simplify_srf_pbes.
// Appends the mandatory closing X_false/X_true equations that every
// srf_pbes produced by pbes2pre_srf ends with, each with the canonical single
// self-referencing summand.
// Note that for the tests here, the initial state is irrelevant, so we just use X_true.
srf_pbes make_srf_pbes(std::vector<srf_equation> equations)
{
  srf_summand false_summand(data::variable_list(), true_(), var("X_false"));
  srf_summand true_summand(data::variable_list(), true_(), var("X_true"));
  equations.emplace_back(fixpoint_symbol::mu(), propositional_variable("X_false", data::variable_list()), std::vector<srf_summand>{false_summand}, false);
  equations.emplace_back(fixpoint_symbol::nu(), propositional_variable("X_true", data::variable_list()), std::vector<srf_summand>{true_summand}, false);
  return srf_pbes(data::data_specification(), equations, var("X_true"));
}

// Unit tests for simplify_srf_pbes.
// Given an already-constructed srf_pbes,
// remove summands whose condition simplifies to true or false.

BOOST_AUTO_TEST_CASE(test_simplify_srf_pbes_drops_dead_disjunct)
{
  // exists <> . false && Y() never holds and must be dropped; the real summand
  // exists <> . true && Z() survives unchanged.
  srf_summand dead(data::variable_list(), false_(), var("Y"));
  srf_summand real(data::variable_list(), true_(), var("Z"));
  srf_equation eqX(fixpoint_symbol::mu(), propositional_variable("X", data::variable_list()), {dead, real}, false);

  srf_pbes result = simplify_srf_pbes(make_srf_pbes({eqX}));

  const srf_equation& X = find_equation(result, "X");
  BOOST_CHECK_EQUAL(X.summands().size(), 1u);
  BOOST_CHECK_EQUAL(X.summands().front().variable().name(), core::identifier_string("Z"));
}

BOOST_AUTO_TEST_CASE(test_simplify_srf_pbes_drops_vacuous_conjunct)
{
  // to_pbes() builds forall params . condition => X (see distribute_not_over_and).
  // A stored condition of false_() makes the guard "false => Y()" vacuously true,
  // regardless of Y, so this conjunct contributes nothing and must be dropped. A stored
  // condition of true_() makes the guard "true => Z()" reduce to the genuine,
  // unconditional dependency Z(), which survives unchanged.
  srf_summand vacuous(data::variable_list(), false_(), var("Y"));
  srf_summand real(data::variable_list(), true_(), var("Z"));
  srf_equation eqX(fixpoint_symbol::mu(), propositional_variable("X", data::variable_list()), {vacuous, real}, true);

  srf_pbes result = simplify_srf_pbes(make_srf_pbes({eqX}));

  const srf_equation& X = find_equation(result, "X");
  BOOST_CHECK_EQUAL(X.summands().size(), 1u);
  BOOST_CHECK_EQUAL(X.summands().front().variable().name(), core::identifier_string("Z"));
}

BOOST_AUTO_TEST_CASE(test_simplify_srf_pbes_disjunctive_equation_collapses_to_true)
{
  // A summand that is unconditionally true (condition true_(), target X_true) makes the
  // whole disjunctive equation true, discarding the other, now irrelevant, summand.
  srf_summand real(data::variable_list(), true_(), var("Y"));
  srf_summand trivial(data::variable_list(), true_(), var("X_true"));
  srf_equation eqX(fixpoint_symbol::mu(), propositional_variable("X", data::variable_list()), {real, trivial}, false);

  srf_pbes result = simplify_srf_pbes(make_srf_pbes({eqX}));

  const srf_equation& X = find_equation(result, "X");
  BOOST_CHECK_EQUAL(X.summands().size(), 1u);
  BOOST_CHECK_EQUAL(X.summands().front().variable().name(), X_true_name(result));
}

BOOST_AUTO_TEST_CASE(test_simplify_srf_pbes_conjunctive_equation_collapses_to_false)
{
  // Dually, a summand whose stored condition true_() negates to false_() (see
  // to_pbes()), targeting X_false, makes the whole conjunctive equation false,
  // discarding the other, now irrelevant, summand.
  srf_summand real(data::variable_list(), true_(), var("Y"));
  srf_summand trivial(data::variable_list(), true_(), var("X_false"));
  srf_equation eqX(fixpoint_symbol::nu(), propositional_variable("X", data::variable_list()), {real, trivial}, true);

  srf_pbes result = simplify_srf_pbes(make_srf_pbes({eqX}));

  const srf_equation& X = find_equation(result, "X");
  BOOST_CHECK_EQUAL(X.summands().size(), 1u);
  BOOST_CHECK_EQUAL(X.summands().front().variable().name(), X_false_name(result));
}

BOOST_AUTO_TEST_CASE(test_simplify_srf_pbes_all_disjuncts_dead_becomes_false)
{
  // If every disjunct is dead, the equation is false overall (the identity for an empty
  // disjunction). simplify_srf_pbes represents this explicitly as a single summand
  // targeting X_false, rather than leaving the equation with zero summands.
  srf_summand dead1(data::variable_list(), false_(), var("Y"));
  srf_summand dead2(data::variable_list(), false_(), var("Z"));
  srf_equation eqX(fixpoint_symbol::mu(), propositional_variable("X", data::variable_list()), {dead1, dead2}, false);

  srf_pbes result = simplify_srf_pbes(make_srf_pbes({eqX}));

  const srf_equation& X = find_equation(result, "X");
  BOOST_CHECK_EQUAL(X.summands().size(), 1u);
  BOOST_CHECK_EQUAL(X.summands().front().variable().name(), X_false_name(result));
}

BOOST_AUTO_TEST_CASE(test_simplify_srf_pbes_all_conjuncts_vacuous_becomes_true)
{
  // Dually, if every conjunct is vacuous, the equation is true overall (the identity
  // for an empty conjunction), represented as a single summand targeting X_true.
  srf_summand vacuous1(data::variable_list(), false_(), var("Y"));
  srf_summand vacuous2(data::variable_list(), false_(), var("Z"));
  srf_equation eqX(fixpoint_symbol::nu(), propositional_variable("X", data::variable_list()), {vacuous1, vacuous2}, true);

  srf_pbes result = simplify_srf_pbes(make_srf_pbes({eqX}));

  const srf_equation& X = find_equation(result, "X");
  BOOST_CHECK_EQUAL(X.summands().size(), 1u);
  BOOST_CHECK_EQUAL(X.summands().front().variable().name(), X_true_name(result));
}

BOOST_AUTO_TEST_CASE(test_simplify_srf_pbes_leaves_ordinary_equations_unchanged)
{
  // An equation with no trivial summands must survive simplify_srf_pbes unchanged, both
  // in summand count and in the (non-constant) condition/target of each summand.
  data::variable b1("b1", data::sort_bool::bool_());
  data::variable b2("b2", data::sort_bool::bool_());
  srf_summand s1(data::variable_list(), atermpp::down_cast<pbes_expression>(b1), var("Y"));
  srf_summand s2(data::variable_list(), atermpp::down_cast<pbes_expression>(b2), var("Z"));
  srf_equation eqX(fixpoint_symbol::nu(), propositional_variable("X", data::variable_list()), {s1, s2}, true);

  srf_pbes result = simplify_srf_pbes(make_srf_pbes({eqX}));

  const srf_equation& X = find_equation(result, "X");
  BOOST_CHECK_EQUAL(X.summands().size(), 2u);
  BOOST_CHECK_EQUAL(X.summands()[0].variable().name(), core::identifier_string("Y"));
  BOOST_CHECK_EQUAL(X.summands()[0].condition(), b1);
  BOOST_CHECK_EQUAL(X.summands()[1].variable().name(), core::identifier_string("Z"));
  BOOST_CHECK_EQUAL(X.summands()[1].condition(), b2);
}

//
// Integration tests for pre_srf2srfpbes: these confirm that substituting true/false for
// counter example variables produces exactly the trivial summand shapes that
// simplify_srf_pbes (tested directly above) recognizes, and that pre_srf2srfpbes indeed
// runs the result through it.
//

BOOST_AUTO_TEST_CASE(test_disjunctive_equation_collapses_to_true)
{
  // Z+ is always true, so X = Z+ || Y is unconditionally true, regardless of Y.
  // The whole equation should collapse to a single summand recursing into X_true,
  // discarding the (now irrelevant) dependency on Y.
  std::string text =
    "pbes\n"
    "mu X = Zpos_0_ || Y;\n"
    "mu Y = false;\n"
    "nu Zpos_0_ = true;\n"
    "init X;\n";
  pbes p = txt2pbes(text);
  srf_pbes srf = to_srf_via_ce(p);

  const srf_equation& eqn = find_equation(srf, "X");
  BOOST_CHECK_EQUAL(eqn.summands().size(), 1u);
  BOOST_CHECK_EQUAL(eqn.summands().front().variable().name(), X_true_name(srf));
  BOOST_CHECK_EQUAL(gauss_elimination(p), 1); // Sanity check: X is indeed semantically true.
}

BOOST_AUTO_TEST_CASE(test_conjunctive_equation_collapses_to_false)
{
  // Z- is always false, so X = Z- && Y is unconditionally false, regardless of Y.
  // The whole equation should collapse to a single summand recursing into X_false,
  // discarding the (now irrelevant) dependency on Y.
  std::string text =
    "pbes\n"
    "nu X = Zneg_0_ && Y;\n"
    "nu Y = true;\n"
    "mu Zneg_0_ = false;\n"
    "init X;\n";
  pbes p = txt2pbes(text);
  srf_pbes srf = to_srf_via_ce(p);

  const srf_equation& eqn = find_equation(srf, "X");
  BOOST_CHECK_EQUAL(eqn.summands().size(), 1u);
  BOOST_CHECK_EQUAL(eqn.summands().front().variable().name(), X_false_name(srf));
  BOOST_CHECK_EQUAL(gauss_elimination(p), 0); // Sanity check: X is indeed semantically false.
}

BOOST_AUTO_TEST_CASE(test_example_stramaglia_paper)
{
  // Example from A. Stramaglia, J.J.A. Keiren, M. Laveaux and T.A.C. Willemse. "Efficient Evidence Generation for Modal mu-Calculus Model Checking". TACAS 2025.
  // Substituting true/false for the Z_a/Z_b/Z_c variables turns
  // the second and fourth disjuncts of X's equation into dead disjuncts, and the guard
  // of Y's third summand into a vacuous conjunct; both should be simplified away.
  std::string text = "pbes\n"
                      "mu X(s: Nat) = (exists n: Nat. (val(s == 1) && val(0 < n && n < 3) && (X(s+n) || Zneg_0_(s, "
                      "s+n)) && Zpos_0_(s, s+n)))\n"
                      "              || Y(s);\n"
                      "nu Y(s: Nat) = val(s == 3) && (Y(s) || Zneg_1_(s,s)) && Zpos_1_(s,s);\n"
                      "nu Zpos_0_(s,s1: Nat) = true;\n"
                      "nu Zpos_1_(s,s1: Nat) = true;\n"
                      "mu Zneg_0_(s,s1: Nat) = false;\n"
                      "mu Zneg_1_(s,s1: Nat) = false;\n"
                      "\n"
                      "init X(0);\n";
  pbes p = txt2pbes(text);

  auto pre_srf = detail::pbes2pre_srf<true>(p, true);
  BOOST_CHECK_EQUAL(find_equation(pre_srf, "X").summands().size(), 3u);

  srf_pbes srf = pre_srf2srfpbes(pre_srf);
  BOOST_CHECK_EQUAL(find_equation(srf, "X").summands().size(), 2u);
  BOOST_CHECK_EQUAL(find_equation(srf, "Y").summands().size(), 2u);
}

BOOST_AUTO_TEST_CASE(test_issue_1949)
{
  // Regression test for https://github.com/mCRL2org/mCRL2/issues/1949: the PBES with
  // counter example information generated a larger structure graph than the one without,
  // because substituting true/false for the counterexample variables left behind trivial
  // summands like "true && X_true" that were never cleaned up. X0's equation below is a
  // disjunction with one top-level disjunct per transition; its first disjunct is the bare
  // counter example variable Zpos_0_c(x, x), which is unconditionally true (see its equation
  // below), so substitution alone makes the whole equation true -- independent of the other
  // three disjuncts, which carry the real (and otherwise irrelevant) recursion into X0.
  //
  // The PBESs below are obtained using lps2pbes on the following specification:
  //
  //   act  a,b,c;
  //   proc P(x: Bool) =
  //          c . P()
  //        + (x) -> a . P(x = false)
  //        + (!x) -> b . P(x = true)
  //        + delta;
  //   init P(false);
  //
  // with formula:
  //
  //   [true*][b]<true*><c>true
  std::string with_ce_text =
    "pbes nu Z(x: Bool) =\n"
    "       (val(!!x) || X0(true) && Zpos_2_b(x, true) || Zneg_2_b(x, true)) && (Z(x) && Zpos_0_c(x, x) || Zneg_0_c(x, x)) && (val(!x) || Z(false) && Zpos_1_a(x, false) || Zneg_1_a(x, false)) && (val(!!x) || Z(true) && Zpos_2_b(x, true) || Zneg_2_b(x, true));\n"
    "     mu X0(x: Bool) =\n"
    "       Zpos_0_c(x, x) || (X0(x) || Zneg_0_c(x, x)) && Zpos_0_c(x, x) || val(x) && (X0(false) || Zneg_1_a(x, false)) && Zpos_1_a(x, false) || val(!x) && (X0(true) || Zneg_2_b(x, true)) && Zpos_2_b(x, true);\n"
    "     nu Zneg_0_c(x,x1: Bool) =\n"
    "       val(false);\n"
    "     nu Zneg_1_a(x,x1: Bool) =\n"
    "       val(false);\n"
    "     nu Zneg_2_b(x,x1: Bool) =\n"
    "       val(false);\n"
    "     nu Zpos_0_c(x,x1: Bool) =\n"
    "       val(true);\n"
    "     nu Zpos_1_a(x,x1: Bool) =\n"
    "       val(true);\n"
    "     nu Zpos_2_b(x,x1: Bool) =\n"
    "       val(true);\n"
    "\n"
    "init Z(false);\n";

  // The equivalent PBES without counter example information: the counter example
  // machinery collapses entirely, leaving X0 = true.
  std::string no_ce_text =
    "pbes nu Z(x: Bool) =\n"
    "       (val(!!x) || X0(true)) && Z(x) && (val(!x) || Z(false)) && (val(!!x) || Z(true));\n"
    "     mu X0(x: Bool) =\n"
    "       val(true);\n"
    "\n"
    "init Z(false);\n";

  pbes p_ce = txt2pbes(with_ce_text);
  pbes p_no_ce = txt2pbes(no_ce_text);

  srf_pbes srf_ce = to_srf_via_ce(p_ce);
  srf_pbes srf_no_ce = detail::pbes2pre_srf<false>(p_no_ce, true);

  // Both must produce the same equations: Z, X0, X_false, X_true (i.e. the counter
  // example bookkeeping equations for Z+/Z- must all have disappeared, and X0 must
  // collapse to a single trivial summand, exactly as in the counter-example-free case).
  BOOST_CHECK_EQUAL(srf_ce.equations().size(), srf_no_ce.equations().size());
  BOOST_CHECK_EQUAL(srf_ce.equations().size(), 4u);

  const srf_equation& X0 = find_equation(srf_ce, "X0");
  BOOST_CHECK(!X0.is_conjunctive());
  BOOST_CHECK_EQUAL(X0.summands().size(), 1u);
  BOOST_CHECK_EQUAL(X0.summands().front().variable().name(), X_true_name(srf_ce));

  const srf_equation& Z = find_equation(srf_ce, "Z");
  BOOST_CHECK(Z.is_conjunctive());
  BOOST_CHECK_EQUAL(Z.summands().size(), find_equation(srf_no_ce, "Z").summands().size());
}

BOOST_AUTO_TEST_CASE(test_issue_1949_simpler_example)
{
  // Regression test for the second, simpler example given in the discussion of
  // https://github.com/mCRL2org/mCRL2/issues/1949. The PBESs
  // below are again the verbatim output of lps2pbes (with -otext, with and without
  // --counter-example) applied to the mCRL2 specification and modal formula from that
  // comment:
  //
  //   act  b,c;
  //   proc P(x: Bool) =
  //          c . P()
  //        + (x) -> b . P(x = false)
  //        + (!x) -> b . P(x = true)
  //        + delta;
  //   init P(false);
  //
  // with formula:
  //
  //   <true*.c>true
  //
  // Without counter example information this reduces (by lps2pbes itself) to the single
  // trivial equation Y = true. With counter example information, Y's equation is a
  // disjunction whose first disjunct is the bare counter example variable Zpos_0_c(x, x),
  // which is unconditionally true, so it should collapse to true here as well.
  std::string with_ce_text =
    "pbes mu Y(x: Bool) =\n"
    "       Zpos_0_c(x, x) || (Y(x) || Zneg_0_c(x, x)) && Zpos_0_c(x, x) || val(x) && (Y(false) || Zneg_2_b(x, false)) && Zpos_2_b(x, false) || val(!x) && (Y(true) || Zneg_2_b(x, true)) && Zpos_2_b(x, true);\n"
    "     nu Zneg_0_c(x,x1: Bool) =\n"
    "       val(false);\n"
    "     nu Zneg_2_b(x,x1: Bool) =\n"
    "       val(false);\n"
    "     nu Zpos_0_c(x,x1: Bool) =\n"
    "       val(true);\n"
    "     nu Zpos_2_b(x,x1: Bool) =\n"
    "       val(true);\n"
    "\n"
    "init Y(false);\n";

  std::string no_ce_text =
    "pbes mu Y(x: Bool) =\n"
    "       val(true);\n"
    "\n"
    "init Y(false);\n";

  pbes p_ce = txt2pbes(with_ce_text);
  pbes p_no_ce = txt2pbes(no_ce_text);

  srf_pbes srf_ce = to_srf_via_ce(p_ce);
  srf_pbes srf_no_ce = detail::pbes2pre_srf<false>(p_no_ce, true);

  // Both must produce the same equations: Y, X_false, X_true, with Y collapsed to a
  // single trivial summand, exactly as in the counter-example-free case.
  BOOST_CHECK_EQUAL(srf_ce.equations().size(), srf_no_ce.equations().size());
  BOOST_CHECK_EQUAL(srf_ce.equations().size(), 3u);

  const srf_equation& Y = find_equation(srf_ce, "Y");
  BOOST_CHECK(!Y.is_conjunctive());
  BOOST_CHECK_EQUAL(Y.summands().size(), 1u);
  BOOST_CHECK_EQUAL(Y.summands().front().variable().name(), X_true_name(srf_ce));
}

BOOST_AUTO_TEST_CASE(test_pre_srf2srfpbes_after_unify_parameters)
{
  // pre_srf2srfpbes must keep collapsing trivial
  // equations even when unify_parameters has already attached a
  // non-empty, unified parameter list to X_true/X_false.
  std::string with_ce_text =
    "pbes mu Y(x: Bool) =\n"
    "       Zpos_0_c(x, x) || (Y(x) || Zneg_0_c(x, x)) && Zpos_0_c(x, x) || val(x) && (Y(false) || Zneg_2_b(x, false)) && Zpos_2_b(x, false) || val(!x) && (Y(true) || Zneg_2_b(x, true)) && Zpos_2_b(x, true);\n"
    "     nu Zneg_0_c(x,x1: Bool) =\n"
    "       val(false);\n"
    "     nu Zneg_2_b(x,x1: Bool) =\n"
    "       val(false);\n"
    "     nu Zpos_0_c(x,x1: Bool) =\n"
    "       val(true);\n"
    "     nu Zpos_2_b(x,x1: Bool) =\n"
    "       val(true);\n"
    "\n"
    "init Y(false);\n";

  pbes p = txt2pbes(with_ce_text);

  auto pre_srf = pbes2pre_srf(p, true); // the public wrapper, allow_ce = true by default
  unify_parameters(pre_srf, true, false);
  resolve_summand_variable_name_clashes(pre_srf, pre_srf.equations().front().variable().parameters());

  // unify_parameters has given every propositional variable, including X_true, the
  // process parameter x, so the target of Y's (would-be) trivial summand is now
  // X_true(x), not the bare X_true() seen in the tests above.
  const propositional_variable_instantiation& X_true_instantiation
    = find_equation(pre_srf, "X_true").summands().front().variable();
  BOOST_CHECK_EQUAL(X_true_instantiation.parameters().size(), 1u);

  srf_pbes srf = pre_srf2srfpbes(pre_srf);

  const srf_equation& Y = find_equation(srf, "Y");
  BOOST_CHECK_EQUAL(Y.summands().size(), 1u);
  BOOST_CHECK_EQUAL(Y.summands().front().variable().name(), X_true_name(srf));
  BOOST_CHECK_EQUAL(Y.summands().front().variable().parameters().size(), 1u); // Still carries x, unlike the tests above.
}
