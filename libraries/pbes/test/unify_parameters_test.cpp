// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file unify_parameters_test.cpp
/// \brief Test for unify_parameters, which rewrites every propositional variable
/// instantiation of a PBES to carry a common, unified parameter list.

#define BOOST_TEST_MODULE unify_parameters_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/txt2pbes.h"
#include "mcrl2/pbes/unify_parameters.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system;

const pbes_equation& find_equation(const pbes& p, const std::string& name)
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

BOOST_AUTO_TEST_CASE(test_has_unified_parameters)
{
  BOOST_CHECK(has_unified_parameters(txt2pbes(
    "pbes\n"
    "mu X(n: Nat) = Y(n);\n"
    "nu Y(n: Nat) = X(n);\n"
    "init X(0);\n")));

  BOOST_CHECK(!has_unified_parameters(txt2pbes(
    "pbes\n"
    "mu X(n: Nat) = Y;\n"
    "nu Y = X(0);\n"
    "init X(0);\n")));

  // A single equation trivially has unified parameters.
  BOOST_CHECK(has_unified_parameters(txt2pbes(
    "pbes\n"
    "mu X(n: Nat) = X(n);\n"
    "init X(0);\n")));
}

BOOST_AUTO_TEST_CASE(test_unify_parameters_reset_false_reuses_callers_variable)
{
  // X and Y have different arity: unify_parameters has something to do.
  std::string text =
    "pbes\n"
    "mu X(n: Nat, b: Bool) = Y(n) || Zpos_0_;\n"
    "nu Y(n: Nat) = X(n+1, true);\n"
    "nu Zpos_0_ = true;\n"
    "init X(0, false);\n";

  pbes p = txt2pbes(text);
  BOOST_CHECK_EQUAL(find_equation(p, "X").variable().parameters().size(), 2u);
  BOOST_CHECK_EQUAL(find_equation(p, "Y").variable().parameters().size(), 1u);

  // With ignore_ce_equations=true and reset=false: X and Y end up with the same
  // (unified) parameter list; the counter example equation Zpos_0_ is left untouched
  // (0 parameters); and the newly introduced parameter of Y is filled in at the call
  // site inside X's equation by reusing X's own variable b (rather than a fresh value).
  unify_parameters(p, true, false);
  BOOST_CHECK(p.is_well_typed());
  BOOST_CHECK_EQUAL(find_equation(p, "X").variable().parameters().size(), 2u);
  BOOST_CHECK_EQUAL(find_equation(p, "Y").variable().parameters(), find_equation(p, "X").variable().parameters());
  BOOST_CHECK_EQUAL(find_equation(p, "Zpos_0_").variable().parameters().size(), 0u);
  BOOST_CHECK_EQUAL(pbes_system::pp(find_equation(p, "X").formula()), "Y(n, b) || Zpos_0_");
}

BOOST_AUTO_TEST_CASE(test_unify_parameters_reset_true_uses_default_value)
{
  // Same PBES as above, but now with reset=true: the newly introduced parameter is
  // filled in with a fresh default value for its sort, not the caller's variable b.
  std::string text =
    "pbes\n"
    "mu X(n: Nat, b: Bool) = Y(n) || Zpos_0_;\n"
    "nu Y(n: Nat) = X(n+1, true);\n"
    "nu Zpos_0_ = true;\n"
    "init X(0, false);\n";

  pbes p = txt2pbes(text);
  unify_parameters(p, true, true);
  BOOST_CHECK(p.is_well_typed());
  BOOST_CHECK_EQUAL(pbes_system::pp(find_equation(p, "X").formula()), "Y(n, false) || Zpos_0_");
}

BOOST_AUTO_TEST_CASE(test_unify_parameters_merges_parameters_with_the_same_name_and_sort)
{
  // X and Y share a parameter n : Nat (same name, same sort). unify_parameters must
  // recognize these as the same parameter, and not duplicate it: the unified parameter
  // list is [n, b, c] (size 3), not [n, b, n, c] (size 4).
  std::string text =
    "pbes\n"
    "mu X(n: Nat, b: Bool) = Y(n, true) || val(b);\n"
    "nu Y(n: Nat, c: Bool) = X(n, c);\n"
    "init X(0, false);\n";

  pbes p = txt2pbes(text);
  unify_parameters(p, false, false);
  BOOST_CHECK(p.is_well_typed());
  BOOST_CHECK_EQUAL(find_equation(p, "X").variable().parameters().size(), 3u);
  BOOST_CHECK_EQUAL(find_equation(p, "X").variable().parameters(), find_equation(p, "Y").variable().parameters());
  BOOST_CHECK(has_unified_parameters(p));
}

BOOST_AUTO_TEST_CASE(test_unify_parameters_yields_has_unified_parameters_without_counterexample_equations)
{
  // has_unified_parameters looks at every equation, without any counter example
  // exception. So once unify_parameters has run on a PBES that has no counter example
  // equations to begin with, has_unified_parameters must become true.
  std::string text =
    "pbes\n"
    "mu X(n: Nat, b: Bool) = Y(n);\n"
    "nu Y(n: Nat) = X(n, true);\n"
    "init X(0, false);\n";

  pbes p = txt2pbes(text);
  BOOST_CHECK(!has_unified_parameters(p));
  unify_parameters(p, true, false);
  BOOST_CHECK(has_unified_parameters(p));
}

BOOST_AUTO_TEST_CASE(test_unify_parameters_counterexample_equation_left_hand_side_is_never_updated)
{
  // Regardless of ignore_ce_equations, a counter example equation's own left-hand side
  // parameter list is never rewritten (only ordinary equations' left-hand sides are).
  std::string text =
    "pbes\n"
    "mu X(n: Nat) = Y(n) || Zpos_0_(n);\n"
    "nu Y(n: Nat) = X(n+1);\n"
    "nu Zpos_0_(m: Nat) = true;\n"
    "init X(0);\n";

  for (bool ignore_ce_equations: {true, false})
  {
    pbes p = txt2pbes(text);
    unify_parameters(p, ignore_ce_equations, false);
    BOOST_CHECK(p.is_well_typed());
    BOOST_CHECK_EQUAL(find_equation(p, "Zpos_0_").variable().parameters().size(), 1u); // still just m
  }
}

BOOST_AUTO_TEST_CASE(test_unify_parameters_ignore_ce_equations_controls_the_union)
{
  // With ignore_ce_equations=true, the counter example equation Zpos_0_'s own parameter
  // m is excluded when computing the unified parameter list, so X and Y only unify to
  // their shared parameter n. With ignore_ce_equations=false, m is pulled into the
  // union as well, so X and Y (which never mention m) gain it as an extra parameter.
  std::string text =
    "pbes\n"
    "mu X(n: Nat) = Y(n) || Zpos_0_(n);\n"
    "nu Y(n: Nat) = X(n+1);\n"
    "nu Zpos_0_(m: Nat) = true;\n"
    "init X(0);\n";

  {
    pbes p = txt2pbes(text);
    unify_parameters(p, true, false);
    BOOST_CHECK_EQUAL(find_equation(p, "X").variable().parameters().size(), 1u);
  }
  {
    pbes p = txt2pbes(text);
    unify_parameters(p, false, false);
    BOOST_CHECK_EQUAL(find_equation(p, "X").variable().parameters().size(), 2u);
  }
}

BOOST_AUTO_TEST_CASE(test_unify_parameters_srf_pbes_counterexample_summand_targets_always_use_defaults)
{
  // The SRF overload of unify_parameters (operating on a pre_srf_pbes rather than a
  // plain pbes) treats references to X_true/X_false from within a counter example
  // equation's own summand specially: unlike ordinary equations, these are always
  // padded with fresh default values, regardless of the reset flag passed in -- since a
  // counter example equation's own (never-unified) parameters give no value to reuse.
  std::string text =
    "pbes\n"
    "mu X(n: Nat, b: Bool) = Y(n) || Zpos_0_;\n"
    "nu Y(n: Nat) = X(n+1, true);\n"
    "nu Zpos_0_ = true;\n"
    "init X(0, false);\n";

  propositional_variable_instantiation zpos_0_target_reset_false;
  propositional_variable_instantiation zpos_0_target_reset_true;
  propositional_variable_instantiation y_target_reset_false;
  propositional_variable_instantiation y_target_reset_true;

  {
    pbes p = txt2pbes(text);
    auto pre_srf = pbes2pre_srf(p, true);
    unify_parameters(pre_srf, true, false);
    zpos_0_target_reset_false = find_equation(pre_srf, "Zpos_0_").summands().front().variable();
    y_target_reset_false = find_equation(pre_srf, "X").summands().front().variable();
  }
  {
    pbes p = txt2pbes(text);
    auto pre_srf = pbes2pre_srf(p, true);
    unify_parameters(pre_srf, true, true);
    zpos_0_target_reset_true = find_equation(pre_srf, "Zpos_0_").summands().front().variable();
    y_target_reset_true = find_equation(pre_srf, "X").summands().front().variable();
  }

  // Both refer to X_true, now carrying the unified [n, b] parameters.
  BOOST_CHECK_EQUAL(zpos_0_target_reset_false.parameters().size(), 2u);
  BOOST_CHECK_EQUAL(zpos_0_target_reset_true.parameters().size(), 2u);

  // The counter example equation's summand target is identical regardless of reset,
  // since it is always padded with defaults.
  BOOST_CHECK_EQUAL(zpos_0_target_reset_false, zpos_0_target_reset_true);

  // In contrast, X's own (ordinary) summand target Y(n) is missing parameter b, and
  // padding it does respect the reset flag: reused variable versus default value.
  BOOST_CHECK_EQUAL(pbes_system::pp(y_target_reset_false), "Y(n, b)");
  BOOST_CHECK_EQUAL(pbes_system::pp(y_target_reset_true), "Y(n, false)");
  BOOST_CHECK(y_target_reset_false != y_target_reset_true);
}
