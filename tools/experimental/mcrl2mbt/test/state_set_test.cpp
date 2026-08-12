// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE state_set_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/data/rewriter.h"
#include "mcrl2/lps/explorer_options.h"
#include "mcrl2/lps/parse.h"

#include "../io_classifier.cpp" // NOLINT
#include "../state_set.cpp" // NOLINT

using namespace mcrl2;

// Helper: construct a state_set from an inline LPS spec.
static state_set make_state_set(const std::string& spec_text)
{
  lps::specification spec = lps::parse_linear_process_specification(spec_text);
  data::rewriter rewr(spec.data());
  lps::explorer_options opts;
  return state_set(spec, opts, rewr);
}

// Helper: build a classifier from two lists of action names.
static io_classifier make_classifier(const std::vector<std::string>& inputs, const std::vector<std::string>& outputs)
{
  // Write a temp file and load it.
  char tmpl[] = "/tmp/mbt_ss_test_XXXXXX";
  const int fd = mkstemp(tmpl);
  BOOST_REQUIRE(fd >= 0);
  std::string content = "input\n";
  for (const auto& a: inputs)
  {
    content += "  " + a + "\n";
  }
  content += "output\n";
  for (const auto& a: outputs)
  {
    content += "  " + a + "\n";
  }
  ::write(fd, content.data(), content.size());
  ::close(fd);
  io_classifier cls;
  cls.load(tmpl);
  std::remove(tmpl);
  return cls;
}

// LPS: single output action looping.
static const std::string SPEC_OUTPUT_LOOP = "act b;\n"
                                            "proc P = b.P;\n"
                                            "init P;\n";

// LPS: single input action looping.
static const std::string SPEC_INPUT_LOOP = "act a;\n"
                                           "proc P = a.P;\n"
                                           "init P;\n";

// LPS: tau leading to a state with action a.
// P(false) --tau--> P(true) --a--> P(true)
static const std::string SPEC_TAU_THEN_A = "act a;\n"
                                           "proc P(reached: Bool) = (!reached) -> tau.P(true)\n"
                                           "                      + reached    -> a.P(true);\n"
                                           "init P(false);\n";

// LPS: both input and output actions from same state.
static const std::string SPEC_IO = "act inp, outp;\n"
                                   "proc P = inp.P + outp.P;\n"
                                   "init P;\n";

// ─────────────────────────────────────────────────────────────────────────────

BOOST_AUTO_TEST_CASE(reset_to_initial_is_singleton)
{
  state_set ss = make_state_set(SPEC_INPUT_LOOP);
  ss.reset_to_initial();
  BOOST_CHECK_EQUAL(ss.size(), 1u);
}

BOOST_AUTO_TEST_CASE(is_action_enabled_true)
{
  state_set ss = make_state_set(SPEC_INPUT_LOOP);
  ss.reset_to_initial();
  BOOST_CHECK(ss.is_action_enabled("a", 0));
}

BOOST_AUTO_TEST_CASE(is_action_enabled_false_wrong_label)
{
  state_set ss = make_state_set(SPEC_INPUT_LOOP);
  ss.reset_to_initial();
  BOOST_CHECK(!ss.is_action_enabled("z", 0));
}

BOOST_AUTO_TEST_CASE(is_action_enabled_via_tau)
{
  state_set ss = make_state_set(SPEC_TAU_THEN_A);
  ss.reset_to_initial();
  // With depth 0: no tau closure — 'a' is not enabled from P(false).
  BOOST_CHECK(!ss.is_action_enabled("a", 0));
  // With depth 1: follows the tau to P(true) — 'a' becomes enabled.
  BOOST_CHECK(ss.is_action_enabled("a", 1));
}

BOOST_AUTO_TEST_CASE(apply_action_enabled_returns_true)
{
  state_set ss = make_state_set(SPEC_INPUT_LOOP);
  ss.reset_to_initial();
  BOOST_CHECK(ss.apply_action("a", 0));
  BOOST_CHECK_EQUAL(ss.size(), 1u); // still one state (loop)
}

BOOST_AUTO_TEST_CASE(apply_action_disabled_returns_false)
{
  state_set ss = make_state_set(SPEC_INPUT_LOOP);
  ss.reset_to_initial();
  BOOST_CHECK(!ss.apply_action("z", 0));
  // S must remain unchanged.
  BOOST_CHECK_EQUAL(ss.size(), 1u);
}

BOOST_AUTO_TEST_CASE(get_enabled_classifies_input_and_output)
{
  state_set ss = make_state_set(SPEC_IO);
  ss.reset_to_initial();
  const auto cls = make_classifier({"inp"}, {"outp"});
  const auto result = ss.get_enabled(0, cls);
  BOOST_REQUIRE_EQUAL(result.inputs.size(), 1u);
  BOOST_REQUIRE_EQUAL(result.outputs.size(), 1u);
  BOOST_CHECK_EQUAL(result.inputs[0], "inp");
  BOOST_CHECK_EQUAL(result.outputs[0], "outp");
  BOOST_CHECK(!result.quiescence); // output present — not quiescent
}

BOOST_AUTO_TEST_CASE(get_enabled_quiescent_when_only_inputs)
{
  state_set ss = make_state_set(SPEC_INPUT_LOOP);
  ss.reset_to_initial();
  const auto cls = make_classifier({"a"}, {});
  const auto result = ss.get_enabled(0, cls);
  BOOST_CHECK(result.quiescence);
  BOOST_REQUIRE_EQUAL(result.inputs.size(), 1u);
  BOOST_CHECK(result.outputs.empty());
}

BOOST_AUTO_TEST_CASE(refine_to_quiescent_false_when_output_present)
{
  state_set ss = make_state_set(SPEC_OUTPUT_LOOP);
  ss.reset_to_initial();
  const auto cls = make_classifier({}, {"b"});
  BOOST_CHECK(!ss.refine_to_quiescent(0, cls));
}

BOOST_AUTO_TEST_CASE(refine_to_quiescent_true_when_only_inputs)
{
  state_set ss = make_state_set(SPEC_INPUT_LOOP);
  ss.reset_to_initial();
  const auto cls = make_classifier({"a"}, {});
  BOOST_CHECK(ss.refine_to_quiescent(0, cls));
  BOOST_CHECK_EQUAL(ss.size(), 1u);
}

BOOST_AUTO_TEST_CASE(tau_depth_zero_does_not_expand)
{
  // With depth 0, apply_action("a", 0) from P(false) fails because 'a' is
  // only reachable after a tau step.
  state_set ss = make_state_set(SPEC_TAU_THEN_A);
  ss.reset_to_initial();
  BOOST_CHECK(!ss.apply_action("a", 0));
}

BOOST_AUTO_TEST_CASE(tau_depth_one_follows_tau_edge)
{
  state_set ss = make_state_set(SPEC_TAU_THEN_A);
  ss.reset_to_initial();
  // With depth 1, apply_action succeeds by first following the tau.
  BOOST_CHECK(ss.apply_action("a", 1));
}
