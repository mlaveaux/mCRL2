// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file rename_test.cpp
/// \brief Add your file description here.

#include <boost/core/no_exceptions_support.hpp>
#define MCRL2_STATE_FORMULA_BUILDER_DEBUG

#define BOOST_TEST_MODULE state_formula
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/lps/linearise.h"
#include "mcrl2/modal_formula/count_fixpoints.h"
#include "mcrl2/modal_formula/maximal_closed_subformula.h"
#include "mcrl2/modal_formula/parse.h"
#include "mcrl2/modal_formula/preprocess_state_formula.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

using namespace mcrl2;
using namespace mcrl2::core;
using namespace mcrl2::lps;
using namespace mcrl2::lps::detail;
using namespace mcrl2::state_formulas;

BOOST_AUTO_TEST_CASE(test_rename)
{
  const std::string SPECIFICATION =
    "act a:Nat;                              \n"
    "                                        \n"
    "map smaller: Nat#Nat -> Bool;           \n"
    "                                        \n"
    "var x,y : Nat;                          \n"
    "                                        \n"
    "eqn smaller(x,y) = x < y;               \n"
    "                                        \n"
    "proc P(n:Nat) = sum m: Nat. a(m). P(m); \n"
    "                                        \n"
    "init P(0);                              \n"
    ;

  using mcrl2::state_formulas::pp;
  specification lpsspec = remove_stochastic_operators(linearise(SPECIFICATION));

  state_formula formula = parse_state_formula("(mu X. X) && (mu X. X)", lpsspec, false);

  data::set_identifier_generator generator;
  generator.add_identifiers(lps::find_identifiers(lpsspec));
  formula = rename_predicate_variables(formula, generator);

  std::cout << "pp(formula) == " << pp(formula) << std::endl;
  BOOST_CHECK(pp(formula) == "(mu X1. X1) && (mu X. X)" || pp(formula) == "(mu X. X) && (mu X1. X1)");

  generator = data::set_identifier_generator();
  generator.add_identifiers(lps::find_identifiers(lpsspec));
  parse_state_formula_options options;
  options.check_monotonicity = false;
  options.resolve_name_clashes = false;
  formula = parse_state_formula("mu X. mu X. X", lpsspec, false, options);
  std::cout << "formula: " << pp(formula) << std::endl;
  formula = rename_predicate_variables(formula, generator);
  std::cout << "formula: " << pp(formula) << std::endl;
  BOOST_CHECK_EQUAL(pp(formula), "mu X. mu X1. X1");
}

BOOST_AUTO_TEST_CASE(test_normalize)
{
  state_formula x = state_formulas::variable(identifier_string("X"), data::data_expression_list());
  state_formula y = state_formulas::variable(identifier_string("Y"), data::data_expression_list());
  state_formula f;
  state_formula f1;
  state_formula f2;

  f = imp(not_(x), y);
  f1 = normalize(f);
  f2 = or_(x, y);
  std::cout << "f  = " << state_formulas::pp(f) << std::endl;
  std::cout << "f1 = " << state_formulas::pp(f1) << std::endl;
  std::cout << "f2 = " << state_formulas::pp(f2) << std::endl;
  BOOST_CHECK_EQUAL(f1, f2);

  f  = not_(and_(not_(x), not_(y)));
  f1 = normalize(f);
  f2 = or_(x, y);
  std::cout << "f  = " << state_formulas::pp(f) << std::endl;
  std::cout << "f1 = " << state_formulas::pp(f1) << std::endl;
  std::cout << "f2 = " << state_formulas::pp(f2) << std::endl;
  BOOST_CHECK_EQUAL(f1, f2);
}

BOOST_AUTO_TEST_CASE(test_type_checking)
{
  specification lpsspec = remove_stochastic_operators(linearise(
                            "sort CPU = struct p1;"
                            "sort CPUs = Set(CPU);"
                            "init delta;"));

  state_formula formula = parse_state_formula("nu X (P : CPUs = {p1}) . val(P != {})", lpsspec, false);
}

BOOST_AUTO_TEST_CASE(test_type_checking_conversion_of_arguments)
{
  specification lpsspec = remove_stochastic_operators(linearise(
                            "sort B = struct d;"
                            "act a: List(B);"
                            "init a([d]);"
                          ));

  state_formula formula = parse_state_formula("<a([d])>true", lpsspec, false);

  BOOST_CHECK(is_may(formula));
  const state_formulas::may& f = atermpp::down_cast<state_formulas::may>(formula);
  BOOST_CHECK(is_regular_formula(f.formula()));
}


static inline
state_formula negate_variable(const variable& x)
{
  return state_formulas::not_(x);
}

BOOST_AUTO_TEST_CASE(test_not)
{
  data::data_expression_list args;
  variable v(core::identifier_string("v"), args);
  state_formula s = not_(v);
  BOOST_CHECK(is_not(s));

  state_formula t = negate_variable(v);
  BOOST_CHECK_EQUAL(s, t);

  // The following is expected to trigger an assertion failure
  // aterm a = v;
  // state_formula t = not_(a);
}

// test case supplied by Jan Friso, 4-1-2011
BOOST_AUTO_TEST_CASE(test_parse)
{
  std::string spec_text =
    "act a:Nat; \n"
    "init a(1); \n"
    ;

  std::string formula_text = "<a(1)>true";

  lps::specification lpsspec = remove_stochastic_operators(lps::linearise(spec_text));
  state_formulas::state_formula f = state_formulas::parse_state_formula(formula_text, lpsspec, false);

  std::cerr << "--- f ---\n" << state_formulas::pp(f) << std::endl;
  std::set<core::identifier_string> ids = state_formulas::find_identifiers(f);
  BOOST_CHECK(ids.find(core::identifier_string("1")) == ids.end());
#ifndef MCRL2_ENABLE_MACHINENUMBERS
  BOOST_CHECK(ids.find(core::identifier_string("@c1")) != ids.end());
#endif
}

BOOST_AUTO_TEST_CASE(test_count_fixpoints)
{
  state_formula formula;
  specification lpsspec;

  formula = parse_state_formula("(mu X. X) && (mu X. X)", lpsspec, false);
  BOOST_CHECK_EQUAL(count_fixpoints(formula), 2u);

  formula = parse_state_formula("exists b:Bool. (mu X. X) || forall b:Bool. (nu X. mu Y. (X || Y))", lpsspec, false);
  BOOST_CHECK_EQUAL(count_fixpoints(formula), 3u);
}

// Test case for bug #1094.
// This is expected to fail, due to the occurrence of an unbound variable Y
// in the formula.
BOOST_AUTO_TEST_CASE(test_1094)
{
  const std::string SPEC =
      "act  a,b,c,d;\n"
      "\n"
      "proc P(s3_X: Pos) =\n"
      "       (s3_X == 1) ->\n"
      "         a .\n"
      "         P(s3_X = 2)\n"
      "     + (s3_X == 2) ->\n"
      "         b .\n"
      "         P(s3_X = 3)\n"
      "     + (s3_X == 3) ->\n"
      "         a .\n"
      "         P(s3_X = 4)\n"
      "     + (s3_X == 3) ->\n"
      "         a .\n"
      "         P(s3_X = 5)\n"
      "     + (s3_X == 4) ->\n"
      "         c .\n"
      "         P(s3_X = 1)\n"
      "     + (s3_X == 5) ->\n"
      "         d .\n"
      "         P(s3_X = 1)\n"
      "     + delta;\n"
      "\n"
      "init P(1);\n"
    ;
  specification s(parse_linear_process_specification(SPEC));

  const std::string FORMULA = "[true*]([b] nu X. mu X.( [!c]X && [!c]Y))";

  BOOST_CHECK_THROW(parse_state_formula(FORMULA, s, false), mcrl2::runtime_error);
}

inline
state_formula sigma(const state_formula& x)
{
  variable X("X", data::data_expression_list());
  return x == X ? false_() : x;
}

BOOST_AUTO_TEST_CASE(test_replace_state_formulas)
{
  specification lpsspec;
  state_formula f = parse_state_formula("(mu X. X) && (mu X. X)", lpsspec, false);
  state_formula result = replace_state_formulas(f, sigma);
  state_formula expected_result = parse_state_formula("(mu X. false) && (mu X. false)", lpsspec, false);
  if (!(result == expected_result))
  {
    std::cout << "error: " << state_formulas::pp(result) << " != " << state_formulas::pp(expected_result) << std::endl;
  }
  BOOST_CHECK(result == expected_result);
}

BOOST_AUTO_TEST_CASE(test_find_state_variables)
{
  specification lpsspec;

  state_formula f = parse_state_formula("(mu X. nu Y. true && mu Z. X && Z)", lpsspec, false);
  std::set<state_formulas::variable> v = state_formulas::find_state_variables(f);
  BOOST_CHECK(v.size() == 2);

  f = parse_state_formula("mu X. nu Y. (true && mu Z. (X && Y || Z))", lpsspec, false);
  v = find_state_variables(f);
  BOOST_CHECK(v.size() == 3);

  state_formulas::variable X("X", data::data_expression_list());
  state_formulas::variable Y("Y", data::data_expression_list());
  state_formulas::variable Z("Z", data::data_expression_list());
  state_formula phi = state_formulas::and_(X, Y);
  v = find_state_variables(phi);
  BOOST_CHECK(v.size() == 2);
  v = find_free_state_variables(phi);
  BOOST_CHECK(v.size() == 2);
  state_formula psi = state_formulas::mu("X", data::assignment_list(), phi);
  v = find_state_variables(psi);
  BOOST_CHECK(v.size() == 2);
  v = find_free_state_variables(psi);
  BOOST_CHECK(v.size() == 1);
}

inline
bool contains(const std::set<state_formulas::state_formula>& v, const std::string& s)
{
  for (const state_formulas::state_formula& f: v)
  {
    if (state_formulas::pp(f) == s)
    {
      return true;
    }
  }
  return false;
}

BOOST_AUTO_TEST_CASE(test_maximal_closed_subformulas)
{
  mcrl2::log::logger::set_reporting_level(mcrl2::log::debug);

  specification lpsspec;
  state_formula f = parse_state_formula("(mu X. nu Y. true && mu Z. X && Z)", lpsspec, false);
  std::set<state_formulas::state_formula> v = maximal_closed_subformulas(f);
  BOOST_CHECK(v.size() == 1);

  f = parse_state_formula("exists b: Bool. forall c: Bool. val(b) && (val(c) || true) && false", lpsspec, false);
  v = maximal_closed_subformulas(f);
  BOOST_CHECK(v.size() == 1);

  state_formula g = exists(f).body();
  std::cout << "g = " << state_formulas::pp(g) << std::endl;
  v = maximal_closed_subformulas(g);
  for (const state_formula& f: v)
  {
    std::cout << "element " << f << std::endl;
  }
  BOOST_CHECK(v.size() == 2);
  BOOST_CHECK(contains(v, "true"));
  BOOST_CHECK(contains(v, "false"));

  state_formula h = forall(g).body();
  v = maximal_closed_subformulas(h);
  std::cout << "h = " << state_formulas::pp(h) << std::endl;
  for (const state_formula& f: v)
  {
    std::cout << "element " << f << std::endl;
  }
  BOOST_CHECK(v.size() == 2);
  BOOST_CHECK(contains(v, "true"));
  BOOST_CHECK(contains(v, "false"));
}

void test_has_unscoped_modal_formulas(const std::string& text, lps::specification& lpsspec, bool expected_result)
{
  state_formula x = parse_state_formula(text, lpsspec, false);
  bool result = state_formulas::detail::has_unscoped_modal_formulas(x);
  if (result != expected_result)
  {
    std::cout << "--- has_unscoped_modal_formulas failure ---" << std::endl;
    std::cout << "formula = " << text << std::endl;
    std::cout << "result  = " << result << std::endl;
  }
  BOOST_CHECK_EQUAL(result, expected_result);
}

BOOST_AUTO_TEST_CASE(has_unscoped_modal_formulas_test)
{
  const std::string text =
    "act a, b, c;                  \n"
    "                              \n"
    "proc P = delta;               \n"
    "                              \n"
    "init P;                       \n"
    ;
  lps::specification lpsspec = lps::parse_linear_process_specification(text);

  test_has_unscoped_modal_formulas("[a]true", lpsspec, true);
  test_has_unscoped_modal_formulas("mu X. X", lpsspec, false);
  test_has_unscoped_modal_formulas("true => [a]true", lpsspec, true);
  test_has_unscoped_modal_formulas("true => mu X. [a]<b>true", lpsspec, false);
}

void test_preprocess_nested_modal_operators(const std::string& text, lps::specification& lpsspec, const std::string& expected_result_text)
{
  state_formula x = parse_state_formula(text, lpsspec, false);
  state_formula expected_result = parse_state_formula(expected_result_text, lpsspec, false);
  state_formula result = state_formulas::preprocess_nested_modal_operators(x);
  if (result != expected_result)
  {
    std::cout << "--- preprocess_nested_modal_operators failure ---" << std::endl;
    std::cout << "formula          = " << text << std::endl;
    std::cout << "result           = " << state_formulas::pp(result) << std::endl;
    std::cout << "expected result  = " << expected_result_text << std::endl;
  }
  BOOST_CHECK_EQUAL(result, expected_result);
}

BOOST_AUTO_TEST_CASE(preprocess_nested_modal_operators_test)
{
  const std::string text =
    "act a, b, c;                  \n"
    "                              \n"
    "proc P = delta;               \n"
    "                              \n"
    "init P;                       \n"
    ;
  lps::specification lpsspec = lps::parse_linear_process_specification(text);

  test_preprocess_nested_modal_operators("[a]true", lpsspec, "[a]true");
  test_preprocess_nested_modal_operators("[a]<b>true", lpsspec, "[a]mu X. <b>true");
  test_preprocess_nested_modal_operators("true && <a><a>true", lpsspec, "true && <a>mu X. <a>true");
  test_preprocess_nested_modal_operators("true => mu X. [a]<b>true", lpsspec, "true => mu X. [a]mu X1. <b>true");
}

BOOST_AUTO_TEST_CASE(parse_modal_formula_test)
{
  const std::string text =
    "act a, b, c;                  \n"
    "                              \n"
    "proc P = delta;               \n"
    "                              \n"
    "init P;                       \n"
    ;
  lps::specification lpsspec = lps::parse_linear_process_specification(text);

  const std::string formula_text = "nu X(n: Nat = 0, b: Bool = false). X(n, b)";
  state_formula x = parse_state_formula(formula_text, lpsspec, false);
  std::string s = state_formulas::pp(x);
  BOOST_CHECK_EQUAL(s, formula_text);
}

BOOST_AUTO_TEST_CASE(check_parameter_name_clashes_test)
{
  const std::string text =
    "act a, b, c;                  \n"
    "                              \n"
    "proc P = delta;               \n"
    "                              \n"
    "init P;                       \n"
    ;
  lps::specification lpsspec = lps::parse_linear_process_specification(text);

  BOOST_CHECK(has_data_variable_name_clashes(parse_state_formula("nu X(n: Nat = 0). mu Y(n: Nat = 1). val(n == 0)", lpsspec, false)));
  BOOST_CHECK(!has_data_variable_name_clashes(parse_state_formula("nu X(n: Nat = 0). mu Y(m: Nat = 1). val(n == 0)", lpsspec, false)));
}

BOOST_AUTO_TEST_CASE(parse_state_formula_specification_test)
{
  const std::string text =
    "sort S;                       \n"
    "cons s0: S;                   \n"
    "act a: S;                     \n"
    "form forall s: S. [a(s)]true; \n"
    ;
  state_formula_specification x = parse_state_formula_specification(text, false);
  std::string text1 = state_formulas::pp(x);
  std::cout << text1 << std::endl;
  BOOST_CHECK_EQUAL(utilities::remove_whitespace(text), utilities::remove_whitespace(text1));
}

BOOST_AUTO_TEST_CASE(parse_state_formula_with_supremum_specification_test)
{
  const std::string text =
    "sort S;                       \n"
    "cons s0: S;                   \n"
    "act a: S;                     \n"
    "form sup s: S. [a(s)]true; \n"
    ;
  state_formula_specification x = parse_state_formula_specification(text, true);
  std::string text1 = state_formulas::pp(x);
  std::cout << text1 << std::endl;
  BOOST_CHECK_EQUAL(utilities::remove_whitespace(text), utilities::remove_whitespace(text1));
}

BOOST_AUTO_TEST_CASE(parse_state_formula_with_sum_specification_test)
{
  const std::string text =
    "sort S;                       \n"
    "cons s0: S;                   \n"
    "act a: S;                     \n"
    "form sum s: S. val(1 / 2) * [a(s)]false; \n"
    ;
  state_formula_specification x = parse_state_formula_specification(text, true);
  std::string text1 = state_formulas::pp(x);
  std::cout << text1 << std::endl;
  BOOST_CHECK_EQUAL(utilities::remove_whitespace(text), utilities::remove_whitespace(text1));
}

BOOST_AUTO_TEST_CASE(parse_state_formula_with_const_multiply_alt_specification_test)
{
  const std::string text =
    "sort S;                       \n"
    "cons s0: S;                   \n"
    "act a: S;                     \n"
    "form inf s: S. [a(s)]false * val(1/5); \n"
    ;
  state_formula_specification x = parse_state_formula_specification(text, true);
  std::string text1 = state_formulas::pp(x);
  std::cout << text1 << std::endl;
  BOOST_CHECK_EQUAL(utilities::remove_whitespace(text), utilities::remove_whitespace(text1));
}

BOOST_AUTO_TEST_CASE(parse_state_formula_with_infimum_specification_test)
{
  const std::string text =
    "form inf n: Nat. val(n); \n"
    ;
  state_formula_specification x = parse_state_formula_specification(text, true);
  std::string text1 = state_formulas::pp(x);
  std::cout << text1 << std::endl;
  BOOST_CHECK_EQUAL(utilities::remove_whitespace(text), utilities::remove_whitespace(text1));
}

BOOST_AUTO_TEST_CASE(parse_state_formula_precedence_1801)
{
  const std::string text =
    "<true> exists j: Int. val(j == 0) && val(j == 0)"
    ;

  // This should fail, the right j is not bound by the exist, because the <true> has precedence over the conjunction.
  BOOST_CHECK_THROW(parse_state_formula_specification(text, false), mcrl2::runtime_error);
}