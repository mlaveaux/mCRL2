// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file parser_test.cpp
/// \brief Regression test for parsing a data specification.

#define BOOST_TEST_MODULE parser_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/data/parse.h"

using namespace mcrl2;

void parser_test()
{
  const std::string text(
    "sort S;\n"
    "cons s:S;\n"
    "map f:S -> List(S);\n"
  );

  data::data_specification spec(data::parse_data_specification(text));

  std::cerr << "number of sorts " << spec.sorts().size() << "\n";
  std::cerr << "number of functions " << spec.mappings().size() << "\n";
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  BOOST_CHECK(spec.sorts().size() == 8); // Bool, Pos, S, List(S), S->List(S), Nat, @NatNatPair.
  BOOST_CHECK(spec.mappings().size() == 204);
#else
  BOOST_CHECK(spec.sorts().size() == 7); // Bool, Pos, S, List(S), S->List(S), Nat, @NatPair.
  BOOST_CHECK(spec.mappings().size() == 114);
#endif
  BOOST_CHECK(spec.constructors(data::basic_sort("S")).size() == 1);

  BOOST_CHECK(data::parse_data_expression("2") == data::sort_pos::pos(2));
  BOOST_CHECK(data::parse_data_expression("0") == data::sort_nat::nat(0));
  BOOST_CHECK(data::parse_sort_expression("Nat") == data::sort_nat::nat());
}

// This test triggers a sort normalization problem.
void test_user_defined_sort()
{
  using namespace data;

  std::string text = "sort D = struct d1 | d2;\n";
  data_specification data_spec = parse_data_specification(text);
  sort_expression s = parse_sort_expression("D", data_spec);
}

void test_whr()
{
  using namespace data;
  data_expression x = parse_data_expression("exists n: Pos . n == 0 whr n = 1 end");
}

// returns true if parsing succeeded
bool parse_sort(const data::data_specification& data_spec, const std::string& text)
{
  try
  {
    data::sort_expression s = data::parse_sort_expression(text, data_spec);
  }
  catch (mcrl2::runtime_error&)
  {
    return false;
  }
  return true;
}

// if expected_result is false, parsing is expected to fail
void test_sort_expression(const data::data_specification& data_spec, const std::string& text, bool expected_result)
{
  bool result = parse_sort(data_spec, text);
  if (result != expected_result)
  {
    std::cout << "ERROR: parsing the sort expression '" << text << "' produced the result " << std::boolalpha << result << std::endl;
  }
  BOOST_CHECK(result == expected_result);
}

void test_ticket_1267()
{
  data::data_specification data_spec = data::parse_data_specification("sort A;");
  test_sort_expression(data_spec, "A -> A # A", false);
  test_sort_expression(data_spec, "(A -> A) # A", false);
  test_sort_expression(data_spec, "A # A", false);
  test_sort_expression(data_spec, "(A # A)", false);

  test_sort_expression(data_spec, "(A # A) -> A", true);
  test_sort_expression(data_spec, "A -> ((A # A) -> (A -> A))", true);
  test_sort_expression(data_spec, "A", true);
  test_sort_expression(data_spec, "A # A -> A", true);
  test_sort_expression(data_spec, "A -> A # A -> A -> A", true);
  test_sort_expression(data_spec, "A -> (A -> A)", true);
  test_sort_expression(data_spec, "(A -> A) -> A", true);
  test_sort_expression(data_spec, "A # A -> (A -> A)", true);
  test_sort_expression(data_spec, "A # (A -> A) -> (A -> A)", true);
  test_sort_expression(data_spec, "(A -> A) # A -> (A -> A)", true);
  test_sort_expression(data_spec, "A # A # A -> A", true);
  test_sort_expression(data_spec, "A -> A -> A -> A", true);
}

BOOST_AUTO_TEST_CASE(test_main)
{
  parser_test();
  test_user_defined_sort();
  test_whr();
  test_ticket_1267();
}
