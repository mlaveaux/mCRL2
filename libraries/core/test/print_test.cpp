// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file print_test.cpp
/// \brief Add your file description here.

#define BOOST_TEST_MODULE print_test

#include <boost/test/included/unit_test.hpp>
#include "mcrl2/core/print.h"

using namespace mcrl2;

void test_print()
{
  // print directly using core::detail::printer
  std::ostringstream out1;
  core::detail::apply_printer<core::detail::printer> printer(out1, true);
  core::identifier_string s("abc");
  printer.apply(s);
  BOOST_CHECK(out1.str() == "abc");

  // print using the one argument function print
  std::ostringstream out2;
  core::stream_printer()(s, out2, true);
  BOOST_CHECK(out2.str() == "abc");

  // print using the two argument function print
  BOOST_CHECK(core::pp(s) == "abc");

}

BOOST_AUTO_TEST_CASE(test_main)
{
  test_print();
}
