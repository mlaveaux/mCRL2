// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbes_remove_counterexample_info_test.cpp
/// \brief Tests for the recognition of counter example variables.

#define BOOST_TEST_MODULE pbes_remove_counterexample_info_test
#include <boost/test/included/unit_test.hpp>
#include "mcrl2/pbes/detail/pbes_remove_counterexample_info.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system::detail;

/// The names that lps2pbes -c generates for the counter example equations.
static const std::array<const char*, 6> counter_example_names = {
  "Zpos_0_r1", "Zpos_12_c6", "Zneg_0_r1", "Zneg_12_c6", "Zpos_3_s4", "Zneg_3_s4"};

/// Names that must not be mistaken for counter example equations.
static const std::array<const char*, 8> other_names
  = {"X", "Y", "Z", "Zpos", "Zneg", "ZposX", "X_Zpos_0_r1", "Zprime_0_r1"};

BOOST_AUTO_TEST_CASE(test_is_counter_example_name)
{
  for (const char* name: counter_example_names)
  {
    BOOST_CHECK(is_counter_example_name(core::identifier_string(name)));
  }
  for (const char* name: other_names)
  {
    BOOST_CHECK(!is_counter_example_name(core::identifier_string(name)));
  }
}

/// Names that carry the prefix but not the digit block that the regular expression requires. A
/// prefix test on Zpos_/Zneg_ would accept these, which is why pbessolvesymbolic uses the regular
/// expression rather than such a test: a user is free to introduce a variable with that name.
BOOST_AUTO_TEST_CASE(test_prefix_alone_is_not_enough)
{
  for (const char* name: {"Zpos_x", "Zneg_", "Zpos_", "Zpos_a_b"})
  {
    BOOST_CHECK(!is_counter_example_name(core::identifier_string(name)));
  }
}

BOOST_AUTO_TEST_CASE(test_short_names)
{
  for (const char* name: {"", "Z", "Zp", "Zpo", "Zpos"})
  {
    BOOST_CHECK(!is_counter_example_name(core::identifier_string(name)));
  }
}

/// is_counter_example_positive distinguishes the two kinds.
BOOST_AUTO_TEST_CASE(test_is_counter_example_positive)
{
  BOOST_CHECK(is_counter_example_positive(core::identifier_string("Zpos_0_r1")));
  BOOST_CHECK(!is_counter_example_positive(core::identifier_string("Zneg_0_r1")));
  BOOST_CHECK(!is_counter_example_positive(core::identifier_string("X")));
}
