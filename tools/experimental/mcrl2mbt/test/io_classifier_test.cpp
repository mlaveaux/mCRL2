// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE io_classifier_test
#include <boost/test/included/unit_test.hpp>

#include <sstream>
#include <unistd.h>

// Pull the implementation directly (no separate library target needed).
#include "../io_classifier.cpp" // NOLINT

using namespace mcrl2;

BOOST_AUTO_TEST_CASE(classify_basic)
{
  const std::string content = "input\n"
                              "  button_press\n"
                              "  pin_enter\n"
                              "output\n"
                              "  display_message\n"
                              "  eject_card\n";

  io_classifier cls;
  std::istringstream input(content);
  cls.read(input);

  BOOST_CHECK(cls.classify("button_press") == io_classifier::action_type::input);
  BOOST_CHECK(cls.classify("pin_enter") == io_classifier::action_type::input);
  BOOST_CHECK(cls.classify("display_message") == io_classifier::action_type::output);
  BOOST_CHECK(cls.classify("eject_card") == io_classifier::action_type::output);
  BOOST_CHECK(cls.classify("unknown_action") == io_classifier::action_type::internal);
}

BOOST_AUTO_TEST_CASE(classify_comments_and_blank_lines)
{
  const std::string content = "# This is a classification file\n"
                              "input\n"
                              "  a # inline comment\n"
                              "\n"
                              "output\n"
                              "  b\n";

  std::istringstream input2(content);
  io_classifier cls;
  cls.read(input2);

  BOOST_CHECK(cls.classify("a") == io_classifier::action_type::input);
  BOOST_CHECK(cls.classify("b") == io_classifier::action_type::output);
}

BOOST_AUTO_TEST_CASE(classify_action_before_section_throws)
{
  const std::string content = "stray_action\ninput\n  a\n";
  std::istringstream input3(content);
  io_classifier cls;
  BOOST_CHECK_THROW(cls.read(input3), mcrl2::runtime_error);
}

BOOST_AUTO_TEST_CASE(classify_empty_file)
{
  const std::string content = "";
  std::istringstream input4(content);
  io_classifier cls;
  BOOST_CHECK_NO_THROW(cls.read(input4));
  BOOST_CHECK(cls.classify("anything") == io_classifier::action_type::internal);
}

BOOST_AUTO_TEST_CASE(is_input_is_output_helpers)
{
  const std::string content = "input\n  i\noutput\n  o\n";
  std::istringstream input5(content);
  io_classifier cls;
  cls.read(input5);

  BOOST_CHECK(cls.is_input("i"));
  BOOST_CHECK(!cls.is_input("o"));
  BOOST_CHECK(cls.is_output("o"));
  BOOST_CHECK(!cls.is_output("i"));
}
