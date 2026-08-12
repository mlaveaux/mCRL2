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

// Write a temporary file containing the given text and return its path.
static std::string write_temp(const std::string& content)
{
  char tmpl[] = "/tmp/mbt_test_XXXXXX";
  const int fd = mkstemp(tmpl);
  BOOST_REQUIRE(fd >= 0);
  const std::string path(tmpl);
  ::write(fd, content.data(), content.size());
  ::close(fd);
  return path;
}

#include <fstream>
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

  const std::string path = write_temp(content);
  io_classifier cls;
  cls.load(path);
  std::remove(path.c_str());

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

  const std::string path = write_temp(content);
  io_classifier cls;
  cls.load(path);
  std::remove(path.c_str());

  BOOST_CHECK(cls.classify("a") == io_classifier::action_type::input);
  BOOST_CHECK(cls.classify("b") == io_classifier::action_type::output);
}

BOOST_AUTO_TEST_CASE(classify_action_before_section_throws)
{
  const std::string content = "stray_action\ninput\n  a\n";
  const std::string path = write_temp(content);
  io_classifier cls;
  BOOST_CHECK_THROW(cls.load(path), mcrl2::runtime_error);
  std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(classify_empty_file)
{
  const std::string path = write_temp("");
  io_classifier cls;
  BOOST_CHECK_NO_THROW(cls.load(path));
  BOOST_CHECK(cls.classify("anything") == io_classifier::action_type::internal);
  std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(is_input_is_output_helpers)
{
  const std::string content = "input\n  i\noutput\n  o\n";
  const std::string path = write_temp(content);
  io_classifier cls;
  cls.load(path);
  std::remove(path.c_str());

  BOOST_CHECK(cls.is_input("i"));
  BOOST_CHECK(!cls.is_input("o"));
  BOOST_CHECK(cls.is_output("o"));
  BOOST_CHECK(!cls.is_output("i"));
}
