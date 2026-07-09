// Author(s): mCRL2 developers
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file review_regression_test.cpp
/// \brief Regression tests for defects found during the atermpp code review.
///
/// Each test below guards against the reintroduction of a specific, previously
/// confirmed defect in the atermpp library headers/sources.

#define BOOST_TEST_MODULE review_regression_test
#include <boost/test/included/unit_test.hpp>

#include <list>
#include <type_traits>

#include "mcrl2/atermpp/aterm_int.h"
#include "mcrl2/atermpp/aterm_io.h"
#include "mcrl2/atermpp/aterm_list.h"

using namespace atermpp;

// Regression: term_appl_iterator::operator--(int) used to return *this (the
// already-decremented iterator) instead of the saved pre-decrement value,
// violating the postfix decrement contract.
BOOST_AUTO_TEST_CASE(test_term_appl_iterator_postfix_decrement)
{
  const aterm t = read_term_from_string("f(a, b, c)");
  const aterm a = read_term_from_string("a");
  const aterm b = read_term_from_string("b");
  const aterm c = read_term_from_string("c");

  // Position an iterator at the last argument (index 2, i.e. c).
  auto it = t.begin();
  ++it;
  ++it;
  BOOST_CHECK(*it == c);

  // Postfix decrement must return the value *before* decrementing.
  const auto previous = it--;
  BOOST_CHECK(*previous == c); // pre-decrement value
  BOOST_CHECK(*it == b);       // iterator itself moved back

  // Sanity check that postfix increment behaves symmetrically.
  const auto previous2 = it++;
  BOOST_CHECK(*previous2 == b);
  BOOST_CHECK(*it == c);
  (void)a;
}

// Regression: the filtering forward constructor of term_list constrained on
// !random_access, which also matched bidirectional iterators. That made the
// call ambiguous with the bidirectional (backward) overload for e.g.
// std::list iterators. It must now compile and select the backward overload.
BOOST_AUTO_TEST_CASE(test_term_list_bidirectional_filter_constructor)
{
  const aterm a = read_term_from_string("a");
  const aterm b = read_term_from_string("b");
  const aterm c = read_term_from_string("c");
  const std::list<aterm> source = { a, b, c };

  auto convert = [](const aterm& x) { return x; };
  auto keep = [](const aterm&) { return true; };

  const term_list<aterm> result(source.begin(), source.end(), convert, keep);
  BOOST_CHECK_EQUAL(result.size(), 3u);
  BOOST_CHECK(result.front() == a);
}

// Regression: the free function make_term_list had the same non-exclusive
// constraint on its filtering forward overload, causing ambiguity for
// bidirectional iterators.
BOOST_AUTO_TEST_CASE(test_make_term_list_bidirectional_filter)
{
  const aterm a = read_term_from_string("a");
  const aterm b = read_term_from_string("b");
  const std::list<aterm> source = { a, b };

  auto convert = [](const aterm& x) { return x; };
  auto keep = [](const aterm&) { return true; };

  term_list<aterm> result;
  make_term_list(result, source.begin(), source.end(), convert, keep);
  BOOST_CHECK_EQUAL(result.size(), 2u);
  BOOST_CHECK(result.front() == a);
}

// Regression: aterm_stream_state is a scope guard that restores stream state in
// its destructor. Copying/moving it would restore the state more than once, so
// the copy and move operations must be deleted.
BOOST_AUTO_TEST_CASE(test_aterm_stream_state_is_non_copyable)
{
  static_assert(!std::is_copy_constructible_v<aterm_stream_state>,
      "aterm_stream_state must not be copy constructible");
  static_assert(!std::is_move_constructible_v<aterm_stream_state>,
      "aterm_stream_state must not be move constructible");
  static_assert(!std::is_copy_assignable_v<aterm_stream_state>,
      "aterm_stream_state must not be copy assignable");
  static_assert(!std::is_move_assignable_v<aterm_stream_state>,
      "aterm_stream_state must not be move assignable");
  BOOST_CHECK(true);
}

// Regression: aterm_int(const aterm&) asserted "type_is_int() || !defined()".
// Because type_is_int() dereferences the internal pointer, it crashed on a
// default-constructed (undefined) aterm before the short-circuit could help.
// The precondition is now expressed as "!defined() || type_is_int()".
BOOST_AUTO_TEST_CASE(test_aterm_int_from_undefined_aterm)
{
  const aterm_int i{aterm()};
  BOOST_CHECK(!i.defined());
}

// Regression: parse_aterm_int wrote the terminating '\0' one past the end of a
// fixed 32-byte buffer when the digit sequence filled the buffer, corrupting
// the stack. Parsing an overly long integer must no longer crash.
BOOST_AUTO_TEST_CASE(test_parse_long_integer_does_not_overflow)
{
  const std::string many_digits(40, '9');
  BOOST_CHECK_NO_THROW(read_int_from_string(many_digits));
}

// Regression sanity check: list parsing must continue to work after tightening
// the identifier terminator set in parse_unquoted_string (which previously
// checked ']' twice instead of also checking '[').
BOOST_AUTO_TEST_CASE(test_list_round_trip)
{
  const aterm original = read_term_from_string("f([1, 2], 3)");
  std::ostringstream out;
  write_term_to_text_stream(original, out);
  std::istringstream in(out.str());
  aterm restored;
  read_term_from_text_stream(in, restored);
  BOOST_CHECK(original == restored);
}
