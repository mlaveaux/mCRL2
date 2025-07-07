// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file concepts.cpp
/// \brief Demonstrates the usage of the algorithms library.

#include "mcrl2/atermpp/function_symbol.h"
#define BOOST_TEST_MODULE aterm_concepts_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/atermpp/concepts.h"
#include "mcrl2/atermpp/aterm.h"

using namespace atermpp;

class derived : aterm {

};

class derived_with_extra : aterm {
  [[maybe_unused]]
  std::size_t index = 0;
};

static_assert(IsATerm<derived>, "This class remains an aterm since it derives from it and has no additional fields.");
static_assert(!IsATerm<derived_with_extra>, "This class has extra fields that are not allowed");


std::size_t not_term_converter(const atermpp::unprotected_aterm& term) {
  return term.function().arity();
}


aterm wrong_argument_convert(const std::size_t index) {
  return aterm(function_symbol("test", index));
}

aterm term_converter(const aterm& input) {
  return *input.begin();
}

// static_assert(!IsTermConverter<decltype(not_term_converter)>, "This function is a not a term converter since the argument is not a term");
// static_assert(!IsTermConverter<decltype(wrong_argument_convert)>, "This function has the wrong arguments");
static_assert(IsTermConverter<decltype(term_converter)>, "This function has the wrong arguments");

BOOST_AUTO_TEST_CASE(test_static)
{

}