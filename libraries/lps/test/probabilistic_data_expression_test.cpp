// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file probabilistic_data_expression_test.cpp
/// \brief Test whether calculations with fractions go well. 

#define BOOST_TEST_MODULE probabilistic_data_expression_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/lps/probabilistic_data_expression.h"

using namespace mcrl2;
using namespace mcrl2::data;
using namespace mcrl2::data::sort_real;
using namespace mcrl2::data::sort_nat;
using namespace mcrl2::data::sort_pos;
using namespace mcrl2::data::sort_int;
using namespace mcrl2::lps;

BOOST_AUTO_TEST_CASE(test_very_few_digits)
{
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  probabilistic_data_expression r1(real_("79","842"));
  probabilistic_data_expression r2(real_("6398","87493"));
#else
  // Without machine numbers, this test was too slow.
  probabilistic_data_expression r1(real_("3","4"));
  probabilistic_data_expression r2(real_("1","5"));
#endif
  probabilistic_data_expression zero;
  BOOST_CHECK(r1>zero);
  BOOST_CHECK(r1+r2>zero);
  BOOST_CHECK(r1-r2>zero);
  BOOST_CHECK(r2>zero);
  BOOST_CHECK(r2<r1);
  BOOST_CHECK(r2<r1+r2);
}

#ifdef MCRL2_ENABLE_MACHINENUMBERS
BOOST_AUTO_TEST_CASE(test_few_digits)
{
  probabilistic_data_expression r1(real_("7979","84298"));
  probabilistic_data_expression r2(real_("63984264","8749379236"));
  probabilistic_data_expression zero;
  BOOST_CHECK(r1>zero);
  BOOST_CHECK(r1+r2>zero);
  BOOST_CHECK(r1-r2>zero);
  BOOST_CHECK(r2>zero);
  BOOST_CHECK(r2<r1);
  BOOST_CHECK(r2<r1+r2);
}

BOOST_AUTO_TEST_CASE(test_many_digits)
{
  probabilistic_data_expression r1(real_("798797978978797979","1348769876316786867868684298"));
  probabilistic_data_expression r2(real_("883298364493924693","9379236492843693629438693824"));
  probabilistic_data_expression zero;
  BOOST_CHECK(r1>zero);
  BOOST_CHECK(r1+r2>zero);
  BOOST_CHECK(r1-r2>zero);
  BOOST_CHECK(r2-r1<zero);
  BOOST_CHECK(r2>zero);
  BOOST_CHECK(r2<r1);
  BOOST_CHECK(r2<r1+r2);
}

BOOST_AUTO_TEST_CASE(test_digits_many_vs_few)
{
  probabilistic_data_expression r1(real_("7979","2989999999999999999999999999"));
  probabilistic_data_expression r2(real_("134","93824"));
  probabilistic_data_expression zero;
  BOOST_CHECK(r1>zero);
  BOOST_CHECK(r1+r2>zero);
  BOOST_CHECK(r1-r2<zero);
  BOOST_CHECK(r2>zero);
  BOOST_CHECK(r2>r1);
  BOOST_CHECK(r2<r1+r2);
}

BOOST_AUTO_TEST_CASE(test_many_digits_few_versus_many)
{
  probabilistic_data_expression r1(real_("7979","898798"));
  probabilistic_data_expression r2(real_("3298364493924693","379236492843693629438693824"));
  probabilistic_data_expression zero;
  BOOST_CHECK(r1>zero);
  BOOST_CHECK(r1+r2>zero);
  BOOST_CHECK(r1-r2>zero);
  BOOST_CHECK(r2>zero);
  BOOST_CHECK(r2<r1);
  BOOST_CHECK(r2<r1+r2);
}
#endif

  
