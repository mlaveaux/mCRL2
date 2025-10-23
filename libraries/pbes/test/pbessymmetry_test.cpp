// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#define BOOST_TEST_MODULE pbessymmetry_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/pbes/pbes_symmetry.h"

using namespace mcrl2::pbes_system;

BOOST_AUTO_TEST_CASE(test_permutation)
{
    std::vector<std::size_t> indices = {3, 7, 10};
    std::vector<permutation> perms = permutation_group(indices);
    for (const auto& p: perms)
    {
       std::cout << p << std::endl;
    }

    BOOST_CHECK_EQUAL(perms.size(), 6);
}