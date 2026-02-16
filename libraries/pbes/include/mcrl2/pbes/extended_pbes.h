// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/find_equalities.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_EXTENDED_PBES_H
#define MCRL2_PBES_EXTENDED_PBES_H

#include "mcrl2/pbes/pbes.h"
#include "mcrl2/lps/io.h"

namespace mcrl2::pbes_system {

struct extended_pbes {
    pbes transformed_pbes;

    /// The original PBES that contains counter example information.
    pbes original_pbes;
    /// The original LPS used to generate the PBES.
    lps::stochastic_specification original_lps;
};

inline
std::string print(const extended_pbes& x, bool precedence_aware)
{
    return pp(x.transformed_pbes, precedence_aware);
}

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_EXTENDED_PBES_H