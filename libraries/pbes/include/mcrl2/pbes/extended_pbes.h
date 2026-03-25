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

#include "mcrl2/atermpp/aterm_list.h"
#include <mcrl2/atermpp/aterm_io.h>
#include <string>
#include "mcrl2/core/detail/function_symbols.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/lps/io.h"

namespace mcrl2::pbes_system {

struct extended_pbes {
    /// The core PBES.
    pbes transformed_pbes;

    /// The original PBES that contains counter example information.
    pbes original_pbes;
    /// The original LPS used to generate the PBES.
    lps::stochastic_specification original_lps;

    /// Additional information from transformations on the core PBES.
    atermpp::aterm_list transformations; // TODO we should protect this such that only the right tool can edit its respective transformation info.
};

/// \brief Gets the transformation done by `pbesparelm` and converts it into a mapping that can be used by `pbessolve`.
/// \param transformations aterm list of transformations on an extended pbes
/// \returns a mapping from equation variable names X to indices of redundant parameters
inline
std::unordered_map<std::string, std::set<int>> parelm_info(const atermpp::aterm_list& transformations)
{
    std::unordered_map<std::string, std::set<int>> R = {};
    for(const auto & t : transformations)
    {
        // identify aterm containing parelm-related info by function symbol
        if (t.function() == core::detail::function_symbol_PBESParelmRemoved())
        {
            for (const auto &l: t)
            {
                for (const auto &i: l)
                {
                    if (!i.empty() && i.type_is_list())
                    {
                        auto j = down_cast<atermpp::aterm_list>(i);
                        std::string var_name = pp(j.front()); // name of equation variable
                        auto k = down_cast<atermpp::aterm_list>(j.tail()); // indices of redundant parameters
                        for (auto &m: k)
                        {
                            R[var_name].insert(down_cast<atermpp::aterm_int>(m).value());
                        }
                    }
                }
            }
            break;
        }
    }
    return R;
}

inline
std::string print(const extended_pbes& x, bool precedence_aware)
{
    return pp(x.transformed_pbes, precedence_aware);
}

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_EXTENDED_PBES_H