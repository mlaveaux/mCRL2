// Author(s): Maurice Laveaux and Menno Bartels
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL_PBES_PBES_QUOTIENT_H
#define MCRL_PBES_PBES_QUOTIENT_H

#include "mcrl2/pbes/pbes_symmetry.h"
#include <boost/process.hpp>

namespace mcrl2::pbes_system {

class pbes_quotient
{
public:
    pbes_quotient(const mcrl2::pbes_system::permutation& pi, std::size_t num_variables)
    {
        boost::asio::io_context ctx;
        boost::process::async_pipe input_pipe(ctx);
        boost::process::async_pipe output_pipe(ctx);
        boost::process::child gap_process(ctx.get_executor(), "gap", 
                        boost::process::std_in < input_pipe,
                        boost::process::std_out > output_pipe);

        // Set the group in gap
        std::stringstream gap_input;
        gap_input << "grp := Group([";

        // Convert permutation to cycle notation
        std::vector<bool> visited(num_variables, false);
        bool first_cycle = true;

        for (size_t i = 0; i < num_variables; ++i) {
            if (!visited[i] && pi[i] != i) {
                if (!first_cycle) {
                    gap_input << ",";
                }
                gap_input << "(";
                
                size_t current = i;
                bool first_element = true;
                do {
                    if (!first_element) {
                        gap_input << ",";
                    }
                    gap_input << (current + 1); // GAP uses 1-based indexing
                    visited[current] = true;
                    current = pi[current];
                    first_element = false;
                } while (current != i);
                
                gap_input << ")";
                first_cycle = false;
            }
        }

        gap_input << "]);\n";

        // Write to GAP process
        boost::asio::write(input_pipe, boost::asio::buffer(gap_input.str()));




    }

private:
};

} // namespace mcrl2::pbes_system

#endif // MCRL_PBES_PBES_QUOTIENT_H