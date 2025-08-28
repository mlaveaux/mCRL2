// Author(s): Maurice Laveaux and Menno Bartels
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL_PBES_PBES_SYMMETRY_H
#define MCRL_PBES_PBES_SYMMETRY_H

#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/stategraph.h"
#include "mcrl2/pbes/tools/pbesstategraph_options.h"

namespace mcrl2::pbes_system
{

/// Uses the stategraph algorithm to extract control flow graphs from a given
/// PBES.
class cliques_algorithm: private detail::stategraph_local_algorithm
{
    using super = detail::stategraph_local_algorithm;

public:
    cliques_algorithm(const pbes& input)
        : super(input, pbesstategraph_options {
            .print_influence_graph = true
        })
    {}

    void run() override
    {        
        // Compute the control flow graphs first.
        start_timer("compute_local_control_flow_graphs");
        compute_local_control_flow_graphs();
        finish_timer("compute_local_control_flow_graphs");
        print_local_control_flow_graphs();

        for (int i = 0; i < m_local_control_flow_graphs.size(); i++)
        {
            // For every other control flow graph check if it is compatible.
            std::vector<detail::local_control_flow_graph> I = {m_local_control_flow_graphs[i]};
            for (int j = 0; j < m_local_control_flow_graphs.size(); j++)
            {
                if (i < j)
                {
                    // Property is symmetrical.
                    if (compatible(m_local_control_flow_graphs[i], m_local_control_flow_graphs[j]))
                    {
                        I.emplace_back(m_local_control_flow_graphs[j]);
                    }
                }
            }

            

        }
    }

    bool compatible(const detail::local_control_flow_graph& c, const detail::local_control_flow_graph& c_prime)
    {
        if (c.vertices.size() != c_prime.vertices.size())
        {
            // If V_c != V_C' return false
            return false;
        }      

        for(const auto& s : c.vertices)
        {
            for (const auto& s_prime: c.vertices)
            {

            }

        }

        return true;        
    }

};

/// Contains all the implementation of the PBES symmetry algorithm, based on the article by Bartels et al.
class pbes_symmetry
{
public:
    pbes_symmetry(const pbes& input, const data::rewriter& rewr)
    {
        cliques(input);
    }

private:
    static void cliques(const pbes& input)
    {
        cliques_algorithm algorithm(input);
        algorithm.run();
    }
};

} // namespace mcrl2::pbes_system

#endif // MCRL_PBES_PBES_SYMMETRY_H