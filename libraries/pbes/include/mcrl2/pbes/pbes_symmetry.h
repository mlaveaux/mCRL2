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
#include "mcrl2/utilities/logger.h"

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
        // Does too much, but okay.
        super::run();

        for (auto i = m_local_control_flow_graphs.begin(); i != m_local_control_flow_graphs.end(); ++i)
        {
            mCRL2log(log::verbose) << "--- computed local control flow graph " << (i - m_local_control_flow_graphs.begin()) << "\n" << *i << std::endl;
        }

        std::vector<std::vector<std::size_t>> cal_I = cliques();
    }

    std::vector<std::vector<std::size_t>> cliques()
    {
        std::vector<std::vector<std::size_t>> cal_I;
        for (int i = 0; i < m_local_control_flow_graphs.size(); i++)
        {
            // For every other control flow graph check if it is compatible.
            std::vector<std::size_t> I = {static_cast<unsigned long>(i)};
            for (int j = 0; j < m_local_control_flow_graphs.size(); j++)
            {
                if (i < j)
                {
                    // Property is symmetrical.
                    if (compatible(i, j))
                    {
                        I.emplace_back(j);
                    }
                }
            }

            mCRL2log(log::verbose) << "--- control flow graph in clique \n";
            for (const auto& graph : I)
            {
                mCRL2log(log::verbose) << graph << std::endl;
            }

            cal_I.emplace_back(I);
        }
    
        return cal_I;
    }

    /// Checks whether two control flow graphs are compatible according to Algorithm 4.
    bool compatible(int i, int j)
    {
        const detail::local_control_flow_graph& c = m_local_control_flow_graphs[i];
        const detail::local_control_flow_graph& c_prime = m_local_control_flow_graphs[j];
        mCRL2log(log::debug) << "Checking compatible(" << i << ", " << j << ")" << std::endl;


        if (c.vertices != c_prime.vertices)
        {
            // If V_c != V_C' return false
            mCRL2log(log::debug) << "Vertex sets don't match" << std::endl;
            return false;
        }      

        for(const auto& s : c.vertices)
        {
            // There exist t such that s and t match according to the definitions in the paper.
            for (const auto& s_prime: c_prime.vertices)
            {
                // X(v) in c and X(v) in c_prime.
                if (s.value() == s_prime.value() && s.variable() == s_prime.variable())
                {
                    mCRL2log(log::debug) << "Comparing vertices s = " << s << " and s'=  " << s_prime << std::endl;
                    for (const auto& [t, e] : s.outgoing_edges())
                    {
                        for (const auto& [t_prime, f] : s_prime.outgoing_edges())
                        {
                            // sizes
                            if (e.size() == f.size())
                            {
                                mCRL2log(log::debug) << "Found different number of edges " << e.size() << " and " << f.size() << std::endl;
                                return false;
                            }

                            // #
                            for (const auto& eqn_X:  m_pbes.equations())
                            {
                                for (const auto& Y : eqn_X.predicate_variables())
                                {

                                }
                            }

                        }
                    }
                }

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
        srf_pbes srf = pbes2srf(input);
        cliques(srf.to_pbes());
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