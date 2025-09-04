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

#include "mcrl2/core/detail/print_utility.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/stategraph.h"
#include "mcrl2/pbes/tools/pbesstategraph_options.h"
#include "mcrl2/pbes/unify_parameters.h"
#include "mcrl2/utilities/logger.h"
#include <cstddef>

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

    /// Given a clique return all relevant data parameters.
    std::set<std::size_t> data_parameters(const std::vector<size_t>& clique)
    {
        std::set<std::size_t> data_parameters;
        for (const auto& i : clique)
        {        
            const detail::local_control_flow_graph& c = m_local_control_flow_graphs[i];
            for(const auto& s : c.vertices)
            {
                // Compute the data parameters
                // Get the changed by, used for and used in
                auto it = s.outgoing_edges().find(&s);

                if (it != s.outgoing_edges().end())
                {
                    for (const std::size_t& label : it->second)
                    {
                        for (const auto& variable : m_pbes.equations().at(label).predicate_variables())
                        {
                            if (variable.name() == s.name())
                            {
                                data_parameters.insert(variable.changed().begin(), variable.changed().end());
                                data_parameters.insert(variable.used().begin(), variable.used().end());
                            }
                        }
                    }
                }
            }
        }
        
        return data_parameters;
    }

    /// Computes the set of candidates we can derive from a single clique
    void clique_candidates(const std::vector<size_t>& clique)
    {
        auto D = data_parameters(clique);

        // auto clique_permutations = permutation_group(clique);
    }

    /// Determine the cliques of the control flow graphs.
    std::vector<std::vector<std::size_t>> cliques()
    {
        std::vector<std::vector<std::size_t>> cal_I;
        for (int i = 0; i < m_local_control_flow_graphs.size(); i++)
        {
            if (std::any_of(cal_I.begin(), cal_I.end(), 
                [i](const auto& clique) {
                    return std::find(clique.begin(), clique.end(), i) != clique.end();
                }))
            {
                // Skip every graph that already belongs to a clique.
                continue;
            }

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

    /// Computes the sizes(c, s, s')
    std::set<std::pair<size_t, size_t>> sizes(const detail::local_control_flow_graph&, const detail::local_control_flow_graph_vertex& s, const detail::local_control_flow_graph_vertex& s_prime)
    {
        // Get the changed by, used for and used in
        auto it = s.outgoing_edges().find(&s_prime);

        std::set<std::pair<size_t, size_t>> result;
        if (it != s.outgoing_edges().end())
        {
            // Find the corresponding equation
            for (const auto& equation: m_pbes.equations())
            {
                if (equation.variable().name() == s.name())
                {
                    for (const std::size_t& label : it->second)
                    {
                        // Compute the sizes.
                        const auto& variable = equation.predicate_variables().at(label);
                        result.insert(std::make_pair(variable.changed().size(), variable.used().size()));
                    }
                }
            }
        }

        return result;
    }

    /// Checks whether two control flow graphs are compatible according to Algorithm 4.
    bool compatible(int i, int j)
    {
        const detail::local_control_flow_graph& c = m_local_control_flow_graphs[i];
        const detail::local_control_flow_graph& c_prime = m_local_control_flow_graphs[j];
        mCRL2log(log::debug) << "Checking compatible(" << i << ", " << j << ")" << std::endl;

        if (!vertex_sets_compatible(c, c_prime))
        {
            // If V_c != V_C' return false
            mCRL2log(log::debug) << "Vertex sets don't match" << std::endl;
            return false;
        }

        // Note that this algorithm is slightly different than the pseudocode, because the graphs in the implementation are over different (compatible) vertex sets.
        for(const auto& s : c.vertices)
        {
            // There exist t such that s and t match according to the definitions in the paper.
            for (const auto& s_c_prime: c_prime.vertices)
            {
                // X(v) in c and X(v) in c_prime.
                if (s.value() == s_c_prime.value() && s.name() == s_c_prime.name())
                {
                    for (const auto& s_prime: c.vertices)
                    {   
                        // There exist t such that s and t match according to the definitions in the paper.
                        for (const auto& s_prime_c_prime: c_prime.vertices)
                        {
                            // Y(v) in c and Y(v) in c_prime.
                            if (s_prime.value() == s_prime_c_prime.value() && s_prime.name() == s_prime_c_prime.name())
                            {
                                mCRL2log(log::trace) << "Comparing vertices s = " << s << " and s'= " << s_prime << std::endl;
                                auto it = s.outgoing_edges().find(&s_prime);
                                auto it_c_prime = s_c_prime.outgoing_edges().find(&s_prime_c_prime);
                                
                                if ((it == s.outgoing_edges().end()) != (it_c_prime == s_c_prime.outgoing_edges().end()))
                                {
                                    mCRL2log(log::trace) << "Found different number of edges " << s << " and " << s_prime << std::endl;
                                    return false;
                                }

                                if (it != s.outgoing_edges().end() && 
                                    it_c_prime != s.outgoing_edges().end() && 
                                    it->second.size() != it_c_prime->second.size())
                                {
                                    mCRL2log(log::trace) << "Found different number of edges " << it->second.size() << " and " << it_c_prime->second.size() << std::endl;
                                    return false;
                                }

                                if (sizes(c, s, s_prime) != sizes(c_prime, s_c_prime, s_prime_c_prime))
                                {
                                    mCRL2log(log::debug) << "Found different sizes " << core::detail::print_container(sizes(c, s, s_prime))
                                         << " and " << core::detail::print_container(sizes(c_prime, s_c_prime, s_prime_c_prime)) << std::endl;
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }

        return true;        
    }

    /// Checks whether two control flow graphs have compatible vertex sets, meaning that the PVI and values of the vertices match.
    bool vertex_sets_compatible(const detail::local_control_flow_graph& c, const detail::local_control_flow_graph& c_prime)
    {
        if (c.vertices.size() != c_prime.vertices.size())
        {
            return false;
        }
        
        for (const auto& vertex : c.vertices)
        {
            if (!std::any_of(c_prime.vertices.begin(), c_prime.vertices.end(), 
                [&vertex](const auto& vertex_prime) {
                    return vertex.name() == vertex_prime.name() && vertex.value() == vertex_prime.value();
                }))
            {
                return false;
            }
        }
        
        return true;
    }

};

/// Contains all the implementation of the PBES symmetry algorithm, based on the article by Bartels et al.
class pbes_symmetry
{
public:
    pbes_symmetry(const pbes& input, const data::rewriter&)
    {
        srf_pbes srf = pbes2srf(input);

        mCRL2log(mcrl2::log::debug) << srf.to_pbes() << std::endl;

        unify_parameters(srf, false, false);

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