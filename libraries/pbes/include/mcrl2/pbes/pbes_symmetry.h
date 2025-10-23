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
#include "mcrl2/pbes/detail/stategraph_algorithm.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/stategraph.h"
#include "mcrl2/pbes/tools/pbesstategraph_options.h"
#include "mcrl2/pbes/unify_parameters.h"
#include "mcrl2/utilities/logger.h"
#include <cstddef>
#include <unordered_map>
#include <unordered_set>

namespace mcrl2::pbes_system
{

/// A representation of a permutation.
class permutation
{
public:
  permutation() = default;

  permutation(const std::unordered_map<std::size_t, std::size_t>& mapping)
    : m_mapping(mapping)
  {}

  std::unordered_map<std::size_t, std::size_t> mapping() const
  {
    return m_mapping;
  }

  std::size_t operator[](std::size_t i) const {
    auto it = m_mapping.find(i);
    if (it != m_mapping.end())
    {
      return it->second;
    }
    return i;
  }

  // Applies the permutation to a set of indices.
  std::set<std::size_t> permute(const std::set<std::size_t>& s) const
  {
    std::set<std::size_t> result;
    for (const auto& i: s)
    {
      result.insert((*this)[i]);
    }
    return result;
  }

  /// Returns the concatenation of this permutation with another permutation.
  permutation concat(const permutation& other) const
  {
    std::unordered_map<std::size_t, std::size_t> new_mapping(m_mapping.size());

    for (const auto& [key, value]: m_mapping)
    {
      new_mapping[key] = other[value];
    }

    for (const auto& [key, value]: other.m_mapping)
    {
      assert (m_mapping.find(key) == m_mapping.end());
      new_mapping[key] = value;
    }

    return permutation(new_mapping);
  }

private:
  std::unordered_map<std::size_t, std::size_t> m_mapping;
};

/// Returns all the permutations for the given indices.
inline std::vector<permutation> permutation_group(const std::vector<std::size_t>& indices)
{
  std::vector<permutation> result;

  // Recursive function to generate permutations using backtracking.
  std::function<void(std::unordered_map<std::size_t, std::size_t>&, std::vector<std::size_t>&, std::size_t)> generate_permutations
    = [&](std::unordered_map<std::size_t, std::size_t>& mapping, std::vector<std::size_t>& available, std::size_t pos)
  {
    if (pos == indices.size())
    {
      result.emplace_back(mapping);
      return;
    }

    for (std::size_t i = 0; i < available.size(); ++i)
    {
      std::size_t target = available[i];
      mapping[indices[pos]] = target;
      
      available.erase(available.begin() + i);
      generate_permutations(mapping, available, pos + 1);
      available.insert(available.begin() + i, target); // backtrack
      
      mapping.erase(indices[pos]);
    }
  };

  std::unordered_map<std::size_t, std::size_t> mapping;
  std::vector<std::size_t> available_indices = indices;
  generate_permutations(mapping, available_indices, 0);

  return result;
}

/// Prints the permutation in cycle notation.
inline std::ostream& operator<<(std::ostream& out, const permutation& p)
{
  out << "[";
  bool first = true;
  for (const auto& [key, value]: p.mapping())
  {
    if (!first)
    {
      out << ", ";
    }

    out << key << " -> " << value;
    first = false;
  }
  out << "]";

  return out;
}

/// Uses the stategraph algorithm to extract control flow graphs from a given
/// PBES.
class cliques_algorithm : private detail::stategraph_local_algorithm
{
  using super = detail::stategraph_local_algorithm;

public:
  cliques_algorithm(const pbes& input)
    : super(input, pbesstategraph_options{.print_influence_graph = true})
  {}

  void run() override
  {
    // We explicitly ignore the virtual call to run in the base class
    detail::stategraph_algorithm::stategraph_algorithm::run(); // NOLINT(bugprone-parent-virtual-call)

    compute_local_control_flow_graphs();

    for (auto i = m_local_control_flow_graphs.begin(); i != m_local_control_flow_graphs.end(); ++i)
    {
      mCRL2log(log::verbose) << "--- computed local control flow graph " << (i - m_local_control_flow_graphs.begin())
        << "\n"
        << *i << std::endl;
    }

    std::vector<std::vector<std::size_t>> cal_I = cliques();

    for (const auto& clique: cal_I)
    {
      auto candidates = clique_candidates(clique);
    }
  }

  /// Takes as input a clique of compatible control flow parameters and return
  /// the set of all data parameters that somehow play a role for any of these
  /// parameters.
  std::set<std::size_t> data_parameters(const std::vector<size_t>& clique)
  {
    std::set<std::size_t> data_parameters;
    for (const auto& i: clique)
    {
      const detail::local_control_flow_graph& c = m_local_control_flow_graphs[i];
      for (const auto& s: c.vertices)
      {
        // Compute the data parameters
        // Get the changed by, used for and used in
        auto it = s.outgoing_edges().find(&s);

        if (it != s.outgoing_edges().end())
        {
          for (const std::size_t& label: it->second)
          {
            for (const auto& variable: m_pbes.equations().at(label).predicate_variables())
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
  std::vector<permutation> clique_candidates(const std::vector<size_t>& I)
  {
    auto D = data_parameters(I);

    std::vector<permutation> result;
    for (const auto& alpha: permutation_group(I))
    {
      for (const auto& beta: permutation_group(std::vector<std::size_t>(D.begin(), D.end())))
      {
        mCRL2log(log::verbose) << "Trying candidate: " << alpha << " and " << beta << std::endl;

        permutation pi = alpha.concat(beta);
        if (complies(pi, I))
        {
          result.emplace_back(std::move(pi));
        }
      }
    }

    
    mCRL2log(log::verbose) << "--- compliant permutations \n";
    for (const auto& p: result)
    {
       mCRL2log(log::verbose) << p << std::endl;
    }

    return result;
  }

  /// Determine the cliques of the control flow graphs.
  std::vector<std::vector<std::size_t>> cliques()
  {
    std::vector<std::vector<std::size_t>> cal_I;
    for (int i = 0; i < m_local_control_flow_graphs.size(); i++)
    {
      if (std::any_of(cal_I.begin(),
            cal_I.end(),
            [i](const auto& clique) { return std::find(clique.begin(), clique.end(), i) != clique.end(); }))
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

      if (I.size() > 1)
      {
        mCRL2log(log::verbose) << "--- control flow graph in clique \n";
        for (const auto& graph: I)
        {
            mCRL2log(log::verbose) << graph << std::endl;
        }
        cal_I.emplace_back(I);
      }
    }

    return cal_I;
  }

  /// Returns true iff all vertices in I comply with the permutation pi.
  bool complies(const permutation& pi, const std::vector<std::size_t>& I)
  {
    return std::all_of(I.begin(), I.end(), [&](std::size_t c) { return complies(pi, c); });
  }

  /// Takes a permutation and a control flow parameter and returns true or
  /// false depending on whether the permutation complies with the control
  /// flow parameter according to Definition
  bool complies(const permutation& pi, std::size_t c) const
  {
    const detail::local_control_flow_graph& graph = m_local_control_flow_graphs[c];
    const detail::local_control_flow_graph& other_graph = m_local_control_flow_graphs[pi[c]];

    // TODO: Is this equivalent to the bijection check in the paper.
    for (const auto& s: graph.vertices)
    {
      for (const auto& s_prime: other_graph.vertices)
      {
        if (s.value() == s_prime.value() && s.name() == s_prime.name())
        {
          // s == s'
          for (const auto& [to, labels]: s.outgoing_edges())
          {
            for (const auto& [to_prime, labels_prime]: s_prime.outgoing_edges())
            {
              if (to->value() == to_prime->value() && to->name() == to_prime->name())
              {
                // t == t'
                // Find the corresponding equation
                for (const auto& equation: m_pbes.equations())
                {
                  if (equation.variable().name() == s.name())
                  {
                    // For each i find a corresponding j.
                    std::set<std::size_t> remaining_j = labels_prime;
                    for (const std::size_t& i: labels)
                    {
                      const auto& variable = equation.predicate_variables().at(i);
                      std::optional<std::size_t> matching_j;
                      for (const std::size_t& j: remaining_j)
                      {
                        const auto& variable_prime = equation.predicate_variables().at(j);
                        if (pi.permute(variable.changed()) == variable_prime.changed()
                            && pi.permute(variable.used()) == variable_prime.used())
                        {
                          matching_j = j;
                          break;
                        }
                      }

                      if (matching_j)
                      {
                        // Found a matching j for i.
                        remaining_j.erase(*matching_j);
                      }
                    }

                    if (!remaining_j.empty())
                    {
                      mCRL2log(log::debug) << "No matching found for edge from " << s << " to " << *to
                                           << " under permutation " << pi << std::endl;
                      return false;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    return true;
  }

  /// Computes the sizes(c, s, s')
  std::set<std::pair<size_t, size_t>> sizes(const detail::local_control_flow_graph&,
    const detail::local_control_flow_graph_vertex& s,
    const detail::local_control_flow_graph_vertex& s_prime) const
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
          for (const std::size_t& label: it->second)
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
  bool compatible(int i, int j) const
  {
    const detail::local_control_flow_graph& c = m_local_control_flow_graphs[i];
    const detail::local_control_flow_graph& c_prime = m_local_control_flow_graphs[j];
    mCRL2log(log::trace) << "Checking compatible(" << i << ", " << j << ")" << std::endl;

    if (!vertex_sets_compatible(c, c_prime))
    {
      // If V_c != V_C' return false
      mCRL2log(log::trace) << "Vertex sets don't match" << std::endl;
      return false;
    }

    // Note that this algorithm is slightly different than the pseudocode, because the graphs in the implementation are
    // over different (compatible) vertex sets.
    for (const auto& s: c.vertices)
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

                if (it != s.outgoing_edges().end() && it_c_prime != s.outgoing_edges().end()
                    && it->second.size() != it_c_prime->second.size())
                {
                  mCRL2log(log::trace) << "Found different number of edges " << it->second.size() << " and "
                                       << it_c_prime->second.size() << std::endl;
                  return false;
                }

                if (sizes(c, s, s_prime) != sizes(c_prime, s_c_prime, s_prime_c_prime))
                {
                  mCRL2log(log::trace) << "Found different sizes "
                                       << core::detail::print_container(sizes(c, s, s_prime)) << " and "
                                       << core::detail::print_container(sizes(c_prime, s_c_prime, s_prime_c_prime))
                                       << std::endl;
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

  /// Checks whether two control flow graphs have compatible vertex sets, meaning that the PVI and values of the
  /// vertices match.
  bool vertex_sets_compatible(const detail::local_control_flow_graph& c,
    const detail::local_control_flow_graph& c_prime) const
  {
    if (c.vertices.size() != c_prime.vertices.size())
    {
      mCRL2log(log::trace) << "Different number of vertices: " << c.vertices.size() << " and "
                           << c_prime.vertices.size() << std::endl;
      return false;
    }

    for (const auto& vertex: c.vertices)
    {
      if (!std::any_of(c_prime.vertices.begin(),
            c_prime.vertices.end(),
            [&vertex](const auto& vertex_prime)
            { return vertex.name() == vertex_prime.name() && vertex.value() == vertex_prime.value(); }))
      {
        mCRL2log(log::trace) << "Vertex " << vertex << " does not occur in the right hand side control flow graph"
                             << std::endl;
        return false;
      }
    }

    for (const auto& vertex_prime: c_prime.vertices)
    {
      if (!std::any_of(c.vertices.begin(),
            c.vertices.end(),
            [&vertex_prime](const auto& vertex)
            { return vertex.name() == vertex_prime.name() && vertex.value() == vertex_prime.value(); }))
      {
        mCRL2log(log::trace) << "Vertex " << vertex_prime << " does not occur in the left hand side control flow graph"
                             << std::endl;
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