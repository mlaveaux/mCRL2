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

#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/core/detail/print_utility.h"
#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/replace.h"
#include "mcrl2/data/substitutions/mutable_map_substitution.h"
#include "mcrl2/pbes/detail/stategraph_algorithm.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/pbes_expression.h"
#include "mcrl2/pbes/replace.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/stategraph.h"
#include "mcrl2/pbes/tools/pbeschain.h"
#include "mcrl2/pbes/tools/pbesstategraph_options.h"
#include "mcrl2/pbes/unify_parameters.h"
#include "mcrl2/utilities/logger.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/container/flat_map.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <ranges>
#include <unordered_map>
#include <unordered_set>

namespace mcrl2::pbes_system
{

/// A representation of a permutation.
class permutation
{
public:
  permutation() = default;

  permutation(const boost::container::flat_map<std::size_t, std::size_t>& mapping)
    : m_mapping(mapping)
  {}

  boost::container::flat_map<std::size_t, std::size_t> mapping() const { return m_mapping; }

  std::size_t operator[](std::size_t i) const
  {
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
    boost::container::flat_map<std::size_t, std::size_t> new_mapping;

    for (const auto& [key, value]: m_mapping)
    {
      new_mapping[key] = other[value];
    }

    for (const auto& [key, value]: other.m_mapping)
    {
      assert(m_mapping.find(key) == m_mapping.end());
      new_mapping[key] = value;
    }

    return permutation(new_mapping);
  }

  bool operator==(const permutation& other) const { return m_mapping == other.m_mapping; }

private:
  boost::container::flat_map<std::size_t, std::size_t> m_mapping;
};

/// Iterator that generates all permutations of a given set of indices
class permutation_range
{
public:
  permutation_range(const std::vector<std::size_t>& indices)
    : m_indices(indices)
  {
    std::sort(m_indices.begin(), m_indices.end());
    m_current_permutation = m_indices;
    next_permutation(); // Skip the identity permutation
  }

  permutation_range begin() const { return *this; }

  permutation_range end() const
  {
    permutation_range iter(*this);
    iter.m_finished = true;
    return iter;
  }

  permutation operator*() const
  {
    boost::container::flat_map<std::size_t, std::size_t> mapping;
    for (std::size_t i = 0; i < m_indices.size(); ++i)
    {
      mapping[m_indices[i]] = m_current_permutation[i];
    }
    return permutation(mapping);
  }

  permutation_range& operator++()
  {
    if (!next_permutation())
    {
      m_finished = true;
    }
    return *this;
  }

  bool operator!=(const permutation_range& other) const { return m_finished != other.m_finished; }

private:
  bool next_permutation()
  {
    // Find the largest index k such that a[k] < a[k + 1]
    int k = -1;
    for (int i = m_current_permutation.size() - 2; i >= 0; --i)
    {
      if (m_current_permutation[i] < m_current_permutation[i + 1])
      {
        k = i;
        break;
      }
    }

    if (k == -1)
    {
      return false; // No next permutation
    }

    // Find the largest index l greater than k such that a[k] < a[l]
    int l = -1;
    for (int i = m_current_permutation.size() - 1; i > k; --i)
    {
      if (m_current_permutation[k] < m_current_permutation[i])
      {
        l = i;
        break;
      }
    }

    // Swap a[k] and a[l]
    std::swap(m_current_permutation[k], m_current_permutation[l]);

    // Reverse the suffix starting at a[k + 1]
    std::reverse(m_current_permutation.begin() + k + 1, m_current_permutation.end());

    return true;
  }

  std::vector<std::size_t> m_indices;
  std::vector<std::size_t> m_current_permutation;
  bool m_finished = false;
};

/// Returns all the permutations for the given indices.
inline permutation_range permutation_group(const std::vector<std::size_t>& indices)
{
  return permutation_range(indices);
}

/// Prints the permutation as a mapping
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

/// Combines the candidates derived from two different cliques.
inline std::vector<std::pair<permutation, permutation>> combine(
  const std::vector<std::pair<permutation, permutation>>& I_1,
  const std::vector<std::pair<permutation, permutation>>& I_2)
{
  std::vector<std::pair<permutation, permutation>> result;
  for (const auto& [alpha_1, beta_1]: I_1)
  {
    for (const auto& [alpha_2, beta_2]: I_2)
    {
      if (beta_1 == beta_2)
      {
        result.emplace_back(alpha_1.concat(alpha_2), beta_1);
      }
    }
  }

  return result;
}

/// Apply the given permutation to an expression
inline pbes_expression
apply_permutation(const pbes_expression& expr, const std::vector<data::variable>& parameters, const permutation& pi)
{
  data::mutable_map_substitution<> sigma;
  for (std::size_t i = 0; i < parameters.size(); ++i)
  {
    sigma[parameters[i]] = parameters[pi[i]];
  }

  auto result = pbes_system::replace_variables(expr, sigma);

  result = replace_propositional_variables(result,
    [sigma, pi, parameters](const pbes_system::propositional_variable_instantiation& x) -> pbes_system::pbes_expression
    {
      std::vector<data::data_expression> new_parameters(x.parameters().size());
      for (std::size_t i = 0; i < x.parameters().size(); ++i)
      {
        new_parameters[i]
          = data::data_expression(*std::next(x.parameters().begin(), pi[i]));
      }
      return propositional_variable_instantiation(x.name(), data::data_expression_list(new_parameters));
    });

  mCRL2log(log::debug) << "pi(phi): \n" << expr << "\n" << result << std::endl;
  return result;
}

/// Fold is only available in C++23 so we provide a simple implementation here.
template<typename T, typename BinaryOperation>
inline T fold_left(const std::vector<T>& vec, BinaryOperation op)
{
  if (vec.empty())
  {
    throw std::invalid_argument("fold_left: input vector is empty");
  }
  T result = vec[0];
  for (std::size_t i = 1; i < vec.size(); ++i)
  {
    result = op(result, vec[i]);
  }
  return result;
}

inline std::size_t variable_index(const detail::local_control_flow_graph& c)
{
  for (const auto& vertex: c.vertices)
  {
    return vertex.index();
  }

  throw std::runtime_error("No vertices in control flow graph");
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
    std::vector<std::vector<std::pair<permutation, permutation>>> candidates;

    for (const auto& clique: cal_I)
    {
      candidates.emplace_back(clique_candidates(clique));
    }

    if (candidates.empty())
    {
      mCRL2log(log::info) << "No symmetry candidates found!" << std::endl;
      return;
    }

    m_result = fold_left(candidates, combine);
  }

  const std::vector<std::pair<permutation, permutation>>& result() const 
  {
    return m_result;
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
        for (const auto& [to, labels]: s.outgoing_edges())
        {
          for (const auto& equation: m_pbes.equations())
          {
            if (equation.variable().name() == s.name())
            {
              for (const auto& label: labels)
              {
                const auto& variable = equation.predicate_variables().at(label);
                data_parameters.insert(variable.changed().begin(), variable.changed().end());
                data_parameters.insert(variable.used().begin(), variable.used().end());
              }
            }
          }
        }
      }
    }

    // Remove the control flow parameters from the data parameters.
    for (const auto& i: clique)
    {
      // Every vertex should have the same index.
      const detail::local_control_flow_graph& c = m_local_control_flow_graphs[i];
      data_parameters.erase(variable_index(c));
    }

    mCRL2log(log::verbose) << "--- data parameters for clique \n";
    for (const auto& dp: data_parameters)
    {
      mCRL2log(log::verbose) << dp << std::endl;
    }

    return data_parameters;
  }

  /// Computes the set of candidates we can derive from a single clique
  std::vector<std::pair<permutation, permutation>> clique_candidates(const std::vector<size_t>& I)
  {
    auto D = data_parameters(I);

    std::vector<std::size_t> parameter_indices;
    for (const auto& i: I)
    {
      parameter_indices.emplace_back(variable_index(m_local_control_flow_graphs[i]));
    }

    std::vector<std::pair<permutation, permutation>> result;
    for (const auto& alpha: permutation_group(parameter_indices))
    {
      for (const auto& beta: permutation_group(std::vector<std::size_t>(D.begin(), D.end())))
      {
        mCRL2log(log::verbose) << "Trying candidate: " << alpha << " and " << beta << std::endl;

        permutation pi = alpha.concat(beta);
        if (complies(pi, I))
        {
          result.emplace_back(std::move(alpha), std::move(beta));
        }
      }
    }

    mCRL2log(log::verbose) << "--- compliant permutations \n";
    for (const auto& p: result)
    {
      mCRL2log(log::verbose) << p.first << ", " << p.second << std::endl;
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
        mCRL2log(log::verbose) << "--- control flow graphs in clique \n";
        for (const auto& graph: I)
        {
          mCRL2log(log::verbose) << graph << " variable index: " << variable_index(m_local_control_flow_graphs[graph]) << std::endl;
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
    const detail::local_control_flow_graph& graph = m_local_control_flow_graphs.at(c);
    
    std::size_t other_c = 0;
    for (std::size_t i = 0; i < m_local_control_flow_graphs.size(); ++i)
    {
      if (variable_index(m_local_control_flow_graphs.at(i)) == pi[variable_index(m_local_control_flow_graphs[c])])
      {
        other_c = i;
        break;
      }
    }

    const detail::local_control_flow_graph& other_graph = m_local_control_flow_graphs.at(other_c);
    

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
                mCRL2log(log::trace) << "Matching edges from " << s << " to " << *to << " and " << s_prime << " to " << *to_prime << std::endl;

                // t == t'
                // Find the corresponding equation
                bool found_match = false;
                for (const auto& equation: m_pbes.equations())
                {
                  if (equation.variable().name() == s.name())
                  {
                    mCRL2log(log::trace) << "Checking equation " << equation.variable().name() << std::endl;

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
                        mCRL2log(log::trace) << "Removing match " << *matching_j << std::endl;
                        remaining_j.erase(*matching_j);
                      }
                    }

                    if (remaining_j.empty())
                    {
                      found_match = true;
                      break;
                    }
                  }
                }

                if (!found_match)
                {
                  mCRL2log(log::verbose) << "No matching found for edge from " << s << " to " << *to << std::endl;
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

private:
  std::vector<std::pair<permutation, permutation>> m_result; 
};

/// Contains all the implementation of the PBES symmetry algorithm, based on the article by Bartels et al.
class pbes_symmetry
{
public:
  pbes_symmetry(const pbes& input, const data::rewriter&)
  {
    srf = pbes2srf(input);

    mCRL2log(mcrl2::log::debug) << srf.to_pbes() << std::endl;

    unify_parameters(srf, false, false);

    // cliques()
    pbes srf_input = srf.to_pbes();
    cliques_algorithm algorithm(srf_input);
    algorithm.run();
    
    if (!input.equations().empty())
    {
      // After unification, all equations have the same parameters.
      auto parameters = input.equations()[0].variable().parameters();
      m_parameters = std::vector<data::variable>(parameters.begin(), parameters.end());
    }

    for (const auto& option: algorithm.result())
    {
      permutation permutation = option.first.concat(option.second);
      if (symcheck(permutation))
      {
        mCRL2log(log::info) << "Found valid symmetry: " << permutation << std::endl;
      }
    }
  }

private:
  static void cliques(const pbes& input)
  {
  }  

  /// Performs the syntactic check defined as symcheck in the paper.
  bool symcheck(const permutation& pi)
  {
    for (const auto& equation: srf.equations())
    {
      for (const auto& summand : equation.summands())
      {        
        mCRL2log(log::trace) << "Checking summand " << summand << " of equation " << equation << std::endl;
        bool matched = false;
        for (const auto& other_equation: srf.equations())
        {
          for (const auto& other_summand : other_equation.summands())
          {     
            mCRL2log(log::trace) << "Against summand " << other_summand << " of equation " << other_equation << std::endl;
            if (equation.variable().name() == other_equation.variable().name()
                && apply_permutation(summand.condition(), m_parameters, pi) == other_summand.condition()
                && apply_permutation(summand.variable(), m_parameters, pi) == other_summand.variable())
            {
              matched = true;
              break;
            }
          }

          if (matched)
          {
            break;
          }
        }

        if (!matched)
        {
          mCRL2log(log::debug) << "No match for equation " << equation << std::endl;
          return false;
        }
      }        
    }

    return true;
  }

  srf_pbes srf;
  std::vector<data::variable> m_parameters;
};

} // namespace mcrl2::pbes_system

#endif // MCRL_PBES_PBES_SYMMETRY_H