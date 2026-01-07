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

#include "mcrl2/data/data_expression.h"
#include "mcrl2/pbes/detail/instantiate_global_variables.h"
#include "mcrl2/pbes/detail/cartesian_product.h"
#include "mcrl2/pbes/detail/fold_left.h"
#include "mcrl2/pbes/detail/permutation.h"
#include "mcrl2/pbes/detail/stategraph_algorithm.h"
#include "mcrl2/pbes/detail/stategraph_local_algorithm.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/unify_parameters.h"
#include "mcrl2/utilities/logger.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/container/flat_map.hpp>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <type_traits>

namespace mcrl2::pbes_system
{

/// Combines the candidates derived from two different cliques.
template<typename Range>
  requires(
    std::ranges::range<Range>
    && std::is_same_v<typename std::ranges::range_value_t<Range>, std::pair<detail::permutation, detail::permutation>>)
inline std::ranges::range auto candidate_combine(Range I_1, Range I_2)
{
  // It seems that structured bindings cannot be used here.
  auto view = detail::cartesian_product(I_1, I_2);
  return view
         | std::views::filter([](const std::pair<std::pair<detail::permutation, detail::permutation>,
                                std::pair<detail::permutation, detail::permutation>>& pair)
           { return pair.first.second == pair.second.second; })
         | std::ranges::views::transform(
           [](const std::pair<std::pair<detail::permutation, detail::permutation>,
             std::pair<detail::permutation, detail::permutation>>& pair)
             -> std::pair<detail::permutation, detail::permutation>
           {
             const auto& [alpha_1, beta_1] = pair.first;
             const auto& [alpha_2, beta_2] = pair.second;
             return std::make_pair(alpha_1.concat(alpha_2), beta_1);
           });
}

/// Returns the index of the variable of this control flow graph.
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
  cliques_algorithm(const pbes& input, std::vector<data::variable> parameters)
    : super(input, pbesstategraph_options{.print_influence_graph = true}),
      parameters(std::move(parameters))
  {}

  void run() override
  {
    // We explicitly ignore the virtual call to run in the base class
    detail::stategraph_algorithm::stategraph_algorithm::run(); // NOLINT(bugprone-parent-virtual-call)

    compute_local_control_flow_graphs();

    for (decltype(m_local_control_flow_graphs)::iterator i = m_local_control_flow_graphs.begin();
      i != m_local_control_flow_graphs.end();
      ++i)
    {
      mCRL2log(log::verbose) << "--- computed local control flow graph " << (i - m_local_control_flow_graphs.begin())
                             << " --- \n"
                             << *i << std::endl;
    }
  }

  /// Computes the set of candidates we can derive from a single clique
  /// TODO: How to specify requires on the return type?
  std::ranges::range auto clique_candidates(const std::vector<size_t>& I, const std::vector<size_t>& all_control_parameters) const
  {
    std::set<std::size_t> D = data_parameters(I, all_control_parameters);

    std::vector<std::size_t> parameter_indices;
    for (const std::size_t& i: I)
    {
      parameter_indices.emplace_back(variable_index(m_local_control_flow_graphs[i]));
    }

    // Determine parameters with the same sort
    std::vector<std::vector<std::size_t>> same_sort_parameters;
    std::size_t index = 0;
    for (const data::variable& param : parameters) 
    {
      bool found = false;
      for (std::vector<std::size_t>& group : same_sort_parameters) 
      {
        if (param.sort() == parameters[group[0]].sort()) 
        {
          group.emplace_back(index);
          found = true;
          break;
        }
      }

      if (!found) 
      {
        same_sort_parameters.push_back({index});
      } 

      index++;
    }

    for (const std::vector<std::size_t>& group : same_sort_parameters) 
    {
      mCRL2log(log::verbose) << "--- parameters with same sort --- \n";
      for (const std::size_t& param_index : group) 
      {
        mCRL2log(log::verbose) << param_index << " : " << parameters[param_index].sort() << std::endl;
      }
    }



    return detail::cartesian_product(detail::permutation_group(parameter_indices),
             detail::permutation_group(std::vector<std::size_t>(D.begin(), D.end())))
           | std::ranges::views::transform(
             [this, I, all_control_parameters](const std::pair<detail::permutation, detail::permutation>& pair)
               -> std::optional<std::pair<detail::permutation, detail::permutation>>
             {
               const auto& [alpha, beta] = pair;

               detail::permutation pi = alpha.concat(beta);
               mCRL2log(log::debug) << "Trying candidate: " << alpha << " and " << beta << std::endl;
               if (complies(pi, I, all_control_parameters))
               {
                 mCRL2log(log::verbose) << "Compliant permutation: " << pi << std::endl;
                 return std::make_pair(alpha, beta);
               }

               return std::nullopt;
             })
           | std::ranges::views::filter(
             [](const std::optional<std::pair<detail::permutation, detail::permutation>>& b) -> bool
             { return b.has_value(); })
           | std::ranges::views::transform([](std::optional<std::pair<detail::permutation, detail::permutation>> b)
                                             -> std::pair<detail::permutation, detail::permutation> { return *b; });
  }

  /// Takes as input a set of parameters and removes the control flow parameters.
  std::set<std::size_t> remove_cfgs(const std::set<size_t>& parameters, const std::vector<size_t>& all_control_parameters) const
  {
    std::set<size_t> data_parameters = parameters;

    // Remove the control flow parameters from the data parameters.
    for (const size_t& i: all_control_parameters)
    {
      // Every vertex should have the same index.
      const detail::local_control_flow_graph& c = m_local_control_flow_graphs[i];
      data_parameters.erase(variable_index(c));
    }

    return data_parameters;
  }

  /// Takes as input a clique of compatible control flow parameters and return
  /// the set of all data parameters that somehow play a role for any of these
  /// parameters.
  std::set<std::size_t> data_parameters(const std::vector<size_t>& clique, const std::vector<size_t>& all_control_parameters) const
  {
    std::set<std::size_t> data_parameters;
    for (const std::size_t& i: clique)
    {
      const detail::local_control_flow_graph& c = m_local_control_flow_graphs[i];
      for (const detail::local_control_flow_graph_vertex& s: c.vertices)
      {
        // Compute the data parameters
        // Get the changed by, used for and used in
        for (const auto& [to, labels]: s.outgoing_edges())
        {
          for (const detail::stategraph_equation& equation: m_pbes.equations())
          {
            if (equation.variable().name() == s.name())
            {
              for (const unsigned long& label: labels)
              {
                const detail::predicate_variable& variable = equation.predicate_variables().at(label);
                data_parameters.insert(variable.changed().begin(), variable.changed().end());
                data_parameters.insert(variable.used().begin(), variable.used().end());
              }
            }
          }
        }
      }
    }

    data_parameters =  remove_cfgs(data_parameters, all_control_parameters);
    mCRL2log(log::verbose) << "--- data parameters for clique --- \n";
    for (const size_t& parameter: data_parameters)
    {
      mCRL2log(log::verbose) << parameter << std::endl;
    }

    return data_parameters;
  }

  /// Determine the cliques of the control flow graphs.
  std::vector<std::vector<std::size_t>> cliques()
  {
    std::vector<std::vector<std::size_t>> cal_I;
    for (int i = 0; i < m_local_control_flow_graphs.size(); i++)
    {
      if (std::any_of(cal_I.begin(),
            cal_I.end(),
            [i](const std::vector<unsigned long>& clique)
            { return std::find(clique.begin(), clique.end(), i) != clique.end(); }))
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
        for (const unsigned long& graph: I)
        {
          mCRL2log(log::verbose) << graph << " variable index: " << variable_index(m_local_control_flow_graphs[graph])
                                 << std::endl;
        }
        cal_I.emplace_back(I);
      }
    }

    return cal_I;
  }

  /// Returns true iff all vertices in I comply with the detail::permutation pi.
  bool complies(const detail::permutation& pi, const std::vector<std::size_t>& I, const std::vector<size_t>& all_control_parameters) const
  {
    return std::all_of(I.begin(), I.end(), [&](std::size_t c) { return complies(pi, c, all_control_parameters); });
  }

  /// Takes a detail::permutation and a control flow parameter and returns true or
  /// false depending on whether the detail::permutation complies with the control
  /// flow parameter according to Definition
  bool complies(const detail::permutation& pi, std::size_t c, const std::vector<size_t>& all_control_parameters) const
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
    for (const detail::local_control_flow_graph_vertex& s: graph.vertices)
    {
      for (const detail::local_control_flow_graph_vertex& s_prime: other_graph.vertices)
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
                mCRL2log(log::trace) << "Matching edges from " << s << " to " << *to << " and " << s_prime << " to "
                                     << *to_prime << std::endl;

                // t == t'
                // Find the corresponding equation
                bool found_match = false;
                for (const detail::stategraph_equation& equation: m_pbes.equations())
                {
                  if (equation.variable().name() == s.name())
                  {
                    mCRL2log(log::trace) << "Checking equation " << equation.variable().name() << std::endl;

                    // For each i find a corresponding j.
                    std::set<std::size_t> remaining_j = labels_prime;
                    for (const std::size_t& i: labels)
                    {
                      const detail::predicate_variable& variable = equation.predicate_variables().at(i);
                      mCRL2log(log::trace) << "Checking summand " << variable << std::endl;

                      std::optional<std::size_t> matching_j;
                      for (const std::size_t& j: remaining_j)
                      {
                        const detail::predicate_variable& variable_prime = equation.predicate_variables().at(j);
                        mCRL2log(log::trace) << "Against summand " << variable_prime << std::endl;

                        auto changed = remove_cfgs(variable.changed(), all_control_parameters);
                        auto used = remove_cfgs(variable.used(), all_control_parameters);
                        auto changed_prime = remove_cfgs(variable_prime.changed(), all_control_parameters);
                        auto used_prime = remove_cfgs(variable_prime.used(), all_control_parameters);

                        mCRL2log(log::trace) << "Changed: " << core::detail::print_set(changed) << ", used: "
                                             << core::detail::print_set(used) << std::endl;
                        mCRL2log(log::trace) << "Permuted changed: " << core::detail::print_set(pi.permute(changed_prime)) << ", permuted used: "
                                             << core::detail::print_set(pi.permute(used_prime)) << std::endl;

                        if (pi.permute(used_prime) == used
                            && pi.permute(changed_prime) == changed)
                        {
                          matching_j = j;
                          break;
                        }
                      }

                      if (matching_j)
                      {
                        // Found a matching j for i.
                        mCRL2log(log::trace) << "Matching " << i << " to " << *matching_j << std::endl;
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
                  mCRL2log(log::debug) << "No matching found for edge from " << s << " to " << *to << std::endl;
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
      for (const detail::stategraph_equation& equation: m_pbes.equations())
      {
        if (equation.variable().name() == s.name())
        {
          for (const std::size_t& label: it->second)
          {
            // Compute the sizes.
            const detail::predicate_variable& variable = equation.predicate_variables().at(label);
            // TODO: Remove cfgs?
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
    for (const detail::local_control_flow_graph_vertex& s: c.vertices)
    {
      // There exist t such that s and t match according to the definitions in the paper.
      for (const detail::local_control_flow_graph_vertex& s_c_prime: c_prime.vertices)
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

    for (const detail::local_control_flow_graph_vertex& vertex: c.vertices)
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

    for (const detail::local_control_flow_graph_vertex& vertex_prime: c_prime.vertices)
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
  std::vector<data::variable> parameters;
};

/// Contains all the implementation of the PBES symmetry algorithm, based on the article by Bartels et al.
class pbes_symmetry
{
public:

  pbes_symmetry(const pbes& input)
  {
    // This has to be done consistently with the LPS for the counter examples.
    pbes pbes = input;
    data::mutable_map_substitution<> sigma = pbes_system::detail::instantiate_global_variables(pbes);
    pbes_system::detail::replace_global_variables(pbes, sigma);
    srf = pbes2srf(pbes);

    unify_parameters(srf, false, true);

    mCRL2log(mcrl2::log::debug) << srf.to_pbes() << std::endl;
    
    if (!srf.equations().empty())
    {
      // After unification, all equations have the same parameters.
      data::variable_list list = srf.equations()[0].variable().parameters();
      parameters = std::vector<data::variable>(list.begin(), list.end());
    }
  }

  void run()
  {
    // cliques()
    pbes srf_input = srf.to_pbes();
    cliques_algorithm algorithm(srf_input, parameters);
    algorithm.run();

    std::vector<size_t> all_control_parameters;
    auto cliques = algorithm.cliques();
    for (const auto& clique: cliques)
    {
      for (const auto& c: clique)
      {
        all_control_parameters.emplace_back(c);
      }
    }

    std::vector<std::vector<std::pair<detail::permutation, detail::permutation>>> candidates;
    for (const auto& clique: cliques)
    {
      // std::ranges::to<std::vector>) is not a thing yet.
      std::vector<std::pair<detail::permutation, detail::permutation>> clique_candidates;
      for (const std::pair<detail::permutation, detail::permutation>& candidate: algorithm.clique_candidates(clique, all_control_parameters))
      {
        clique_candidates.emplace_back(candidate);
      }

      candidates.emplace_back(std::move(clique_candidates));
    }

    for (const auto& result:
      detail::fold_left<std::vector<std::pair<detail::permutation, detail::permutation>>>(candidates,
        [](const auto& acc, const auto& x)
        {
          // std::ranges::to<std::vector>) is not a thing yet.
          std::vector<std::pair<detail::permutation, detail::permutation>> clique_candidates;
          for (const std::pair<detail::permutation, detail::permutation>& candidate: candidate_combine(acc, x))
          {
            clique_candidates.emplace_back(candidate);
          }

          return clique_candidates;
        }))
    {
      detail::permutation permutation = result.first.concat(result.second);
      std::cout << "Checking permutation: " << permutation << std::endl;
      if (symcheck(permutation))
      {
        std::cout << "Found symmetry: " << permutation << std::endl;

        // Stop after finding the first symmetry.
        return;
      }
    }
  }

  /// Checks whether a given permutation is a symmetry for the PBES.
  bool check_permutation(const detail::permutation& pi)
  {
    return symcheck(pi);
  }

  /// Returns the PBES after preprocessing it for symmetry reduction.
  pbes pbesspec() const
  {
    return srf.to_pbes();
  }

private:
  /// Performs the syntactic check defined as symcheck in the paper.
  bool symcheck(const detail::permutation& pi)
  {
    for (const auto& equation: srf.equations())
    {
      mCRL2log(log::trace) << "Checking equation " << equation << std::endl;
      for (const auto& summand: equation.summands())
      {
        mCRL2log(log::trace) << "Summand " << summand << std::endl;

        bool matched = false;
        for (const auto& other_equation: srf.equations())
        {
          mCRL2log(log::trace) << "Against equation " << other_equation << std::endl;
          for (const auto& other_summand: other_equation.summands())
          {
            mCRL2log(log::trace) << "and summand " << other_summand << std::endl;
            if (equation.variable().name() == other_equation.variable().name()
                && detail::apply_permutation(summand.condition(), parameters, pi) == other_summand.condition()
                && detail::apply_permutation(summand.variable(), parameters, pi) == other_summand.variable())
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

  std::vector<data::variable> parameters;
  srf_pbes srf;
};

} // namespace mcrl2::pbes_system

#endif // MCRL_PBES_PBES_SYMMETRY_H