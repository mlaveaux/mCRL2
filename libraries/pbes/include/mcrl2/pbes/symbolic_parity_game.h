// Author(s): Maurice Laveaux and Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_PBES_SYMBOLIC_PARITY_GAME_H
#define MCRL2_PBES_SYMBOLIC_PARITY_GAME_H

#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/pbes/pbes_equation_index.h"
#include "mcrl2/symbolic/alternative_relprod.h"
#include "mcrl2/symbolic/data_index.h"
#include "mcrl2/symbolic/print.h"
#include "mcrl2/symbolic/symbolic_reachability.h"
#include "mcrl2/utilities/text_utility.h"
#include "mcrl2/utilities/stopwatch.h"

namespace mcrl2 {

namespace pbes_system {

using sylvan::ldds::ldd;

namespace detail
{

inline
std::string print_pbes_info(const srf_pbes& pbesspec)
{
  std::ostringstream out;
  pbes_equation_index equation_index(pbesspec);
  for (const auto& equation: pbesspec.equations())
  {
    const auto& name = equation.variable().name();
    out << name << " rank = " << equation_index.rank(name) << " decoration = " << (equation.is_conjunctive() ? "conjunctive" : "disjunctive") << std::endl;
  }
  return out.str();
}

// print the subgraph U of V
template <typename SummandGroup>
std::string print_graph(
  const ldd& U,
  const ldd& V,
  const std::vector<SummandGroup>& R,
  const std::vector<symbolic::data_expression_index>& data_index,
  const ldd& V0, // disjunctive nodes
  const std::map<std::size_t, ldd>& rank_map // maps rank to the corresponding set of nodes
)
{
  using namespace sylvan::ldds;
  using utilities::detail::contains;

  auto rank = [&](const ldd& u)
  {
    for (const auto& [r, U]: rank_map)
    {
      if (includes(U, u))
      {
        return r;
      }
    }
    throw mcrl2::runtime_error("print_graph: could not find a rank");
  };

  auto index = [](const std::vector<ldd>& v, const ldd& x)
  {
    auto i = std::find(v.begin(), v.end(), x);
    if (i == v.end())
    {
      throw mcrl2::runtime_error("print_graph: index error");
    }
    return i - v.begin();
  };

  auto values = [](const ldd& X)
  {
    std::vector<ldd> result;
    auto X_elements = ldd_solutions(X);
    for (const auto& x: X_elements)
    {
      result.push_back(cube(x));
    }
    return std::make_pair(result, X_elements);
  };

  auto succ = [&](const ldd& U)
  {
    ldd result = empty_set();
    for (std::size_t i = 0; i < R.size(); i++)
    {
      result = union_(result, alternative_relprod(U, R[i]));
    }
    return result;
  };

  auto [U_values, U_solutions] = values(U);
  auto [V_values, V_solutions] = values(V);

  std::vector<std::string> text(U_values.size());

  for (std::size_t i = 0; i < U_values.size(); i++)
  {
    ldd u = U_values[i];
    std::size_t u_index = index(V_values, u);

    ldd W = succ(u);
    auto [W_values, W_solutions] = values(W);
    std::vector<std::uint32_t> u_successors;
    for (const ldd& w: W_values)
    {
      if (contains(U_values, w))
      {
        u_successors.push_back(index(V_values, w));
      }
    }
    text[i] = std::to_string(u_index) + " " + symbolic::print_state(data_index, U_solutions[i]) + ", decoration = " + (includes(V0, u) ? "disjunctive" : "conjunctive") + ", rank = " + std::to_string(rank(u)) + ", successors = " + core::detail::print_list(u_successors);
  }
  return utilities::string_join(text, "\n");
}

// print the indices of U (subset of V)
std::string print_nodes(const ldd& U, const ldd& V)
{
  using namespace sylvan::ldds;

  auto index = [](const std::vector<ldd>& v, const ldd& x)
  {
    auto i = std::find(v.begin(), v.end(), x);
    if (i == v.end())
    {
      throw mcrl2::runtime_error("print_graph: index error");
    }
    return i - v.begin();
  };

  auto values = [](const ldd& X)
  {
    std::vector<ldd> result;
    auto X_elements = ldd_solutions(X);
    for (const auto& x: X_elements)
    {
      result.push_back(cube(x));
    }
    return std::make_pair(result, X_elements);
  };

  auto [V_values, V_solutions] = values(V);
  auto [U_values, W0_solutions] = values(U);

  std::vector<std::size_t> u;
  for (const ldd& x: U_values)
  {
    u.push_back(index(V_values, x));
  }
  return core::detail::print_set(u);
}

/// \brief maps proposition variable ldd values to (rank, is_disjunctive)
std::map<std::size_t, std::pair<std::size_t, bool>> compute_equation_info(const pbes_system::srf_pbes& pbes, const std::vector<symbolic::data_expression_index>& data_index)
{
  pbes_system::pbes_equation_index equation_index(pbes);

  // map propositional variable names to the corresponding ldd value
  std::map<core::identifier_string, std::uint32_t> propvar_index;
  for (const data::data_expression& X: data_index[0])
  {
    const auto& X_ = atermpp::down_cast<data::function_symbol>(X);
    std::uint32_t i = propvar_index.size();
    propvar_index[X_.name()] = i;
  }

  // maps ldd values to (rank, is_disjunctive)
  std::map<std::size_t, std::pair<std::size_t, bool>> equation_info;
  for (const auto& equation: pbes.equations())
  {
    const core::identifier_string& name = equation.variable().name();
    std::size_t rank = equation_index.rank(name);
    bool is_disjunctive = !equation.is_conjunctive();
    auto i = propvar_index.find(name);
    if (i != propvar_index.end())
    {
      std::uint32_t ldd_value = i->second;
      equation_info[ldd_value] = { rank, is_disjunctive };
    }
  }

  return equation_info;
}

} // namespace detail

/// \brief This class represents a symbolic (incomplete) parity game with sinks.
///        Many functions have a parameter V that restricts operations to that set of vertices.
class symbolic_parity_game
{
  protected:
    ldd m_V[2]; // m_V[0] is the set of even nodes, m_V[1] is the set of odd nodes
    const std::vector<symbolic::summand_group> m_summand_groups;
    std::map<std::size_t, ldd> m_rank_map;
    bool m_no_relprod = false;
    bool m_chaining = false;

    const std::vector<symbolic::data_expression_index>& m_data_index; // for debugging only
    ldd m_all_nodes; // for debugging only

  public:

    /// \brief Determine a symbolic parity game from the given pbes, transition groups and index.
    /// \param V the set of reachable vertices.
    symbolic_parity_game(
      const srf_pbes& pbes,
      const std::vector<symbolic::summand_group> summand_groups,
      const std::vector<symbolic::data_expression_index>& data_index,
      const ldd& V,
      bool no_relprod,
      bool chaining
    )
      : m_summand_groups(summand_groups), m_no_relprod(no_relprod), m_chaining(chaining), m_data_index(data_index), m_all_nodes(V)
    {
      using namespace sylvan::ldds;
      using utilities::detail::contains;

      // Determine priority and owner from the given pbes.
      auto equation_info = detail::compute_equation_info(pbes, data_index);

      m_V[0] = empty_set();
      m_V[1] = empty_set();

      // determine the rank and owner of all states.
      for (const auto& [value, p]: equation_info)
      {
        auto rank = p.first;
        auto is_disjunctive = p.second;
        ldd X = fix_first_element(V, value);

        auto j = m_rank_map.find(rank);
        if (j == m_rank_map.end())
        {
          m_rank_map[rank] = X;
        }
        else
        {
          j->second = union_(j->second, X);
        }

        if (is_disjunctive)
        {
          m_V[0] = union_(m_V[0], X);
        }
        else
        {
          m_V[1] = union_(m_V[1], X);
        }
      }
    }
    
    /// \brief Determine a symbolic parity game from the given pbes, transition groups and index.
    /// \param V the set of reachable vertices.
    symbolic_parity_game(
      const std::vector<symbolic::summand_group>& summand_groups,
      const std::vector<symbolic::data_expression_index>& data_index,
      const ldd& V,
      const ldd& Veven,
      const std::vector<ldd>& prio,
      bool no_relprod,
      bool chaining
    )
      : m_summand_groups(summand_groups), m_no_relprod(no_relprod), m_chaining(chaining), m_data_index(data_index), m_all_nodes(V)
    {
      m_V[0] = Veven;
      m_V[1] = minus(V, Veven);

      std::size_t i = 0;
      for (const auto& p : prio)
      {
        m_rank_map[i] = p;
        ++i;
      }
    }

    /// \returns Prints basic parity game information such as number of vertices per priority and per owners.
    void print_information()
    {
      mCRL2log(log::verbose) << "--- parity game information ---" << std::endl;
      for (const auto&[rank, Vrank] : m_rank_map)
      {
        mCRL2log(log::verbose) << "priority " << rank << ": there are " << satcount(Vrank) << " vertices.\n";
      }

      mCRL2log(log::verbose) << "there are " << satcount(m_V[0]) << " even vertices and " << satcount(m_V[1]) << " odd vertices.\n";
    }

    /// \returns A string representing the given vertex set in human readable form.
    std::string print_nodes(const ldd& V) const
    {
      return detail::print_nodes(V, m_all_nodes);
    }

    /// \returns A string representing the graph restricted to V.
    std::string print_graph(const ldd& V) const
    {
      return detail::print_graph(V, m_all_nodes, m_summand_groups, m_data_index, m_V[0], m_rank_map);
    }

    /// \brief Compute the attractor set for U assuming that sinks(V) = emptyset.
    /// \param alpha the current player
    /// \param V is the set of states
    /// \param Vplayer a partitioning of the nodes into the sets of even nodes V[0] and odd V[1].
    /// \param I is a set of incomplete vertices.
    /// \param T A set of states such that iteration stops when one of them occurs in the attractor.
    ldd safe_attractor(const ldd& U,
      std::size_t alpha,
      const ldd& V,
      const std::array<const ldd, 2>& Vplayer,
      const ldd& I = sylvan::ldds::empty_set(),
      const ldd& T = sylvan::ldds::empty_set()) const
    {
      stopwatch attractor_watch;
      mCRL2log(log::debug) << "start attractor set computation\n";

      using namespace sylvan::ldds;

      std::size_t iter = 0;
      ldd Z = U;
      ldd todo = U;

      ldd Zoutside = minus(V, Z);
      while (todo != empty_set())
      {
        // Terminate early when a vertex in T was found.
        if (intersect(T, Z) != empty_set() )
        {
          return Z;
        }

        mCRL2log(log::debug1) << "todo = " << print_nodes(todo) << std::endl;
        mCRL2log(log::debug1) << "Zoutside = " << print_nodes(Zoutside) << std::endl;
        stopwatch iter_start;

        //todo = minus(safe_control_predecessors_impl(alpha, todo, V, Zoutside, V, Vplayer, I), Z);
        todo = minus(safe_control_predecessors_impl2(alpha, todo, Zoutside, Vplayer, I), Z);
        Z = union_(Z, todo);
        Zoutside = minus(Zoutside, todo);

        mCRL2log(log::debug) << "attractor set iteration " << iter << " (time = " << std::setprecision(2) << std::fixed << iter_start.seconds() << "s)" << std::endl;

        ++iter;
      }

      mCRL2log(log::debug) << "finished attractor set computation (time = " << std::setprecision(2) << std::fixed << attractor_watch.seconds() << "s)" << std::endl;
      return Z;
    }

    /// \brief Compute the monotone attractor set for U assuming that sinks(V) = emptyset.
    /// \param alpha the current player.
    /// \param c priority.
    /// \param V the set of all vertices.
    /// \param Vplayer a partitioning of the nodes into the sets of even nodes V[0] and odd V[1].
    ldd safe_monotone_attractor(const ldd& U,
       std::size_t alpha,
       std::size_t c,
       const ldd& V,
       const std::array<const ldd, 2>& Vplayer,
       const ldd& I = sylvan::ldds::empty_set(),
       const ldd& T = sylvan::ldds::empty_set()) const
    {
      using namespace sylvan::ldds;

      stopwatch attractor_watch;
      mCRL2log(log::debug) << "start monotone attractor set computation\n";

      using namespace sylvan::ldds;

      // Compute the set of states with at least priority c.
      ldd Vc = empty_set();
      for (const auto&[rank, Vrank] : m_rank_map)
      {
        if (rank >= c)
        {
          Vc = union_(Vc, Vrank);
        }
      }

      // Vertices of player alpha and priority c
      std::size_t iter = 0;
      ldd Z = empty_set();
      ldd todo = U; // union of U and X

      ldd Zoutside = V; // V minus Z
      while (todo != empty_set())
      {
        // Terminate early when a vertex in T was found.
        if (intersect(T, Z) != empty_set() )
        {
          return Z;
        }

        mCRL2log(log::debug1) << "todo = " << print_nodes(todo) << std::endl;
        mCRL2log(log::debug1) << "Zoutside = " << print_nodes(Zoutside) << std::endl;
        stopwatch iter_start;

        todo = intersect(Vc, minus(safe_control_predecessors_impl(alpha, union_(todo, U), V, minus(Zoutside, U), Vc, Vplayer, I), Z));
        Z = union_(Z, todo);
        Zoutside = minus(Zoutside, todo);

        mCRL2log(log::debug) << "monotone attractor set iteration " << iter << " (time = " << std::setprecision(2) << std::fixed << iter_start.seconds() << "s)" << std::endl;

        ++iter;
      }

      mCRL2log(log::debug) << "finished monotone attractor set computation (time = " << std::setprecision(2) << std::fixed << attractor_watch.seconds() << "s)" << std::endl;
      return Z;
    }

    /// \returns (min, Vmin) where min is the minimum rank in V and Vmin is the set of vertices with the minimum rank in V
    std::pair<std::size_t, ldd> get_min_rank(const ldd& V) const
    {
      using namespace sylvan::ldds;

      for (auto i = m_rank_map.begin(); i != m_rank_map.end(); ++i)
      {
        ldd Vmin = intersect(V, i->second);
        if (Vmin != empty_set())
        {
          std::size_t min_rank = i->first;
          return { min_rank, Vmin };
        }
      }

      throw mcrl2::runtime_error("get_min_rank did not find any nodes");
    }

    /// \brief Computes the pair of even and odd vertices.
    std::array<const ldd, 2> players(const ldd& V) const
    {
      return { intersect(V, m_V[0]), intersect(V, m_V[1]) };
    }

    /// \brief Computes the vertices with even parity priority that are not sinks and the same for odd.
    std::array<const ldd, 2> parity(const ldd& V) const
    {
      std::array<ldd, 2> parity;
      for (const auto&[rank, Vrank] : ranks())
      {
        parity[rank % 2] = sylvan::ldds::union_(parity[rank % 2], Vrank);
      }

      ldd Vother = minus(V, sinks(V, V));
      return { intersect(Vother, parity[0]), intersect(Vother, parity[1]) };
    }

    /// \brief Computes all vertices above priority c.
    ldd prio_above(const ldd& V, std::size_t c) const
    {      
      // Compute the set of states with at least priority c.
      ldd Vc = sylvan::ldds::empty_set();
      for (const auto&[rank, Vrank] : m_rank_map)
      {
        if (rank >= c)
        {
          Vc = union_(Vc, Vrank);
        }
      }

      return intersect(V, Vc);
    }

    /// \brief Removes all winning states (and updates winning partition).
    ldd compute_total_graph(const ldd& V, const ldd& I, const ldd& Vsinks, std::array<ldd, 2>& winning) const
    {
      using namespace sylvan::ldds;
      std::array<const ldd, 2> Vplayer = players(V);

      // After removing the deadlock (winning) states the resulting set of states is a total graph.
      mCRL2log(log::debug) << "removing winning regions" << std::endl;
      if (Vsinks != empty_set())
      {
        winning[0] = union_(winning[0], intersect(Vsinks, m_V[1]));
        winning[1] = union_(winning[1], intersect(Vsinks, m_V[0]));
      }

      winning[0] = safe_attractor(winning[0], 0, V, Vplayer, I);
      winning[1] = safe_attractor(winning[1], 1, V, Vplayer, I);

      // After removing the deadlock (winning) states the resulting set of states is a total graph.
      return minus(minus(V, winning[0]), winning[1]);
    }

    /// \brief Computes the set of vertices for which partial solving is safe w.r.t. alpha.
    ldd compute_safe_vertices(
        std::size_t alpha,
        const ldd& V,
        const ldd& I) const
    {
      using namespace sylvan::ldds;

      // Compute the safe sets from the resulting subgraph.
      std::array<const ldd, 2> Vplayer = players(V);
      ldd S = sinks(I, V);
      return minus(V, safe_attractor(union_(intersect(I, Vplayer[1-alpha]), S), 1-alpha, V, Vplayer));
    }

    /// \brief Returns the mapping from priorities (ranks) to vertex sets.
    const std::map<std::size_t, ldd>& ranks() const { return m_rank_map; }

    /// \returns The set { u in U | exists v in V: u ->* v } where with chaining ->* only visits states in V (if U subseteq V)
    ldd predecessors(const ldd& U, const ldd& V) const
    {
      using namespace sylvan::ldds;

      ldd result;
      for (int i = m_summand_groups.size() - 1; i >= 0; --i)
      {
        const symbolic::summand_group& group = m_summand_groups[i];

        stopwatch watch;
        result = union_(result, predecessors(U, V, group));
        mCRL2log(log::debug1) << "added predecessors for group " << i << " out of " << m_summand_groups.size()
                               << " (time = " << std::setprecision(2) << std::fixed << watch.seconds() << "s)\n";
      }

      return result;
    }

    /// \returns A set of vertices { u in U | exists v in V: u ->* v } where ->* only visits intermediate vertices in W (without chaining ->* = ->)
    /// \pre U,W subseteq V.
    ldd predecessors(const ldd& U, const ldd& V, const ldd& W) const
    {
      using namespace sylvan::ldds;

      if (m_chaining)
      {
        ldd P = empty_set();
        ldd todo = V;
        for (int i = m_summand_groups.size() - 1; i >= 0; --i)
        {
          const symbolic::summand_group& group = m_summand_groups[i];

          stopwatch watch;
          ldd todo1 = predecessors(U, todo, group);
          mCRL2log(log::debug1) << "added predecessors for group " << i << " out of " << m_summand_groups.size()
                                 << " (time = " << std::setprecision(2) << std::fixed << watch.seconds() << "s)\n";

          P = union_(P, todo1);
          todo = union_(todo, intersect(todo1, W));
        }

        return P;
      }
      else
      {
        return predecessors(U, V);
      }
    }

    /// \brief Compute the safe control attractor set for U w.r.t. vertices in V.
    ///        The set W is a set of vertices that restrict chaining.
    ldd safe_control_predecessors(std::size_t alpha,
      const ldd& U,
      const ldd& V,
      const ldd& W,
      const std::array<const ldd, 2>& Vplayer,
      const ldd& I = sylvan::ldds::empty_set()) const
    {
      ldd outside = minus(V, U);
      return safe_control_predecessors_impl(alpha, U, V, outside, W, Vplayer, I);
    }

    /// \brief Computes the set of vertices in U subseteq V that are sinks (no outgoing edges into V).
    ldd sinks(const ldd& U, const ldd& V) const
    {
      return minus(U, predecessors(U, V));
    }   

private:
    /// \returns The set { u in U | exists v in V: u -> v }, where -> is described by the given group.
    ldd predecessors(const ldd& U, const ldd& V, const symbolic::summand_group& group) const
    {
      return m_no_relprod ? symbolic::alternative_relprev(V, group, U) : relprev(V, group.L, group.Ir, U);
    }

    /// \brief Compute the safe control attractor set for todo where chaining is restricted to W.
    ///        The set outside should be minus(V, U)
    ldd safe_control_predecessors_impl(std::size_t alpha,
      const ldd& U,
      const ldd& V,
      const ldd& outside,
      const ldd& W,
      const std::array<const ldd, 2>& Vplayer,
      const ldd& I = sylvan::ldds::empty_set()) const
    {
      using namespace sylvan::ldds;

      ldd P = predecessors(V, U, intersect(Vplayer[alpha], W));
      ldd Palpha = intersect(P, Vplayer[alpha]);
      ldd Pforced = minus(intersect(P, Vplayer[1-alpha]), I);

      for (std::size_t i = 0; i < m_summand_groups.size(); ++i)
      {
        const symbolic::summand_group& group = m_summand_groups[i];

        stopwatch watch;
        Pforced = minus(Pforced, predecessors(Pforced, outside, group));

        mCRL2log(log::debug1) << "removed 1 - alpha predecessors for group " << i << " out of " << m_summand_groups.size()
                               << " (time = " << std::setprecision(2) << std::fixed << watch.seconds() << "s)\n";
      }

      return union_(Palpha, Pforced);
    } 
    
    /// \brief Compute the safe control attractor set for U.
    ldd safe_control_predecessors_impl2(std::size_t alpha,
      const ldd& U,
      const ldd& Zoutside,
      const std::array<const ldd, 2>& Vplayer,
      const ldd& I = sylvan::ldds::empty_set()) const
    {
      using namespace sylvan::ldds;

      ldd P = predecessors(Zoutside, U);
      ldd Palpha = intersect(P, Vplayer[alpha]);
      ldd Pforced = minus(intersect(P, Vplayer[1-alpha]), I);

      for (std::size_t i = 0; i < m_summand_groups.size(); ++i)
      {
        const symbolic::summand_group& group = m_summand_groups[i];

        stopwatch watch;
        Pforced = minus(Pforced, predecessors(Pforced, Zoutside, group));

        mCRL2log(log::debug) << "removed 1 - alpha predecessors for group " << i << " out of " << m_summand_groups.size()
                               << " (time = " << std::setprecision(2) << std::fixed << watch.seconds() << "s)\n";
      }

      return union_(Palpha, Pforced);
    }

};

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_SYMBOLIC_PBESSOLVE_H
