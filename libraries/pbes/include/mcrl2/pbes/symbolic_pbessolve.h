// Author(s): Maurice Laveaux and Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_PBES_SYMBOLIC_PBESSOLVE_H
#define MCRL2_PBES_SYMBOLIC_PBESSOLVE_H

#include "sylvan_ldd.hpp"
#ifdef MCRL2_ENABLE_SYLVAN

#include "symbolic_parity_game.h"
#include "mcrl2/utilities/logger.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system {

using sylvan::ldds::ldd;

class symbolic_pbessolve_algorithm
{
  private:
    const symbolic_parity_game& m_G;
    bool m_check_strategy = false;

  public:
    symbolic_pbessolve_algorithm(const symbolic_parity_game& G, bool check_strategy = false) :
      m_G(G),
      m_check_strategy(check_strategy)
    {}

    std::tuple<ldd, ldd, ldd, ldd> zielonka(const ldd& V)
    {
      using namespace sylvan::ldds;

      if (V == empty_set())
      {
        return { empty_set(), empty_set(), empty_set(), empty_set() };
      }

      stopwatch timer;
      mCRL2log(log::debug) << "start zielonka recursion\n";
      auto [m, U] = m_G.get_min_rank(V);

      std::size_t alpha = m % 2; // 0 = disjunctive, 1 = conjunctive

      // Compute the partitioning of V for players 0 (in V[0]) and 1 (in V[1]).
      std::array<const ldd, 2> Vplayer = m_G.players(V);

      ldd W[2];
      ldd strategy[2];
      ldd W_1[2];
      ldd strategy_1[2];

      const auto [A, A_strategy] = m_G.safe_attractor(U, alpha, V, Vplayer);
      mCRL2log(log::trace) << "A = attractor(" << m_G.print_nodes(U) << ", " << m_G.print_nodes(V) << ") = " << m_G.print_nodes(A) << std::endl;
      std::tie(W_1[0], W_1[1], strategy_1[0], strategy_1[1]) = zielonka(minus(V, A));

      // Original Zielonka version
      if (W_1[1 - alpha] == empty_set())
      {
        W[alpha] = union_(A, W_1[alpha]);
        W[1 - alpha] = empty_set();
        strategy[alpha] = union_(union_(A_strategy, strategy_1[alpha]), merge(intersect(U, Vplayer[alpha]), V));
        strategy[1 - alpha] = empty_set();
      }
      else
      {
        const auto [B, B_strategy] = m_G.safe_attractor(W_1[1 - alpha], 1 - alpha, V, Vplayer);
        mCRL2log(log::trace) << "B = attractor(" << m_G.print_nodes(W_1[1 - alpha]) << ", " << m_G.print_nodes(V) << ") = " << m_G.print_nodes(B) << std::endl;
        std::tie(W[0], W[1], strategy[0], strategy[1]) = zielonka(minus(V, B));
        W[1 - alpha] = union_(W[1 - alpha], B);
        strategy[1 - alpha] = union_(union_(strategy_1[1 - alpha], B_strategy), strategy[1 - alpha]);
      }

      mCRL2log(log::debug) << "finished zielonka recursion (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";

      mCRL2log(log::trace) << "\n  --- zielonka solution for ---\n" << m_G.print_graph(V) << std::endl;
      mCRL2log(log::trace) << "W0 = " << m_G.print_nodes(W[0]) << std::endl;
      mCRL2log(log::trace) << "W1 = " << m_G.print_nodes(W[1]) << std::endl;
      mCRL2log(log::trace) << "S0 = " << m_G.print_strategy(strategy[0]) << std::endl;
      mCRL2log(log::trace) << "S1 = " << m_G.print_strategy(strategy[1]) << std::endl;

      assert(union_(W[0], W[1]) == V);
      return { W[0], W[1], strategy[0], strategy[1] };
    }

  public:

    /// \brief Solve the given game restricted to V with Zielonka's recursive algorithm as solver.
    ///        The remaining parameters are sinks, vertices won by even and odd respectively.
    ///        Terminates early when initial_vertex has been solved.
    /// \returns The winner and W0, W1, S0, S1. Where S0 and S1 are the strategies.
    std::tuple<bool, ldd, ldd, ldd, ldd> solve(const ldd& initial_vertex,
        const ldd& V,
        const ldd& Vsinks = sylvan::ldds::empty_set(),
        const ldd& W0 = sylvan::ldds::empty_set(),
        const ldd& W1 = sylvan::ldds::empty_set())
    {
      using namespace sylvan::ldds;
      stopwatch timer;
      
      std::array<ldd, 2> winning = { W0, W1 };
      std::array<ldd, 2> strategy;

      ldd Vtotal = m_G.compute_total_graph(V, empty_set(), Vsinks, winning, strategy);

      if (includes(winning[0], initial_vertex))
      {
        mCRL2log(log::verbose) << "finished solving (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";
        mCRL2log(log::trace) << "W0 = " << m_G.print_nodes(winning[0]) << std::endl;
        mCRL2log(log::trace) << "W1 = " << m_G.print_nodes(winning[1]) << std::endl;
        mCRL2log(log::trace) << "S0 = " << m_G.print_strategy(strategy[0]) << std::endl;
        mCRL2log(log::trace) << "S1 = " << m_G.print_strategy(strategy[1]) << std::endl;
        return std::make_tuple(true, winning[0], winning[1], strategy[0], strategy[1]);
      }
      else if (includes(winning[1], initial_vertex))
      {
        mCRL2log(log::verbose) << "finished solving (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";
        mCRL2log(log::trace) << "W0 = " << m_G.print_nodes(winning[0]) << std::endl;
        mCRL2log(log::trace) << "W1 = " << m_G.print_nodes(winning[1]) << std::endl;
        mCRL2log(log::trace) << "S0 = " << m_G.print_strategy(strategy[0]) << std::endl;
        mCRL2log(log::trace) << "S1 = " << m_G.print_strategy(strategy[1]) << std::endl;
        return std::make_tuple(false, winning[0], winning[1], strategy[0], strategy[1]);
      }

      // If the initial vertex has not yet been won then run the zielonka solver as well.
      mCRL2log(log::trace) << "\n--- apply zielonka to ---\n" << m_G.print_graph(V) << std::endl;
      auto [solved0, solved1, strategy0, strategy1] = zielonka(Vtotal);

      // Ensure that previously solved sets are included.
      solved0 = union_(solved0, winning[0]);
      solved1 = union_(solved1, winning[1]);
      strategy0 = union_(strategy0, strategy[0]);
      strategy1 = union_(strategy1, strategy[1]);

      mCRL2log(log::verbose) << "finished solving (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";
      mCRL2log(log::trace) << "W0 = " << m_G.print_nodes(solved0) << std::endl;
      mCRL2log(log::trace) << "W1 = " << m_G.print_nodes(solved1) << std::endl;
      mCRL2log(log::trace) << "S0 = " << m_G.print_strategy(strategy0) << std::endl;
      mCRL2log(log::trace) << "S1 = " << m_G.print_strategy(strategy1) << std::endl;

      if (includes(solved0, initial_vertex))
      {
        if (m_check_strategy)
        {
          check_strategy(initial_vertex, V, Vsinks, solved0, solved1, false, strategy0);
        }
        return std::make_tuple(true, solved0, solved1, strategy0, strategy1);
      }
      else if (includes(solved1, initial_vertex))
      {
        if (m_check_strategy)
        {
          check_strategy(initial_vertex, V, Vsinks, solved0, solved1, true, strategy1);
        }
        return std::make_tuple(false, solved0, solved1, strategy0, strategy1);
      }
      else
      {
        throw mcrl2::runtime_error("No solution found!");
      }
    }

    /// \brief Solve the given incomplete parity game (m_G, I) restricted to V.
    ///        The remaining parameters are sinks, vertices won by even and odd respectively.
    ///        Terminates early when initial_vertex has been solved.
    std::pair<const ldd, const ldd> partial_solve(const ldd& initial_vertex,
      const ldd& V,
      const ldd& I,
      const ldd& Vsinks = sylvan::ldds::empty_set(),
      const ldd& W0 = sylvan::ldds::empty_set(),
      const ldd& W1 = sylvan::ldds::empty_set())
    {
      // Make the game total.
      using namespace sylvan::ldds;
      std::array<ldd, 2> winning = { W0, W1 };
      std::array<ldd, 2> strategy;
      ldd Vtotal = m_G.compute_total_graph(V, I, Vsinks, winning, strategy);
      if (includes(winning[0], initial_vertex) || includes(winning[1], initial_vertex))
      {
        return { winning[0], winning[1] };
      }

      // Solve with zielonka twice for the safe sets.
      ldd solved0 = union_(std::get<0>(zielonka(m_G.compute_safe_vertices(0, Vtotal, I))), winning[0]);
      if (includes(solved0, initial_vertex))
      {
        return { solved0, winning[1] };
      }
      ldd solved1 = union_(std::get<1>(zielonka(m_G.compute_safe_vertices(1, Vtotal, I))), winning[1]);
      return { solved0, solved1 };
    }

    /// \brief Detect solitair winning cycles for the given incomplete parity game (m_G, I) restricted to V.
    ///        The remaining parameters are sinks, vertices won by even and odd respectively.
    ///        Terminates early when initial_vertex has been solved.
    /// \param safe Whether to use the safe attractor variant (as opposed to computing safe vertices first).
    std::pair<ldd, ldd> detect_solitair_cycles(const ldd& initial_vertex,
      const ldd& V,
      const ldd& I,
      bool safe_variant,
      const ldd& Vsinks,
      const ldd& W0 = sylvan::ldds::empty_set(),
      const ldd& W1 = sylvan::ldds::empty_set())
    {
      using namespace sylvan::ldds;

      // Make the game total and removed winning sets.
      std::array<ldd, 2> winning = { W0, W1 };
      std::array<ldd, 2> strategy;

      ldd Vtotal = m_G.compute_total_graph(V, I, Vsinks, winning, strategy);
      if (includes(winning[0], initial_vertex) || includes(winning[1], initial_vertex))
      {
        return { winning[0], winning[1] };
      }

      mCRL2log(log::trace) << "\n--- apply solitair winning cycle detection to ---\n" << m_G.print_graph(V) << std::endl;

      // Computes two vertex sets of all even priority and odd priority nodes respectively.
      std::array<ldd, 2> parity;
      std::array<const ldd, 2> Vplayer = m_G.players(Vtotal);
      for (const auto&[rank, Vrank] : m_G.ranks())
      {
        parity[rank % 2] = union_(parity[rank % 2], Vrank);
      }

      std::array<ldd, 2> Vsafe;
      if (!safe_variant)
      {
        Vsafe = { m_G.compute_safe_vertices(0, Vtotal, I), m_G.compute_safe_vertices(1, Vtotal, I) };
      }

      for (std::size_t alpha = 0; alpha <= 1; ++alpha)
      {
        // Determine the cycles for this player.
        ldd U = empty_set();
        ldd Unext = intersect(parity[alpha], Vplayer[alpha]);
        if (!safe_variant)
        {
          Unext = intersect(Unext, Vsafe[alpha]);
        }

        mCRL2log(log::debug) << "solitair winning cycle detection for player " << alpha << "\n";

        std::size_t iter = 0;
        while (U != Unext)
        {
          stopwatch timer;
          U = Unext;
          Unext = m_G.predecessors(U, U);

          mCRL2log(log::debug) << "iteration " << iter << " (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";

          ++iter;
        }

        mCRL2log(log::debug) << "found " << std::setw(12) << satcount(U) << " states in cycles for player " << alpha << "\n";

        if (safe_variant)
        {
          winning[alpha] = union_(winning[alpha], m_G.safe_attractor(U, alpha, Vtotal, Vplayer, I).first);
        }
        else
        {
          winning[alpha] = union_(winning[alpha], m_G.safe_attractor(U, alpha, Vsafe[alpha], Vplayer).first);
        }
      }

      mCRL2log(log::trace) << "W0 = " << m_G.print_nodes(winning[0]) << std::endl;
      mCRL2log(log::trace) << "W1 = " << m_G.print_nodes(winning[1]) << std::endl;

      return { winning[0], winning[1] };
    }

    /// \brief Detect forced winning cycles for the given incomplete parity game (m_G, I) restricted to V.
    ///        The remaining parameters are sinks, vertices won by even and odd respectively.
    ///        Terminates early when initial_vertex has been solved.
    /// \param safe Whether to use the safe attractor variant (as opposed to computing safe vertices first).
    std::pair<ldd, ldd> detect_forced_cycles(const ldd& initial_vertex,
      const ldd& V,
      const ldd& I,
      bool safe_variant,
      const ldd& Vsinks,
      const ldd& W0 = sylvan::ldds::empty_set(),
      const ldd& W1 = sylvan::ldds::empty_set())
    {
      using namespace sylvan::ldds;

      // Make the game total and removed winning sets.
      std::array<ldd, 2> winning = { W0, W1 };
      std::array<ldd, 2> strategy;

      ldd Vtotal = m_G.compute_total_graph(V, I, Vsinks, winning, strategy);
      if (includes(winning[0], initial_vertex) || includes(winning[1], initial_vertex))
      {
        return { winning[0], winning[1] };
      }

      mCRL2log(log::trace) << "\n--- apply forced winning cycle detection to ---\n" << m_G.print_graph(V) << std::endl;

      // Computes two vertex sets of all even priority and odd priority nodes respectively.
      std::array<ldd, 2> parity;
      std::array<const ldd, 2> Vplayer = m_G.players(Vtotal);
      for (const auto&[rank, Vrank] : m_G.ranks())
      {
        parity[rank % 2] = union_(parity[rank % 2], Vrank);
      }

      std::array<ldd, 2> Vsafe;
      if (!safe_variant)
      {
        Vsafe = { m_G.compute_safe_vertices(0, Vtotal, I), m_G.compute_safe_vertices(1, Vtotal, I) };
      }

      for (std::size_t alpha = 0; alpha <= 1; ++alpha)
      {
        // Determine the cycles for this player.
        ldd U = empty_set();
        ldd Unext = parity[alpha];
        if (!safe_variant)
        {
          Unext = intersect(Unext, Vsafe[alpha]);
        }

        mCRL2log(log::debug) << "forced winning cycle detection for player " << alpha << "\n";

        std::size_t iter = 0;
        while (U != Unext)
        {
          stopwatch timer;
          U = Unext;
          if (safe_variant)
          {
            Unext = intersect(U, m_G.safe_control_predecessors(alpha, U, Vtotal, U, Vplayer, I));
          }
          else
          {
            Unext = intersect(U, m_G.safe_control_predecessors(alpha, U, Vsafe[alpha], U, Vplayer));
          }

          mCRL2log(log::debug) << "iteration " << iter << " (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";

          ++iter;
        }

        mCRL2log(log::debug) << "found " << std::setw(12) << satcount(U) << " states in cycles for player " << alpha << "\n";

        if (safe_variant)
        {
          winning[alpha] = union_(winning[alpha], m_G.safe_attractor(U, alpha, Vtotal, Vplayer, I).first);
        }
        else
        {
          winning[alpha] = union_(winning[alpha], m_G.safe_attractor(U, alpha, Vsafe[alpha], Vplayer).first);
        }
      }

      mCRL2log(log::trace) << "W0 = " << m_G.print_nodes(winning[0]) << std::endl;
      mCRL2log(log::trace) << "W1 = " << m_G.print_nodes(winning[1]) << std::endl;

      return { winning[0], winning[1] };
    }

    /// \returns Partial solve using the fatal attractors.
    std::pair<ldd, ldd> detect_fatal_attractors(const ldd& initial_vertex,
        const ldd& V,
        const ldd& I,
        bool safe_variant,
        const ldd& Vsinks,
        const ldd& W0 = sylvan::ldds::empty_set(),
        const ldd& W1 = sylvan::ldds::empty_set())
    {
      using namespace sylvan::ldds;
      stopwatch timer;

      std::array<ldd, 2> winning = { W0, W1 };
      std::array<ldd, 2> strategy;

      ldd Vtotal = m_G.compute_total_graph(V, I, Vsinks, winning, strategy);
      if (includes(winning[0], initial_vertex) || includes(winning[1], initial_vertex))
      {
        return { winning[0], winning[1] };
      }

      std::array<const ldd, 2> Vplayer = m_G.players(Vtotal);

      // Compute safe vertices if necessary
      std::array<ldd, 2> Vsafe;
      if (!safe_variant)
      {
        Vsafe = { m_G.compute_safe_vertices(0, Vtotal, I), m_G.compute_safe_vertices(1, Vtotal, I) };
      }

      mCRL2log(log::trace) << "\n--- apply fatal attractor detection to ---\n" << m_G.print_graph(Vtotal) << std::endl;

      // For priorities in descending order
      for (auto it = m_G.ranks().rbegin(); it != m_G.ranks().rend(); it++)
      {
        std::size_t c = it->first;
        std::size_t alpha = c % 2;
        mCRL2log(log::debug) << "fatal attractor detection for priority " << c << "\n";
        ldd X = safe_variant ? it->second : intersect(it->second, Vsafe[alpha]);
        ldd Y = empty_set();

        while (X != empty_set() && X != Y)
        {
          Y = X;
          ldd Z = m_G.safe_monotone_attractor(X, alpha, c, safe_variant ? Vtotal : Vsafe[alpha], Vplayer, safe_variant ? I : empty_set());


          if (includes(Z, X))
          {
            if (safe_variant)
            {
              winning[alpha] = union_(winning[alpha], m_G.safe_attractor(Z, alpha, Vtotal, Vplayer, I).first);
            }
            else
            {
              winning[alpha] = union_(winning[alpha], m_G.safe_attractor(Z, alpha, Vsafe[alpha], Vplayer).first);
            }
            mCRL2log(log::debug) << "found " << std::setw(12) << satcount(Z) << " states in fatal attractors for priority " << c << "\n";
            break;
          }
          else
          {
            X = intersect(X, Z);
          }
        }
      }

      mCRL2log(log::debug) << "finished fatal attractor detection (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";
      mCRL2log(log::trace) << "W0 = " << m_G.print_nodes(winning[0]) << std::endl;
      mCRL2log(log::trace) << "W1 = " << m_G.print_nodes(winning[1]) << std::endl;

      return { winning[0], winning[1] };
    }

    /// Checks whether the computed strategy is indeed a correct certificate for the winning partition.
    ///
    /// Throws an exception when the strategy is invalid.
    void check_strategy(const ldd& initial_vertex,
      const ldd& V,
      const ldd& Vsinks,
      const ldd& W0, 
      const ldd& W1, 
      bool alpha, 
      const ldd& strategy)
    {
      mCRL2log(log::debug) << "Checking the strategy of the solved parity game..." << std::endl;
      symbolic_parity_game new_G = m_G.apply_strategy(alpha, strategy);
      mCRL2log(log::trace) << "Minimal parity game G = " << new_G.print_graph(V) << std::endl;
      

      symbolic_pbessolve_algorithm check(new_G);

      auto[result, W0_prime, W1_prime, S0, S1] = check.solve(initial_vertex, V, Vsinks);      
      if (!(W0 == W0_prime && W1 == W1_prime && result != alpha))
      {
        throw mcrl2::runtime_error("Computed strategy does not match the winning partition");
      }
    }
};

} // namespace mcrl2::pbes_system



#endif // MCRL2_ENABLE_SYLVAN

#endif // MCRL2_PBES_SYMBOLIC_PBESSOLVE_H
