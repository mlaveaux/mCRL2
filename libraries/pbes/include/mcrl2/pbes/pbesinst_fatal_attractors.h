// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesinst_fatal_attractors.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PBESINST_FATAL_ATTRACTORS_H
#define MCRL2_PBES_PBESINST_FATAL_ATTRACTORS_H

#include "mcrl2/pbes/pbessolve_attractors.h"
#include "mcrl2/pbes/simple_structure_graph.h"





namespace mcrl2::pbes_system::detail {

template <typename StructureGraph, typename Compare>
deque_vertex_set attr_min_rank_todo_generic(const StructureGraph& G, const vertex_set& A, const vertex_set& U, std::size_t j, Compare compare)
{
  std::size_t n = G.extent();
  deque_vertex_set todo(n);
  for (auto v: A.vertices())
  {
    for (auto u: G.predecessors(v))
    {
      if (U.contains(u) && compare(G.rank(u), j) && !A.contains(u))
      {
        todo.insert(u);
      }
    }
  }
  return todo;
}

template <typename StructureGraph>
deque_vertex_set attr_min_rank_original_todo(const StructureGraph& G, const vertex_set& A, const vertex_set& U, std::size_t j)
{
  std::size_t n = G.extent();
  deque_vertex_set todo(n);
  for (auto v: A.vertices())
  {
    for (auto u: G.predecessors(v))
    {
      if (U.contains(u) && G.rank(u) >= j)
      {
        todo.insert(u);
      }
    }
  }
  return todo;
}

inline
void insert_in_rank_map(std::map<std::size_t, vertex_set>& U_rank_map, structure_graph::index_type u, std::size_t j, std::size_t n)
{
  auto i = U_rank_map.find(j);
  if (i == U_rank_map.end())
  {
    i = U_rank_map.insert({ j, vertex_set(n) }).first;
  }
  i->second.insert(u);
}

template <typename StructureGraph>
std::map<std::size_t, vertex_set> compute_U_j_map(const StructureGraph& G, const vertex_set& V)
{
  std::size_t n = G.extent();

  // compute U_j_map, such that U_j_map[j] = U_j
  std::map<std::size_t, vertex_set> U_j_map;
  for (structure_graph::index_type u: V.vertices())
  {
    std::size_t j = G.rank(u);
    if (j == data::undefined_index())
    {
      continue;
    }
    auto alpha = j % 2;
    if ((alpha == 0 && G.decoration(u) == structure_graph::decoration_type::d_false) || (alpha == 1 && G.decoration(u) == structure_graph::decoration_type::d_true))
    {
      continue;
    }
    insert_in_rank_map(U_j_map, u, j, n);
  }
  return U_j_map;
}

// Computes an attractor set, by extending A. Only predecessors in U are considered with a rank of at least j.
// alpha = 0: disjunctive
// alpha = 1: conjunctive
// StructureGraph is either structure_graph or simple_structure_graph
template <typename StructureGraph, typename Compare>
vertex_set attr_min_rank_generic(const StructureGraph& G, vertex_set A, std::size_t alpha, const vertex_set& U, std::size_t j, Compare compare)
{
  // put all predecessors of elements in A in todo
  deque_vertex_set todo = attr_min_rank_todo_generic(G, A, U, j, compare);

  while (!todo.is_empty())
  {
    // N.B. Use a breadth first search, to minimize counter examples
    auto u = todo.pop_front();

    if (G.decoration(u) == alpha || includes_successors(G, u, A))
    {
      global_strategy<StructureGraph>(G).set_strategy(u, find_successor_in(G, u, A));

      A.insert(u);

      for (auto v: G.predecessors(u))
      {
        if (U.contains(v) && (compare(G.rank(v), j) || (G.rank(v) == data::undefined_index() && G.decoration(v) <= 1)) && !A.contains(v))
        {
          todo.insert(v);
        }
      }
    }
  }

  return A;
}

// calculation_steps is used to count how many steps are required to calculate a fatal_attractor.
template <typename Compare>
void fatal_attractors_generic(const simple_structure_graph& G,
                              std::array<vertex_set, 2>& S,
                              std::array<strategy_vector, 2>& tau,
                              std::size_t& calculation_steps,
                              std::size_t equation_count,
                              Compare compare
                             )
{
  mCRL2log(log::debug) << "\n  === fatal attractors (equation " << equation_count << ") ===\n" << G << std::endl;
  mCRL2log(log::debug) << "  S0 = " << S[0] << std::endl;
  mCRL2log(log::debug) << "  S1 = " << S[1] << std::endl;

  // count the number of insertions in the sets S[0] and S[1]
  std::size_t insertion_count = 0;
  std::size_t n = G.extent();

  // compute V
  vertex_set V(n);
  for (std::size_t u = 0; u < n; u++)
  {
    calculation_steps++;
    V.insert(u);
  }

  // compute U_j_map, such that U_j_map[j] = U_j
  std::map<std::size_t, vertex_set> U_j_map = compute_U_j_map(G, V);

  S[0] = attr_default_with_tau(G, S[0], 0, tau);
  S[1] = attr_default_with_tau(G, S[1], 1, tau);

  for (auto& p: U_j_map)
  {
    std::size_t j = p.first;
    auto alpha = j % 2;
    mCRL2log(log::debug) << "  --- iteration j = " << j << " ---" << std::endl;

    vertex_set& U_j = p.second;
    U_j = set_minus(U_j, S[1 - alpha]);
    mCRL2log(log::debug) << "  U_" << std::to_string(j) << " = " << U_j << std::endl;
    vertex_set U = set_union(U_j, S[alpha]);
    vertex_set X = detail::attr_min_rank_generic(G, U, alpha, V, j, compare);
    vertex_set Y = set_minus(V, attr_default(G, set_minus(V, X), 1 - alpha));

    while (X != Y)
    {
      calculation_steps++;
      mCRL2log(log::debug) << "  X = " << X << std::endl;
      mCRL2log(log::debug) << "  Y = " << Y << std::endl;
      X = detail::attr_min_rank_generic(G, set_intersection(U, Y), alpha, V, j, compare);
      Y = set_minus(Y, attr_default(G, set_minus(Y, X), 1 - alpha));
    }
    mCRL2log(log::debug) << "  X (final) = " << X << std::endl;

    // set strategy for v \in X \ S[alpha]
    for (structure_graph::index_type v: X.vertices())
    {
      calculation_steps++;
      if (S[alpha].contains(v))
      {
        continue;
      }
      if ((alpha == 0 && G.decoration(v) == structure_graph::d_disjunction) || (alpha == 1 && G.decoration(v) == structure_graph::d_conjunction))
      {
        if (U_j.contains(v))
        {
          auto w = find_successor_in(G, v, Y);
          local_strategy(tau, alpha).set_strategy(v, w);
        }
        else
        {
          auto tau_v = G.find_vertex(v).strategy;
          local_strategy(tau, alpha).set_strategy(v, tau_v);
        }
      }
    }

    // S_alpha := S_alpha U X
    for (structure_graph::index_type x: X.vertices())
    {
      calculation_steps++;
      insertion_count++;
      S[alpha].insert(x);
      mCRL2log(log::debug) << "  insert vertex " << x << " in S" << alpha << std::endl;
    }

    S[alpha] = attr_default_with_tau(G, S[alpha], alpha, tau);
  }
  mCRL2log(log::debug) << "\n  === result of fatal attractors (equation " << equation_count << ") ===" << std::endl;
  mCRL2log(log::debug) << "  S0 = " << S[0] << std::endl;
  mCRL2log(log::debug) << "  S1 = " << S[1] << std::endl;
  mCRL2log(log::debug) << "  tau0 = " << print_strategy_vector(S[0], tau[0]) << std::endl;
  mCRL2log(log::debug) << "  tau1 = " << print_strategy_vector(S[1], tau[1]) << std::endl;
  mCRL2log(log::debug) << "  inserted " << insertion_count << " vertices." << std::endl;
}

inline
void fatal_attractors(const simple_structure_graph& G,
                      std::array<vertex_set, 2>& S,
                      std::array<strategy_vector, 2>& tau,
                      std::size_t& calculation_steps,
                      std::size_t equation_count
                     )
{
  fatal_attractors_generic(G, S, tau, calculation_steps, equation_count, std::greater_equal<>());
}

// calculation_steps returns how many steps were needed to find the loops. 
inline
void find_loops2(const simple_structure_graph& G,
                 std::array<vertex_set, 2>& S,
                 std::array<strategy_vector, 2>& tau,
                 std::size_t& calculation_steps,
                 std::size_t equation_count
)
{
  fatal_attractors_generic(G, S, tau, calculation_steps, equation_count, std::equal_to<>());
}

// Computes an attractor set, by extending A. Only predecessors in U are considered with a rank of at least j.
// alpha = 0: disjunctive
// alpha = 1: conjunctive
// StructureGraph is either structure_graph or simple_structure_graph
// This version is based on the work of Michael Huth, Jim Huan-Pu Kuo, and Nir Piterman.
template <typename StructureGraph>
vertex_set attr_min_rank_original(const StructureGraph& G, vertex_set A, std::size_t alpha, const vertex_set& U, std::size_t j)
{
  // put all predecessors of elements in A in todo
  deque_vertex_set todo = attr_min_rank_original_todo(G, A, U, j);

  vertex_set X(A.extent());
  for (auto u: todo.vertices())
  {
    if (A.contains(u) && (G.decoration(u) == alpha || includes_successors(G, u, A)))
    {
      X.insert(u);
    }
  }

  while (!todo.is_empty())
  {
    // N.B. Use a breadth first search, to minimize counter examples
    auto u = todo.pop_front();

    if (G.decoration(u) == alpha || includes_successors(G, u, lazy_union(A, X)))
    {
      global_strategy<StructureGraph>(G).set_strategy(u, find_successor_in(G, u, lazy_union(A, X)));

      X.insert(u);

      for (auto v: G.predecessors(u))
      {
        if (U.contains(v) && (G.rank(v) >= j || (G.rank(v) == data::undefined_index() && G.decoration(v) <= 1)) && !X.contains(v))
        {
          todo.insert(v);
        }
      }
    }
  }

  return X;
}

inline
void fatal_attractors_original(const simple_structure_graph& G,
                               std::array<vertex_set, 2>& S,
                               std::array<strategy_vector, 2>& tau,
                               std::size_t& calculation_steps,
                               std::size_t equation_count
)
{
  mCRL2log(log::debug) << "\n  === fatal attractors original (equation " << equation_count << ") ===\n" << G << std::endl;
  mCRL2log(log::debug) << "  S0 = " << S[0] << std::endl;
  mCRL2log(log::debug) << "  S1 = " << S[1] << std::endl;

  // count the number of insertions in the sets S[0] and S[1]
  std::size_t insertion_count = 0;
  std::size_t n = G.extent();

  // compute V
  vertex_set V(n);
  calculation_steps=+n;
  for (std::size_t u = 0; u < n; u++)
  {
    V.insert(u);
  }

  // compute U_j_map, such that U_j_map[j] = U_j
  std::map<std::size_t, vertex_set> U_j_map = compute_U_j_map(G, V);

  S[0] = attr_default_with_tau(G, S[0], 0, tau);
  S[1] = attr_default_with_tau(G, S[1], 1, tau);

  for (auto& p: U_j_map)
  {
    std::size_t j = p.first;
    auto alpha = j % 2;
    mCRL2log(log::debug) << "  --- iteration j = " << j << " ---" << std::endl;

    vertex_set& U_j = p.second;
    U_j = set_minus(U_j, S[1 - alpha]);

    vertex_set X(n);

    while (!U_j.is_empty() && U_j != X)
    {
      mCRL2log(log::debug) << "  U_" + std::to_string(j) << " = " << U_j << std::endl;
      X = U_j;
      mCRL2log(log::debug) << "  X = " << X << std::endl;
      mCRL2log(log::debug) << "  X U S_" + std::to_string(alpha) << " = " << set_union(X, S[alpha]) << std::endl;
      vertex_set Y = detail::attr_min_rank_original(G, set_union(X, S[alpha]), alpha, V, j);
      mCRL2log(log::debug) << "  Y = " << Y << std::endl;
      mCRL2log(log::debug) << "  U_" + std::to_string(j) << " is " << (is_subset_of(U_j, Y) ? "a" : "no") << " subset of Y" << std::endl;
      if (is_subset_of(U_j, Y))
      {
        // set strategy for v \in Y \ S[alpha]
        for (structure_graph::index_type v: Y.vertices())
        {
          calculation_steps++;
          if (S[alpha].contains(v))
          {
            continue;
          }
          if ((alpha == 0 && G.decoration(v) == structure_graph::d_disjunction) || (alpha == 1 && G.decoration(v) == structure_graph::d_conjunction))
          {
            if (U_j.contains(v))
            {
              auto w = find_successor_in(G, v, Y);
              local_strategy(tau, alpha).set_strategy(v, w);
            }
            else
            {
              auto tau_v = G.find_vertex(v).strategy;
              local_strategy(tau, alpha).set_strategy(v, tau_v);
            }
          }
        }

        // S_alpha := S_alpha U Y
        for (structure_graph::index_type y: Y.vertices())
        {
          calculation_steps++;
          insertion_count++;
          S[alpha].insert(y);
          mCRL2log(log::debug) << "  insert vertex " << y << " in S" << alpha << std::endl;
        }

        S[alpha] = attr_default_with_tau(G, S[alpha], alpha, tau);
        break;
      }
      else
      {
        U_j = set_intersection(U_j, Y);
      }
    }
  }
  mCRL2log(log::debug) << "\n  === result of fatal attractors original (equation " << equation_count << ") ===" << std::endl;
  mCRL2log(log::debug) << "  S0 = " << S[0] << std::endl;
  mCRL2log(log::debug) << "  S1 = " << S[1] << std::endl;
  mCRL2log(log::debug) << "  tau0 = " << print_strategy_vector(S[0], tau[0]) << std::endl;
  mCRL2log(log::debug) << "  tau1 = " << print_strategy_vector(S[1], tau[1]) << std::endl;
  mCRL2log(log::debug) << "  inserted " << insertion_count << " vertices." << std::endl;
}

} // namespace mcrl2::pbes_system::detail





#endif // MCRL2_PBES_PBESINST_FATAL_ATTRACTORS_H
