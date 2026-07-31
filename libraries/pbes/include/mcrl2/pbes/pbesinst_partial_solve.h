// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesinst_partial_solve.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PBESINST_PARTIAL_SOLVE_H
#define MCRL2_PBES_PBESINST_PARTIAL_SOLVE_H

#include "mcrl2/pbes/pbesinst_lazy.h"
#include "mcrl2/pbes/simple_structure_graph.h"
#include "mcrl2/pbes/solve_structure_graph.h"
#include "mcrl2/pbes/structure_graph_builder.h"





namespace mcrl2::pbes_system::detail {

inline
void partial_solve(structure_graph& G,
                   const pbesinst_lazy_todo& todo,
                   std::array<vertex_set, 2>& S,
                   std::array<strategy_vector, 2>& tau,
                   std::size_t& calculation_steps,
                   std::size_t equation_count,
                   const detail::structure_graph_builder& graph_builder
                  )
{
  mCRL2log(log::debug) << "\n  === partial solve (equation " << equation_count << ") ===\n" << G << std::endl;
  mCRL2log(log::debug) << "  S0 = " << S[0] << std::endl;
  mCRL2log(log::debug) << "  S1 = " << S[1] << std::endl;

  std::size_t N = G.extent();

  // The size of S[0] and S[1] may be bigger than N, because of the resize strategy in vertex_set.
  S[0].truncate(N);
  S[1].truncate(N);

  mCRL2log(log::debug) << "  computing S0 = attr_default_with_tau(G, S0, 0, tau0)" << std::endl;
  S[0] = attr_default_with_tau(G, S[0], 0, tau);
  mCRL2log(log::debug) << "  computing S1 = attr_default_with_tau(G, S1, 1, tau1)" << std::endl;
  S[1] = attr_default_with_tau(G, S[1], 1, tau);

  // Si_todo := Si U todo
  std::array<vertex_set, 2> S_todo = S;

  for (const propositional_variable_instantiation& X: todo.elements())
  {
    calculation_steps++;
    structure_graph::index_type u = graph_builder.find_vertex(X);
    S_todo[0].insert(u);
    S_todo[1].insert(u);
  }

  bool check_strategy = false;
  bool use_toms_optimization = false;
  solve_structure_graph_algorithm algorithm(check_strategy, use_toms_optimization);

  vertex_set W[2] = { vertex_set(N), vertex_set(N) }; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
  std::tie(W[0], W[1]) = algorithm.solve_recursive(G, set_union(S[1], attr_default_no_strategy(G, S_todo[0], 0)), calculation_steps);
  for (structure_graph::index_type v: W[1].vertices())
  {
    calculation_steps++;
    if (S[1].contains(v))
    {
      continue;
    }
    if (G.decoration(v) == structure_graph::d_conjunction)
    {
      auto tau_v = G.decoration(v);
      local_strategy(tau, 1).set_strategy(v, tau_v);
    }
  }
  std::tie(W[0], W[1]) = algorithm.solve_recursive(G, set_union(S[0], attr_default_no_strategy(G, S_todo[1], 1)), calculation_steps);
  for (structure_graph::index_type v: W[0].vertices())
  {
    calculation_steps++;
    if (S[0].contains(v))
    {
      continue;
    }
    if (G.decoration(v) == structure_graph::d_disjunction)
    {
      auto tau_v = G.decoration(v);
      local_strategy(tau, 0).set_strategy(v, tau_v);
    }
  }

  mCRL2log(log::debug) << "\n  === result of partial solve (iteration " << equation_count << ") ===" << std::endl;
  mCRL2log(log::debug) << "  S0 = " << S[0] << std::endl;
  mCRL2log(log::debug) << "  S1 = " << S[1] << std::endl;
  mCRL2log(log::debug) << "  tau0 = " << print_strategy_vector(S[0], tau[0]) << std::endl;
  mCRL2log(log::debug) << "  tau1 = " << print_strategy_vector(S[1], tau[1]) << std::endl;
}

} // namespace mcrl2::pbes_system::detail





#endif // MCRL2_PBES_PBESINST_PARTIAL_SOLVE_H
