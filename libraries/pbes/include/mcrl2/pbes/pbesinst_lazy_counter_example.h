// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesinst_lazy_counter_example.h
/// \brief A lazy algorithm for instantiating a PBES, ported from bes_deprecated.h.

#ifndef MCRL2_PBES_PBESINST_LAZY_COUNTER_EXAMPLE_H
#define MCRL2_PBES_PBESINST_LAZY_COUNTER_EXAMPLE_H

#include "mcrl2/pbes/pbesinst_structure_graph.h"
#include "mcrl2/pbes/pbesinst_structure_graph2.h"
#include "mcrl2/pbes/replace.h"
#include "mcrl2/pbes/rewriters/simplify_quantifiers_rewriter.h"
#include "mcrl2/pbes/structure_graph.h"

#include <regex>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system
{

/// Replaces propositional variables in the expression psi that are irrelevant for the given proof_graph.
static void rewrite_star(pbes_expression& result,
    const fixpoint_symbol& /* symbol */,
    const propositional_variable_instantiation& X,
    const pbes_expression& psi,
    const structure_graph& G,
    bool alpha,
    const std::unordered_map<pbes_expression, structure_graph::index_type>& mapping)
{
  bool changed = false;
  std::smatch match;

  // Now we need to find all reachable X --> Y, following vertices that are not ranked.
  mCRL2log(log::debug) << "X = " << X << ", psi = " << psi << std::endl;

  std::unordered_set<pbes_expression> Ys;

  // If X is won by player alpha, i.e. in the winning set W.
  auto it = mapping.find(X);
  if (it != mapping.end())
  {
    structure_graph::index_type index = it->second;
    std::unordered_set<structure_graph::index_type> todo = {index};
    std::unordered_set<structure_graph::index_type> done;

    while (!todo.empty())
    {
      structure_graph::index_type u = *todo.begin();
      todo.erase(todo.begin());
      done.insert(u);

      if (mapping.count(G.find_vertex(u).formula()) != 0)
      {
        // This vertex is won by alpha
        if (G.strategy(u) != undefined_vertex() &&
            ((!alpha && G.decoration(u) == structure_graph::d_disjunction)
            || (alpha && G.decoration(u) == structure_graph::d_conjunction)))
        {
          // The strategy is defined so only explore the strategy edge.
          auto v = G.strategy(u);
          if (G.rank(v) == data::undefined_index())
          {
            if (!mcrl2::utilities::detail::contains(done, v))
            {
              // explore all outgoing edges that are unranked
              todo.insert(v);
            }
          }
          else if (mapping.count(G.find_vertex(v).formula()) != 0)
          {
            // Insert the outgoing edge, but do not add it to the todo set to stop exploring this vertex.
            Ys.insert(G.find_vertex(v).formula());
          }
        }
        else
        {
          // Explore all edges.
          for (structure_graph::index_type v : G.all_successors(u))
          {
            if (G.rank(v) == data::undefined_index())
            {
              if (!mcrl2::utilities::detail::contains(done, v))
              {
                // explore all outgoing edges that are unranked
                todo.insert(v);
              }
            }
            else if (mapping.count(G.find_vertex(v).formula()) != 0)
            {
              // Insert the outgoing edge, but do not add it to the todo set to stop exploring this vertex.
              Ys.insert(G.find_vertex(v).formula());
            }
          }
        }
      }
      else
      {
        // Explore all edges.
        for (structure_graph::index_type v : G.all_successors(u))
        {
          if (G.rank(v) == data::undefined_index())
          {
            if (!mcrl2::utilities::detail::contains(done, v))
            {
              // explore all outgoing edges that are unranked
              todo.insert(v);
            }
          }
          else if (mapping.count(G.find_vertex(v).formula()) != 0)
          {
            // Insert the outgoing edge, but do not add it to the todo set to stop exploring this vertex.
            Ys.insert(G.find_vertex(v).formula());
          }
        }
      }
    }
  }

  mCRL2log(log::debug) << "Ys := " << core::detail::print_set(Ys) << std::endl;

  replace_propositional_variables(result,
      psi,
      [&](const propositional_variable_instantiation& Y) -> pbes_expression
      {
        if (std::regex_match(static_cast<const std::string&>(Y.name()),
                match,
                mcrl2::pbes_system::detail::positive_or_negative))
        {
          // If Y in L return Y
          mCRL2log(log::debug) << "rewrite_star " << Y << " is counter example equation (in L)" << std::endl;
          return Y;
        }
        else
        {
          if (mcrl2::utilities::detail::contains(Ys, Y))
          {
            mCRL2log(log::debug) << "rewrite_star " << Y << " is reachable" << std::endl;
            return Y;
          }
          else
          {
            changed = true;
            if (alpha == 0)
            {
              // If Y is not reachable, replace it by false
              mCRL2log(log::debug) << "rewrite_star " << Y << " is not reachable, becomes false" << std::endl;
              return false_();
            }
            else
            {
              // If Y is not reachable, replace it by true
              mCRL2log(log::debug) << "rewrite_star " << Y << " is not reachable, becomes true" << std::endl;
              return true_();
            }
          }
        }
      });

  if (changed)
  {
    simplify_rewriter simplify;
    const pbes_expression result1 = result;
    simplify(result, result1);
  }

  mCRL2log(log::debug) << "result = " << psi << std::endl;
}

class pbesinst_counter_example_structure_graph_algorithm : public pbesinst_structure_graph_algorithm
{
public:
  pbesinst_counter_example_structure_graph_algorithm(const pbessolve_options& options,
      const pbes& p,
      const structure_graph& SG,
      bool _alpha,
      const std::unordered_map<pbes_expression, structure_graph::index_type>& _mapping,
      structure_graph& G,
      std::optional<data::rewriter> rewriter = std::nullopt)
      : pbesinst_structure_graph_algorithm(options, p, G, rewriter),
        G(SG),
        alpha(_alpha),
        mapping(_mapping)
  {}

  void rewrite_psi(const std::size_t thread_index,
      pbes_expression& result,
      const fixpoint_symbol& symbol,
      const propositional_variable_instantiation& X,
      const pbes_expression& psi) override
  {
    rewrite_star(result, symbol, X, psi, G, alpha, mapping);
    pbesinst_structure_graph_algorithm::rewrite_psi(thread_index, result, symbol, X, psi);
  }

private:
  const structure_graph& G;
  bool alpha;
  const std::unordered_map<pbes_expression, structure_graph::index_type>& mapping;
};

class pbesinst_counter_example_structure_graph_algorithm2 : public pbesinst_structure_graph_algorithm2
{
public:
  pbesinst_counter_example_structure_graph_algorithm2(const pbessolve_options& options,
      const pbes& p,
      const structure_graph& SG,
      bool _alpha,
      const std::unordered_map<pbes_expression, structure_graph::index_type>& _mapping,
      structure_graph& G,
      std::optional<data::rewriter> rewriter = std::nullopt)
      : pbesinst_structure_graph_algorithm2(options, p, G, rewriter),
        G(SG),
        alpha(_alpha),
        mapping(_mapping)
  {}

  void rewrite_psi(const std::size_t thread_index,
      pbes_expression& result,
      const fixpoint_symbol& symbol,
      const propositional_variable_instantiation& X,
      const pbes_expression& psi) override
  {
    rewrite_star(result, symbol, X, psi, G, alpha, mapping);
    pbesinst_structure_graph_algorithm2::rewrite_psi(thread_index, result, symbol, X, psi);
  }

private:
  const structure_graph& G;
  bool alpha;
  const std::unordered_map<pbes_expression, structure_graph::index_type>& mapping;
};

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_PBESINST_LAZY_COUNTER_EXAMPLE_H