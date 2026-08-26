// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/solve_structure_graph.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_SOLVE_STRUCTURE_GRAPH_H
#define MCRL2_PBES_SOLVE_STRUCTURE_GRAPH_H

#include "mcrl2/atermpp/standard_containers/vector.h"
#include "mcrl2/data/join.h"
#include "mcrl2/lts/lts_algorithm.h"
#include "mcrl2/pbes/pbes_equation_index.h"
#include "mcrl2/pbes/pbessolve_attractors.h"
#include "mcrl2/pbes/detail/pbes_remove_counterexample_info.h"

namespace mcrl2::pbes_system {

namespace detail {

/// \brief Walks the extracted evidence vertex set \a V of a structure graph that was built from
/// the lps2pbes -c encoding, and calls \a yield once for every Zpos_i/Zneg_i vertex it contains,
/// with the summand it originates from and the transition it encodes. Shared between the
/// LPS-evidence and LTS-evidence sinks so they cannot drift apart.
/// \param yield A callback with signature
///        void(std::size_t summand_index, const data::data_expression_vector& source,
///             const lps::multi_action& a, const data::data_expression_vector& target)
///        where \a source and \a target have lpsspec.process().process_parameters().size() elements,
///        positionally corresponding to those process parameters.
template <typename Callback>
void for_each_evidence_transition(structure_graph& G,
                                   const std::set<structure_graph::index_type>& V,
                                   const lps::specification& lpsspec,
                                   const pbes& p,
                                   const pbes_equation_index& p_index,
                                   Callback yield)
{
  std::size_t n = lpsspec.process().process_parameters().size();

  for (structure_graph::index_type vi: V)
  {
    const auto& v = G.find_vertex(vi);
    if (!is_propositional_variable_instantiation(v.formula()))
    {
      continue;
    }

    // The variable Z below should be a reference, but this leads to crashes with the GCC compiler (March 2022).
    // JFG: I think this is a GCC problem, which may resolve itself in due time.
    const auto Z = atermpp::down_cast<propositional_variable_instantiation>(v.formula());
    std::size_t summand_index;
    if (!counter_example_index(Z.name(), summand_index))
    {
      continue;
    }

    if (summand_index >= lpsspec.process().action_summands().size())
    {
      throw mcrl2::runtime_error("Counter-example cannot be reconstructed from this LPS. Did you supply the correct file?");
    }

    // The parameters are [from] + [action_parameters] + [to]
    const lps::action_summand& summand = lpsspec.process().action_summands().at(summand_index);
    std::size_t equation_index = p_index.index(Z.name());
    const pbes_equation& eqn = p.equations().at(equation_index);
    const data::variable_list& d = eqn.variable().parameters();
    std::size_t m = d.size() - 2 * n;

    const data::data_expression_list& e = Z.parameters();
    data::data_expression_vector e1(e.begin(), e.end());

    data::data_expression_vector source(e1.begin(), e1.begin() + static_cast<std::ptrdiff_t>(n));
    data::data_expression_vector target(e1.begin() + static_cast<std::ptrdiff_t>(n + m), e1.begin() + static_cast<std::ptrdiff_t>(2 * n + m));

    process::action_vector actions;
    std::size_t index = 0;
    for (const process::action& a: summand.multi_action().actions())
    {
      if (index > e1.size() || index + a.arguments().size() > e1.size())
      {
        throw mcrl2::runtime_error("Invalid parameter index");
      }

      process::action a1(a.label(), data::data_expression_list(e1.begin() + static_cast<std::ptrdiff_t>(n + index), e1.begin() + static_cast<std::ptrdiff_t>(n + index + a.arguments().size())));
      actions.push_back(a1);
      index = index + a.arguments().size();
    }

    lps::multi_action a(process::action_list(actions.begin(), actions.end()), summand.multi_action().time());
    yield(summand_index, source, a, target);
  }
}

} // namespace detail

inline
std::tuple<std::size_t, std::size_t, vertex_set> get_minmax_rank(const structure_graph& G)
{
  std::size_t min_rank = (std::numeric_limits<std::size_t>::max)();
  std::size_t max_rank = 0;
  std::vector<structure_graph::index_type> M; // vertices with minimal rank
  std::size_t N = G.all_vertices().size();

  for (std::size_t vi = 0; vi < N; vi++)
  {
    if (!G.contains(vi))
    {
      continue;
    }
    const auto& v = G.find_vertex(vi);
    if (v.rank <= min_rank)
    {
      if (v.rank < min_rank)
      {
        M.clear();
        min_rank = v.rank;
      }
      M.push_back(vi);
    }
    if (v.rank > max_rank)
    {
      max_rank = v.rank;
    }
  }
  return std::make_tuple(min_rank, max_rank, vertex_set(N, M.begin(), M.end()));
}

class solve_structure_graph_algorithm
{
  protected:
    // do a sanity check on the computed strategy
    bool check_strategy = false;

    bool use_toms_optimization = false;

    // find a successor of u
    static structure_graph::index_type succ(const structure_graph& G, structure_graph::index_type u)
    {
      for (structure_graph::index_type v: G.successors(u))
      {
        return v;
      }
      return undefined_vertex();
    }

    // find a successor of u in U, or a random one if no successor in U exists
    static inline
    structure_graph::index_type succ(const structure_graph& G, structure_graph::index_type u, const vertex_set& U)
    {
      auto result = undefined_vertex();
      for (structure_graph::index_type v: G.successors(u))
      {
        if (U.contains(v))
        {
          return v;
        }
        result = v;
      }
      return result;
    }

  public:
    // computes solve_recursive(G \ A)
    inline
    std::pair<vertex_set, vertex_set> solve_recursive(structure_graph& G, const vertex_set& A, std::size_t& calculation_steps)
    {
      calculation_steps++;
      auto exclude = G.exclude() | A.include();
      std::swap(G.exclude(), exclude);
      auto result = solve_recursive(G, calculation_steps);
      std::swap(G.exclude(), exclude);
      return result;
    }

  protected:
    // pre: G does not contain nodes with decoration true or false.
    //
    // N.B. If use_toms_optimization is true, then the oomputed strategy may be incorrect.
    // So this flag should only be used to compute the solution.
    inline
    std::pair<vertex_set, vertex_set> solve_recursive(structure_graph& G, std::size_t& calculation_steps)
    {
      calculation_steps++;
      mCRL2log(log::debug) << "\n  --- solve_recursive input ---\n" << G << std::endl;
      std::size_t N = G.extent();

      if (G.is_empty())
      {
        return { vertex_set(N), vertex_set(N) };
      }

      auto q = get_minmax_rank(G);
      std::size_t m = std::get<0>(q);
      const vertex_set& U = std::get<2>(q);

      std::size_t alpha = m % 2; // 0 = disjunctive, 1 = conjunctive

      // set strategy
      for (structure_graph::index_type ui: U.vertices())
      {
        calculation_steps++;
        const auto& u = G.find_vertex(ui);
        if (u.decoration == alpha)
        {
          // auto v = succ(G, ui); // N.B. this may lead to a wrong strategy!
          auto v = succ(G, ui, U);
          if (v != undefined_vertex())
          {
            global_strategy<structure_graph>(G).set_strategy(ui, v);
          }
        }
      }

      vertex_set W[2]   = { vertex_set(N), vertex_set(N) }; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
      vertex_set W_1[2]; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

      vertex_set A = attr_default(G, U, alpha);
      std::tie(W_1[0], W_1[1]) = solve_recursive(G, A, calculation_steps);

      if (use_toms_optimization)
      {
        // More efficient than Zielonka, because some recursive calls are skipped.
        // As a consequence, the computed strategy may be wrong.
        vertex_set B = attr_default(G, W_1[1 - alpha], 1 - alpha);
        if (W_1[1 - alpha].size() == B.size())
        {
          W[alpha] = set_union(A, W_1[alpha]);
          W[1 - alpha] = B;
        }
        else
        {
          std::tie(W[0], W[1]) = solve_recursive(G, B, calculation_steps);
          W[1 - alpha] = set_union(W[1 - alpha], B);
        }
      }
      else
      {
         // Original Zielonka version
         if (W_1[1 - alpha].is_empty())
         {
           W[alpha] = set_union(A, W_1[alpha]);
           W[1 - alpha].clear();
         }
         else
         {
           vertex_set B = attr_default(G, W_1[1 - alpha], 1 - alpha);
           std::tie(W[0], W[1]) = solve_recursive(G, B, calculation_steps);
           W[1 - alpha] = set_union(W[1 - alpha], B);
         }
      }

      mCRL2log(log::debug) << "\n  --- solution for solve_recursive input ---\n" << G;
      mCRL2log(log::debug) << "   W0 = " << W[0] << std::endl;
      mCRL2log(log::debug) << "   W1 = " << W[1] << std::endl;
      assert(W[0].size() + W[1].size() + G.exclude().count() == N);
      return { W[0], W[1] };
    }

    // Handles nodes with decoration true or false.
    inline
    std::pair<vertex_set, vertex_set> solve_recursive_extended(structure_graph& G, std::size_t& calculation_steps)
    {
      mCRL2log(log::debug) << "\n  --- solve_recursive_extended input ---\n" << G << std::endl;
      calculation_steps++;

      std::size_t N = G.extent();
      vertex_set Vconj(N);
      vertex_set Vdisj(N);

      // find vertices Vconj with decoration false and Vdisj with decoration true
      for (std::size_t vi = 0; vi < N; vi++)
      {
        calculation_steps++;
        if (!G.contains(vi))
        {
          continue;
        }
        const auto& v = G.find_vertex(vi);
        if (v.decoration == structure_graph::d_false)
        {
          Vconj.insert(vi);
        }
        else if (v.decoration == structure_graph::d_true)
        {
          Vdisj.insert(vi);
        }
      }

      // extend Vconj and Vdisj
      if (!Vconj.is_empty())
      {
        Vconj = attr_default(G, Vconj, 1);
      }
      if (!Vdisj.is_empty())
      {
        Vdisj = attr_default(G, Vdisj, 0);
      }

      // default case
      if (Vconj.is_empty() && Vdisj.is_empty())
      {
        return solve_recursive(G, calculation_steps);
      }
      else
      {
        vertex_set Wconj(N);
        vertex_set Wdisj(N);
        vertex_set Vunion = set_union(Vconj, Vdisj);
        std::tie(Wdisj, Wconj) = solve_recursive(G, Vunion, calculation_steps);
        return std::make_pair(set_union(Wdisj, Vdisj), set_union(Wconj, Vconj));
      }
    }

    static void insert_edge(structure_graph::vertex_vector& V, structure_graph::index_type ui, structure_graph::index_type vi)
    {
      using utilities::detail::contains;
      structure_graph::vertex& u = V[ui];
      structure_graph::vertex& v = V[vi];
      if (!contains(u.successors, vi))
      {
        u.successors.push_back(vi);
        v.predecessors.push_back(ui);
      }
    }

    void check_solve_recursive_solution(const structure_graph& G, bool is_disjunctive, const vertex_set& Wdisj, const vertex_set& Wconj, std::size_t& calculation_steps)
    {
      using utilities::detail::contains;

      mCRL2log(log::debug) << "\n--- CHECK STRATEGY ---" << std::endl;
      log_vertex_set(G, Wconj, "Wconj");
      log_vertex_set(G, Wdisj, "Wdisj");

      calculation_steps++;
      using vertex = structure_graph::vertex;
      structure_graph::index_type init = G.initial_vertex();

      // V contains the vertices of G, but not the edges
      structure_graph::vertex_vector V = G.all_vertices();
      for (vertex& v: V)
      {
        calculation_steps++;
        v.successors.clear();
        v.predecessors.clear();
      }

      std::set<structure_graph::index_type> todo = { init };
      std::set<structure_graph::index_type> done;

      while (!todo.empty())
      {
        calculation_steps++;
        structure_graph::index_type u = *todo.begin();
        todo.erase(todo.begin());
        done.insert(u);
        if ((is_disjunctive && G.decoration(u) == structure_graph::d_disjunction) || (!is_disjunctive && G.decoration(u) == structure_graph::d_conjunction))
        {
          // explore only the strategy edge
          structure_graph::index_type v = G.strategy(u);
          insert_edge(V, u, v);
          if (v != undefined_vertex() && !contains(done, v))
          {
            todo.insert(v);
          }
        }
        else
        {
          // explore all outgoing edges
          for (structure_graph::index_type v: G.successors(u))
          {
            calculation_steps++;
            insert_edge(V, u, v);
            if (!contains(done, v))
            {
              todo.insert(v);
            }
          }
        }
      }

      vertex_set Wconj1;
      vertex_set Wdisj1;

      structure_graph Gcopy(V, G.initial_vertex(), G.exclude());
      std::tie(Wdisj1, Wconj1) = solve_recursive_extended(Gcopy, calculation_steps);
      bool is_disjunctive1;
      if (Wdisj1.contains(G.initial_vertex()))
      {
        is_disjunctive1 = true;
      }
      else if (Wconj1.contains(G.initial_vertex()))
      {
        is_disjunctive1 = false;
      }
      else
      {
        throw mcrl2::runtime_error("No solution found!!!");
      }
      if (is_disjunctive != is_disjunctive1)
      {
        log_vertex_set(Gcopy, Wconj1, "Wconj1");
        log_vertex_set(Gcopy, Wdisj1, "Wdisj1");
        throw mcrl2::runtime_error("check_solve_recursive_solution failed!");
      }
    }

  public:
    explicit solve_structure_graph_algorithm(bool check_strategy_ = false, bool use_toms_optimization_ = false)
      : check_strategy(check_strategy_),
        use_toms_optimization(use_toms_optimization_)
    {}

    /// Returns the winning player (alpha)
    inline
    bool solve(structure_graph& G, std::size_t& calculation_steps)
    {
      calculation_steps++;
      auto W = solve_partitions(G, calculation_steps);

      bool is_disjunctive;
      if (W.first.contains(G.initial_vertex()))
      {
        is_disjunctive = true;
      }
      else
      {
        is_disjunctive = false;
      }

      return is_disjunctive;
    }

    inline
    bool solve(structure_graph& G)
    {
      std::size_t calculation_steps=0;
      return solve(G, calculation_steps);
    }

    /// Returns the winning partition
    inline
    std::pair<vertex_set, vertex_set> solve_partitions(structure_graph& G, std::size_t& calculation_steps)
    {
      mCRL2log(log::verbose) << "Solving parity game..." << std::endl;
      mCRL2log(log::debug) << G << std::endl;
      assert(G.extent() > 0);
      assert(G.is_defined());
      calculation_steps++;
      auto W = solve_recursive_extended(G,calculation_steps);
      bool is_disjunctive;
      if (W.first.contains(G.initial_vertex()))
      {
        is_disjunctive = true;
      }
      else if (W.second.contains(G.initial_vertex()))
      {
        is_disjunctive = false;
      }
      else
      {
        throw mcrl2::runtime_error("No solution found!!!");
      }

      if (check_strategy)
      {
        check_solve_recursive_solution(G, is_disjunctive, W.first, W.second, calculation_steps);
      }

      mCRL2log(log::debug) << "\nSolved structure graph " << std::endl;
      mCRL2log(log::debug) << G << std::endl;
      return W;
    }

  protected:
    // Solves G, and extracts the minimal vertex set that forms the evidence for the computed
    // solution. Shared preamble of solve_with_counter_example / solve_with_evidence_lts in the
    // subclasses below.
    inline
    std::pair<bool, std::set<structure_graph::index_type>> extract_evidence(structure_graph& G)
    {
      mCRL2log(log::verbose) << "Solving parity game..." << std::endl;
      vertex_set Wconj;
      vertex_set Wdisj;
      std::size_t calculation_steps = 0;
      std::tie(Wdisj, Wconj) = solve_recursive_extended(G, calculation_steps);
      structure_graph::index_type init = G.initial_vertex();

      mCRL2log(log::verbose) << "Extracting evidence..." << std::endl;
      std::set<structure_graph::index_type> W = extract_minimal_structure_graph(G, init, Wdisj, Wconj);
      return { Wdisj.contains(init), W };
    }
};

class lps_solve_structure_graph_algorithm: public solve_structure_graph_algorithm
{
  protected:
    static lps::specification create_counter_example_lps(structure_graph& G, const std::set<structure_graph::index_type>& V, const lps::specification& lpsspec, const pbes& p, const pbes_equation_index& p_index)
    {
      try {
        lps::specification result = lpsspec;
        result.process().action_summands().clear();
        result.process().deadlock_summands().clear();
        auto& action_summands = result.process().action_summands();
        const data::variable_list& process_parameters = lpsspec.process().process_parameters();

        detail::for_each_evidence_transition(G, V, lpsspec, p, p_index,
          [&](std::size_t summand_index, const data::data_expression_vector& source, const lps::multi_action& a, const data::data_expression_vector& target)
          {
            lps::action_summand summand = lpsspec.process().action_summands().at(summand_index);

            data::data_expression_vector condition;
            data::assignment_vector next_state_assignments;
            std::size_t i = 0;
            for (const data::variable& d_i: process_parameters)
            {
              condition.push_back(data::equal_to(d_i, source.at(i)));
              next_state_assignments.emplace_back(d_i, target.at(i));
              i++;
            }

            summand.summation_variables() = data::variable_list();
            summand.condition() = data::join_and(condition.begin(), condition.end());
            summand.multi_action() = a;
            summand.assignments() = data::assignment_list(next_state_assignments.begin(), next_state_assignments.end());

            action_summands.push_back(summand);
          });

        if (!check_well_typedness(result))
        {
          throw mcrl2::runtime_error("The counter example LPS is not well typed, either wrong file provided or an internal error occurred.");
        }

        return result;
      }
      catch (const std::exception& e)
      {
        throw mcrl2::runtime_error(std::string("Counter-example cannot be reconstructed, either wrong file provided or an internal error occurred. ") + e.what());
      }
    }

  public:
    lps_solve_structure_graph_algorithm() = default;

    /// \brief Solve a pbes for some equation, while constructing a counter example or wittness based on the accompanying linear process.
    /// \param G       A structure graph.
    /// \param lpsspec The original LPS that was used to create the PBES.
    /// \param p       The pbes to be solved.
    /// \param p_index The index of the pbes equation to be solved.
    /// \return A boolean indicating the solution and a linear process that represents the counter example.
    std::pair<bool, lps::specification> solve_with_counter_example(structure_graph& G, const lps::specification& lpsspec, const pbes& p, const pbes_equation_index& p_index)
    {
      if (!lpsspec.global_variables().empty())
      {
        throw mcrl2::runtime_error("solve_with_counter_example requires an LPS without global variables.");
      }
      if (!p.global_variables().empty())
      {
        throw mcrl2::runtime_error("solve_with_counter_example requires a PBES without global variables.");
      }

      auto [is_disjunctive, W] = extract_evidence(G);
      return { is_disjunctive, create_counter_example_lps(G, W, lpsspec, p, p_index) };
    }
};

class lts_solve_structure_graph_algorithm: public solve_structure_graph_algorithm
{
  protected:
    // Removes all transitions from ltsspec, except the ones in transition_indices.
    // After that, the unreachable parts of the LTS are removed.
    static inline
    void filter_transitions(lts::lts_lts_t& ltsspec, const std::set<std::size_t>& transition_indices)
    {
      // remove transitions
      const auto& lts_transitions = ltsspec.get_transitions();
      std::vector<lts::transition> transitions;
      for (std::size_t i: transition_indices)
      {
        if (i >= lts_transitions.size())
        {
          throw mcrl2::runtime_error("Counter-example cannot be reconstructed from this LTS. Did you supply the correct file?");
        }
        transitions.push_back(lts_transitions.at(i));
      }
      ltsspec.get_transitions() = transitions;

      // remove unreachable states
      lts::reachability_check(ltsspec, true);
    }

    // modifies ltsspec
    static inline
    void create_counter_example_lts(structure_graph& G, const std::set<structure_graph::index_type>& V, lts::lts_lts_t& ltsspec)
    {
      try
      {
        std::set<std::size_t> transition_indices;
        for (structure_graph::index_type vi: V)
        {
          const auto& v = G.find_vertex(vi);
          if (is_propositional_variable_instantiation(v.formula()))
          {
            const propositional_variable_instantiation& Z = atermpp::down_cast<propositional_variable_instantiation>(v.formula());
            std::size_t transition_index;
            if (detail::counter_example_index(Z.name(), transition_index))
            {
              transition_indices.insert(transition_index);
            }
          }
        }
        filter_transitions(ltsspec, transition_indices);
      }
      catch (const std::exception& e)
      {
        throw mcrl2::runtime_error(std::string("Counter-example cannot be reconstructed, either wrong file provided or an internal error occurred. ") + e.what());
      }
    }

  public:
    lts_solve_structure_graph_algorithm() = default;

    /// \brief Solve a boolean equation system while generating a counter example.
    /// \param G       A structure graph.
    /// \param ltsspec The original LTS that was used to create the PBES.
    inline
    bool solve_with_counter_example(structure_graph& G, lts::lts_lts_t& ltsspec)
    {
      auto [is_disjunctive, W] = extract_evidence(G);
      create_counter_example_lts(G, W, ltsspec);
      return is_disjunctive;
    }
};

inline
bool solve_structure_graph(structure_graph& G, std::size_t& calculation_steps, bool check_strategy = false)
{
  bool use_toms_optimization = !check_strategy;
  solve_structure_graph_algorithm algorithm(check_strategy, use_toms_optimization);
  return algorithm.solve(G, calculation_steps);
}

inline
bool solve_structure_graph(structure_graph& G, bool check_strategy = false)
{
  std::size_t calculation_steps=0;
  return solve_structure_graph(G, calculation_steps, check_strategy);
}

/// Returns a mapping from PBES variable instantations to vertices in the structure graph for vertices won by player alpha.
inline
std::pair<bool, std::unordered_map<pbes_expression, structure_graph::index_type>> solve_structure_graph_winning_mapping(structure_graph& G, bool check_strategy = false)
{
  bool use_toms_optimization = !check_strategy;
  solve_structure_graph_algorithm algorithm(check_strategy, use_toms_optimization);
  std::size_t calculation_steps=0;
  auto W = algorithm.solve_partitions(G, calculation_steps);

  bool is_disjunctive;
  if (W.first.contains(G.initial_vertex()))
  {
    is_disjunctive = true;
  }
  else
  {
    is_disjunctive = false;
  }

  // Make a mapping from the formula to the index it belongs to.
  std::unordered_map<pbes_expression, structure_graph::index_type> mapping;
  for (structure_graph::index_type index : W.first.vertices()) {
    mapping.insert(std::make_pair(G.find_vertex(index).formula(), index));
  }

  for (structure_graph::index_type index : W.second.vertices()) {
    mapping.insert(std::make_pair(G.find_vertex(index).formula(), index));
  }

  return { is_disjunctive, mapping };
}

inline
std::pair<bool, lps::specification> solve_structure_graph_with_counter_example(structure_graph& G, const lps::specification& lpsspec, const pbes& p, const pbes_equation_index& p_index)
{
  lps_solve_structure_graph_algorithm algorithm;
  return algorithm.solve_with_counter_example(G, lpsspec, p, p_index);
}

/// \brief Solve this pbes_system using a structure graph generating a counter example.
/// \param G       The structure graph.
/// \param ltsspec The original LTS that was used to create the PBES.
inline
bool solve_structure_graph_with_counter_example(structure_graph& G, lts::lts_lts_t& ltsspec)
{
  lts_solve_structure_graph_algorithm algorithm;
  return algorithm.solve_with_counter_example(G, ltsspec);
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_SOLVE_STRUCTURE_GRAPH_H
