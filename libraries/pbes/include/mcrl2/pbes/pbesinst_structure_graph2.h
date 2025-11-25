// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesinst_structure_graph2.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PBESINST_STRUCTURE_GRAPH2_H
#define MCRL2_PBES_PBESINST_STRUCTURE_GRAPH2_H

#include "mcrl2/atermpp/standard_containers/deque.h"
#include "mcrl2/atermpp/standard_containers/indexed_set.h"
#include "mcrl2/atermpp/standard_containers/vector.h"
#include "mcrl2/pbes/pbesinst_fatal_attractors.h"
#include "mcrl2/pbes/pbesinst_find_loops.h"
#include "mcrl2/pbes/pbesinst_partial_solve.h"
#include "mcrl2/pbes/pbesinst_structure_graph.h"
#include "mcrl2/pbes/pbessolve_options.h"
#include "mcrl2/utilities/stopwatch.h"



namespace mcrl2::pbes_system {

namespace detail {

class computation_guard
{
  protected:
     std::size_t m_count = 64;

  public:
    explicit computation_guard(std::size_t initial_count = 64)
      : m_count(initial_count)
    {}

    bool operator()(std::size_t count)
    {
      bool result = count >= m_count;
      while (m_count <= count)
      {
        m_count *= 2;
      }
      return result;
    }
};

class periodic_guard
{
  protected:
    std::size_t count = 0;
    std::size_t regeneration_period = 100;

  public:
    bool operator()(std::size_t period)
    {
      if (++count == regeneration_period)
      {
        count = 0;
        regeneration_period = period;
        return true;
      }
      return false;
    }
};

} // namespace detail

/// \brief Adds an optimization to pbesinst_structure_graph.
class pbesinst_structure_graph_algorithm2: public pbesinst_structure_graph_algorithm
{
  protected:
    std::array<vertex_set, 2> S;
    std::array<strategy_vector, 2> tau;
    std::array<detail::computation_guard, 2> S_guard;

    atermpp::vector<pbes_expression> b; // to store the result of the Rplus computation
    detail::computation_guard find_loops_guard;
    detail::computation_guard fatal_attractors_guard;
    detail::periodic_guard reset_guard;

    template<typename T>
    pbes_expression expr(const T& x) const
    {
      return static_cast<const pbes_expression&>(x);
    }

    struct Rplus_traverser: public pbes_expression_traverser<Rplus_traverser>
    {
      using super = pbes_expression_traverser<Rplus_traverser>;
      using super::enter;
      using super::leave;
      using super::apply;

      struct stack_element
      {
        atermpp::detail::reference_aterm<pbes_expression> b;
        atermpp::detail::reference_aterm<pbes_expression> f;
        atermpp::detail::reference_aterm<pbes_expression> g0;
        atermpp::detail::reference_aterm<pbes_expression> g1;

        stack_element(
          pbes_expression b_,
          pbes_expression f_,
          pbes_expression g0_,
          pbes_expression g1_
        )
         : b(std::move(b_)), f(std::move(f_)), g0(std::move(g0_)), g1(std::move(g1_))
        {}

        void mark(atermpp::term_mark_stack& todo) const
        {
          mark_term(b, todo);
          mark_term(f, todo);
          mark_term(g0, todo);
          mark_term(g1, todo);
        }
      };

      std::array<vertex_set, 2>& S;

      detail::structure_graph_builder& graph_builder;
      // TODO: replace with aterm container
      atermpp::vector<stack_element> stack;

      Rplus_traverser(std::array<vertex_set, 2>& S_, detail::structure_graph_builder& graph_builder_)
       : S(S_), graph_builder(graph_builder_)
      {}

      void push(const stack_element& elem)
      {
        stack.push_back(elem);
      }

      stack_element pop()
      {
        auto result = stack.back();
        stack.pop_back();
        return result;
      }

      // Return the top element of result_stack
      stack_element& top()
      {
        return stack.back();
      }

      // Return the top element of result_stack
      const stack_element& top() const
      {
        return stack.back();
      }

      // TODO: use a heuristic for the smallest term
      static bool less(const pbes_expression& /* x1 */, const pbes_expression& /* x2 */)
      {
        return true;
      }

      void leave(const data::data_expression& x)
      {
        if (is_true(x) || is_false(x))
        {
          stack.emplace_back(x, x, true_(), false_());
        }
        else
        {
         throw mcrl2::runtime_error("Fail to evaluate the expression " + pbes_system::pp(x) + " as it should be equal to true or false.");
        }
      }

      void leave(const propositional_variable_instantiation& x)
      {
        auto u = graph_builder.find_vertex(x);
        if (u == undefined_vertex())
        {
          // if x is not yet in the graph, then it certainly isn't in S[0] or S[1]
          stack.emplace_back(data::undefined_data_expression(), x, true_(), false_());
        }
        else if (S[0].contains(u))
        {
          stack.emplace_back(true_(), x, x, false_());
        }
        else if (S[1].contains(u))
        {
          stack.emplace_back(false_(), x, true_(), x);
        }
        else
        {
          stack.emplace_back(data::undefined_data_expression(), x, true_(), false_());
        }
      }

      void leave(const and_& /* x */)
      {
        stack_element elem2 = pop();
        stack_element& elem1 = top();
        auto& b_1 = elem1.b;
        auto& f1_prime = elem1.f;
        auto& g0_1 = elem1.g0;
        auto& g1_1 = elem1.g1;
        auto& b_2 = elem2.b;
        auto& f2_prime = elem2.f;
        auto& g0_2 = elem2.g0;
        auto& g1_2 = elem2.g1;

        // put the result in (b1, f1_prime, g0, g1)
        if (is_true(b_1) && is_true(b_2))
        {
          b_1 = true_();
          f1_prime = and_(f1_prime, f2_prime);
          g0_1 = and_(g0_1, g0_2);
          g1_1 = false_();
        }
        else if (is_false(b_1) && !is_false(b_2))
        {
          b_1 = false_();
          // f1_prime = f1_prime;
          g0_1 = true_();
          // g1_1 = g1_1;
        }
        else if (!is_false(b_1) && is_false(b_2))
        {
          b_1 = false_();
          f1_prime = f2_prime;
          g0_1 = true_();
          g1_1 = g1_2;
        }
        else if (is_false(b_1) && is_false(b_2))
        {
          if (less(f1_prime, f2_prime))
          {
            b_1 = false_();
            // f1_prime = f1_prime;
            g0_1 = true_();
            // g1_1 = g1_1;
          }
          else
          {
            b_1 = false_();
            f1_prime = f2_prime;
            g0_1 = true_();
            g1_1 = g1_2;
          }
        }
        else // if (b1 == data::undefined_data_expression() && b2 == data::undefined_data_expression())
        {
          b_1 = data::undefined_data_expression();
          f1_prime = and_(f1_prime, f2_prime);
          g0_1 = true_();
          g1_1 = false_();
        }
      }

      void leave(const or_& /* x */)
      {
        stack_element elem2 = pop();
        stack_element& elem1 = top();
        auto& b_1 = elem1.b;
        auto& f1_prime = elem1.f;
        auto& g0_1 = elem1.g0;
        auto& g1_1 = elem1.g1;
        auto& b_2 = elem2.b;
        auto& f2_prime = elem2.f;
        auto& g0_2 = elem2.g0;
        auto& g1_2 = elem2.g1;

        // put the result in (b1, f1_prime, g0, g1)
        if (is_false(b_1) && is_false(b_2))
        {
          b_1 = false_();
          f1_prime = or_(f1_prime, f2_prime);
          g0_1 = true_();
          g1_1 = or_(g1_1, g1_2);
        }
        else if (is_true(b_1) && !is_true(b_2))
        {
          b_1 = true_();
          // f1_prime = f1_prime;
          // g0_1 = g0_1;
          g1_1 = false_();
        }
        else if (!is_true(b_1) && is_true(b_2))
        {
          b_1 = true_();
          f1_prime = f2_prime;
          g0_1 = g0_2;
          g1_1 = false_();
        }
        else if (is_true(b_1) && is_true(b_2))
        {
          if (less(f1_prime, f2_prime))
          {
            b_1 = true_();
            // f1_prime = f1_prime;
            // g0_1 = g0_1;
            g1_1 = false_();
          }
          else
          {
            b_1 = true_();
            f1_prime = f2_prime;
            g0_1 = g0_2;
            g1_1 = false_();
          }
        }
        else // if (b1 == data::undefined_data_expression() && b2 == data::undefined_data_expression())
        {
          b_1 = data::undefined_data_expression();
          f1_prime = or_(f1_prime, f2_prime);
          g0_1 = true_();
          g1_1 = false_();
        }
      }

      void enter(const imp& x)
      {
        throw mcrl2::runtime_error("Fail to evaluate the expression " + pbes_system::pp(x) + " as the routine Rplus does not expect an implication, which is an internal error.");
      }

      void leave(const imp& x)
      {
        enter(x);  // Print an error message, although this should never be reached..
      }

      void enter(const exists& x)
      {
        throw mcrl2::runtime_error("Fail to evaluate the expression " + pbes_system::pp(x) + " as enumeration of this exists is not possible.");
      }

      void leave(const exists& x)
      {
        enter(x); // Print an error message. 
      }

      void enter(const forall& x)
      {
        throw mcrl2::runtime_error("Fail to evaluate the expression " + pbes_system::pp(x) + " as enumeration of this forall is not possible.");
      }

      void leave(const forall& x)
      {
        enter(x); // Print an error. 
      }
    };

    Rplus_traverser::stack_element Rplus(const pbes_expression& x)
    {
      Rplus_traverser f(S, m_graph_builder);
      f.apply(x);
      return f.top();
    }

    bool solution_found(const propositional_variable_instantiation& init) const override
    {
      auto u = m_graph_builder.find_vertex(init);
      return S[0].contains(u) || S[1].contains(u);
    }

    // Returns true if all nodes in the todo list are undefined (i.e. have not been processed yet)
    bool todo_has_only_undefined_nodes() const
    {
      /* for (const propositional_variable_instantiation& X: todo.all_elements())  all_elements does not seem to work. Therefore split below. 
      {
        structure_graph::index_type u = m_graph_builder.find_vertex(X);
        const structure_graph::vertex& u_ = m_graph_builder.vertex(u);
        if (u_.is_defined())
        {
          return false;
        }
      } */
      for (const propositional_variable_instantiation& X: todo.elements())
      {
        structure_graph::index_type u = m_graph_builder.find_vertex(X);
        const structure_graph::vertex& u_ = m_graph_builder.vertex(u);
        if (u_.is_defined())
        {
          return false;
        }
      }
      for (const propositional_variable_instantiation& X: todo.irrelevant_elements())
      {
        structure_graph::index_type u = m_graph_builder.find_vertex(X);
        const structure_graph::vertex& u_ = m_graph_builder.vertex(u);
        if (u_.is_defined())
        {
          return false;
        }
      }
      return true;
    }

    void prune_todo_list(
      const propositional_variable_instantiation& init,
      pbesinst_lazy_todo& todo,
      std::size_t regeneration_period
    )
    {
      using utilities::detail::contains;

      if (!reset_guard(regeneration_period) && !m_options.aggressive && !todo.elements().empty())
      {
        return;
      }

      simple_structure_graph G(m_graph_builder.vertices());
      atermpp::deque<pbes_expression> todo1{init};
      atermpp::indexed_set<pbes_expression> done1;
      atermpp::indexed_set<propositional_variable_instantiation> new_todo;

      pbes_expression X;
      while (!todo1.empty())
      {
        using utilities::detail::contains;

        X = todo1.front();
        todo1.pop_front();
        done1.insert(X);
        auto u = m_graph_builder.find_vertex(X);
        const auto& u_ = m_graph_builder.vertex(u);

        if (u_.decoration == structure_graph::d_none && u_.successors.empty())
        {
          assert(is_propositional_variable_instantiation(u_.formula()));
          new_todo.insert(atermpp::down_cast<propositional_variable_instantiation>(u_.formula()));
        }
        else
        {
          if (!S[0].contains(u) && !S[1].contains(u))
          {
            // todo' := todo' U (succ(u) \ done')
            for (auto v: G.successors(u))
            {
              const auto& v_ = m_graph_builder.vertex(v);
              const auto& Y = v_.formula();
              if (contains(done1, Y))
              {
                continue;
              }
              todo1.push_back(Y);
            }
          }
        }
      }

      // new_todo_list := new_todo \cap (todo U irrelevant)
      // N.B. An attempt is made to preserve the order of the current todo list, to not
      // disturb breadth first and depth first search.
      atermpp::deque<propositional_variable_instantiation> new_todo_list;
      for (const propositional_variable_instantiation& X: todo.irrelevant_elements())
      {
        if (contains(new_todo, X))
        {
          new_todo_list.push_back(X);
        }
      }
      for (const propositional_variable_instantiation& X: todo.elements())
      {
        if (contains(new_todo, X))
        {
          new_todo_list.push_back(X);
        }
      }
      todo.set_todo(new_todo_list);
      assert(todo_has_only_undefined_nodes());
    };

    bool strategies_are_set_in_solved_nodes() const
    {
      simple_structure_graph G(m_graph_builder.vertices());
      for (structure_graph::index_type u: S[0].vertices())
      {
        // if (G.decoration(u) == structure_graph::d_disjunction && G.strategy(u) == undefined_vertex())
        if (G.decoration(u) == structure_graph::d_disjunction && tau[0][u] == undefined_vertex())
        {
          mCRL2log(log::debug) << "Error: no strategy has been set for disjunctive node " << u << " in S0." << std::endl;
          mCRL2log(log::debug) << G << std::endl;
          mCRL2log(log::debug) << "S0 = " << S[0] << std::endl;
          mCRL2log(log::debug) << "S1 = " << S[1] << std::endl;
          return false;
        }
      }
      for (structure_graph::index_type u: S[1].vertices())
      {
        // if (G.decoration(u) == structure_graph::d_conjunction && G.strategy(u) == undefined_vertex())
        if (G.decoration(u) == structure_graph::d_conjunction && tau[1][u] == undefined_vertex())
        {
          mCRL2log(log::debug) << "Error: no strategy has been set for conjunctive node " << u << " in S1." << std::endl;
          mCRL2log(log::debug) << G << std::endl;
          mCRL2log(log::debug) << "S0 = " << S[0] << std::endl;
          mCRL2log(log::debug) << "S1 = " << S[1] << std::endl;
          return false;
        }
      }
      return true;
    }

  public:
    using super = pbesinst_structure_graph_algorithm;

    pbesinst_structure_graph_algorithm2(
      const pbessolve_options& options,
      const pbes& p,
      structure_graph& G,
      pbes_quotient& quotient,
      std::optional<data::rewriter> rewriter = std::nullopt
    )
      : pbesinst_structure_graph_algorithm(options, p, G, quotient, rewriter),
        b(options.number_of_threads+1), find_loops_guard(2), fatal_attractors_guard(2)
    {}

    // Optimization 2 is implemented by overriding the function rewrite_psi.
    void rewrite_psi(const std::size_t thread_index,
                     pbes_expression& result,
                     const fixpoint_symbol& symbol,
                     const propositional_variable_instantiation& X,
                     const pbes_expression& psi
                    ) override
    {
      super::rewrite_psi(thread_index, result, symbol, X, psi);
      auto rplus_result = Rplus(result);
      b[thread_index] = rplus_result.b;
      if (is_true(rplus_result.b))
      {
        result = rplus_result.g0;
        return;
      }
      else if (is_false(rplus_result.b))
      {
        result = rplus_result.g1;
        return;
      }
      result = rplus_result.f;
    }

    void on_report_equation(const std::size_t thread_index,
                            const propositional_variable_instantiation& X,
                            const pbes_expression& psi, std::size_t k
                           ) override
    {
      super::on_report_equation(thread_index, X, psi, k);

      // The structure graph has just been extended, so S[0] and S[1] need to be resized.
      S[0].resize(m_graph_builder.extent());
      S[1].resize(m_graph_builder.extent());

      auto u = m_graph_builder.find_vertex(X);
      if (is_true(b[thread_index]))
      {
        S[0].insert(u);
      }
      else if (is_false(b[thread_index]))
      {
        S[1].insert(u);
      }
    }


    void on_discovered_elements(const std::set<propositional_variable_instantiation>& elements) override
    {
      using utilities::detail::contains;
      stopwatch timer;

      bool report = false;
      if (m_options.optimization == partial_solve_strategy::propagate_solved_equations_using_attractor)
      {
        if (S_guard[0](S[0].size()))
        {
          simple_structure_graph G(m_graph_builder.vertices());
          S[0] = attr_default_with_tau(G, S[0], 0, tau);
        }
        if (S_guard[1](S[1].size()))
        {
          simple_structure_graph G(m_graph_builder.vertices());
          S[1] = attr_default_with_tau(G, S[1], 1, tau);
        }
        assert(strategies_are_set_in_solved_nodes());
      }
      else if (m_options.optimization == partial_solve_strategy::detect_winning_loops_using_fatal_attractor && (m_options.aggressive || find_loops_guard(m_iteration_count)))
      {
        mCRL2log(log::verbose) << "start partial solving\n"; report = true;

        simple_structure_graph G(m_graph_builder.vertices());
        detail::find_loops2(G, S, tau, m_iteration_count); // modifies S[0] and S[1]
        assert(strategies_are_set_in_solved_nodes());

      }
      else if ((partial_solve_strategy::solve_subgames_using_fatal_attractor_local <= m_options.optimization && m_options.optimization <= partial_solve_strategy::solve_subgames_using_solver) && (m_options.aggressive || fatal_attractors_guard(m_iteration_count)))
      {
        mCRL2log(log::verbose) << "start partial solving\n"; report = true;

        simple_structure_graph G(m_graph_builder.vertices());
        if (m_options.optimization == partial_solve_strategy::solve_subgames_using_fatal_attractor_local)
        {
          detail::fatal_attractors(G, S, tau, m_iteration_count); // modifies S[0] and S[1]
          assert(strategies_are_set_in_solved_nodes());
        }
        else if (m_options.optimization == partial_solve_strategy::solve_subgames_using_fatal_attractor_original)
        {
          detail::fatal_attractors_original(G, S, tau, m_iteration_count); // modifies S[0] and S[1]
          assert(strategies_are_set_in_solved_nodes());
        }
        else if (m_options.optimization == partial_solve_strategy::solve_subgames_using_solver)
        {
          m_graph_builder.finalize();
          detail::partial_solve(m_graph_builder.m_graph, todo, S, tau, m_iteration_count, m_graph_builder); // modifies S[0] and S[1]
          assert(strategies_are_set_in_solved_nodes());
        }
      }
      else if (m_options.optimization == partial_solve_strategy::detect_winning_loops_original && (m_options.aggressive || find_loops_guard(m_iteration_count)))
      {        
        mCRL2log(log::verbose) << "start partial solving\n"; report = true;

        simple_structure_graph G(m_graph_builder.vertices());
        detail::find_loops(G, discovered, todo, S, tau, m_iteration_count, m_graph_builder); // modifies S[0] and S[1]
        assert(strategies_are_set_in_solved_nodes());
      }

      if (report)
      {
        mCRL2log(log::verbose) << "found solution for" << std::setw(12) << S[0].size() + S[1].size() << " BES equations" << std::endl;
        mCRL2log(log::verbose) << "finished partial solving (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s)\n";
      }

      if (m_options.prune_todo_list)
      {
        for (const propositional_variable_instantiation& e: elements)
        {
          todo.irrelevant_elements().erase(e);
        }
        prune_todo_list(init, todo, (discovered.size() - todo.size()) / 2);
      }
    }

    void on_end_while_loop() override
    {
      using  utilities::detail::contains;

      simple_structure_graph G(m_graph_builder.vertices());

      structure_graph::index_type u = m_graph_builder.find_vertex(init);
      assert(strategies_are_set_in_solved_nodes());

      std::set<structure_graph::index_type> V = extract_minimal_structure_graph(G, u, S[0], S[1], tau[0], tau[1]);

      std::size_t n = m_graph_builder.extent();
      vertex_set to_be_removed(n);
      for (std::size_t v = 0; v < n; v++)
      {
        if (!contains(V, v))
        {
          to_be_removed.insert(v);
        }
      }
      m_graph_builder.erase_vertices(to_be_removed);

      mCRL2log(log::debug) << "\nFinal structure graph " << std::endl;
      mCRL2log(log::debug) << G << std::endl;
    }
};

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_PBESINST_STRUCTURE_GRAPH2_H
