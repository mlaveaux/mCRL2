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
#include "mcrl2/atermpp/standard_containers/unordered_set.h"
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
    std::size_t m_count = 0;
    std::size_t m_regeneration_period = 100;

  public:
    periodic_guard()
    {}

    periodic_guard(std::size_t initial_regeneration_period)
      : m_regeneration_period(initial_regeneration_period)
    {}

    bool is_expired() 
    {
      ++m_count;
      if (m_count >= m_regeneration_period)
      {
        m_count=0;
        return true;
      }
      return false;
    }
  
    void set_expiration_steps(const std::size_t p)
    {
      m_regeneration_period = p;
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
    detail::periodic_guard on_the_fly_solve_trigger;
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

      // A reference_aterm_stack_element is a stack element to be used inside an aterm container.
      // Its elements are not protected individually, when used outside a container. In that 
      // case the stack_element below should be used. 
      struct reference_aterm_stack_element
      {
        atermpp::detail::reference_aterm<pbes_expression> b;
        atermpp::detail::reference_aterm<pbes_expression> f;
        atermpp::detail::reference_aterm<pbes_expression> g0;
        atermpp::detail::reference_aterm<pbes_expression> g1;

        reference_aterm_stack_element(
          const pbes_expression& b_,
          const pbes_expression& f_,
          const pbes_expression& g0_,
          const pbes_expression& g1_
        )
         : b(b_), f(f_), g0(g0_), g1(g1_)
        {}

        template<class E1, class E2>
        reference_aterm_stack_element(
          const E1& b_,
          const E2& f_,
          const pbes_expression& g0_,
          const pbes_expression& g1_
        )
         : b(atermpp::down_cast<pbes_expression>(b_)),
           f(atermpp::down_cast<pbes_expression>(f_)), 
           g0(g0_), 
           g1(g1_)
        {}

        void mark(atermpp::term_mark_stack& todo) const
        {
          mark_term(b, todo);
          mark_term(f, todo);
          mark_term(g0, todo);
          mark_term(g1, todo);
        }
      };

      // A stack element is an reference_aterm_stack_element of which the elements are protected against garbage collection. 
      struct stack_element
      {
        pbes_expression b;
        pbes_expression f;
        pbes_expression g0;
        pbes_expression g1;

        stack_element(const reference_aterm_stack_element& s)
         : b(s.b), f(s.f), g0(s.g0), g1(s.g1)
        {}
      };

      std::array<vertex_set, 2>& S;

      detail::structure_graph_builder& graph_builder;
      atermpp::vector<reference_aterm_stack_element> stack;

      Rplus_traverser(std::array<vertex_set, 2>& S_, detail::structure_graph_builder& graph_builder_)
       : S(S_), graph_builder(graph_builder_)
      {}

      void push(const reference_aterm_stack_element& elem)
      {
        stack.push_back(elem);
      }

      stack_element pop()
      {
        stack_element result(stack.back());
        stack.pop_back();
        return result;
      }

      // Return the top element of result_stack
      reference_aterm_stack_element& top()
      {
        return stack.back();
      }

      // Return the top element of result_stack
      const reference_aterm_stack_element& top() const
      {
        return stack.back();
      }

      static bool less(const pbes_expression& x1, const pbes_expression& x2)
      {
        return x1 < x2;
      }

      void leave(const data::data_expression& x)
      {
        if (is_true(x) || is_false(x))
        {
          stack.emplace_back(x, x, true_(), false_());
        }
        else
        {
         throw mcrl2::runtime_error("Fail to evaluate the expression " + data::pp(x) + " as it should be equal to true or false.");
        }
      }

      void leave(const propositional_variable_instantiation& x)
      {
        const structure_graph::index_type u = graph_builder.find_vertex(x);
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
        const stack_element elem2 = pop();
        reference_aterm_stack_element& elem1 = top();
        atermpp::detail::reference_aterm<pbes_expression>& b_1 = elem1.b;
        atermpp::detail::reference_aterm<pbes_expression>& f1_prime = elem1.f;
        atermpp::detail::reference_aterm<pbes_expression>& g0_1 = elem1.g0;
        atermpp::detail::reference_aterm<pbes_expression>& g1_1 = elem1.g1;
        const pbes_expression& b_2 = elem2.b;
        const pbes_expression& f2_prime = elem2.f;
        const pbes_expression& g0_2 = elem2.g0;
        const pbes_expression& g1_2 = elem2.g1;

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
          g0_1 = true_();
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
            g0_1 = true_();
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
        const stack_element elem2 = pop();
        reference_aterm_stack_element& elem1 = top();
        atermpp::detail::reference_aterm<pbes_expression>& b_1 = elem1.b;
        atermpp::detail::reference_aterm<pbes_expression>& f1_prime = elem1.f;
        atermpp::detail::reference_aterm<pbes_expression>& g0_1 = elem1.g0;
        atermpp::detail::reference_aterm<pbes_expression>& g1_1 = elem1.g1;
        const pbes_expression& b_2 = elem2.b;
        const pbes_expression& f2_prime = elem2.f;
        const pbes_expression& g0_2 = elem2.g0;
        const pbes_expression& g1_2 = elem2.g1;

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
      return Rplus_traverser::stack_element(f.top());  // Protection is added explicitly. 
    }

    bool solution_found(const propositional_variable_instantiation& init) const override
    {
      const structure_graph::index_type u = m_graph_builder.find_vertex(init);
      return S[0].contains(u) || S[1].contains(u);
    }

    // Returns true if all nodes in the todo list are undefined (i.e. have not been processed yet)
    bool todo_has_only_undefined_nodes() const
    {
      for (const propositional_variable_instantiation& X: todo.elements())
      {
        const structure_graph::index_type u = m_graph_builder.find_vertex(X);
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
      std::size_t& calculation_steps)
    {
      using utilities::detail::contains;
      global_current_prune_round++;

      std::size_t old_todo_size = todo.elements().size();

      simple_structure_graph G(m_graph_builder.vertices());

      atermpp::deque<pbes_expression> todo1{init};
      // The index of the indexed_set below does not appear to be used. 
      // Maybe the use of an atermpp::unordered_set is a better choice.
      atermpp::indexed_set<pbes_expression> done1;
      std::ignore = done1.insert(init);
      
      atermpp::unordered_set<propositional_variable_instantiation> new_todo;
      pbes_expression X;
      while (!todo1.empty())
      {
        using utilities::detail::contains;

        X = todo1.front();
        todo1.pop_front();
        const structure_graph::index_type u = m_graph_builder.find_vertex(X);
        const structure_graph::vertex& u_ = m_graph_builder.vertex(u);
        calculation_steps++;
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
            for (const structure_graph::index_type& v: G.successors(u))
            {
              calculation_steps++;
              const structure_graph::vertex& v_ = m_graph_builder.vertex(v);
              const pbes_expression& Y = v_.formula();
              if (!contains(done1, Y))
              {
                todo1.emplace_back(Y);
                std::ignore = done1.insert(Y);
              }
            }
          }
        }
      }

      // new_todo_list := new_todo \cap todo 
      // N.B. An attempt is made to preserve the order of the current todo list, to not
      // disturb breadth first and depth first search.
      atermpp::deque<propositional_variable_instantiation> new_todo_list;
      calculation_steps=calculation_steps+todo.elements().size();
      for (const propositional_variable_instantiation& X: todo.elements())
      {
        if (new_todo.contains(X))
        {
          new_todo_list.push_back(X);
          new_todo.erase(X);
        }
      }
      calculation_steps=calculation_steps+new_todo.size();
      for(const propositional_variable_instantiation& X: new_todo)
      {
        if (m_options.exploration_strategy == breadth_first)
        {
          new_todo_list.push_back(X);
        }
        else
        {
          new_todo_list.push_front(X);
        }
      }
      todo.set_todo(new_todo_list);
      assert(todo_has_only_undefined_nodes());
      if (todo.elements().size() == old_todo_size)
      { 
        mCRL2log(log::verbose) << "Pruning of the  todo list had no effect on its size. ";
      }
      else if (todo.elements().size() > old_todo_size)
      { 
        mCRL2log(log::verbose) << "Pruned the todo list. Added " << todo.elements().size() - old_todo_size << " elements. ";
      }
      else
      { 
        mCRL2log(log::verbose) << "Pruned the todo list. Removed " << old_todo_size - todo.elements().size() << " elements. ";
      }
      mCRL2log(log::verbose) << "The todo list has size " << todo.elements().size() << ".\n";
   };

   // Execute a prune_todo_list if m_options.prune_todo_list is set. 
   void prune_todo_list_conditional(
      const propositional_variable_instantiation& init,
      pbesinst_lazy_todo& todo,
      std::size_t& calculation_steps)
    {
      if (m_options.prune_todo_list)
      {
        prune_todo_list(init, todo, calculation_steps);
      }
    }

    // Trigger the prune_todo_list at fixed intervals, provided m_options.prune_todo_list is set, 
    // keeping the workload limited related to other calculations.
    // Trigger always when the todo set is empty. 
    void prune_todo_list_time_triggered(
      const propositional_variable_instantiation& init,
      pbesinst_lazy_todo& todo)
    {
      if (m_options.prune_todo_list&&
          (reset_guard.is_expired() || todo.elements().empty() || m_options.aggressive))
      {
        std::size_t calculation_steps=0;
        prune_todo_list(init, todo, calculation_steps);
        reset_guard.set_expiration_steps(m_options.prune_and_solve_frequently?calculation_steps:calculation_steps*100);
      }
    }


    bool strategies_are_set_in_solved_nodes() const
    {
      simple_structure_graph G(m_graph_builder.vertices());
      for (structure_graph::index_type u: S[0].vertices())
      {
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
      std::optional<data::rewriter> rewriter = std::nullopt
    )
      : pbesinst_structure_graph_algorithm(options, p, G, rewriter),
        b(options.number_of_threads+1), on_the_fly_solve_trigger(2)
    {}

    // Optimization 2 is implemented by overriding the function rewrite_psi.
    void rewrite_psi(const std::size_t thread_index,
                     pbes_expression& result,
                     const fixpoint_symbol& symbol,
                     const propositional_variable_instantiation& X,
                     const pbes_expression& psi
                    ) override
    {
      assert(&result != &psi); // required by super::rewrite_psi
      super::rewrite_psi(thread_index, result, symbol, X, psi);
      const Rplus_traverser::stack_element& rplus_result = Rplus(result);
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

      const structure_graph::index_type u = m_graph_builder.find_vertex(X);
      if (is_true(b[thread_index]))
      {
        S[0].insert(u);
      }
      else if (is_false(b[thread_index]))
      {
        S[1].insert(u);
      }
    }

    void report_found_solutions(stopwatch& timer)
    {
      mCRL2log(log::verbose) << "Found solution for" << std::setw(12) << S[0].size() + S[1].size() << " BES equations." << std::endl;
      mCRL2log(log::verbose) << "Finished partial solving (time = " << std::setprecision(2) << std::fixed << timer.seconds() << "s).\n";
    }

    void on_discovered_elements(const std::set<propositional_variable_instantiation>& elements) override
    {
      using utilities::detail::contains;
      stopwatch timer;

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
      else if (m_options.optimization == partial_solve_strategy::detect_winning_loops_using_fatal_attractor && 
               (m_options.aggressive || on_the_fly_solve_trigger.is_expired()))
      {
        mCRL2log(log::verbose) << "Start partial solving.\n"; 

        std::size_t calculation_steps=0;  // Count how many calculation steps it takes to find loops, and retry this after on_discovered_elements have been called that many times. 
        simple_structure_graph G(m_graph_builder.vertices());
        detail::find_loops2(G, S, tau, calculation_steps, m_iteration_count); // modifies S[0] and S[1]
        assert(strategies_are_set_in_solved_nodes());
        report_found_solutions(timer);
        prune_todo_list_conditional(init, todo, calculation_steps);
        on_the_fly_solve_trigger.set_expiration_steps(m_options.prune_and_solve_frequently?calculation_steps/1000:calculation_steps/10);
      }
      else if ((partial_solve_strategy::solve_subgames_using_fatal_attractor_local <= m_options.optimization && 
                m_options.optimization <= partial_solve_strategy::solve_subgames_using_solver) && 
                (m_options.aggressive || on_the_fly_solve_trigger.is_expired()))
      {
        mCRL2log(log::verbose) << "Start partial solving.\n"; 

        std::size_t calculation_steps=0;  // Count how many calculation steps it takes to find loops, and retry this after on_discovered_elements have been called that many times. 

        simple_structure_graph G(m_graph_builder.vertices());
        if (m_options.optimization == partial_solve_strategy::solve_subgames_using_fatal_attractor_local)
        {
          detail::fatal_attractors(G, S, tau, calculation_steps, m_iteration_count); // modifies S[0] and S[1]
          assert(strategies_are_set_in_solved_nodes());
        }
        else if (m_options.optimization == partial_solve_strategy::solve_subgames_using_fatal_attractor_original)
        {
          detail::fatal_attractors_original(G, S, tau, calculation_steps, m_iteration_count); // modifies S[0] and S[1]
          assert(strategies_are_set_in_solved_nodes());
        }
        else if (m_options.optimization == partial_solve_strategy::solve_subgames_using_solver)
        {
          m_graph_builder.finalize();
          detail::partial_solve(m_graph_builder.m_graph, todo, S, tau, calculation_steps, m_iteration_count, m_graph_builder); // modifies S[0] and S[1]
          assert(strategies_are_set_in_solved_nodes());
        }
        report_found_solutions(timer);
        prune_todo_list_conditional(init, todo, calculation_steps);
        on_the_fly_solve_trigger.set_expiration_steps(m_options.prune_and_solve_frequently?calculation_steps/1000:calculation_steps/10);
      }
      else if (m_options.optimization == partial_solve_strategy::detect_winning_loops_original && 
               (m_options.aggressive || on_the_fly_solve_trigger.is_expired()))
      {
        mCRL2log(log::verbose) << "Start partial solving.\n"; 

        std::size_t calculation_steps=0;  // Count how many calculation steps it takes to find loops, and retry this after on_discovered_elements have been called that many times. 

        simple_structure_graph G(m_graph_builder.vertices());
        detail::find_loops(G, discovered, todo, S, tau, m_iteration_count, m_graph_builder); // modifies S[0] and S[1]
        assert(strategies_are_set_in_solved_nodes());
        prune_todo_list_conditional(init, todo, calculation_steps);
        on_the_fly_solve_trigger.set_expiration_steps(m_options.prune_and_solve_frequently?calculation_steps/1000:calculation_steps/10);
      }

      prune_todo_list_time_triggered(init, todo);
    }

    void on_end_while_loop() override
    {
      using  utilities::detail::contains;

      simple_structure_graph G(m_graph_builder.vertices());

      const structure_graph::index_type u = m_graph_builder.find_vertex(init);
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
