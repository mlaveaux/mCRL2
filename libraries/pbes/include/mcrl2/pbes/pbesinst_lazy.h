// Author(s): Wieger Wesselink 2017-2019
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesinst_lazyh
/// \brief A lazy algorithm for instantiating a PBES, ported from bes_deprecated.h.

#ifndef MCRL2_PBES_PBESINST_LAZY_H
#define MCRL2_PBES_PBESINST_LAZY_H

#include <optional>
#include <string>
#include <thread>
#include <mutex>
#include <functional>
#include <regex>

#include "mcrl2/atermpp/standard_containers/deque.h"
#include "mcrl2/atermpp/standard_containers/indexed_set.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/pbes/detail/bes_equation_limit.h"
#include "mcrl2/pbes/detail/instantiate_global_variables.h"
#include "mcrl2/pbes/pbes_equation_index.h"
#include "mcrl2/pbes/pbes_expression.h"
#include "mcrl2/pbes/pbessolve_options.h"
#include "mcrl2/pbes/remove_equations.h"
#include "mcrl2/pbes/replace_constants_by_variables.h"
#include "mcrl2/pbes/rewriters/enumerate_quantifiers_rewriter.h"
#include "mcrl2/pbes/rewriters/one_point_rule_rewriter.h"
#include "mcrl2/pbes/rewriters/simplify_quantifiers_rewriter.h"
#include "mcrl2/pbes/structure_graph.h"
#include "mcrl2/pbes/transformation_strategy.h"
#include "mcrl2/pbes/transformations.h"
#include "mcrl2/utilities/detail/container_utility.h"


namespace mcrl2::pbes_system
{

// This todo set maintains elements that were removed by the reset procedure.
class pbesinst_lazy_todo
{
  protected:
    std::unordered_set<propositional_variable_instantiation> irrelevant;
    atermpp::deque<propositional_variable_instantiation> todo;

    // checks some invariants on the internal state
    bool check_invariants() const
    {
      using utilities::detail::contains;
      for (const auto& X: irrelevant)
      {
        if (contains(todo, X))
        {
          return false;
        }
      }
      std::unordered_set<propositional_variable_instantiation> tmp(todo.begin(), todo.end());
      return tmp.size() == todo.size();
    }

  public:
    const propositional_variable_instantiation& front() const
    {
      return todo.front();
    }

    const propositional_variable_instantiation& back() const
    {
      return todo.back();
    }

    bool empty() const
    {
      return todo.empty() && irrelevant.empty();
    }

    std::size_t size() const
    {
      return todo.size();
    }

    const atermpp::deque<propositional_variable_instantiation>& elements() const
    {
      return todo;
    }

    const std::unordered_set<propositional_variable_instantiation>& irrelevant_elements() const
    {
      return irrelevant;
    }

    std::unordered_set<propositional_variable_instantiation>& irrelevant_elements()
    {
      return irrelevant;
    }

    void pop_front()
    {
      todo.pop_front();
    }

    void pop_back()
    {
      todo.pop_back();
    }

    void insert(const propositional_variable_instantiation& x)
    {
      irrelevant.erase(x);
      todo.push_back(x);
    }

    template <typename FwdIter, bool ThreadSafe>
    void insert(FwdIter first, 
                FwdIter last, 
                const atermpp::indexed_set<propositional_variable_instantiation, ThreadSafe>& discovered,
                const std::size_t thread_index)
    {
      using utilities::detail::contains;

      for (FwdIter i = first; i != last; ++i)
      {
        auto j = irrelevant.find(*i);
        if (j != irrelevant.end())
        {
          todo.push_back(*j);
          irrelevant.erase(j);
        }
        else if (!contains(discovered, *i, thread_index))
        {
          todo.push_back(*i);
        }
      }
    }

    void set_todo(atermpp::deque<propositional_variable_instantiation>& new_todo)
    {
      using utilities::detail::contains;
      std::size_t size_before = todo.size() + irrelevant.size();
      std::unordered_set<propositional_variable_instantiation> new_irrelevant;
      for (const propositional_variable_instantiation& x: todo)
      {
        if (!contains(new_todo, x))
        {
          new_irrelevant.insert(x);
        }
      } 
      for (const propositional_variable_instantiation& x: irrelevant)
      {
        if (!contains(new_todo, x))
        {
          new_irrelevant.insert(x);
        }
      } 
      std::swap(todo, new_todo);
      std::swap(irrelevant, new_irrelevant);

      std::size_t size_after = todo.size() + irrelevant.size();
      if (size_before != size_after)
      {
        throw mcrl2::runtime_error("sizes do not match in pbesinst_lazy_todo::set_todo");
      }
      assert(check_invariants());
    }
};

inline
std::ostream& operator<<(std::ostream& out, const pbesinst_lazy_todo& todo)
{
  return out << "todo = " << core::detail::print_list(todo.elements()) << " irrelevant = " << core::detail::print_list(todo.irrelevant_elements()) << std::endl;
}

/// \brief A PBES instantiation algorithm that uses a lazy strategy
class pbesinst_lazy_algorithm
{
  protected:
    /// \brief Algorithm options.
    const pbessolve_options& m_options;

    /// \brief Data rewriter.
    data::rewriter datar;

    /// \brief A PBES.
    pbes m_pbes;

    /// \brief A lookup map for PBES equations.
    pbes_equation_index m_equation_index;

    /// \brief The propositional variable instantiations that need to be handled.
    pbesinst_lazy_todo todo;

    /// \brief The propositional variable instantiations that have been discovered (not necessarily handled).
    atermpp::indexed_set<propositional_variable_instantiation, true> discovered;

    /// \brief The initial value (after rewriting).
    propositional_variable_instantiation init;

    // \brief The number of iterations
    std::size_t m_iteration_count = 0;

    // The data structures that must be separate per thread.
    /// \brief The rewriter.
    enumerate_quantifiers_rewriter m_global_R;

    // Mutexes
    utilities::mutex m_todo_access;

    volatile bool m_must_abort = false;

    // \brief Returns a status message about the progress
    virtual std::optional<std::string> status_message(std::size_t equation_count)
    {
      if (equation_count > 0 && equation_count % 1000 == 0)
      {
        std::ostringstream out;
        out << "Generated " << equation_count << " BES equations" << std::endl;
        return out.str();
      }

      return std::nullopt;
    }

    // instantiates global variables
    // simplifies the pbes
    pbes preprocess(const pbes& x) const
    {
      pbes p = x;
      pbes_system::detail::instantiate_global_variables(p);
      pbes_system::one_point_rule_rewriter one_point_rule_rewriter;
      pbes_system::simplify_quantifiers_data_rewriter<mcrl2::data::rewriter> simplify_rewriter(datar);
      for (pbes_equation& eqn: p.equations())
      {
        eqn.formula() = order_quantified_variables(one_point_rule_rewriter(simplify_rewriter(eqn.formula())), p.data());
      }
      return p;
    }

    static void rewrite_true_false(pbes_expression& result,
                                   const fixpoint_symbol& symbol,
                                   const propositional_variable_instantiation& X,
                                   const pbes_expression& psi
                                  )
    {
      bool changed = false;
      replace_propositional_variables(
        result,
        psi,
        [&](const propositional_variable_instantiation& Y) -> pbes_expression
        {
          if (Y == X)
          {
            changed = true;
            if (symbol.is_mu())
            {
              return false_();
            }
            else
            {
              return true_();
            }
          }
          return Y;
        }
      );
      if (changed)
      {
        simplify_rewriter simplify;
        const pbes_expression result1=result;
        simplify(result, result1);
      }
    }

    data::rewriter construct_rewriter(const pbes& pbesspec)
    {
      if (m_options.remove_unused_rewrite_rules)
      {
        return data::rewriter(pbesspec.data(),
                              data::used_data_equation_selector(pbesspec.data(), pbes_system::find_function_symbols(pbesspec), pbesspec.global_variables(), false),
                              m_options.rewrite_strategy);
      }
      else
      {
        return data::rewriter(pbesspec.data(), m_options.rewrite_strategy);
      }
    }

  public:

    /// \brief Constructor.
    /// \param p The pbes used in the exploration algorithm.
    /// \param options Te.
    /// \param search_strategy The search strategy used to explore the pbes, typically depth or breadth first.
    /// \param optimization An indication of the optimisation level.
    explicit pbesinst_lazy_algorithm(
      const pbessolve_options& options,
      const pbes& p,
      std::optional<data::rewriter> rewriter = std::nullopt
    )
     : m_options(options),
       datar(rewriter.has_value() ? rewriter.value() : construct_rewriter(p)),
       m_pbes(preprocess(p)),
       m_equation_index(p),
       discovered(m_options.number_of_threads),
       m_global_R(datar, p.data())
    { }

    virtual ~pbesinst_lazy_algorithm() = default;

    /// \brief Reports BES equations that are produced by the algorithm.
    /// This function is called for every BES equation X = psi with rank k that is produced. By default it does nothing.
    virtual void on_report_equation(const std::size_t /* thread_index */,
                                    const propositional_variable_instantiation& /* X */,
                                    const pbes_expression& /* psi */, std::size_t /* k */
                                   )
    { }

    /// \brief This function is called when new elements are added to discovered.
    virtual void on_discovered_elements(const std::set<propositional_variable_instantiation>& /* elements */)
    { }

    /// \brief This function is called right after the while loop is finished.
    virtual void on_end_while_loop()
    { }

    void next_todo(propositional_variable_instantiation& result)
    {
      if (m_options.exploration_strategy == breadth_first)
      {
        result = todo.front();
        todo.pop_front();
      }
      else
      {
        result = todo.back();
        todo.pop_back();
      }
    }

    const fixpoint_symbol& symbol(std::size_t i) const
    {
      return m_pbes.equations()[i].symbol();
    }

    // rewrite the right hand side of the equation X = psi
    virtual void rewrite_psi(const std::size_t /* thread_index */,
                             pbes_expression& result,
                             const fixpoint_symbol& symbol,
                             const propositional_variable_instantiation& X,
                             const pbes_expression& psi
                            )
    {  
      if (m_options.optimization >= partial_solve_strategy::remove_self_loops)
      {
        rewrite_true_false(result, symbol, X, psi);
      }
    }

    virtual bool solution_found(const propositional_variable_instantiation& /* init */) const
    {
      return false;
    }

    virtual void run_thread(const std::size_t thread_index,
                            pbesinst_lazy_todo& todo,
                            std::atomic<std::size_t>& number_of_active_processes,
                            data::mutable_indexed_substitution<> sigma,
                            enumerate_quantifiers_rewriter R
                           )
    {
      using utilities::detail::contains;

      if (m_options.number_of_threads > 1)
      {
        mCRL2log(log::debug) << "Start thread " << thread_index << ".\n";
      }
      R.thread_initialise();

      propositional_variable_instantiation X_e;
      pbes_expression psi_e;

      while (number_of_active_processes > 0)
      {
        m_todo_access.lock();
        while (!todo.elements().empty() && !m_must_abort)
        {
          ++m_iteration_count;

          if (std::optional<std::string> message = status_message(m_iteration_count))
          {
            mCRL2log(log::status) << *message;
          }

          detail::check_bes_equation_limit(m_iteration_count);

          next_todo(X_e);
          m_todo_access.unlock();

          std::size_t index = m_equation_index.index(X_e.name());
          const pbes_equation& eqn = m_pbes.equations()[index];
          const auto& phi = eqn.formula();
          data::add_assignments(sigma, eqn.variable().parameters(), X_e.parameters());
          R(psi_e, phi, sigma);
          R.clear_identifier_generator();
          data::remove_assignments(sigma, eqn.variable().parameters());

          // optional step
          m_todo_access.lock();
          rewrite_psi(thread_index, psi_e, eqn.symbol(), X_e, psi_e);
          m_todo_access.unlock();

          std::set<propositional_variable_instantiation> occ = find_propositional_variable_instantiations(psi_e);

          // report the generated equation
          std::size_t k = m_equation_index.rank(X_e.name());
          m_todo_access.lock();
          mCRL2log(log::debug) << "generated equation " << X_e << " = " << psi_e
                               << " with rank " << k << std::endl;
          on_report_equation(thread_index, X_e, psi_e, k);
          todo.insert(occ.begin(), occ.end(), discovered, thread_index);
          for (const auto& i : occ)
          {
            discovered.insert(i, thread_index);
          }
          on_discovered_elements(occ);

          if (solution_found(init))
          {
            break;
          }
        }
        m_todo_access.unlock();

        // Check whether all processes are ready. If so the
        // number_of_active_processes becomes 0. Otherwise, this thread becomes
        // active again, and tries to see whether the todo buffer is not empty,
        // to take up more work.
        number_of_active_processes--;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        if (number_of_active_processes > 0)
        {
          number_of_active_processes++;
        }
      }

      if (m_options.number_of_threads > 1)
      {
        mCRL2log(log::debug) << "Stop thread " << thread_index << ".\n";
      }
    }

    /// \brief Runs the algorithm. The result is obtained by calling the function \p get_result.
    virtual void run()
    {
      m_iteration_count = 0;

      const std::size_t number_of_threads = m_options.number_of_threads;
      const std::size_t initialisation_thread_index = (number_of_threads==1?0:1);
      std::atomic<std::size_t> number_of_active_processes = number_of_threads;
      std::vector<std::thread> threads;

      data::mutable_indexed_substitution<> sigma;
      if (m_options.replace_constants_by_variables)
      {
        pbes_system::replace_constants_by_variables(m_pbes, datar, sigma);
      }

      init = atermpp::down_cast<propositional_variable_instantiation>(m_global_R(m_pbes.initial_state(), sigma));
      todo.insert(init);
      discovered.insert(init, initialisation_thread_index);

      if (number_of_threads>1)
      {
        threads.reserve(number_of_threads);
        for (std::size_t i = 1; i <= number_of_threads; ++i)
        {
          std::thread tr([&, i](){
            run_thread(i,
                       todo,
                       number_of_active_processes,
                       sigma.clone(),
                       m_global_R.clone()
                      );
          });
          threads.push_back(std::move(tr));
        }

        for (std::size_t i = 1; i <= number_of_threads; ++i)
        {
          threads[i-1].join();
        }
      }
      else 
      {
        // There is only one thread. Run the process in the main thread, without cloning sigma or the rewriter.
        const std::size_t single_thread_index=0;
        run_thread(single_thread_index,
                   todo,
                   number_of_active_processes,
                   sigma,
                   m_global_R
                  );
      }
      on_end_while_loop();

      mCRL2log(log::verbose) << "Generated " << m_iteration_count << " BES equations" << std::endl;
    }

    const pbes_equation_index& equation_index() const
    {
      return m_equation_index;
    }

    enumerate_quantifiers_rewriter& rewriter()
    {
      return m_global_R;
    }

    const data::rewriter& data_rewriter() const
    {
      return datar;
    };
};

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_PBESINST_LAZY_H
