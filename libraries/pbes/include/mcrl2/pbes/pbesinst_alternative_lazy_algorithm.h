// Author(s): Xiao Qi
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesinst_alternative_lazy_algorithm.h
/// \brief An alternative lazy algorithm for instantiating a PBES, ported from
///        bes_deprecated.h.

#ifndef MCRL2_PBES_PBESINST_ALTERNATIVE_LAZY_ALGORITHM_H
#define MCRL2_PBES_PBESINST_ALTERNATIVE_LAZY_ALGORITHM_H

#include "mcrl2/pbes/remove_level.h"
#include "mcrl2/pbes/detail/check_well_formed_bes.h"
#include "mcrl2/pbes/pbesinst_algorithm.h"
#include "mcrl2/pbes/search_strategy.h"
#include "mcrl2/pbes/transformation_strategy.h"



namespace mcrl2::pbes_system
{

namespace detail
{
  // The following function is a helper function to allow to create m_pv_renaming outside
  // the class such that the class becomes a lightweight object.

inline std::unordered_map<propositional_variable_instantiation, propositional_variable_instantiation>
create_pv_renaming(std::vector<std::vector<propositional_variable_instantiation>>& instantiations,
    bool short_renaming_scheme)
{
  std::size_t index = 0;
  std::unordered_map<propositional_variable_instantiation, propositional_variable_instantiation> pv_renaming;
  for (const std::vector<propositional_variable_instantiation>& vec: instantiations)
  {
    for (const propositional_variable_instantiation& inst: vec)
    {
      if (short_renaming_scheme)
      {
        std::stringstream ss;
        ss << "X" << index;
        pv_renaming[inst] = propositional_variable_instantiation(ss.str(), data::data_expression_list());
      }
      else
      {
        pv_renaming[inst] = pbesinst_rename()(inst);
      }
      index++;
    }
  }
  return pv_renaming;
  }

  class rename_pbesinst_consecutively
  {
    protected:
      const std::unordered_map<propositional_variable_instantiation,propositional_variable_instantiation>& m_pv_renaming;

    public:
      rename_pbesinst_consecutively(const std::unordered_map<propositional_variable_instantiation,propositional_variable_instantiation>& pv_renaming)
       :  m_pv_renaming(pv_renaming)
      {}


      propositional_variable_instantiation operator()(const propositional_variable_instantiation& v) const
      {
        assert(m_pv_renaming.count(v)>0);
        return m_pv_renaming.at(v);
      }
  };
} // end namespace detail


inline mcrl2::pbes_system::pbes_expression pbes_expression_order_quantified_variables(
              const mcrl2::pbes_system::pbes_expression& p, const mcrl2::data::data_specification& data_spec)
{
  using namespace mcrl2;
  using namespace mcrl2::pbes_system;

  if (is_pbes_and(p))
  {
    const and_& pa=atermpp::down_cast<and_>(p);
    return and_(pbes_expression_order_quantified_variables(pa.left(),data_spec),
                pbes_expression_order_quantified_variables(pa.right(),data_spec));
  }
  else if (is_pbes_or(p))
  {
    const or_& po=atermpp::down_cast<or_>(p);
    return or_(pbes_expression_order_quantified_variables(po.left(),data_spec),
               pbes_expression_order_quantified_variables(po.right(),data_spec));
  }
  else if (is_pbes_imp(p))
  {
    const imp& pi=atermpp::down_cast<imp>(p);
    return imp(pbes_expression_order_quantified_variables(pi.left(),data_spec),
               pbes_expression_order_quantified_variables(pi.right(),data_spec));
  }
  else if (is_pbes_not(p))
  {
    return not_(pbes_expression_order_quantified_variables(atermpp::down_cast<not_>(p).operand(),data_spec));
  }
  else if (is_pbes_forall(p))
  {
    const forall& pf=atermpp::down_cast<forall>(p);
    const pbes_expression expr = pbes_expression_order_quantified_variables(pf.body(),data_spec);
    return optimized_forall(mcrl2::data::order_variables_to_optimise_enumeration(pf.variables(),data_spec),expr);
  }
  else if (is_pbes_exists(p))
  {
    const exists& pe=atermpp::down_cast<exists>(p);
    const pbes_expression expr = pbes_expression_order_quantified_variables(pe.body(),data_spec);
    return optimized_exists(mcrl2::data::order_variables_to_optimise_enumeration(pe.variables(),data_spec),expr);
  }
  else
  {
    return p;
  }
}

/// \brief An alternative lazy algorithm for instantiating a PBES, ported from
///         bes_deprecated.h.
class pbesinst_alternative_lazy_algorithm
{
  protected:
    const data::data_specification& m_data_spec;

    /// \brief Data rewriter.
    const data::rewriter& m_datar;

    /// \brief The rewriter.
    enumerate_quantifiers_rewriter R;

    /// \brief Initial value for regeneration_period.
    static const std::size_t regeneration_count_init = 100;

    /// The maximum size that the todo buffer is allowed to have.
    const std::size_t m_maximum_todo_size;

    /// The variable m_approximimate_true indicates whether
    /// boolean variables that cannot be dealt with as the todo buffer
    /// would otherwise exceeds m_maximum_todo_size are set to true or
    /// to false
    const bool m_approximate_true;

    /// When the todo buffer is limited, due to m_maximum_todo_size, then the variable below counts how
    /// many elements are dropped out of the todo buffer.
    std::size_t m_elements_not_stored_in_todo_buffer = 0;

    /// Indicate to which extent explored bes equations that turn out not to reachable can be thrown away.
    /// Values are: none, some or all.
    const mcrl2::pbes_system::remove_level m_erase_unused_bes_variables;

    /// \brief Propositional variable instantiations that need to be handled.
    std::deque<propositional_variable_instantiation> todo;

    /// \brief The content of todo as a set.
    std::unordered_set<propositional_variable_instantiation> todo_set;

    /// \brief Map a variable instantiation to a set of other variable
    ///        instantiations on whose right hand sides it appears.
    std::unordered_map<propositional_variable_instantiation, std::unordered_set<propositional_variable_instantiation> > occurrence;

    /// \brief Map a variable instantiation to its right hand side.
    std::unordered_map<propositional_variable_instantiation, pbes_expression> equation;

    /// \brief Map a variable instantiations to its right hand side
    ///        when the latter is trivial (either true or false).
    std::unordered_map<propositional_variable_instantiation, pbes_expression> trivial;

    /// \brief instantiations[i] contains all instantiations of the variable
    ///        of the i-th equation in the PBES.
    std::vector<std::vector<propositional_variable_instantiation> > instantiations;

    /// \brief symbols[i] contains the fixedpoint symbol of the i-th equation
    ///        in the PBES.
    std::vector<fixpoint_symbol> symbols;

    /// \brief ranks[i] contains the rank of the i-th equation in the PBES.
    std::vector<std::size_t> ranks;

    /// \brief The initial value.
    propositional_variable_instantiation init;

    /// \brief A lookup map for PBES equations.
    std::unordered_map<core::identifier_string, std::size_t> equation_index;

    /// \brief The search strategy to use when exploring the state space.
    search_strategy m_search_strategy;

    /// \brief Transformation strategy.
    transformation_strategy m_transformation_strategy;

    /// \brief Prints a log message for every 1000-th equation
    void print_equation_count(const std::size_t nr_of_processed_variables,
                              const std::size_t nr_of_generated_variables,
                              const std::size_t todo_size) const
    {
      static time_t last_log_time = time(nullptr) - 1;
      time_t new_log_time=0;
      if (time(&new_log_time) > last_log_time)
      {
        last_log_time = new_log_time;
        mCRL2log(mcrl2::log::status) << "Processed " << nr_of_processed_variables <<
                       " and generated " << nr_of_generated_variables <<
                       " boolean variables";
        if (m_maximum_todo_size != std::numeric_limits<std::size_t>::max())
        {
          mCRL2log(mcrl2::log::status) << " with a todo buffer of size " << todo_size << ".    \n";
        }
        else
        {
          mCRL2log(mcrl2::log::status) << ".   \n"; // Extra spaces, as the log messages are shown repeatedly on top of each other and can shorten, leaving residual letter on the output. 
        }
      }
    }

  public:

    /// \brief Constructor.
    /// \param data_spec A data specification.
    /// \param datar A data rewriter.
    /// \param search_strategy A search strategy (e.g. breadth_first).
    /// \param transformation_strategy A strategy to transform the PBES into a BES.
    /// \param erase_unused_bes_variables An indicator how often unreachable variables must be garbage collected.
    /// \param maximum_todo_size An indication of the maximal size of the stack with explorable
    ///                          BES variables. If the stack exceeds this size, the ignored variable is set to
    ///                          true or false, depending on the parameter \p approximate_true.
    /// \param approximate_true If true BES variables that are not investigated are set to false. If false
    ///                         these variables are set to true.
    pbesinst_alternative_lazy_algorithm(const data::data_specification& data_spec,
        const data::rewriter& datar,
        search_strategy search_strategy = breadth_first,
        transformation_strategy transformation_strategy = lazy,
        const mcrl2::pbes_system::remove_level erase_unused_bes_variables = mcrl2::pbes_system::none,
        const std::size_t maximum_todo_size = std::numeric_limits<std::size_t>::max(),
        const bool approximate_true = true)
        : m_data_spec(data_spec),
          m_datar(datar),
          R(datar, data_spec),
          m_maximum_todo_size(maximum_todo_size),
          m_approximate_true(approximate_true),

          m_erase_unused_bes_variables(erase_unused_bes_variables),
          m_search_strategy(search_strategy),
          m_transformation_strategy(transformation_strategy)
    {
      // Initialize the random generator, with an arbitrary seed, depending on a new time.
      time_t t=time(nullptr);
      for (; t == time(nullptr);)
      {
        ; // Wait until time changes.
      }
      srand((unsigned)time(nullptr));

      if (m_search_strategy == breadth_first_short)
      {
        m_search_strategy = breadth_first;
      }
      else if (m_search_strategy == depth_first_short)
      {
        m_search_strategy = depth_first;
      }
    }

    inline propositional_variable_instantiation next_todo()
    {
      assert(todo.size()==todo_set.size());
      if (m_search_strategy == breadth_first)
      {
        const propositional_variable_instantiation X_e = todo.front();
        todo.pop_front();
        todo_set.erase(X_e);
        return X_e;
      }
      else
      {
        const propositional_variable_instantiation X_e = todo.back();
        todo.pop_back();
        todo_set.erase(X_e);
        return X_e;
      }
    }

    inline void add_todo(const propositional_variable_instantiation& X)
    {
      if (todo.size()<m_maximum_todo_size)  // If there is no limit on todo, m_maximimum_todo_size is equal to npos.
      {
        todo.push_back(X);
        todo_set.insert(X); // XXXXX
      }
      else
      {
        ++m_elements_not_stored_in_todo_buffer;
        if ((rand() % (todo.size() + m_elements_not_stored_in_todo_buffer)) < todo.size())
        {
          std::size_t index = rand() % (todo.size());
          const propositional_variable_instantiation Y=todo[index];
          equation[Y]=(m_approximate_true?false_():true_());
          trivial[Y]=equation[Y];
          instantiations[equation_index[Y.name()]].push_back(Y);

          todo_set.erase(Y);
          todo[index]=X;
          todo_set.insert(X);
        }
        else
        {
          equation[X]=(m_approximate_true?false_():true_());
          instantiations[equation_index[X.name()]].push_back(X);
          trivial[X]=equation[X];
        }
      }
    }

    std::size_t get_rank(const propositional_variable_instantiation& X)
    {
      return ranks[equation_index[X.name()]];
    }

    template <bool is_mu>
    bool find_loop_rec(
        const pbes_expression& expr,
        propositional_variable_instantiation X,
        std::size_t rank,
        std::unordered_map<propositional_variable_instantiation, bool>& visited)
    {
      if (is_false(expr) || is_true(expr))
      {
        return false;
      }
      if (is_propositional_variable_instantiation(expr))
      {
        const propositional_variable_instantiation& Y = atermpp::vertical_cast<propositional_variable_instantiation>(expr);
        if (Y == X)
        {
          return true;
        }
        else if (get_rank(Y) != rank)
        {
          return false;
        }
        if (visited.count(Y))
        {
          return visited[Y];
        }
        if (equation.count(Y) == 0)
        {
          return false;
        }
        visited[Y] = false;
        bool b = find_loop_rec<is_mu>(equation[Y], X, rank, visited);
        visited[Y] = b;
        return b;
      }

      if (is_mu)
      {
        if (is_and(expr))
        {
          const and_& expra=atermpp::down_cast<and_>(expr);
          return find_loop_rec<is_mu>(expra.left(), X, rank, visited) ||
                 find_loop_rec<is_mu>(expra.right(), X, rank, visited);
        }
        if (is_or(expr))
        {
          const or_& expro=atermpp::down_cast<or_>(expr);
          return find_loop_rec<is_mu>(expro.left(), X, rank, visited) &&
                 find_loop_rec<is_mu>(expro.right(), X, rank, visited);
        }
      }
      else
      {
        if (is_and(expr))
        {
          const and_& expra=atermpp::down_cast<and_>(expr);
          return find_loop_rec<is_mu>(expra.left(), X, rank, visited) &&
                 find_loop_rec<is_mu>(expra.right(), X, rank, visited);
        }
        if (is_or(expr))
        {
          const or_& expro=atermpp::down_cast<or_>(expr);
          return find_loop_rec<is_mu>(expro.left(), X, rank, visited) ||
                 find_loop_rec<is_mu>(expro.right(), X, rank, visited);
        }
      }
      return false;
    }

    template <bool is_mu>
    bool find_loop(pbes_expression expr, propositional_variable_instantiation X)
    {
      std::unordered_map<propositional_variable_instantiation, bool> visited;
      return find_loop_rec<is_mu>(expr, X, get_rank(X), visited);
    }

    void regenerate_states()
    {
      if (m_erase_unused_bes_variables==pbes_system::none)
      {
        // Nothing will be thrown away. Therefore, it makes no sense to rebuild
        // the data structures on the basis of reachable bes variables.
        return;
      }

      // Create a set of reachable propositional_variable_instantiations
      // and use that to clean up the set of equations.
      std::unordered_set<propositional_variable_instantiation> reachable;

      todo.clear();
      todo_set.clear();
      occurrence.clear();
      for(std::vector<propositional_variable_instantiation>& vec: instantiations)
      {
        vec.clear();
      }
      assert(instantiations.size()==0 || instantiations[0].size()==0);


      std::stack<pbes_expression> stack;
      stack.push(init);

      while (!stack.empty())
      {
        const pbes_expression expr = stack.top();
        stack.pop();

        if (is_propositional_variable_instantiation(expr))
        {
          const propositional_variable_instantiation& X = atermpp::vertical_cast<propositional_variable_instantiation>(expr);
          if (reachable.count(X) == 0)
          {
            if (equation.count(X)>0)
            {
              stack.push(equation[X]);
              reachable.insert(X);
            }
            else
            {
              add_todo(X);
            }
          }
        }
        else if (is_and(expr))
        {
          const and_& expra=atermpp::down_cast<and_>(expr);
          stack.push(expra.left());
          stack.push(expra.right());
        }
        else if (is_or(expr))
        {
          const or_& expro=atermpp::down_cast<or_>(expr);
          stack.push(expro.left());
          stack.push(expro.right());
        }
      }
      // erase non reachable equations.
      std::unordered_map<propositional_variable_instantiation, pbes_expression> new_equations;
      for (const auto& i: equation)
      {
        // Insert the new equation if it is reachable, or if it equal to true or false and m_erase_unused_bes_variables is set to some.
        if (reachable.count(i.first) > 0
            || (m_erase_unused_bes_variables == pbes_system::some && (is_true(i.second) || is_false(i.second))))
        {
          new_equations.insert(i);
          std::size_t index = equation_index[i.first.name()];
          instantiations[index].push_back(i.first);
          std::set<propositional_variable_instantiation> rhs_variables
            = find_propositional_variable_instantiations(i.second);
          for (const propositional_variable_instantiation& v: rhs_variables)
          {
            occurrence[v].insert(i.first);
          }

        }
      }
      equation.swap(new_equations);
    }

    // The function below simplifies a boolean_expression, given the knowledge that some propositional variables in trivial
    // are known to be true or false. The idea is that variables that are redundant can be removed. If p = p1 && p2, and p1 is
    // false, then p2 can be removed, as its value does not influence the rewrite system.
    // The result of the function is a pair, with the simplified expression as first term, and the expression that is rewritten under the
    // simplifications in trivial as the second term.
    using pbes_expression_pair = std::pair<pbes_expression, pbes_expression>;
    pbes_expression_pair simplify_pbes_expression(const pbes_expression& p, const std::unordered_map<propositional_variable_instantiation, pbes_expression>& trivial)
    {
      if (is_propositional_variable_instantiation(p))
      {
        const std::unordered_map<propositional_variable_instantiation, pbes_expression>::const_iterator i=
                                            trivial.find(atermpp::down_cast<propositional_variable_instantiation>(p));
        if (i!=trivial.end() && (is_true(i->second) || is_false(i->second)))
        {
          return pbes_expression_pair(p,i->second);
        }
        return pbes_expression_pair(p,p);
      }
      else if (is_true(p)||is_false(p))
      {
        return pbes_expression_pair(p,p);
      }
      else if (is_and(p))
      {
        const and_& pa=atermpp::down_cast<and_>(p);
        const pbes_expression_pair lhs=simplify_pbes_expression(pa.left(),trivial);
        const pbes_expression_pair rhs=simplify_pbes_expression(pa.right(),trivial);
        if (is_false(lhs.second))
        {
          return lhs;
        }
        if (is_false(rhs.second))
        {
          return rhs;
        }
        if (is_true(lhs.second))
        {
          return pbes_expression_pair(and_(lhs.first,rhs.first),rhs.second);
        }
        if (is_true(rhs.second))
        {
          return pbes_expression_pair(and_(lhs.first,rhs.first),lhs.second);
        }
        return pbes_expression_pair(and_(lhs.first,rhs.first),and_(lhs.first,rhs.first));
      }
      assert(is_or(p));
      const or_& po=atermpp::down_cast<or_>(p);
      const pbes_expression_pair lhs=simplify_pbes_expression(po.left(),trivial);
      const pbes_expression_pair rhs=simplify_pbes_expression(po.right(),trivial);
      if (is_true(lhs.second))
      {
        return lhs;
      }
      if (is_true(rhs.second))
      {
        return rhs;
      }
      if (is_false(lhs.second))
      {
        return pbes_expression_pair(or_(lhs.first,rhs.first),rhs.second);
      }
      if (is_false(rhs.second))
      {
        return pbes_expression_pair(or_(lhs.first,rhs.first),lhs.second);;
      }
      return pbes_expression_pair(or_(lhs.first,rhs.first),or_(lhs.first,rhs.first));
    }



    /// \brief Runs the algorithm. The result is obtained by calling the function \p get_result.
    /// \param p A PBES
    void run(pbes& p)
    {
      using utilities::detail::pick_element;
      using utilities::detail::contains;

      std::size_t regeneration_count=regeneration_count_init;
      pbes_system::detail::instantiate_global_variables(p);

      std::vector<pbes_equation>& pbes_equations = p.equations();

      // simplify all right hand sides of p
      //
      // NOTE: This is not just an optimization. There are certain PBES
      // equations for which applying enumerate_quantifiers_rewriter directly
      // won't terminate, like:
      //
      // forall m: Nat . exists k: Nat . val(m == k)
      pbes_system::one_point_rule_rewriter one_point_rule_rewriter;
      pbes_system::simplify_quantifiers_data_rewriter<mcrl2::data::rewriter> simplify_rewriter(m_datar);
      for (pbes_equation& eq: pbes_equations)  // & is important as we change the equation.
      {
        eq.formula() = pbes_expression_order_quantified_variables(one_point_rule_rewriter(simplify_rewriter(eq.formula())), m_data_spec);
      }

      // initialize equation_index, instantiations, symbols and ranks
      std::size_t eqn_index = 0;
      ranks.resize(pbes_equations.size());
      instantiations.resize(pbes_equations.size());
      for (const pbes_equation& eqn: pbes_equations)
      {
        equation_index[eqn.variable().name()] = eqn_index;
        symbols.push_back(eqn.symbol());
        if (eqn_index > 0)
        {
          ranks[eqn_index] = ranks[eqn_index-1] + (symbols[eqn_index] == symbols[eqn_index-1] ? 0 : 1);
        }
        ++eqn_index;
      }

      init = atermpp::down_cast<propositional_variable_instantiation>(R(p.initial_state()));
      add_todo(init);
      while (!todo.empty())
      {
        const propositional_variable_instantiation X_e = next_todo();
        std::size_t index = equation_index[X_e.name()];
        instantiations[index].push_back(X_e);

        const pbes_equation& eqn = pbes_equations[index];
        data::rewriter::substitution_type sigma;
        make_pbesinst_substitution(eqn.variable().parameters(), X_e.parameters(), sigma);
        const pbes_expression& phi = eqn.formula();
        pbes_expression psi_e = R(phi, sigma);
        R.clear_identifier_generator();
        try
        {
          check_whether_argument_is_a_well_formed_bes(psi_e);
        }
        catch (mcrl2::runtime_error& e)
        {
          throw mcrl2::runtime_error("Generated boolean equation system is not well formed.\n" + std::string(e.what()));
        }
        pbes_expression rewritten_psi_e;

        if (m_transformation_strategy >= optimize)
        {
          // Substitute all trivial variable instantiations by their values
          pbes_expression_pair p=simplify_pbes_expression(psi_e,trivial);
          psi_e=p.first;
          rewritten_psi_e=p.second;
        }
        else
        {
          rewritten_psi_e=psi_e;
        }
        // Store the result
        equation[X_e] = psi_e;

        if (m_transformation_strategy >= on_the_fly_with_fixed_points)
        {
          // Find mu or nu loop
          if (eqn.symbol() == fixpoint_symbol::mu())
          {
            if (find_loop<true>(psi_e, X_e))
            {
              rewritten_psi_e = false_();
            }
          }
          else
          {
            if (find_loop<false>(psi_e, X_e))
            {
              rewritten_psi_e = true_();
            }
          }
        }


        // Add all variable instantiations in psi_e to todo and generated,
        // and augment the occurrence sets
        std::set<propositional_variable_instantiation> psi_variables = find_propositional_variable_instantiations(psi_e);
        for (const propositional_variable_instantiation& v: psi_variables)
        {
          if (todo_set.count(v) == 0 && equation.count(v) == 0)
          {
            add_todo(v);
          }
          occurrence[v].insert(X_e);
        }

        if (m_transformation_strategy >= optimize && (is_true(rewritten_psi_e) || is_false(rewritten_psi_e)))
        {
          trivial[X_e] = rewritten_psi_e;
          if (m_transformation_strategy >= on_the_fly)
          {
            // Substitute X_e to its value in all its occurrences, and
            // substitute all other variables to their values that are found
            // to be either true or false in all their occurrences.
            std::stack<propositional_variable_instantiation> new_trivials;
            new_trivials.push(X_e);
            while (!new_trivials.empty())
            {
              const propositional_variable_instantiation X = new_trivials.top();
              new_trivials.pop();

              const std::unordered_set<propositional_variable_instantiation> oc = occurrence[X];
              std::unordered_map<propositional_variable_instantiation, pbes_expression> trivial_X;
              trivial_X[X] = psi_e;
              for (const propositional_variable_instantiation& Y: oc)
              {
                pbes_expression_pair p=simplify_pbes_expression(equation[Y],trivial);
                equation[Y]=p.first;
                const pbes_expression f=p.second;
                if (is_true(f) || is_false(f))
                {
                  trivial[Y] = f;
                  new_trivials.push(Y);
                }
              }
              occurrence.erase(X);
            }
          }
        }

        if (m_transformation_strategy >= on_the_fly)
        {
          if (--regeneration_count == 0 || trivial.count(init)>0 )
          {
            regeneration_count = equation.size() / 2;
            regenerate_states();
          }
        }

        print_equation_count(equation.size(), equation.size()+todo.size(), todo.size()); // Print the number of equations every second in verbose mode.
        detail::check_bes_equation_limit(equation.size());
      }
      // Remove unnessary equations.
      regenerate_states();
    }

    /// \brief Returns the computed bes in pbes format
    /// \return The computed bes in pbes format
    pbes get_result(bool short_rename_scheme=true)
    {
      mCRL2log(log::verbose) << "Generated " << equation.size() << " BES equations in total, generating BES" << std::endl;
      pbes result;
      std::size_t index = 0;
      const std::unordered_map<propositional_variable_instantiation,propositional_variable_instantiation>
            pv_renaming = detail::create_pv_renaming(instantiations,short_rename_scheme);
      detail::rename_pbesinst_consecutively renamer(pv_renaming);
      for (const std::vector<propositional_variable_instantiation>& vec: instantiations)
      {
        const fixpoint_symbol symbol = symbols[index++];
        for (const propositional_variable_instantiation& X_e: vec)
        {
          const propositional_variable lhs = propositional_variable(renamer(X_e).name(), data::variable_list());
          const pbes_expression rhs = replace_propositional_variables(equation[X_e], renamer);
          result.equations().emplace_back(symbol, lhs, rhs);
          mCRL2log(log::debug) << "BESEquation: " << atermpp::aterm(symbol) << " " << lhs << " = " << rhs << std::endl;

        }
      }

      result.initial_state() = renamer(init);
      return result;
    }

    enumerate_quantifiers_rewriter& rewriter()
    {
      return R;
    }
};

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_PBESINST_ALTERNATIVE_LAZY_ALGORITHM_H
