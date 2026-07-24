// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/pres2res.h
/// \brief A straightforward algorithm that gets pres, and instantiates it to a res. 

#ifndef MCRL2_PRES_PRES2RES_H
#define MCRL2_PRES_PRES2RES_H

#include "mcrl2/atermpp/standard_containers/vector.h"
#include "mcrl2/atermpp/standard_containers/indexed_set.h"
#include "mcrl2/atermpp/function_symbol_generator.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/pbes/pbes_equation_index.h"
#include "mcrl2/pres/pressolve_options.h"
#include "mcrl2/pres/rewriters/enumerate_quantifiers_rewriter.h"
#include "mcrl2/pres/is_res.h"

namespace mcrl2::pres_system
{

class variable_replace_builder: public pres_expression_builder <variable_replace_builder>
{
    const atermpp::indexed_set<propositional_variable_instantiation>& m_stored_variables;
    const std::vector<propositional_variable_instantiation>& m_new_pres_variables;

  public:
    using super = pres_expression_builder<variable_replace_builder>;
    using super::apply;

    variable_replace_builder(const atermpp::indexed_set<propositional_variable_instantiation>& stored_variables,
                             const std::vector<propositional_variable_instantiation>& new_pres_variables)
      :  m_stored_variables(stored_variables),
         m_new_pres_variables(new_pres_variables)
    {}

    template <class T>
    void apply(T& result, const propositional_variable_instantiation& x)
    {
      std::size_t index = m_stored_variables.index(x);
      assert(m_new_pres_variables.size()>index);
      result = propositional_variable_instantiation(m_new_pres_variables.at(index));
    }
};

/// \brief A simple algorithm that takes a pres, instantiates it into a res, without parameterized 
///        variables, and sum, infimum and supremum operators.
/// \details It does not attempt to optimize the solution, by partially solving it. It only applies some
///          simple rewriting, and small ad hoc optimisations. 
/// The result will be put in the structure graph that is passed in the constructor.
class pres2res_algorithm
{
  protected:
    pressolve_options m_options;
    pres m_input_pres;
    atermpp::indexed_set<propositional_variable_instantiation> m_stored_variables;
    pbes_system::pbes_equation_index m_equation_index;
    std::size_t m_next_res_variable_index=0;
    data::mutable_indexed_substitution<> m_sigma;
    atermpp::function_symbol_generator fresh_identifier_generator;

    enumerate_quantifiers_rewriter m_R;   // The rewriter.

  public:
    pres2res_algorithm(
         const pressolve_options& options,
         const pres& input_pres,
         enumerate_quantifiers_rewriter& R)
     : m_options(options),
       m_input_pres(input_pres),
       m_equation_index(input_pres),
       fresh_identifier_generator("X"),
       m_R(R) 
    {}

    pres run()
    {
      std::size_t max_rank=m_equation_index.max_rank();
      std::vector<atermpp::vector<pres_equation>> generated_equations(max_rank+1);
      std::ignore = m_stored_variables.insert(atermpp::down_cast<propositional_variable_instantiation>(m_R(m_input_pres.initial_state())));
      assert(m_stored_variables.index(atermpp::down_cast<propositional_variable_instantiation>(m_R(m_input_pres.initial_state())))==0);
      std::vector<propositional_variable_instantiation> new_pres_variables;
      while (m_next_res_variable_index<m_stored_variables.size())
      {
        const propositional_variable_instantiation& X_e=m_stored_variables.at(m_next_res_variable_index);
        m_next_res_variable_index++;
        new_pres_variables.emplace_back(fresh_identifier_generator(), data::data_expression_list());

        std::size_t index = m_equation_index.index(X_e.name());
        const pres_equation& eqn = m_input_pres.equations()[index];
        const pres_expression& phi = eqn.formula();
        data::add_assignments(m_sigma, eqn.variable().parameters(), X_e.parameters());
        pres_expression psi_e;
        m_R(psi_e, phi, m_sigma);
        std::string error_message;
        if (!is_res(psi_e, error_message))
        {
          throw mcrl2::runtime_error("Generation of RES failed. " + error_message);
        }
        m_R.clear_identifier_generator();
        data::remove_assignments(m_sigma, eqn.variable().parameters());

        std::set<propositional_variable_instantiation> occ = find_propositional_variable_instantiations(psi_e);

        for (const propositional_variable_instantiation& v: occ)
        {
          std::ignore = m_stored_variables.insert(v);
        }

        // save the generated equation
        std::size_t k = m_equation_index.rank(X_e.name());
        
        mCRL2log(log::debug) << "generated equation " << (k%2==0?pres_system::fixpoint_symbol::nu():pres_system::fixpoint_symbol::mu()) << " " 
                             << new_pres_variables[m_stored_variables.index(X_e)].name() << " = " << psi_e
                             << " with rank " << k << " for " << X_e << "\n";

        assert(generated_equations.size()>k);
        generated_equations[k].push_back(pres_equation((k%2==0?pres_system::fixpoint_symbol::nu():pres_system::fixpoint_symbol::mu()),
                                                        pbes_system::propositional_variable(new_pres_variables[m_stored_variables.index(X_e)].name(),
                                                                                            data::variable_list()), 
                                                        psi_e));
        if (generated_equations[k].size() % 1000 == 0)
        {
          mCRL2log(log::verbose) << "Generated " << (generated_equations[k].size() / 1000)
                                 << "k RES variables with rank " << k << ".\n";
        }
      }
      
      variable_replace_builder variable_replacer(m_stored_variables, new_pres_variables);
      std::vector<pres_equation> resulting_equations;
      pres_expression result;
      mCRL2log(log::debug) << "generated equation in final form.\n";
      for(atermpp::vector<pres_equation>& eqns: generated_equations)
      {
        for(const pres_equation& eqn: eqns)
        {
          variable_replacer.apply(result, eqn.formula());
          resulting_equations.emplace_back(eqn.symbol(), eqn.variable(), result);
          mCRL2log(log::debug) <<  resulting_equations.back() << "\n";
        }
        eqns=atermpp::vector<pres_equation>(); // clear the equations.
      } 

      return pres(m_input_pres.data(), 
                  resulting_equations, 
                  propositional_variable_instantiation(new_pres_variables.front().name(),data::data_expression_list()));
    }
};

} // namespace mcrl2::pres_system

#endif // MCRL2_PRES_PRES2RES_H
