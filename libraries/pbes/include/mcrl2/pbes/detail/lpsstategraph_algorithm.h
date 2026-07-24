// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/detail/lpsstategraph_algorithm.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_DETAIL_LPSSTATEGRAPH_ALGORITHM_H
#define MCRL2_PBES_DETAIL_LPSSTATEGRAPH_ALGORITHM_H

#include "mcrl2/lps/detail/instantiate_global_variables.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/is_well_typed.h"
#include "mcrl2/pbes/detail/stategraph_global_reset_variables.h"
#include "mcrl2/pbes/detail/stategraph_local_reset_variables.h"

namespace mcrl2::pbes_system::detail {

class lpsstategraph_algorithm;
void lps_reset_variables(lpsstategraph_algorithm& algorithm,
                         const pbes_expression& x,
                         const stategraph_equation& eq_X,
                         const data::variable_list& process_parameters,
                         std::vector<lps::action_summand>& summands
);

class lpsstategraph_algorithm: public local_reset_variables_algorithm
{
  protected:
    friend struct lps_reset_traverser;
    using super = local_reset_variables_algorithm;

    const lps::specification& m_original_lps;
    lps::specification m_transformed_lps; // will contain the result of the computation

    pbes construct_stategraph_pbes(const lps::specification& lpsspec)
    {
      pbes pbesspec;
      data::data_specification dataspec = lpsspec.data();

      // create 'used' functions for every sort in the free variables of actions
      data::set_identifier_generator generator;
      core::identifier_string used = generator("used");
      std::set<data::sort_expression> action_variable_sorts;
      for (const lps::action_summand& summand: lpsspec.process().action_summands())
      {
        for (const data::variable& v: lps::find_free_variables(summand.multi_action()))
        {
          action_variable_sorts.insert(v.sort());
        }
      }
      std::map<data::sort_expression, data::function_symbol> functions;
      for (const data::sort_expression& s: action_variable_sorts)
      {
        data::function_symbol f(used, data::function_sort({ s }, data::sort_bool::bool_()));
        functions[s] = f;
        dataspec.add_mapping(f);
      }

      // the PBES variable
      core::identifier_string X("X");

      // create a PBES equation that contains a conjunct for each action summand
      std::vector<pbes_expression> conjuncts;
      for (const lps::action_summand& summand: lpsspec.process().action_summands())
      {
        const auto& ei = summand.summation_variables();
        const auto& ci = summand.condition();
        data::data_expression_list gi = summand.next_state(lpsspec.process().process_parameters());

        data::data_expression x = ci;
        for (const data::variable& v: lps::find_free_variables(summand.multi_action()))
        {
          x = data::and_(x, data::application(functions[v.sort()], v ));
        }
        propositional_variable_instantiation Xi(X, gi);
        conjuncts.push_back(pbes_system::optimized_forall(ei, imp(atermpp::down_cast<pbes_expression>(x), Xi)));
      }
      for (const lps::deadlock_summand& summand: lpsspec.process().deadlock_summands())
      {
        const auto& ei = summand.summation_variables();
        const auto& ci = summand.condition();
        data::data_expression_list gi = data::make_data_expression_list(lpsspec.process().process_parameters());
        propositional_variable_instantiation Xi(X, gi);
        conjuncts.push_back(pbes_system::optimized_forall(ei, imp(atermpp::down_cast<pbes_expression>(ci), Xi)));
      }

      // N.B. It is essential that the order in which the conjuncts are traversed in a PBES matches the order of the corresponding summands.
      pbes_expression rhs = join_and(conjuncts.begin(), conjuncts.end());
      pbes_equation eqn(fixpoint_symbol::nu(), propositional_variable(X, lpsspec.process().process_parameters()), rhs);

      pbesspec.data() = dataspec;
      pbesspec.initial_state() = propositional_variable_instantiation(X, lpsspec.initial_process().expressions());
      pbesspec.equations().push_back(eqn);
      return pbesspec;
    }

    // Applies resetting of variables to the original PBES p.
    void reset_variables_to_original(lps::specification& lpsspec)
    {
      mCRL2log(log::debug) << "=== resetting variables to the original LPS ---" << std::endl;

      stategraph_equation& eqn = m_pbes.equations().front();
      lps_reset_variables(*this, eqn.formula(), eqn, lpsspec.process().process_parameters(), lpsspec.process().action_summands());
    }

  public:
    lpsstategraph_algorithm(const lps::specification& lpsspec, const pbesstategraph_options& options)
      : local_reset_variables_algorithm(construct_stategraph_pbes(lpsspec), options), m_original_lps(lpsspec)
    {}

    void execute_postprocessing() override
    {
      m_transformed_lps = m_original_lps;
      compute_occurring_data_parameters();
      start_timer("reset_variables_to_original");
      reset_variables_to_original(m_transformed_lps);
      finish_timer("reset_variables_to_original");
    }

    const lps::specification& original_lps() const
    {
      return m_original_lps;
    }

    const lps::specification& result() const
    {
      return m_transformed_lps;
    }
};

inline
data::assignment_list make_assignments(const data::variable_list& d, const data::data_expression_list& e)
{
  std::vector<data::assignment> result;
  auto di = d.begin();
  auto ei = e.begin();
  for (; di != d.end(); ++di, ++ei)
  {
    if (*di != *ei)
    {
      result.emplace_back(*di, *ei);
    }
  }
  return data::assignment_list(result.begin(), result.end());
}

struct lps_reset_traverser: public pbes_expression_traverser<lps_reset_traverser>
{
  using super = pbes_expression_traverser<lps_reset_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  lpsstategraph_algorithm& algorithm;
  const stategraph_equation& eq_X;
  const data::variable_list& process_parameters;
  std::vector<lps::action_summand>& summands;
  std::size_t i = 0; // the index of a summand

  lps_reset_traverser(lpsstategraph_algorithm& algorithm_,
                      const stategraph_equation& eq_X_,
                      const data::variable_list& process_parameters_,
                      std::vector<lps::action_summand>& summands_
                      )
    : algorithm(algorithm_), eq_X(eq_X_), process_parameters(process_parameters_), summands(summands_)
  {}

  void leave(const pbes_system::propositional_variable_instantiation& x)
  {
    data::data_expression_list g_i = algorithm.reset_variable_parameters(x, eq_X, i);
    if (i < summands.size())
    {
      summands[i].assignments() = make_assignments(process_parameters, g_i);
      i++;
    }
  }
};

inline
void lps_reset_variables(lpsstategraph_algorithm& algorithm,
                         const pbes_expression& x,
                         const stategraph_equation& eq_X,
                         const data::variable_list& process_parameters,
                         std::vector<lps::action_summand>& summands
                        )
{
  lps_reset_traverser f(algorithm, eq_X, process_parameters, summands);
  f.apply(x);
}

} // namespace mcrl2::pbes_system::detail

#endif // MCRL2_PBES_DETAIL_LPSSTATEGRAPH_ALGORITHM_H
