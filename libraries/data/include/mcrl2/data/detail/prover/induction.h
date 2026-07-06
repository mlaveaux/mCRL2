// Author(s): Luc Engelen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/detail/prover/induction.h
/// \brief Proving with induction on lists

#ifndef MCRL2_DATA_DETAIL_PROVER_INDUCTION_H
#define MCRL2_DATA_DETAIL_PROVER_INDUCTION_H

#include "mcrl2/data/list.h"
#include "mcrl2/data/replace.h"
#include "mcrl2/data/representative_generator.h"
#include "mcrl2/data/substitutions/mutable_map_substitution.h"

namespace mcrl2::data::detail
{
/// The class Induction generates statements corresponding to

class Induction
{
  private:
    set_identifier_generator fresh_identifier_generator;

    /// \brief The number of variables used during the last application of induction.
    std::size_t f_count = 0UL;

    /// \brief An expression of sort Bool in mCRL2 format.
    data_expression f_formula;

    /// \brief
    variable_vector f_list_variables;

    /// \brief
    variable_vector recurse_expression_for_lists(const data_expression& a_expression) const
    {
      variable_vector result;
      for(const variable& v: find_free_variables(a_expression))
      {
        if(sort_list::is_list(v.sort()))
        {
          result.push_back(v);
        }
      }
      return result;
    }

    /// \brief
    sort_expression get_sort_of_list_elements(const variable& a_list_variable) const
    {
      assert(sort_list::is_list(a_list_variable.sort()));
      return atermpp::down_cast<container_sort>(a_list_variable.sort()).element_sort();
    }

    /// \brief
    variable get_fresh_dummy(const sort_expression& a_sort)
    {
      return variable(fresh_identifier_generator("dummy$"),a_sort);
    }


    //TODO check if this special case can be integrated into create_clauses
    data_expression apply_induction_one(const core::identifier_string& fresh_name) const
    {
      const variable v_induction_variable = f_list_variables.front();
      assert(sort_list::is_list(v_induction_variable.sort()));

      const sort_expression v_dummy_sort = get_sort_of_list_elements(v_induction_variable);
      const variable v_dummy_variable(fresh_name, v_dummy_sort);

      mutable_map_substitution<> v_substitution1;
      v_substitution1[v_induction_variable] = sort_list::empty(v_dummy_sort);
      const data_expression v_base_case = replace_variables_capture_avoiding(f_formula, v_substitution1);

      mutable_map_substitution<> v_substitution2;
      v_substitution2[v_induction_variable] = sort_list::cons_(v_dummy_sort, v_dummy_variable, v_induction_variable);
      const data_expression v_induction_step = sort_bool::implies(f_formula, replace_variables_capture_avoiding(f_formula, v_substitution2));

      return sort_bool::and_(v_base_case, v_induction_step);
    }

    /// \brief
    data_expression create_hypotheses(
                        const data_expression& a_hypothesis,
                        variable_list a_list_of_variables,
                        variable_list a_list_of_dummies) const
    {
      if (a_list_of_variables.empty())
      {
        return sort_bool::true_();
      }
      else
      {
        data_expression v_clause = a_hypothesis;
        if (a_list_of_variables.size() > 1)
        {
          for (const variable& v_variable: a_list_of_variables)
          {
            const variable v_dummy(a_list_of_dummies.front());
            a_list_of_dummies.pop_front();

            mutable_map_substitution<> v_substitution;
            v_substitution[v_variable] = sort_list::cons_(v_dummy.sort(), v_dummy, v_variable);
            v_clause = sort_bool::and_(v_clause, replace_variables_capture_avoiding(a_hypothesis, v_substitution));
          }
        }

        return v_clause;
      }
    }

    /// \brief
    data_expression_list create_clauses(const data_expression& a_formula,
                             const data_expression& a_hypothesis,
                             const std::size_t a_variable_number,
                             const std::size_t a_number_of_variables,
                             const variable_list& a_list_of_variables,
                             const variable_list& a_list_of_dummies)
    {
      const variable v_variable = f_list_variables[a_variable_number];
      variable_list v_list_of_variables = a_list_of_variables;
      v_list_of_variables.push_front(v_variable);
      const sort_expression v_dummy_sort = get_sort_of_list_elements(v_variable);
      const variable v_dummy = get_fresh_dummy(v_dummy_sort);
      variable_list v_list_of_dummies = a_list_of_dummies;
      v_list_of_dummies.push_front(v_dummy);

      mutable_map_substitution<> v_substitution1;
      v_substitution1[v_variable] = sort_list::cons_(v_dummy_sort, v_dummy, v_variable);
      const data_expression v_formula_1 = replace_variables_capture_avoiding(a_formula, v_substitution1);

      mutable_map_substitution<> v_substitution2;
      assert(sort_list::is_list(v_variable.sort()));
      v_substitution2[v_variable] = sort_list::empty(v_dummy_sort);
      const data_expression v_formula_2 = replace_variables_capture_avoiding(a_formula, v_substitution2);
      const data_expression v_hypothesis = replace_variables_capture_avoiding(a_hypothesis, v_substitution2);

      if (a_variable_number < a_number_of_variables - 1)
      {
        const data_expression_list v_list_1 = create_clauses(v_formula_1, a_hypothesis, a_variable_number + 1, a_number_of_variables, v_list_of_variables, v_list_of_dummies);
        const data_expression_list v_list_2 = create_clauses(v_formula_2, v_hypothesis, a_variable_number + 1, a_number_of_variables, a_list_of_variables, a_list_of_dummies);
        return v_list_1 + v_list_2;
      }
      else
      {
        data_expression v_hypotheses_1 = create_hypotheses(a_hypothesis, v_list_of_variables, v_list_of_dummies);
        data_expression v_hypotheses_2 = create_hypotheses(v_hypothesis, a_list_of_variables, a_list_of_dummies);
        data_expression_list result;
        result.push_front(sort_bool::implies(v_hypotheses_2, v_formula_2));
        result.push_front(sort_bool::implies(v_hypotheses_1, v_formula_1));
        return result;
      }
    }

  public:
    void initialize(const data_expression& a_formula)
    {
      f_formula = a_formula;
      f_list_variables = recurse_expression_for_lists(a_formula);
      f_count = 0;
    }

    bool can_apply_induction() const
    {
      return f_list_variables.size() != f_count;
    }

    /// \pre can_apply_induction()
    data_expression apply_induction()
    {
      assert(can_apply_induction());
      data_expression v_result;

      f_count++;
      if (f_count == 1)
      {
        mCRL2log(log::verbose) << "Induction on one variable." << std::endl;
        v_result = apply_induction_one(fresh_identifier_generator("dummy$"));
      }
      else
      {
        mCRL2log(log::verbose) << "Induction on " << f_count << " variables." << std::endl;
        data_expression_list v_list_of_clauses = create_clauses(f_formula, f_formula, 0, f_count, variable_list(), variable_list());
        v_result = v_list_of_clauses.front();
        v_list_of_clauses.pop_front();
        while (!v_list_of_clauses.empty())
        {
          data_expression v_clause(v_list_of_clauses.front());
          v_list_of_clauses.pop_front();
          v_result = sort_bool::and_(v_result, v_clause);
        }
      }
      return v_result;
    }
};

} // namespace mcrl2::data::detail

#endif
