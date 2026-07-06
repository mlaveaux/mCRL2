// Author(s): Luc Engelen, Djurre van der Wal
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/confluence_checker.h
/// \brief Checks confluence of tau-summands in a linear process specification.

// Interface to class Confluence_Checker
// file: confluence_checker.h

#ifndef MCRL2_LPS_CONFLUENCE_CHECKER_H
#define MCRL2_LPS_CONFLUENCE_CHECKER_H

#include "mcrl2/lps/disjointness_checker.h"
#include "mcrl2/lps/invariant_checker.h"
#include <iomanip>


/** \brief A class that takes a linear process specification and checks all tau-summands of that LPS for confluence.
    \brief The tau actions of all confluent tau-summands are renamed to ctau.

     Given an LPS,

       P(d: D) = ...
               + sum ei: Ei. ci(d, ei) -> ai(fi(d, ei)) . P(gi(d, ei))
               + ...
               + sum ej: Ej. cj(d, ej) -> tau . P(gj(d, ej))
               + ...;


    tau-summand j is confluent with summand i if the following condition holds for all d: D, for all ei: Ei and for all
     ej: Ej:

       (inv(d) /\ ci(d, ei) /\ cj(d, ej))
       =>
       (
        ci(gj(d, ej), ei) /\
        cj(gi(d, ei), ej) /\
        fi(d, ei) == fi(gj(d, ej), ei) /\
        gi(gj(d, ej), ei) == gj(gi(d, ei), ej)
       )

    where inv() is the invariant specified using the parameter a_invariant of the function
    Confluence_Checker::check_confluence_and_mark. In case ai is also a tau-action, the formula above can be weakened to
    the following:

       (inv(d) /\ ci(d, ei) /\ cj(d, ej))
       =>
       (
        gi(d, ei) == gj(d, ej) \/
        (
         ci(gj(d, ej), ei) /\
         cj(gi(d, ei), ej) /\
         gi(gj(d, ej), ei) == gj(gi(d, ei), ej)
        )
       )

    The class Confluence_Checker can determine whether two summands are confluent in three ways and will indicate which
    of the methods was used while proving confluence. The three ways of determining confluence are as follows:

         If summand number 1 has been proven confluent with summand number 2, summand number 2 is obviously confluent
    with summand number 1. This method of checking confluence is called checking confluence by symmetry. If two summands
    are confluent by symmetry, the class Confluence_Checker indicates this by printing a dot ('.').

         Another way of checking the confluence of two summands is determining whether the two summands are
    syntactically disjoint. Two summands are syntactically disjoint if the following holds:
    - The set of variables used by one summand is disjoint from the set of variables changed by the other summand and
      vice versa.
    - The set of variables changed by one summand is disjoint from the set of variables changed by the other summand.
    If two summands are confluent because of syntactic disjointness, the class Confluence_Checker indicates this by
    printing a colon (':').

         The most time consuming way of checking the confluence of two summands is generating the confluence condition
    and then checking if this condition is a tautology using a prover for expressions of sort Bool. If two summands are
    proven confluent using the prover, the class Confluence_Checker indicates this by printing a plus sign ('+'). If the
    parameter a_generate_invariants is set to true, the class Confluence_Checker will try to prove that the reduced
    confluence condition is an invariant of the LPS, in case the confluence condition is not a tautology. If the reduced
    confluence condition is indeed an invariant, the two summands are proven confluent. The class Confluence_Checker
    indicates this by printing an 'i'.

    The class Confluence_Checker uses an instance of the class BDD_Prover, an instance of the class Disjointness_Checker
    and an instance of the class Invariant_Checker to determine which tau-summands of an mCRL2 LPS are confluent.
    Confluent tau-summands will be marked by renaming their tau-actions to ctau. The constructor
    Confluence_Checker::Confluence_Checker initializes the BDD based prover with the parameters a_rewrite_strategy,
    a_time_limit, a_path_eliminator, a_solver_type, a_apply_induction and a_lps. The parameter a_rewrite_strategy
    specifies which rewrite strategy is used by the prover's rewriter. It can be set to either
    GS_REWR_JITTY or GS_REWR_JITTYC. The parameter a_time_limit specifies the maximum amount of time in
    seconds to be spent by the prover on proving a single expression. If a_time_limit is set to 0, no time limit will be
    enforced. The parameter a_path_eliminator specifies whether or not path elimination is applied. When path
    elimination is applied, the prover uses an SMT solver to remove inconsistent paths from BDDs. The parameter
    a_solver_type specifies which SMT solver is used for path elimination. Either the SMT solver ario
    (http://www.eecs.umich.edu/~ario/) or cvc-lite (http://www.cs.nyu.edu/acsys/cvcl/) can be used. To use one of these
    solvers, the directory containing the corresponding executable must be in the path. If the parameter
    a_path_eliminator is set to false, the parameter a_solver_type is ignored. The parameter a_apply_induction indicates
    whether or not induction on list will be applied.

    The parameter a_dot_file_name specifies whether a file in dot format of the resulting BDD is saved each time the
    prover cannot determine whether an expression of sort Bool is a contradiction or a tautology. If the parameter is
    set to 0, no .dot files are saved. If a string is passed as parameter a_dot_file_name, this string will be used as
    the prefix of the filenames. An instance of the class BDD2Dot is used to save these files in dot format.

    If the parameter a_counter_example is set to true, a so called counter example is printed to stderr each time the
    prover indicates that two summands are not confluent. A counter example is a valuation for which the confluence
    condition to be proven does not hold.

    If the parameter a_check_all is set to true, the confluence of the tau-summands regarding all other summands will be
    checked. If the parameter is set to false, Confluence_Checker continues with the next tau-summand as soon as a
    summand is encountered that is not confluent with the current tau-summand.

    If the parameter a_generate_invariants is set, an invariant checker is used to check if the reduced confluence
    condition is an invariant of the LPS passed as parameter a_lps. If the reduced confluence condition is an invariant,
    the two summands are confluent.

    The function Confluence_Checker::check_confluence_and_mark returns an LPS with all tau-actions of confluent
    tau-summands renamed to ctau, unless the parameter a_no_marking is set to true. In case the parameter a_no_marking
    was set to true, the confluent tau-summands will not be marked, only the results of the confluence checking will be
    displayed.

    If there already is an action named ctau present in the LPS passed as parameter a_lps, an error will be reported. */

namespace mcrl2::lps::detail
{

/**
 * \brief Creates an identifier for the for the ctau action
 **/
inline process::action_label make_ctau_act_id()
{
  static atermpp::aterm ctau_act_id = atermpp::aterm(core::detail::function_symbol_ActId(), atermpp::aterm(atermpp::function_symbol("ctau", 0)), atermpp::aterm_list());

  assert(atermpp::detail::address(ctau_act_id));

  return process::action_label(ctau_act_id);
}

/**
 * \brief Creates the ctau action
 **/
inline process::action make_ctau_action()
{
  static atermpp::aterm ctau_action = atermpp::aterm(core::detail::function_symbol_Action(), make_ctau_act_id(), atermpp::aterm_list());

  assert(atermpp::detail::address(ctau_action));

  return process::action(ctau_action);
}

template <typename Specification>
class Confluence_Checker
{
  using process_type = typename Specification::process_type;
  using action_summand_type = typename process_type::action_summand_type;
  using action_summand_vector_type = std::vector<action_summand_type>;

private:
  /// \brief Class that can check if two summands are disjoint.
  Disjointness_Checker f_disjointness_checker;

  /// \brief Class that checks if an invariant holds for an LPS.
  Invariant_Checker<Specification> f_invariant_checker;

  /// \brief BDD based prover.
  data::detail::BDD_Prover f_bdd_prover;

  /// \brief Class that prints BDDs in dot format.
  data::detail::BDD2Dot f_bdd2dot;

  /// \brief A linear process specification.
  Specification& f_lps;

  /// \brief Flag indicating whether or not the process of checking the confluence of a summand stops when
  /// \brief a summand is encountered that is not confluent with the tau summand at hand.
  bool f_check_all;

  /// \brief Do not rewrite summands with sum operators.
  bool f_no_sums;

  /// \brief Confluence types for which the tool should check.
  std::string f_conditions;

  /// \brief Flag indicating whether or not counter examples are printed.
  bool f_counter_example;

  /// \brief The prefix for the names of the files written in dot format.
  std::string f_dot_file_name;

  /// \brief Flag indicating whether or not invariants are generated and checked each time a
  /// \brief summand is encountered that is not confluent with the tau summand at hand.
  bool f_generate_invariants;

  /// \brief The number of summands of the current LPS.
  std::size_t f_number_of_summands = 0UL;

  /// \brief An integer array, storing intermediate results per summand.
  std::vector<std::size_t> f_intermediate;

  /// \brief Identifier generator to allow variables to be uniquely renamed.
  data::set_identifier_generator f_set_identifier_generator;

  /// \brief Writes a dot file of the BDD created when checking the confluence of summands a_summand_number_1 and
  /// a_summand_number_2.
  void save_dot_file(std::size_t a_summand_number_1, std::size_t a_summand_number_2);

  /// \brief Outputs a path in the BDD corresponding to the condition at hand that leads to a node labelled false.
  void print_counter_example();

  /// \brief Checks the confluence of summand a_summand_1 and a_summand_2
  bool check_summands(const data::data_expression& a_invariant,
    action_summand_type a_summand_1,
    std::size_t a_summand_number_1,
    action_summand_type a_summand_2,
    std::size_t a_summand_number_2,
    char a_condition_type);

  /// \brief Checks and updates the confluence of summand a_summand concerning all other tau-summands.
  void check_confluence_and_mark_summand(action_summand_type& a_summand,
    std::size_t a_summand_number,
    const data::data_expression& a_invariant,
    char a_condition_type,
    bool& a_is_marked);

  // Returns a modified instance of a summand in which summation variables are uniquely renamed.
  void uniquely_rename_summation_variables(action_summand_type& summand);

public:
  /// \brief Constructor that initializes Confluence_Checker::f_lps, Confluence_Checker::f_bdd_prover,
  /// \brief Confluence_Checker::f_generate_invariants and Confluence_Checker::f_dot_file_name.
  /// precondition: the argument passed as parameter a_lps is a valid mCRL2 LPS
  /// precondition: the argument passed as parameter a_time_limit is greater than or equal to 0. If the argument is
  /// equal to 0, no time limit will be enforced
  Confluence_Checker(Specification& a_lps,
      data::rewriter::strategy a_rewrite_strategy = data::jitty,
      int a_time_limit = 0,
      bool a_path_eliminator = false,
      data::detail::smt_solver_type a_solver_type = data::detail::solver_type_cvc,
      bool a_apply_induction = false,
      bool a_check_all = false,
      bool a_no_sums = false,
      std::string a_conditions = "c",
      bool a_counter_example = false,
      bool a_generate_invariants = false,
      std::string const& a_dot_file_name = std::string());

  /// \brief Check the confluence of the LPS Confluence_Checker::f_lps.
  /// precondition: the argument passed as parameter a_invariant is an expression of sort Bool in internal mCRL2 format
  /// precondition: the argument passed as parameter a_summand_number corresponds with a summand of the LPS for which
  /// confluence must be checked (lowest summand has number 1). If this number is 0 confluence for all summands is
  /// checked.
  void check_confluence_and_mark(const data::data_expression& a_invariant, std::size_t a_summand_number);
};

// Auxiliary functions ----------------------------------------------------------------------------

static inline
data::mutable_map_substitution<> get_substitutions_from_assignments(const data::assignment_list& a_assignments)
{
  data::mutable_map_substitution<> v_substitutions;

  for (const data::assignment& a_assignment : a_assignments)
  {
    v_substitutions[a_assignment.lhs()]=a_assignment.rhs();
  }
  return v_substitutions;
}

// ----------------------------------------------------------------------------------------------
static inline
data::data_expression get_subst_equation_from_assignments(
  const data::variable_list& a_variables,
  data::assignment_list a_assignments_1,
  data::assignment_list a_assignments_2,
  data::mutable_map_substitution<>& a_substitutions_1,
  data::mutable_map_substitution<>& a_substitutions_2)
{
  data::data_expression v_result = data::sort_bool::true_();

  const data::assignment_list v_assignment_1;
  const data::assignment_list v_assignment_2;
  data::variable v_variable_1;
  data::variable v_variable_2;
  data::data_expression v_expression_1;
  data::data_expression v_expression_2;
  bool v_next_1 = true;
  bool v_next_2 = true;

  for (data::variable_list::const_iterator i=a_variables.begin(); i!=a_variables.end();)
  {
    data::variable v_variable = *i;
    ++i;

    if (!a_assignments_1.empty() && v_next_1)
    {
      const data::assignment v_assignment_1 = a_assignments_1.front();
      a_assignments_1.pop_front();
      v_variable_1 = v_assignment_1.lhs();
      v_expression_1 = v_assignment_1.rhs();
      v_expression_1 = data::replace_variables_capture_avoiding(v_expression_1, a_substitutions_2);
    }
    if (!a_assignments_2.empty() && v_next_2)
    {
      const data::assignment v_assignment_2 = a_assignments_2.front();
      a_assignments_2.pop_front();
      v_variable_2 = v_assignment_2.lhs();
      v_expression_2 = v_assignment_2.rhs();
      v_expression_2 = data::replace_variables_capture_avoiding(v_expression_2, a_substitutions_1);
    }
    while (v_variable != v_variable_1 && v_variable != v_variable_2 && i!=a_variables.end())
    {
      v_variable = *i;
      ++i;
    }
    if (v_variable_1 == v_variable_2)
    {
      v_result = data::sort_bool::and_(v_result, equal_to(v_expression_1, v_expression_2));
      v_next_1 = true;
      v_next_2 = true;
    }
    else if (v_variable == v_variable_1)
    {
      data::data_expression expr = data::replace_variables_capture_avoiding(data::data_expression(v_variable_1), a_substitutions_1);
      v_result = data::sort_bool::and_(data::data_expression(v_result), equal_to(v_expression_1, expr));
      v_next_1 = true;
      v_next_2 = false;
    }
    else if (v_variable == v_variable_2)
    {
      data::data_expression expr = data::replace_variables_capture_avoiding(data::data_expression(v_variable_2), a_substitutions_2);
      v_result = data::sort_bool::and_(data::data_expression(v_result), equal_to(data::data_expression(v_expression_2), expr));
      v_next_1 = false;
      v_next_2 = true;
    }
  }
  return v_result;
}

// ----------------------------------------------------------------------------------------------

static inline
data::data_expression get_equation_from_assignments(
  const data::variable_list& a_variables,
  data::assignment_list a_assignments_1,
  data::assignment_list a_assignments_2)
{
  data::data_expression v_result = data::sort_bool::true_();

  for (const data::variable& v_variable : a_variables)
  {
    if (!a_assignments_1.empty() && v_variable == a_assignments_1.front().lhs())
    {
      if (!a_assignments_2.empty() && v_variable == a_assignments_2.front().lhs())
      {
        // Create a condition from the assigments from both lists.
        v_result = data::sort_bool::and_(v_result, equal_to(a_assignments_1.front().rhs(), a_assignments_2.front().rhs()));
        a_assignments_2.pop_front();
      }
      else
      {
        // Create a condition from first assigment only.
        v_result = data::sort_bool::and_(v_result, equal_to(a_assignments_1.front().rhs(), v_variable));
      }
      a_assignments_1.pop_front();
    }
    else
    {
      if (!a_assignments_2.empty() && v_variable == a_assignments_2.front().lhs())
      {
        // Create a condition from the second assigments only.
        v_result = data::sort_bool::and_(v_result, equal_to(v_variable, a_assignments_2.front().rhs()));
        a_assignments_2.pop_front();
      }
    }
  }

  assert(a_assignments_1.empty()); // If this is not the case, the assignments do not have the
  assert(a_assignments_2.empty()); // same order as the list of variables. This means that some equations
                                   // have not been generated.
  return v_result;
}

// ----------------------------------------------------------------------------------------------
inline
data::data_expression get_subst_equation_from_actions(
  const process::action_list& a_actions,
  data::mutable_map_substitution<>& a_substitutions)
{
  data::data_expression v_result = data::sort_bool::true_();

  for (const process::action& a_action: a_actions)
  {
    const data::data_expression_list v_expressions = a_action.arguments();
    for (const data::data_expression& v_expression : v_expressions)
    {
      const data::data_expression v_subst_expression = data::replace_variables_capture_avoiding(v_expression, a_substitutions);
      v_result = data::sort_bool::and_(data::data_expression(v_result), equal_to(v_expression, v_subst_expression));
    }
  }
  return v_result;
}

// ----------------------------------------------------------------------------------------------

static inline data::assignment_list get_full_assignment_list(
  data::assignment_list a_assignment_list,
  const data::variable_list& a_variables)
{
  data::assignment_list v_new_list;

  for (const data::variable& v_variable: a_variables)
  {
    bool v_assignment_added = false;

    if (!a_assignment_list.empty())
    {
      const data::assignment v_assignment = a_assignment_list.front();
      a_assignment_list.pop_front();

      if (v_assignment.lhs() == v_variable)
      {
        v_new_list.push_front(v_assignment);
        v_assignment_added = true;
      }
    }

    if (!v_assignment_added)
    {
      v_new_list.push_front(data::assignment(v_variable, v_variable));
    }
  }

  return atermpp::reverse(v_new_list);
}

// ----------------------------------------------------------------------------------------------

template <typename ActionSummand>
data::data_expression get_confluence_condition(
  const data::data_expression& a_invariant,
  const ActionSummand& a_summand_1,
  const ActionSummand& a_summand_2,
  const data::variable_list& a_variables,
  const char condition_type)
{
  assert(a_summand_1.is_tau());

  if (condition_type == 'c' || condition_type == 'C') //Commutative confluence.
  {
    const data::data_expression& v_condition_1 = a_summand_1.condition();
    const data::assignment_list& v_assignments_1 = a_summand_1.assignments();

    data::mutable_map_substitution<> v_substitutions_1 = get_substitutions_from_assignments(v_assignments_1);
    const data::data_expression& v_condition_2 = a_summand_2.condition();
    const data::data_expression v_lhs = data::sort_bool::and_(data::sort_bool::and_(v_condition_1, v_condition_2), a_invariant);
    const data::assignment_list& v_assignments_2 = a_summand_2.assignments();

    data::mutable_map_substitution<> v_substitutions_2 = get_substitutions_from_assignments(v_assignments_2);
    const data::data_expression v_subst_condition_1 = data::replace_variables_capture_avoiding(v_condition_1, v_substitutions_2);
    const data::data_expression v_subst_condition_2 = data::replace_variables_capture_avoiding(v_condition_2, v_substitutions_1);

    const data::data_expression v_subst_equation = get_subst_equation_from_assignments(a_variables, v_assignments_1, v_assignments_2, v_substitutions_1, v_substitutions_2);

    const process::action_list v_actions = a_summand_2.multi_action().actions();
    data::data_expression v_rhs;

    if (v_actions.empty())
    {
      // tau-summand
      const data::data_expression v_equation = get_equation_from_assignments(a_variables, v_assignments_1, v_assignments_2);
      v_rhs = data::sort_bool::and_(v_subst_condition_1, v_subst_condition_2);
      v_rhs = data::sort_bool::and_(v_rhs, v_subst_equation);
      v_rhs = data::sort_bool::or_(v_equation, v_rhs);
    }
    else
    {
      // non-tau-summand
      const data::data_expression v_actions_equation = get_subst_equation_from_actions(v_actions, v_substitutions_1);
      v_rhs = data::sort_bool::and_(v_subst_condition_1, v_subst_condition_2);
      v_rhs = data::sort_bool::and_(v_rhs, v_actions_equation);
      v_rhs = data::sort_bool::and_(v_rhs, data::data_expression(v_subst_equation));
    }
    return data::sort_bool::implies(v_lhs, v_rhs);
  }
  else if (condition_type == 'T') //Triangular confluence
  {
    const data::data_expression& v_condition_1 = a_summand_1.condition();
    const data::data_expression& v_condition_2 = a_summand_2.condition();

    const data::assignment_list& v_assignments_1 = a_summand_1.assignments();
    const data::assignment_list v_assignments_2 = get_full_assignment_list(a_summand_2.assignments(), a_variables);

    const process::action_list v_actions = a_summand_2.multi_action().actions();

    data::mutable_map_substitution<> v_substitutions_1 = get_substitutions_from_assignments(v_assignments_1);
    data::mutable_map_substitution<> v_substitutions_2;

    const data::data_expression v_subst_condition_2 = data::replace_variables_capture_avoiding(v_condition_2, v_substitutions_1);
    const data::data_expression v_subst_equation = get_subst_equation_from_assignments(a_variables, v_assignments_2, v_assignments_2, v_substitutions_1, v_substitutions_2);

    data::data_expression v_lhs = data::sort_bool::and_(data::sort_bool::and_(v_condition_1, v_condition_2), a_invariant);
    data::data_expression v_rhs = data::sort_bool::and_(v_subst_condition_2, v_subst_equation);

    if (v_actions.empty())
    {
      // tau-summand
      const data::data_expression v_equation = get_equation_from_assignments(a_variables, v_assignments_1, v_assignments_2);
      v_rhs = data::sort_bool::or_(v_rhs, v_equation);
    }
    else
    {
      // non-tau-summand
      const data::data_expression v_actions_equation = get_subst_equation_from_actions(v_actions, v_substitutions_1);
      v_rhs = data::sort_bool::and_(v_rhs, v_actions_equation);
    }
    return data::sort_bool::implies(v_lhs, v_rhs);
  }
  else if (condition_type == 'Z') //Trivial confluence
  {
    const data::data_expression& v_condition_1 = a_summand_1.condition();
    const data::data_expression& v_condition_2 = a_summand_2.condition();

    const data::assignment_list& v_assignments_1 = a_summand_1.assignments();
    const data::assignment_list& v_assignments_2 = a_summand_2.assignments();

    data::mutable_map_substitution<> v_substitutions_1;
    data::mutable_map_substitution<> v_substitutions_2;

    const process::action_list v_actions = a_summand_2.multi_action().actions();

    const data::data_expression& v_lhs = v_condition_1;
    data::data_expression v_rhs;

    if (v_actions.empty())
    {
      // tau-summand
      data::data_expression v_subst_equation = get_subst_equation_from_assignments(a_variables, v_assignments_1, v_assignments_2, v_substitutions_1, v_substitutions_2);
      v_rhs = data::sort_bool::or_(data::sort_bool::not_(v_condition_2), v_subst_equation);
    }
    else
    {
      // non-tau-summand
      v_rhs = data::sort_bool::not_(v_condition_2);
    }
    return data::sort_bool::implies(v_lhs, v_rhs);
  }

  return data::sort_bool::false_();
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
bool has_ctau_action(const Specification& a_lps)
{
  const process::action_label_list& v_action_specification = a_lps.action_labels();
  return std::find(v_action_specification.begin(),v_action_specification.end(),make_ctau_act_id())!=v_action_specification.end();
}


// Class Confluence_Checker -----------------------------------------------------------------------
// Class Confluence_Checker - Functions declared private ----------------------------------------

template <typename Specification>
void Confluence_Checker<Specification>::save_dot_file(std::size_t a_summand_number_1, std::size_t a_summand_number_2)
{
  if (!f_dot_file_name.empty())
  {
    f_bdd2dot.output_bdd(f_bdd_prover.get_bdd(), f_dot_file_name + "-" + std::to_string(a_summand_number_1) + "-" + std::to_string(a_summand_number_2) + ".dot");
  }
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
void Confluence_Checker<Specification>::print_counter_example()
{
  if (f_counter_example)
  {
    const data::data_expression v_counter_example(f_bdd_prover.get_counter_example());
    mCRL2log(log::info) << "  Counter example: " << v_counter_example << "\n";
  }
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
void Confluence_Checker<Specification>::uniquely_rename_summation_variables(
  action_summand_type& summand)
{
  data::mutable_map_substitution<> v_substitutions;
  std::set<data::variable> v_substitution_variables;

  data::variable_list summation_variables = summand.summation_variables();
  data::assignment_list assignments = summand.assignments();
  data::data_expression condition = summand.condition();
  process::action_list actions = summand.multi_action().actions();

  //Compose the substitution map.
  //Simultaneously, create a list of the new summation variables.
  data::variable_list new_summation_variables = data::variable_list();

  for (const data::variable& summation_variable : summation_variables)
  {
    core::identifier_string new_name = f_set_identifier_generator(summation_variable.name());

    data::variable renamed_variable = data::variable(new_name, summation_variable.sort());
    new_summation_variables.push_front(renamed_variable);
    v_substitutions[summation_variable] = renamed_variable;
  }

  new_summation_variables = atermpp::reverse(new_summation_variables);
  v_substitution_variables = data::substitution_variables(v_substitutions);

  //Create the new condition.
  data::data_expression new_condition = data::replace_variables(summand.condition(), v_substitutions);

  //Create the new (multi-)action.
  process::action_list new_actions = process::action_list();

  for (const process::action& i : actions)
  {
    data::data_expression_list new_arguments = data::data_expression_list();

    for (const data::data_expression& argument : i.arguments())
    {
      new_arguments.push_front(data::replace_variables(argument, v_substitutions));
    }

    new_actions.push_front(process::action(i.label(), atermpp::reverse(new_arguments)));
  }

  multi_action new_multi_action = multi_action(atermpp::reverse(new_actions));

  //Create the new assignments (essentially the summand's next-state function).
  data::assignment_list new_assignments = data::assignment_list();

  for (const data::assignment& i: assignments)
  {
    data::data_expression rhs = data::replace_variables(i.rhs(), v_substitutions);
    new_assignments.push_front(data::assignment(i.lhs(), rhs));
  }

  new_assignments = atermpp::reverse(new_assignments);

  //Create a new summand from the components.
  summand.summation_variables() = new_summation_variables;
  summand.condition() = new_condition;
  summand.assignments() = new_assignments;
  summand.multi_action() = new_multi_action;
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
bool Confluence_Checker<Specification>::check_summands(
  const data::data_expression& a_invariant,
  const action_summand_type a_summand_1,
  const std::size_t a_summand_number_1,
  const action_summand_type a_summand_2,
  const std::size_t a_summand_number_2,
  const char a_condition_type)
{
  assert(a_summand_1.is_tau());

  const data::variable_list v_variables = f_lps.process().process_parameters();
  bool v_is_confluent = true;

  if ((a_condition_type == 'c' || a_condition_type == 'd') && f_disjointness_checker.disjoint(a_summand_number_1, a_summand_number_2))
  {
    mCRL2log(log::info) << ":";
  }
  else
  {
    action_summand_type tagged = a_summand_2;

    if (!f_no_sums)
    {
      uniquely_rename_summation_variables(tagged);
    }

    const data::data_expression v_condition = get_confluence_condition(a_invariant, a_summand_1, tagged, v_variables, a_condition_type);
    f_bdd_prover.set_formula(v_condition);
    if (f_bdd_prover.is_tautology() == data::detail::answer_yes)
    {
      mCRL2log(log::info) << "+";
    }
    else
    {
      if (f_generate_invariants)
      {
        const data::data_expression v_new_invariant(f_bdd_prover.get_bdd());
        mCRL2log(log::verbose) << "\nChecking invariant: " << data::pp(v_new_invariant) << "\n";
        if (f_invariant_checker.check_invariant(v_new_invariant))
        {
          mCRL2log(log::verbose) << "Invariant holds" << std::endl;
          mCRL2log(log::info) << "i";
        }
        else
        {
          mCRL2log(log::verbose) << "Invariant doesn't hold" << std::endl;
          v_is_confluent = false;
          if (f_check_all)
          {
            mCRL2log(log::info) << "-";
          }
          else
          {
            mCRL2log(log::info) << "Not confluent with summand " << a_summand_number_2 << ".";
          }
          print_counter_example();
          save_dot_file(a_summand_number_1, a_summand_number_2);
        }
      }
      else
      {
        v_is_confluent = false;
        if (f_check_all)
        {
          mCRL2log(log::info) << "-";
        }
        else
        {
          mCRL2log(log::info) << "Not confluent with summand " << a_summand_number_2 << ".";
        }
        print_counter_example();
        save_dot_file(a_summand_number_1, a_summand_number_2);
      }
    }
  }
  return v_is_confluent;
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
void Confluence_Checker<Specification>::check_confluence_and_mark_summand(
  action_summand_type& a_summand,
  const std::size_t a_summand_number,
  const data::data_expression& a_invariant,
  const char a_condition_type,
  bool& a_is_marked)
{
  using action_summand_type = typename Specification::process_type::action_summand_type;
  assert(a_summand.is_tau());
  std::vector<action_summand_type>& v_summands = f_lps.process().action_summands();
  std::size_t v_summand_number = 1;
  bool v_is_confluent = true;

  // Add here that the sum variables of a_summand must be empty otherwise
  // the confluence of the summand must be checked with respect to itself,
  // which requires quantification. Otherwise tau.a+tau.b will be designated
  // tau-confluent, if linearised with summand clustering.

  if (f_no_sums)
  {
    const data::variable_list a_summand_sum_variables=a_summand.summation_variables();
    if (!a_summand_sum_variables.empty())
    {
      v_is_confluent = false;
      mCRL2log(log::info) << "Summand " << a_summand_number << " is not proven confluent because it contains a sum operator.";
    }
  }

  for (typename std::vector<action_summand_type>::const_iterator i=v_summands.begin(); i!=v_summands.end() && (v_is_confluent || f_check_all); ++i)
  {
    const action_summand_type v_summand = *i;

    if (v_summand_number < a_summand_number)
    {
      // Check the cache
      if (f_intermediate[v_summand_number] > a_summand_number)
      {
        mCRL2log(log::info) << ".";
      }
      else
      {
        if (f_intermediate[v_summand_number] == a_summand_number)
        {
          if (f_check_all)
          {
            mCRL2log(log::info) << "-";
          }
          else
          {
            mCRL2log(log::info) << "Not confluent with summand " << v_summand_number << ".";
          }
          v_is_confluent = false;
        }
        else
        {
          v_is_confluent &= check_summands(a_invariant, a_summand, a_summand_number, v_summand, v_summand_number, a_condition_type);
        }
      }
    }
    else
    {
      v_is_confluent &= check_summands(a_invariant, a_summand, a_summand_number, v_summand, v_summand_number, a_condition_type);
    }
    if (v_is_confluent || f_check_all)
    {
      // Only increase number if we will continue
      v_summand_number++;
    }
  }

  if (!f_check_all)
  {
    f_intermediate[a_summand_number] = v_summand_number;
  }

  if (v_is_confluent)
  {
    mCRL2log(log::info) << "Confluent with all summands.";
    a_is_marked = true;
    a_summand.multi_action() = multi_action(make_ctau_action());
  }
}

// Class Confluence_Checker - Functions declared public -----------------------------------------

template <typename Specification>
Confluence_Checker<Specification>::Confluence_Checker(
  Specification& a_lps,
  data::rewriter::strategy a_rewrite_strategy,
  int a_time_limit,
  bool a_path_eliminator,
  data::detail::smt_solver_type a_solver_type,
  bool a_apply_induction,
  bool a_check_all,
  bool a_no_sums,
  std::string a_conditions,
  bool a_counter_example,
  bool a_generate_invariants,
  std::string const& a_dot_file_name):
  f_disjointness_checker(a_lps.process()),
  f_invariant_checker(a_lps, a_rewrite_strategy, a_time_limit, a_path_eliminator, a_solver_type, false, false, 0),
  f_bdd_prover(a_lps.data(), data::used_data_equation_selector(a_lps.data()), a_rewrite_strategy,
                     a_time_limit, a_path_eliminator, a_solver_type, a_apply_induction),
  f_lps(a_lps),
  f_check_all(a_check_all),
  f_no_sums(a_no_sums),
  f_conditions(a_conditions),
  f_counter_example(a_counter_example),
  f_dot_file_name(a_dot_file_name),
  f_generate_invariants(a_generate_invariants)
{
  if (has_ctau_action(a_lps))
  {
    throw mcrl2::runtime_error("An action named \'ctau\' already exists.\n");
  }

  std::string v_conditions = std::string(f_conditions);

  while (v_conditions.length() > 0)
  {
    char v_ct = v_conditions.at(0);
    v_conditions = v_conditions.substr(1);

    if (v_ct != 'c' && v_ct != 'd' && v_ct != 'C' && v_ct != 'T' && v_ct != 'Z')
    {
      std::string msg = std::string("Unknown confluence type: ");
      msg += v_ct;
      msg += '.';

      throw mcrl2::runtime_error(msg);
    }
  }
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
void Confluence_Checker<Specification>::check_confluence_and_mark(const data::data_expression& a_invariant, const std::size_t a_summand_number)
{
  using action_summand_type = typename Specification::process_type::action_summand_type;
  typename Specification::process_type& v_process_equation = f_lps.process();
  std::vector<action_summand_type>& v_summands = v_process_equation.action_summands();
  bool v_is_marked = false;

  std::size_t v_summand_number = 1;
  std::set<std::size_t> v_unmarked_summands;
  std::set<std::size_t> v_marked_summands;

  for (const action_summand_type& s: v_summands)
  {
    if ((a_summand_number == v_summand_number) || (a_summand_number == 0))
    {
      if (s.is_tau())
      {
        v_unmarked_summands.insert(v_summand_number);
      }
    }
    v_summand_number++;
  }

  f_set_identifier_generator.clear_context();
  f_set_identifier_generator.add_identifiers(find_identifiers(f_lps));
  f_set_identifier_generator.add_identifiers(data::function_and_mapping_identifiers(f_lps.data()));

  f_number_of_summands = v_summands.size();
  std::string v_conditions = std::string(f_conditions);

  while (v_conditions.length() > 0)
  {
    f_intermediate = std::vector<std::size_t>(f_number_of_summands + 2, 0);
    char v_condition_type = v_conditions.at(0);
    v_conditions = v_conditions.substr(1);
    v_summand_number = 1;

    for (action_summand_type& s: v_summands)
    {
      std::set<std::size_t>::iterator it = v_unmarked_summands.find(v_summand_number);

      if (it != v_unmarked_summands.end())
      {
        bool summand_is_marked = false;

        mCRL2log(log::info) << "summand " << std::setw(3) << v_summand_number << " of " << v_summands.size() << " (condition = " << v_condition_type << "): ";
        check_confluence_and_mark_summand(s, v_summand_number, a_invariant, v_condition_type, summand_is_marked);
        mCRL2log(log::info) << std::endl;

        if (summand_is_marked)
        {
          v_marked_summands.insert(v_summand_number);
          v_unmarked_summands.erase(it);
          v_is_marked = true;
        }
      }
      v_summand_number++;
    }
  }

  if (v_is_marked && !has_ctau_action(f_lps))
  {
    f_lps.action_labels().push_front(make_ctau_act_id());
  }

  mCRL2log(log::info) << v_marked_summands.size() << " of " << (v_marked_summands.size() + v_unmarked_summands.size()) <<
                         " tau summands were found to be confluent" << std::endl;

  f_intermediate = std::vector<std::size_t>();
}

} // namespace mcrl2::lps::detail

#endif
