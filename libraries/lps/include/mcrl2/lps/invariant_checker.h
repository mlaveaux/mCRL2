// Author(s): Luc Engelen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/invariant_checker.h
/// \brief Add your file description here.

// Interface to class Invariant_Checker

#ifndef MCRL2_LPS_INVARIANT_CHECKER_H
#define MCRL2_LPS_INVARIANT_CHECKER_H

#include <limits>

#include "mcrl2/data/detail/prover/bdd_prover.h"
#include "mcrl2/data/detail/prover/bdd2dot.h"
#include "mcrl2/lps/stochastic_specification.h"

/// The class Invariant_Checker is initialized with an LPS using the constructor Invariant_Checker::Invariant_Checker.
/// After initialization, the function Invariant_Checker::check_invariant can be called any number of times to check
/// whether an expression of sort Bool passed as argument a_invariant is an invariant of this LPS. A new instance of the
/// class Invariant_Checker has to be created for each new LPS that has to be checked.
///
/// The class Invariant_Checker uses an instance of the class BDD_Prover to check whether a formula is a valid invariant
/// of an mCRL2 LPS. The constructor Invariant_Checker::Invariant_Checker initializes the BDD based prover with the
/// parameters a_rewrite_strategy, a_time_limit, a_path_eliminator, a_solver_type and the data specification of the LPS
/// passed as parameter a_lps. The parameter a_rewrite_strategy specifies which rewrite strategy is used by the prover's
/// rewriter. It can be set to either GS_REWR_JITTY or GS_REWR_JITTYC. The parameter
/// a_time_limit specifies the maximum amount of time in seconds to be spent by the prover on proving a single expression.
/// If a_time_limit is set to 0, no time limit will be enforced. The parameter a_path_eliminator specifies whether or not
/// path elimination is applied. When path elimination is applied, the prover uses an SMT solver to remove inconsistent
/// paths from BDDs. The parameter a_solver_type specifies which SMT solver is used for path elimination. Either the SMT
/// solver ario (http://www.eecs.umich.edu/~ario/) or cvc-lite (http://www.cs.nyu.edu/acsys/cvcl/) can be used. To use one
/// of these solvers, the directory containing the corresponding executable must be in the path. If the parameter
/// a_path_eliminator is set to false, the parameter a_solver_type is ignored. The parameter a_apply_induction indicates
/// whether or not induction on list will be applied.
///
/// The parameter a_dot_file_name specifies whether a file in dot format of the resulting BDD is saved each time the
/// prover cannot determine whether an expression is a contradiction or a tautology. If the parameter is set to 0, no .dot
/// files are saved. If a string is passed as parameter a_dot_file_name, this string will be used as the prefix of the
/// filenames. An instance of the class BDD2Dot is used to save these files in dot format.
///
/// If the parameter a_counter_example is set to true, a so called counter example will be printed to stderr each time a
/// summand is encountered that violates the invariant. A counter example is a valuation for which the expression to be
/// proven does not hold. If the parameter a_all_violations is set to true, the invariant checker will not stop as soon as
/// a violation of the invariant is found, but will report all violations instead.
///
/// Given an LPS,
///
///    P(d: D) = ...
///            + sum ei: Ei. ci(d, ei) -> ai(fi(d, ei)) . P(gi(d, ei))
///            + ...;
///
/// an instance of the class Invariant_Checker will generate a formula of the form
///
///    inv(d) /\ ci(d, ei) => inv(gi(d, ei))
///
/// for each of the summands, where inv() is the expression passed as parameter a_invariant. If this expression passed as
/// parameter a_invariant holds for the initial state and all the generated formulas are tautologies according to the
/// prover, it is an invariant.

namespace mcrl2::lps::detail
{

template <typename Specification>
class Invariant_Checker
{
  using process_type = typename Specification::process_type;
  using action_summand_type = typename process_type::action_summand_type;
  using action_summand_vector_type = std::vector<action_summand_type>;

private:
  const Specification& f_spec;
  data::detail::BDD_Prover f_bdd_prover;
  data::detail::BDD2Dot f_bdd2dot;
  process_initializer f_init;
  action_summand_vector_type f_summands;
  bool f_counter_example;
  bool f_all_violations;
  std::string f_dot_file_name;
  void print_counter_example();
  void save_dot_file(std::size_t a_summand_number);
  bool check_init(const data::data_expression& a_invariant);
  bool check_summand(const data::data_expression& a_invariant,
    const action_summand_type& a_summand,
    std::size_t a_summand_number);
  bool check_summands(const data::data_expression& a_invariant);
public:

  /// precondition: the argument passed as parameter a_lps is a valid mCRL2 LPS
  /// precondition: the argument passed as parameter a_time_limit is greater than or equal to 0. If the argument is
  /// equal to 0, no time limit will be enforced
  Invariant_Checker(const Specification& a_lps,
      data::rewriter::strategy a_rewrite_strategy = data::jitty,
      int a_time_limit = 0,
      bool a_path_eliminator = false,
      data::detail::smt_solver_type a_solver_type = data::detail::solver_type_cvc,
      bool a_apply_induction = false,
      bool a_counter_example = false,
      bool a_all_violations = false,
      const std::string& a_dot_file_name = std::string());

  /// precondition: the argument passed as parameter a_invariant is a valid expression in internal mCRL2 format
  bool check_invariant(const data::data_expression& a_invariant);
};

// Class Invariant_Checker ------------------------------------------------------------------------
// Class Invariant_Checker - Functions declared private -----------------------------------------

template <typename Specification>
void Invariant_Checker<Specification>::print_counter_example()
{
  if (f_counter_example)
  {
    data::data_expression v_counter_example(f_bdd_prover.get_counter_example());
    assert(v_counter_example.defined());
    mCRL2log(log::info) << "  Counter example: " << data::pp(v_counter_example) << "\n";
  }
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
void Invariant_Checker<Specification>::save_dot_file(std::size_t a_summand_number)
{
  if (!f_dot_file_name.empty())
  {
    std::string v_file_name=f_dot_file_name;

    if (a_summand_number == std::numeric_limits<std::size_t>::max()) // Sentinel for the initial state.
    {
      v_file_name += "-init.dot";
    }
    else
    {
      v_file_name +=  "-" + std::to_string(a_summand_number) + ".dot";
    }
    f_bdd2dot.output_bdd(f_bdd_prover.get_bdd(), v_file_name);
  }
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
bool Invariant_Checker<Specification>::check_init(const data::data_expression& a_invariant)
{
  data::mutable_map_substitution<> v_substitutions;
  const data::variable_list& d = f_spec.process().process_parameters();
  const data::data_expression_list& e = f_init.expressions();
  auto di = d.begin();
  auto ei = e.begin();
  for (; di != d.end(); ++di, ++ei)
  {
    v_substitutions[*di] = *ei;
  }

  data::data_expression b_invariant = data::replace_variables_capture_avoiding(a_invariant, v_substitutions);
  f_bdd_prover.set_formula(b_invariant);
  if (f_bdd_prover.is_tautology() == data::detail::answer_yes)
  {
    return true;
  }
  else
  {
    if (f_bdd_prover.is_contradiction() != data::detail::answer_yes)
    {
      print_counter_example();
      save_dot_file(std::numeric_limits<std::size_t>::max());
    }
    return false;
  }
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
bool Invariant_Checker<Specification>::check_summand(
  const data::data_expression& a_invariant,
  const action_summand_type& a_summand,
  const std::size_t a_summand_number)
{
  using namespace data::sort_bool;
  const data::data_expression& v_condition = a_summand.condition();

  data::mutable_map_substitution<> v_substitutions;

  for (const data::assignment& a: a_summand.assignments())
  {
    v_substitutions[a.lhs()] = a.rhs();
  }

  const data::data_expression v_subst_invariant = data::replace_variables_capture_avoiding(a_invariant, v_substitutions);

  const data::data_expression v_formula = implies(and_(a_invariant, v_condition), v_subst_invariant);
  f_bdd_prover.set_formula(v_formula);
  if (f_bdd_prover.is_tautology() == data::detail::answer_yes)
  {
    mCRL2log(log::verbose) << "The invariant holds for summand " << a_summand_number << "." << std::endl;
    return true;
  }
  else
  {
    mCRL2log(log::info) << "The invariant does not hold for summand " << a_summand_number << std::endl;
    if (f_bdd_prover.is_contradiction() != data::detail::answer_yes)
    {
      print_counter_example();
      save_dot_file(a_summand_number);
    }
    return false;
  }
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
bool Invariant_Checker<Specification>::check_summands(const data::data_expression& a_invariant)
{
  bool v_result = true;
  std::size_t v_summand_number = 1;

  for (auto i = f_summands.begin(); i != f_summands.end() && (f_all_violations || v_result); ++i)
  {
    v_result = check_summand(a_invariant, *i, v_summand_number) && v_result;
    v_summand_number++;
  }
  return v_result;
}

// Class Invariant_Checker<Specification> - Functions declared public --------------------------------------------

template <typename Specification>
Invariant_Checker<Specification>::Invariant_Checker(
  const Specification& a_lps,
  data::rewriter::strategy a_rewrite_strategy, int a_time_limit, bool a_path_eliminator, data::detail::smt_solver_type a_solver_type,
  bool a_apply_induction, bool a_counter_example, bool a_all_violations, std::string const& a_dot_file_name
):
  f_spec(a_lps),
  f_bdd_prover(a_lps.data(), data::used_data_equation_selector(a_lps.data()), a_rewrite_strategy, a_time_limit, a_path_eliminator, a_solver_type, a_apply_induction)
{
  f_init = a_lps.initial_process();
  f_summands = a_lps.process().action_summands();
  f_counter_example = a_counter_example;
  f_all_violations = a_all_violations;
  f_dot_file_name = a_dot_file_name;
}

// --------------------------------------------------------------------------------------------

template <typename Specification>
bool Invariant_Checker<Specification>::check_invariant(const data::data_expression& a_invariant)
{
  bool v_result = true;

  if (check_init(a_invariant))
  {
    mCRL2log(log::verbose) << "The invariant holds for the initial state." << std::endl;
  }
  else
  {
    mCRL2log(log::info) << "The invariant does not hold for the initial state." << std::endl;
    v_result = false;
  }
  if ((f_all_violations || v_result))
  {
    if (check_summands(a_invariant))
    {
      mCRL2log(log::verbose) << "The invariant holds for all summands." << std::endl;
    }
    else
    {
      mCRL2log(log::info) << "The invariant does not hold for all summands." << std::endl;
      v_result = false;
    }
  }
  if (v_result)
  {
    mCRL2log(log::info) << "The invariant holds for this LPS." << std::endl;
  }
  else
  {
    mCRL2log(log::info) << "The invariant does not hold for this LPS." << std::endl;
  }

  return v_result;
}

} // namespace mcrl2::lps::detail

#endif
