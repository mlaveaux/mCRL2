// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file modal_formula.cpp
/// \brief

#include "mcrl2/modal_formula/algorithms.h"
#include "mcrl2/modal_formula/is_timed.h"
#include "mcrl2/modal_formula/normalize.h"
#include "mcrl2/modal_formula/parse.h"
#include "mcrl2/modal_formula/parse_impl.h"
#include "mcrl2/modal_formula/print.h"
#include "mcrl2/modal_formula/replace.h"

namespace mcrl2
{

namespace action_formulas
{

//--- start generated action_formulas overloads ---//
std::string pp(const action_formulas::action_formula& x, bool arg0) { return action_formulas::pp< action_formulas::action_formula >(x, arg0); }
std::string pp(const action_formulas::and_& x, bool arg0) { return action_formulas::pp< action_formulas::and_ >(x, arg0); }
std::string pp(const action_formulas::at& x, bool arg0) { return action_formulas::pp< action_formulas::at >(x, arg0); }
std::string pp(const action_formulas::exists& x, bool arg0) { return action_formulas::pp< action_formulas::exists >(x, arg0); }
std::string pp(const action_formulas::false_& x, bool arg0) { return action_formulas::pp< action_formulas::false_ >(x, arg0); }
std::string pp(const action_formulas::forall& x, bool arg0) { return action_formulas::pp< action_formulas::forall >(x, arg0); }
std::string pp(const action_formulas::imp& x, bool arg0) { return action_formulas::pp< action_formulas::imp >(x, arg0); }
std::string pp(const action_formulas::multi_action& x, bool arg0) { return action_formulas::pp< action_formulas::multi_action >(x, arg0); }
std::string pp(const action_formulas::not_& x, bool arg0) { return action_formulas::pp< action_formulas::not_ >(x, arg0); }
std::string pp(const action_formulas::or_& x, bool arg0) { return action_formulas::pp< action_formulas::or_ >(x, arg0); }
std::string pp(const action_formulas::true_& x, bool arg0) { return action_formulas::pp< action_formulas::true_ >(x, arg0); }
std::set<data::variable> find_all_variables(const action_formulas::action_formula& x) { return action_formulas::find_all_variables< action_formulas::action_formula >(x); }
//--- end generated action_formulas overloads ---//

namespace detail {

action_formula parse_action_formula(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("ActFrm");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  action_formula result = action_formula_actions(p).parse_ActFrm(node);
  return result;
}

} // namespace detail

} // namespace action_formulas

namespace regular_formulas
{

//--- start generated regular_formulas overloads ---//
std::string pp(const regular_formulas::alt& x, bool arg0) { return regular_formulas::pp< regular_formulas::alt >(x, arg0); }
std::string pp(const regular_formulas::regular_formula& x, bool arg0) { return regular_formulas::pp< regular_formulas::regular_formula >(x, arg0); }
std::string pp(const regular_formulas::seq& x, bool arg0) { return regular_formulas::pp< regular_formulas::seq >(x, arg0); }
std::string pp(const regular_formulas::trans& x, bool arg0) { return regular_formulas::pp< regular_formulas::trans >(x, arg0); }
std::string pp(const regular_formulas::trans_or_nil& x, bool arg0) { return regular_formulas::pp< regular_formulas::trans_or_nil >(x, arg0); }
std::string pp(const regular_formulas::untyped_regular_formula& x, bool arg0) { return regular_formulas::pp< regular_formulas::untyped_regular_formula >(x, arg0); }
//--- end generated regular_formulas overloads ---//

namespace detail
{

regular_formula parse_regular_formula(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("RegFrm");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  regular_formula result = regular_formula_actions(p).parse_RegFrm(node);
  return result;
}

} // namespace detail

} // namespace regular_formulas

namespace state_formulas
{

//--- start generated state_formulas overloads ---//
std::string pp(const state_formulas::and_& x, bool arg0) { return state_formulas::pp< state_formulas::and_ >(x, arg0); }
std::string pp(const state_formulas::const_multiply& x, bool arg0) { return state_formulas::pp< state_formulas::const_multiply >(x, arg0); }
std::string pp(const state_formulas::const_multiply_alt& x, bool arg0) { return state_formulas::pp< state_formulas::const_multiply_alt >(x, arg0); }
std::string pp(const state_formulas::delay& x, bool arg0) { return state_formulas::pp< state_formulas::delay >(x, arg0); }
std::string pp(const state_formulas::delay_timed& x, bool arg0) { return state_formulas::pp< state_formulas::delay_timed >(x, arg0); }
std::string pp(const state_formulas::exists& x, bool arg0) { return state_formulas::pp< state_formulas::exists >(x, arg0); }
std::string pp(const state_formulas::false_& x, bool arg0) { return state_formulas::pp< state_formulas::false_ >(x, arg0); }
std::string pp(const state_formulas::forall& x, bool arg0) { return state_formulas::pp< state_formulas::forall >(x, arg0); }
std::string pp(const state_formulas::imp& x, bool arg0) { return state_formulas::pp< state_formulas::imp >(x, arg0); }
std::string pp(const state_formulas::infimum& x, bool arg0) { return state_formulas::pp< state_formulas::infimum >(x, arg0); }
std::string pp(const state_formulas::may& x, bool arg0) { return state_formulas::pp< state_formulas::may >(x, arg0); }
std::string pp(const state_formulas::minus& x, bool arg0) { return state_formulas::pp< state_formulas::minus >(x, arg0); }
std::string pp(const state_formulas::mu& x, bool arg0) { return state_formulas::pp< state_formulas::mu >(x, arg0); }
std::string pp(const state_formulas::must& x, bool arg0) { return state_formulas::pp< state_formulas::must >(x, arg0); }
std::string pp(const state_formulas::not_& x, bool arg0) { return state_formulas::pp< state_formulas::not_ >(x, arg0); }
std::string pp(const state_formulas::nu& x, bool arg0) { return state_formulas::pp< state_formulas::nu >(x, arg0); }
std::string pp(const state_formulas::or_& x, bool arg0) { return state_formulas::pp< state_formulas::or_ >(x, arg0); }
std::string pp(const state_formulas::plus& x, bool arg0) { return state_formulas::pp< state_formulas::plus >(x, arg0); }
std::string pp(const state_formulas::state_formula& x, bool arg0) { return state_formulas::pp< state_formulas::state_formula >(x, arg0); }
std::string pp(const state_formulas::state_formula_specification& x, bool arg0) { return state_formulas::pp< state_formulas::state_formula_specification >(x, arg0); }
std::string pp(const state_formulas::sum& x, bool arg0) { return state_formulas::pp< state_formulas::sum >(x, arg0); }
std::string pp(const state_formulas::supremum& x, bool arg0) { return state_formulas::pp< state_formulas::supremum >(x, arg0); }
std::string pp(const state_formulas::true_& x, bool arg0) { return state_formulas::pp< state_formulas::true_ >(x, arg0); }
std::string pp(const state_formulas::variable& x, bool arg0) { return state_formulas::pp< state_formulas::variable >(x, arg0); }
std::string pp(const state_formulas::yaled& x, bool arg0) { return state_formulas::pp< state_formulas::yaled >(x, arg0); }
std::string pp(const state_formulas::yaled_timed& x, bool arg0) { return state_formulas::pp< state_formulas::yaled_timed >(x, arg0); }
state_formulas::state_formula normalize_sorts(const state_formulas::state_formula& x, const data::sort_specification& sortspec) { return state_formulas::normalize_sorts< state_formulas::state_formula >(x, sortspec); }
state_formulas::state_formula translate_user_notation(const state_formulas::state_formula& x) { return state_formulas::translate_user_notation< state_formulas::state_formula >(x); }
std::set<data::sort_expression> find_sort_expressions(const state_formulas::state_formula& x) { return state_formulas::find_sort_expressions< state_formulas::state_formula >(x); }
std::set<data::variable> find_all_variables(const state_formulas::state_formula& x) { return state_formulas::find_all_variables< state_formulas::state_formula >(x); }
std::set<data::variable> find_free_variables(const state_formulas::state_formula& x) { return state_formulas::find_free_variables< state_formulas::state_formula >(x); }
std::set<core::identifier_string> find_identifiers(const state_formulas::state_formula& x) { return state_formulas::find_identifiers< state_formulas::state_formula >(x); }
std::set<process::action_label> find_action_labels(const state_formulas::state_formula& x) { return state_formulas::find_action_labels< state_formulas::state_formula >(x); }
//--- end generated state_formulas overloads ---//

namespace detail {

state_formula parse_state_formula(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("StateFrm");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  state_formula result = state_formula_actions(p).parse_StateFrm(node);
  return result;
}

state_formula_specification parse_state_formula_specification(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("StateFrmSpec");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  core::warn_left_merge_merge(node);

  untyped_state_formula_specification untyped_statespec = state_formula_actions(p).parse_StateFrmSpec(node);
  state_formula_specification result = untyped_statespec.construct_state_formula_specification();
  return result;
}

} // namespace detail

namespace algorithms {

state_formula parse_state_formula(std::istream& in, lps::stochastic_specification& lpsspec, const bool formula_is_quantitative)
{
  return state_formulas::parse_state_formula(in, lpsspec, formula_is_quantitative);
}

state_formula parse_state_formula(const std::string& text, lps::stochastic_specification& lpsspec, const bool formula_is_quantitative)
{
  return state_formulas::parse_state_formula(text, lpsspec, formula_is_quantitative);
}

state_formula_specification parse_state_formula_specification(std::istream& in, const bool formula_is_quantitative)
{
  return state_formulas::parse_state_formula_specification(in, formula_is_quantitative);
}

state_formula_specification parse_state_formula_specification(const std::string& text, const bool formula_is_quantitative)
{
  return state_formulas::parse_state_formula_specification(text, formula_is_quantitative);
}

state_formula_specification parse_state_formula_specification(std::istream& in, lps::stochastic_specification& lpsspec, const bool formula_is_quantitative)
{
  return state_formulas::parse_state_formula_specification(in, lpsspec, formula_is_quantitative);
}

state_formula_specification parse_state_formula_specification(const std::string& text, lps::stochastic_specification& lpsspec, const bool formula_is_quantitative)
{
  return state_formulas::parse_state_formula_specification(text, lpsspec, formula_is_quantitative);
}

bool is_monotonous(const state_formula& f)
{
  return state_formulas::is_monotonous(f);
}

state_formula normalize(const state_formula& x, bool quantitative=false, bool negated=false)
{
  return state_formulas::normalize(x, quantitative, negated);
}

bool is_normalized(const state_formula& x)
{
  return state_formulas::is_normalized(x);
}

bool is_timed(const state_formula& x)
{
  return state_formulas::is_timed(x);
}

std::set<core::identifier_string> find_state_variable_names(const state_formula& x)
{
  return state_formulas::find_state_variable_names(x);
}

} // namespace algorithms

} // namespace state_formulas

} // namespace mcrl2
