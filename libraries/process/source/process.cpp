// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file process.cpp
/// \brief

#include "mcrl2/process/detail/alphabet_push_block.h"
#include "mcrl2/process/parse_impl.h"
#include "mcrl2/process/translate_user_notation.h"
#include "mcrl2/process/remove_equations.h"


namespace mcrl2
{

namespace process
{

//--- start generated process overloads ---//
std::string pp(const process::action_list& x, bool arg0) { return process::pp< process::action_list >(x, arg0); }
std::string pp(const process::action_vector& x, bool arg0) { return process::pp< process::action_vector >(x, arg0); }
std::string pp(const process::action_label_list& x, bool arg0) { return process::pp< process::action_label_list >(x, arg0); }
std::string pp(const process::action_label_vector& x, bool arg0) { return process::pp< process::action_label_vector >(x, arg0); }
std::string pp(const process::process_identifier_list& x, bool arg0) { return process::pp< process::process_identifier_list >(x, arg0); }
std::string pp(const process::process_identifier_vector& x, bool arg0) { return process::pp< process::process_identifier_vector >(x, arg0); }
std::string pp(const process::process_expression_list& x, bool arg0) { return process::pp< process::process_expression_list >(x, arg0); }
std::string pp(const process::process_expression_vector& x, bool arg0) { return process::pp< process::process_expression_vector >(x, arg0); }
std::string pp(const process::process_equation_list& x, bool arg0) { return process::pp< process::process_equation_list >(x, arg0); }
std::string pp(const process::process_equation_vector& x, bool arg0) { return process::pp< process::process_equation_vector >(x, arg0); }
std::string pp(const process::action& x, bool arg0) { return process::pp< process::action >(x, arg0); }
std::string pp(const process::action_label& x, bool arg0) { return process::pp< process::action_label >(x, arg0); }
std::string pp(const process::action_name_multiset& x, bool arg0) { return process::pp< process::action_name_multiset >(x, arg0); }
std::string pp(const process::allow& x, bool arg0) { return process::pp< process::allow >(x, arg0); }
std::string pp(const process::at& x, bool arg0) { return process::pp< process::at >(x, arg0); }
std::string pp(const process::block& x, bool arg0) { return process::pp< process::block >(x, arg0); }
std::string pp(const process::bounded_init& x, bool arg0) { return process::pp< process::bounded_init >(x, arg0); }
std::string pp(const process::choice& x, bool arg0) { return process::pp< process::choice >(x, arg0); }
std::string pp(const process::comm& x, bool arg0) { return process::pp< process::comm >(x, arg0); }
std::string pp(const process::communication_expression& x, bool arg0) { return process::pp< process::communication_expression >(x, arg0); }
std::string pp(const process::delta& x, bool arg0) { return process::pp< process::delta >(x, arg0); }
std::string pp(const process::hide& x, bool arg0) { return process::pp< process::hide >(x, arg0); }
std::string pp(const process::if_then& x, bool arg0) { return process::pp< process::if_then >(x, arg0); }
std::string pp(const process::if_then_else& x, bool arg0) { return process::pp< process::if_then_else >(x, arg0); }
std::string pp(const process::left_merge& x, bool arg0) { return process::pp< process::left_merge >(x, arg0); }
std::string pp(const process::merge& x, bool arg0) { return process::pp< process::merge >(x, arg0); }
std::string pp(const process::process_equation& x, bool arg0) { return process::pp< process::process_equation >(x, arg0); }
std::string pp(const process::process_expression& x, bool arg0) { return process::pp< process::process_expression >(x, arg0); }
std::string pp(const process::process_identifier& x, bool arg0) { return process::pp< process::process_identifier >(x, arg0); }
std::string pp(const process::process_instance& x, bool arg0) { return process::pp< process::process_instance >(x, arg0); }
std::string pp(const process::process_instance_assignment& x, bool arg0) { return process::pp< process::process_instance_assignment >(x, arg0); }
std::string pp(const process::process_specification& x, bool arg0) { return process::pp< process::process_specification >(x, arg0); }
std::string pp(const process::rename& x, bool arg0) { return process::pp< process::rename >(x, arg0); }
std::string pp(const process::rename_expression& x, bool arg0) { return process::pp< process::rename_expression >(x, arg0); }
std::string pp(const process::seq& x, bool arg0) { return process::pp< process::seq >(x, arg0); }
std::string pp(const process::stochastic_operator& x, bool arg0) { return process::pp< process::stochastic_operator >(x, arg0); }
std::string pp(const process::sum& x, bool arg0) { return process::pp< process::sum >(x, arg0); }
std::string pp(const process::sync& x, bool arg0) { return process::pp< process::sync >(x, arg0); }
std::string pp(const process::tau& x, bool arg0) { return process::pp< process::tau >(x, arg0); }
std::string pp(const process::untyped_multi_action& x, bool arg0) { return process::pp< process::untyped_multi_action >(x, arg0); }
std::string pp(const process::untyped_process_assignment& x, bool arg0) { return process::pp< process::untyped_process_assignment >(x, arg0); }
process::action normalize_sorts(const process::action& x, const data::sort_specification& sortspec) { return process::normalize_sorts< process::action >(x, sortspec); }
process::action_label_list normalize_sorts(const process::action_label_list& x, const data::sort_specification& sortspec) { return process::normalize_sorts< process::action_label_list >(x, sortspec); }
void normalize_sorts(process::process_equation_vector& x, const data::sort_specification& sortspec) { process::normalize_sorts< process::process_equation_vector >(x, sortspec); }
void normalize_sorts(process::process_specification& x, const data::sort_specification& /* sortspec */) { process::normalize_sorts< process::process_specification >(x, x.data()); }
process::action translate_user_notation(const process::action& x) { return process::translate_user_notation< process::action >(x); }
process::process_expression translate_user_notation(const process::process_expression& x) { return process::translate_user_notation< process::process_expression >(x); }
void translate_user_notation(process::process_specification& x) { process::translate_user_notation< process::process_specification >(x); }
std::set<data::sort_expression> find_sort_expressions(const process::action_label_list& x) { return process::find_sort_expressions< process::action_label_list >(x); }
std::set<data::sort_expression> find_sort_expressions(const process::process_equation_vector& x) { return process::find_sort_expressions< process::process_equation_vector >(x); }
std::set<data::sort_expression> find_sort_expressions(const process::process_expression& x) { return process::find_sort_expressions< process::process_expression >(x); }
std::set<data::sort_expression> find_sort_expressions(const process::process_specification& x) { return process::find_sort_expressions< process::process_specification >(x); }
std::set<data::variable> find_all_variables(const process::action& x) { return process::find_all_variables< process::action >(x); }
std::set<data::variable> find_free_variables(const process::action& x) { return process::find_free_variables< process::action >(x); }
std::set<data::variable> find_free_variables(const process::process_specification& x) { return process::find_free_variables< process::process_specification >(x); }
std::set<core::identifier_string> find_identifiers(const process::process_specification& x) { return process::find_identifiers< process::process_specification >(x); }
//--- end generated process overloads ---//

void alphabet_reduce(process_specification& procspec, std::size_t duplicate_equation_limit)
{
  mCRL2log(log::verbose) << "applying alphabet reduction..." << std::endl;
  process_expression init = procspec.init();

  // cache the alphabet of pcrl equations and apply alphabet reduction to block({}, init)
  std::vector<process_equation>& equations = procspec.equations();
  std::map<process_identifier, multi_action_name_set> pcrl_equation_cache;
  data::set_identifier_generator id_generator;
  for (process_equation& equation: equations)
  {
    id_generator.add_identifier(equation.identifier().name());
  }
  pcrl_equation_cache = detail::compute_pcrl_equation_cache(equations, init);
  core::identifier_string_list empty_blockset;
  procspec.init() = push_block(empty_blockset, init, equations, id_generator, pcrl_equation_cache);

  // remove duplicate equations
  if (procspec.equations().size() < duplicate_equation_limit)
  {
    mCRL2log(log::debug) << "removing duplicate equations..." << std::endl;
    remove_duplicate_equations(procspec);
    mCRL2log(log::debug) << "removing duplicate equations finished" << std::endl;
  }

  mCRL2log(log::debug) << "alphabet reduction finished" << std::endl;
}

process::action_label_list parse_action_declaration(const std::string& text, const data::data_specification& data_spec)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("ActDecl");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  action_label_vector result;
  detail::action_actions(p).callback_ActDecl(node, result);
  process::action_label_list v(result.begin(), result.end());
  v = process::normalize_sorts(v, data_spec);
  return v;
}

namespace detail {

process_expression parse_process_expression_new(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("ProcExpr");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  core::warn_left_merge_merge(node);
  process_expression result = process_actions(p).parse_ProcExpr(node);
  return result;
}

process_specification parse_process_specification_new(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("mCRL2Spec");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  core::warn_left_merge_merge(node);
  untyped_process_specification untyped_procspec = process_actions(p).parse_mCRL2Spec(node);
  process_specification result = untyped_procspec.construct_process_specification();
  return result;
}

void complete_process_specification(process_specification& x, bool alpha_reduce)
{
  typecheck_process_specification(x);
  process::translate_user_notation(x);
  if (alpha_reduce)
  {
    alphabet_reduce(x, 1000ul);
  }
}

} // namespace detail

} // namespace process

} // namespace mcrl2

