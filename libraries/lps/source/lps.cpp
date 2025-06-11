// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lps.cpp
/// \brief

#include "mcrl2/lps/is_well_typed.h"
#include "mcrl2/lps/normalize_sorts.h"
#include "mcrl2/lps/parse_impl.h"
#include "mcrl2/lps/print.h"
#include "mcrl2/lps/replace.h"
#include "mcrl2/lps/translate_user_notation.h"

namespace mcrl2
{

namespace lps
{

//--- start generated lps overloads ---//
std::string pp(const lps::action_summand& x, bool arg0) { return lps::pp< lps::action_summand >(x, arg0); }
std::string pp(const lps::deadlock& x, bool arg0) { return lps::pp< lps::deadlock >(x, arg0); }
std::string pp(const lps::deadlock_summand& x, bool arg0) { return lps::pp< lps::deadlock_summand >(x, arg0); }
std::string pp(const lps::linear_process& x, bool arg0) { return lps::pp< lps::linear_process >(x, arg0); }
std::string pp(const lps::multi_action& x, bool arg0) { return lps::pp< lps::multi_action >(x, arg0); }
std::string pp(const lps::process_initializer& x, bool arg0) { return lps::pp< lps::process_initializer >(x, arg0); }
std::string pp(const lps::specification& x, bool arg0) { return lps::pp< lps::specification >(x, arg0); }
std::string pp(const lps::stochastic_action_summand& x, bool arg0) { return lps::pp< lps::stochastic_action_summand >(x, arg0); }
std::string pp(const lps::stochastic_distribution& x, bool arg0) { return lps::pp< lps::stochastic_distribution >(x, arg0); }
std::string pp(const lps::stochastic_linear_process& x, bool arg0) { return lps::pp< lps::stochastic_linear_process >(x, arg0); }
std::string pp(const lps::stochastic_process_initializer& x, bool arg0) { return lps::pp< lps::stochastic_process_initializer >(x, arg0); }
std::string pp(const lps::stochastic_specification& x, bool arg0) { return lps::pp< lps::stochastic_specification >(x, arg0); }
lps::multi_action normalize_sorts(const lps::multi_action& x, const data::sort_specification& sortspec) { return lps::normalize_sorts< lps::multi_action >(x, sortspec); }
void normalize_sorts(lps::specification& x, const data::sort_specification& /* sortspec */) { lps::normalize_sorts< lps::specification >(x, x.data()); }
void normalize_sorts(lps::stochastic_specification& x, const data::sort_specification& /* sortspec */) { lps::normalize_sorts< lps::stochastic_specification >(x, x.data()); }
lps::multi_action translate_user_notation(const lps::multi_action& x) { return lps::translate_user_notation< lps::multi_action >(x); }
std::set<data::sort_expression> find_sort_expressions(const lps::specification& x) { return lps::find_sort_expressions< lps::specification >(x); }
std::set<data::sort_expression> find_sort_expressions(const lps::stochastic_specification& x) { return lps::find_sort_expressions< lps::stochastic_specification >(x); }
std::set<data::variable> find_all_variables(const lps::linear_process& x) { return lps::find_all_variables< lps::linear_process >(x); }
std::set<data::variable> find_all_variables(const lps::stochastic_linear_process& x) { return lps::find_all_variables< lps::stochastic_linear_process >(x); }
std::set<data::variable> find_all_variables(const lps::specification& x) { return lps::find_all_variables< lps::specification >(x); }
std::set<data::variable> find_all_variables(const lps::stochastic_specification& x) { return lps::find_all_variables< lps::stochastic_specification >(x); }
std::set<data::variable> find_all_variables(const lps::deadlock& x) { return lps::find_all_variables< lps::deadlock >(x); }
std::set<data::variable> find_all_variables(const lps::multi_action& x) { return lps::find_all_variables< lps::multi_action >(x); }
std::set<data::variable> find_free_variables(const lps::linear_process& x) { return lps::find_free_variables< lps::linear_process >(x); }
std::set<data::variable> find_free_variables(const lps::stochastic_linear_process& x) { return lps::find_free_variables< lps::stochastic_linear_process >(x); }
std::set<data::variable> find_free_variables(const lps::specification& x) { return lps::find_free_variables< lps::specification >(x); }
std::set<data::variable> find_free_variables(const lps::stochastic_specification& x) { return lps::find_free_variables< lps::stochastic_specification >(x); }
std::set<data::variable> find_free_variables(const lps::deadlock& x) { return lps::find_free_variables< lps::deadlock >(x); }
std::set<data::variable> find_free_variables(const lps::multi_action& x) { return lps::find_free_variables< lps::multi_action >(x); }
std::set<data::variable> find_free_variables(const lps::process_initializer& x) { return lps::find_free_variables< lps::process_initializer >(x); }
std::set<data::variable> find_free_variables(const lps::stochastic_process_initializer& x) { return lps::find_free_variables< lps::stochastic_process_initializer >(x); }
std::set<data::function_symbol> find_function_symbols(const lps::specification& x) { return lps::find_function_symbols< lps::specification >(x); }
std::set<data::function_symbol> find_function_symbols(const lps::stochastic_specification& x) { return lps::find_function_symbols< lps::stochastic_specification >(x); }
std::set<core::identifier_string> find_identifiers(const lps::specification& x) { return lps::find_identifiers< lps::specification >(x); }
std::set<core::identifier_string> find_identifiers(const lps::stochastic_specification& x) { return lps::find_identifiers< lps::stochastic_specification >(x); }
std::set<process::action_label> find_action_labels(const lps::linear_process& x) { return lps::find_action_labels< lps::linear_process >(x); }
std::set<process::action_label> find_action_labels(const lps::process_initializer& x) { return lps::find_action_labels< lps::process_initializer >(x); }
std::set<process::action_label> find_action_labels(const lps::specification& x) { return lps::find_action_labels< lps::specification >(x); }
std::set<process::action_label> find_action_labels(const lps::stochastic_specification& x) { return lps::find_action_labels< lps::stochastic_specification >(x); }
//--- end generated lps overloads ---//

data::data_expression_list action_summand::next_state(const data::variable_list& process_parameters) const
{
  return data::replace_variables(atermpp::container_cast<data::data_expression_list>(process_parameters),
                                 data::assignment_sequence_substitution(assignments()));
}

std::string pp_with_summand_numbers(const specification& x)
{
  std::ostringstream out;
  core::detail::apply_printer<lps::detail::printer> printer(out, true);
  printer.print_summand_numbers() = true;
  printer.apply(x);
  return out.str();
}

std::string pp_with_summand_numbers(const stochastic_specification& x)
{
  std::ostringstream out;
  core::detail::apply_printer<lps::detail::printer> printer(out, true);
  printer.print_summand_numbers() = true;
  printer.apply(x);
  return out.str();
}

bool check_well_typedness(const linear_process& x)
{
  return lps::detail::check_well_typedness(x);
}

bool check_well_typedness(const stochastic_linear_process& x)
{
  return lps::detail::check_well_typedness(x);
}

bool check_well_typedness(const specification& x)
{
  return lps::detail::check_well_typedness(x);
}

bool check_well_typedness(const stochastic_specification& x)
{
  return lps::detail::check_well_typedness(x);
}

namespace detail {

process::untyped_multi_action parse_multi_action_new(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("MultAct");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  process::untyped_multi_action result = multi_action_actions(p).parse_MultAct(node);
  return result;
}

multi_action complete_multi_action(process::untyped_multi_action& x, multi_action_type_checker& typechecker, const data::data_specification& data_spec = data::detail::default_specification())
{
  multi_action result = lps::typecheck_multi_action(x, typechecker);
  result = lps::translate_user_notation(result);
  lps::normalize_sorts(result, data_spec);
  return result;
}

multi_action complete_multi_action(process::untyped_multi_action& x, const process::action_label_list& action_decls, const data::data_specification& data_spec = data::detail::default_specification())
{
  multi_action result = lps::typecheck_multi_action(x, data_spec, action_decls);
  result = lps::translate_user_notation(result);
  lps::normalize_sorts(result, data_spec);
  return result;
}

action_rename_specification parse_action_rename_specification_new(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("ActionRenameSpec");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  action_rename_specification result = action_rename_actions(p).parse_ActionRenameSpec(node);
  return result;
}

void complete_action_rename_specification(action_rename_specification& x, const lps::stochastic_specification& spec)
{
  using namespace mcrl2::data;
  x = lps::typecheck_action_rename_specification(x, spec);
  x = action_rename_specification(x.data() + spec.data(), x.action_labels(), x.rules());
  detail::translate_user_notation(x);
}

} // namespace detail

} // namespace lps

} // namespace mcrl2
