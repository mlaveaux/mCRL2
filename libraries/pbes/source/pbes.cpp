// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbes.cpp
/// \brief

#include "mcrl2/pbes/detail/has_propositional_variables.h"
#include "mcrl2/pbes/detail/instantiate_global_variables.h"
#include "mcrl2/pbes/detail/is_well_typed.h"
#include "mcrl2/pbes/detail/occurring_variable_visitor.h"
#include "mcrl2/pbes/is_bes.h"
#include "mcrl2/pbes/parse_impl.h"
#include "mcrl2/pbes/print.h"
#include "mcrl2/pbes/translate_user_notation.h"

namespace mcrl2
{

namespace pbes_system
{

//--- start generated pbes_system overloads ---//
std::string pp(const pbes_system::pbes_equation_vector& x, bool arg0) { return pbes_system::pp< pbes_system::pbes_equation_vector >(x, arg0); }
std::string pp(const pbes_system::pbes_expression_list& x, bool arg0) { return pbes_system::pp< pbes_system::pbes_expression_list >(x, arg0); }
std::string pp(const pbes_system::pbes_expression_vector& x, bool arg0) { return pbes_system::pp< pbes_system::pbes_expression_vector >(x, arg0); }
std::string pp(const pbes_system::propositional_variable_list& x, bool arg0) { return pbes_system::pp< pbes_system::propositional_variable_list >(x, arg0); }
std::string pp(const pbes_system::propositional_variable_vector& x, bool arg0) { return pbes_system::pp< pbes_system::propositional_variable_vector >(x, arg0); }
std::string pp(const pbes_system::propositional_variable_instantiation_list& x, bool arg0) { return pbes_system::pp< pbes_system::propositional_variable_instantiation_list >(x, arg0); }
std::string pp(const pbes_system::propositional_variable_instantiation_vector& x, bool arg0) { return pbes_system::pp< pbes_system::propositional_variable_instantiation_vector >(x, arg0); }
std::string pp(const pbes_system::and_& x, bool arg0) { return pbes_system::pp< pbes_system::and_ >(x, arg0); }
std::string pp(const pbes_system::exists& x, bool arg0) { return pbes_system::pp< pbes_system::exists >(x, arg0); }
std::string pp(const pbes_system::fixpoint_symbol& x, bool arg0) { return pbes_system::pp< pbes_system::fixpoint_symbol >(x, arg0); }
std::string pp(const pbes_system::forall& x, bool arg0) { return pbes_system::pp< pbes_system::forall >(x, arg0); }
std::string pp(const pbes_system::imp& x, bool arg0) { return pbes_system::pp< pbes_system::imp >(x, arg0); }
std::string pp(const pbes_system::not_& x, bool arg0) { return pbes_system::pp< pbes_system::not_ >(x, arg0); }
std::string pp(const pbes_system::or_& x, bool arg0) { return pbes_system::pp< pbes_system::or_ >(x, arg0); }
std::string pp(const pbes_system::pbes& x, bool arg0) { return pbes_system::pp< pbes_system::pbes >(x, arg0); }
std::string pp(const pbes_system::pbes_equation& x, bool arg0) { return pbes_system::pp< pbes_system::pbes_equation >(x, arg0); }
std::string pp(const pbes_system::pbes_expression& x, bool arg0) { return pbes_system::pp< pbes_system::pbes_expression >(x, arg0); }
std::string pp(const pbes_system::propositional_variable& x, bool arg0) { return pbes_system::pp< pbes_system::propositional_variable >(x, arg0); }
std::string pp(const pbes_system::propositional_variable_instantiation& x, bool arg0) { return pbes_system::pp< pbes_system::propositional_variable_instantiation >(x, arg0); }
void normalize_sorts(pbes_system::pbes_equation_vector& x, const data::sort_specification& sortspec) { pbes_system::normalize_sorts< pbes_system::pbes_equation_vector >(x, sortspec); }
void normalize_sorts(pbes_system::pbes& x, const data::sort_specification& /* sortspec */) { pbes_system::normalize_sorts< pbes_system::pbes >(x, x.data()); }
pbes_system::pbes_expression normalize_sorts(const pbes_system::pbes_expression& x, const data::sort_specification& sortspec) { return pbes_system::normalize_sorts< pbes_system::pbes_expression >(x, sortspec); }
void translate_user_notation(pbes_system::pbes& x) { pbes_system::translate_user_notation< pbes_system::pbes >(x); }
pbes_system::pbes_expression translate_user_notation(const pbes_system::pbes_expression& x) { return pbes_system::translate_user_notation< pbes_system::pbes_expression >(x); }
std::set<data::sort_expression> find_sort_expressions(const pbes_system::pbes& x) { return pbes_system::find_sort_expressions< pbes_system::pbes >(x); }
std::set<data::variable> find_all_variables(const pbes_system::pbes& x) { return pbes_system::find_all_variables< pbes_system::pbes >(x); }
std::set<data::variable> find_free_variables(const pbes_system::pbes& x) { return pbes_system::find_free_variables< pbes_system::pbes >(x); }
std::set<data::variable> find_free_variables(const pbes_system::pbes_expression& x) { return pbes_system::find_free_variables< pbes_system::pbes_expression >(x); }
std::set<data::variable> find_free_variables(const pbes_system::pbes_equation& x) { return pbes_system::find_free_variables< pbes_system::pbes_equation >(x); }
std::set<data::function_symbol> find_function_symbols(const pbes_system::pbes& x) { return pbes_system::find_function_symbols< pbes_system::pbes >(x); }
std::set<pbes_system::propositional_variable_instantiation> find_propositional_variable_instantiations(const pbes_system::pbes_expression& x) { return pbes_system::find_propositional_variable_instantiations< pbes_system::pbes_expression >(x); }
std::set<core::identifier_string> find_identifiers(const pbes_system::pbes_expression& x) { return pbes_system::find_identifiers< pbes_system::pbes_expression >(x); }
bool search_variable(const pbes_system::pbes_expression& x, const data::variable& v) { return pbes_system::search_variable< pbes_system::pbes_expression >(x, v); }
//--- end generated pbes_system overloads ---//

namespace algorithms {

void instantiate_global_variables(pbes& p)
{
  pbes_system::detail::instantiate_global_variables(p);
}

bool is_bes(const pbes& x)
{
  return pbes_system::is_bes(x);
}

} // namespace algorithms

bool is_well_typed(const pbes_equation& eqn)
{
  return pbes_system::detail::is_well_typed(eqn);
}

bool is_well_typed_equation(const pbes_equation& eqn,
                            const std::set<data::sort_expression>& declared_sorts,
                            const std::set<data::variable>& declared_global_variables,
                            const data::data_specification& data_spec
                           )
{
  return pbes_system::detail::is_well_typed_equation(eqn, declared_sorts, declared_global_variables, data_spec);
}

bool is_well_typed_pbes(const std::set<data::sort_expression>& declared_sorts,
                        const std::set<data::variable>& declared_global_variables,
                        const std::set<data::variable>& occurring_global_variables,
                        const std::set<propositional_variable>& declared_variables,
                        const std::set<propositional_variable_instantiation>& occ,
                        const propositional_variable_instantiation& init,
                        const data::data_specification& data_spec
                       )
{
  return pbes_system::detail::is_well_typed_pbes(declared_sorts, declared_global_variables, occurring_global_variables, declared_variables, occ, init, data_spec);
}

bool pbes_equation::is_solved() const
{
  return !detail::has_propositional_variables(formula());
}

std::set<propositional_variable_instantiation> pbes::occurring_variable_instantiations() const
{
  std::set<propositional_variable_instantiation> result;
  for (const pbes_equation& eqn: equations())
  {
    detail::occurring_variable_visitor f;
    f.apply(eqn.formula());
    result.insert(f.variables.begin(), f.variables.end());
  }
  return result;
}

namespace detail {

pbes_expression parse_pbes_expression_new(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("PbesExpr");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  pbes_expression result = pbes_actions(p).parse_PbesExpr(node);
  return result;
}

untyped_pbes parse_pbes_new(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("PbesSpec");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  untyped_pbes result = pbes_actions(p).parse_PbesSpec(node);
  return result;
}

void complete_pbes(pbes& x)
{
  typecheck_pbes(x);
  pbes_system::translate_user_notation(x);
  complete_data_specification(x);
}

propositional_variable parse_propositional_variable(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("PropVarDecl");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  return detail::pbes_actions(p).parse_PropVarDecl(node);
}

pbes_expression parse_pbes_expression(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("PbesExpr");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  return detail::pbes_actions(p).parse_PbesExpr(node);
}

} // namespace detail

} // namespace pbes_system

} // namespace mcrl2

