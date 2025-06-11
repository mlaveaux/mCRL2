// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file data.cpp
/// \brief

#include "mcrl2/data/normalize_sorts.h"
#include "mcrl2/data/parse_impl.h"
#include "mcrl2/data/print.h"
#include "mcrl2/data/substitutions/mutable_map_substitution.h"
#include "mcrl2/data/translate_user_notation.h"

namespace mcrl2
{

namespace data
{

//--- start generated data overloads ---//
std::string pp(const data::sort_expression_list& x, bool arg0) { return data::pp< data::sort_expression_list >(x, arg0); }
std::string pp(const data::sort_expression_vector& x, bool arg0) { return data::pp< data::sort_expression_vector >(x, arg0); }
std::string pp(const data::data_expression_list& x, bool arg0) { return data::pp< data::data_expression_list >(x, arg0); }
std::string pp(const data::data_expression_vector& x, bool arg0) { return data::pp< data::data_expression_vector >(x, arg0); }
std::string pp(const data::assignment_list& x, bool arg0) { return data::pp< data::assignment_list >(x, arg0); }
std::string pp(const data::assignment_vector& x, bool arg0) { return data::pp< data::assignment_vector >(x, arg0); }
std::string pp(const data::variable_list& x, bool arg0) { return data::pp< data::variable_list >(x, arg0); }
std::string pp(const data::variable_vector& x, bool arg0) { return data::pp< data::variable_vector >(x, arg0); }
std::string pp(const data::function_symbol_list& x, bool arg0) { return data::pp< data::function_symbol_list >(x, arg0); }
std::string pp(const data::function_symbol_vector& x, bool arg0) { return data::pp< data::function_symbol_vector >(x, arg0); }
std::string pp(const data::structured_sort_constructor_list& x, bool arg0) { return data::pp< data::structured_sort_constructor_list >(x, arg0); }
std::string pp(const data::structured_sort_constructor_vector& x, bool arg0) { return data::pp< data::structured_sort_constructor_vector >(x, arg0); }
std::string pp(const data::data_equation_list& x, bool arg0) { return data::pp< data::data_equation_list >(x, arg0); }
std::string pp(const data::data_equation_vector& x, bool arg0) { return data::pp< data::data_equation_vector >(x, arg0); }
std::string pp(const data::abstraction& x, bool arg0) { return data::pp< data::abstraction >(x, arg0); }
std::string pp(const data::alias& x, bool arg0) { return data::pp< data::alias >(x, arg0); }
std::string pp(const data::application& x, bool arg0) { return data::pp< data::application >(x, arg0); }
std::string pp(const data::assignment& x, bool arg0) { return data::pp< data::assignment >(x, arg0); }
std::string pp(const data::assignment_expression& x, bool arg0) { return data::pp< data::assignment_expression >(x, arg0); }
std::string pp(const data::bag_comprehension& x, bool arg0) { return data::pp< data::bag_comprehension >(x, arg0); }
std::string pp(const data::bag_comprehension_binder& x, bool arg0) { return data::pp< data::bag_comprehension_binder >(x, arg0); }
std::string pp(const data::bag_container& x, bool arg0) { return data::pp< data::bag_container >(x, arg0); }
std::string pp(const data::basic_sort& x, bool arg0) { return data::pp< data::basic_sort >(x, arg0); }
std::string pp(const data::binder_type& x, bool arg0) { return data::pp< data::binder_type >(x, arg0); }
std::string pp(const data::container_sort& x, bool arg0) { return data::pp< data::container_sort >(x, arg0); }
std::string pp(const data::container_type& x, bool arg0) { return data::pp< data::container_type >(x, arg0); }
std::string pp(const data::data_equation& x, bool arg0) { return data::pp< data::data_equation >(x, arg0); }
std::string pp(const data::data_expression& x, bool arg0) { return data::pp< data::data_expression >(x, arg0); }
std::string pp(const data::data_specification& x, bool arg0) { return data::pp< data::data_specification >(x, arg0); }
std::string pp(const data::exists& x, bool arg0) { return data::pp< data::exists >(x, arg0); }
std::string pp(const data::exists_binder& x, bool arg0) { return data::pp< data::exists_binder >(x, arg0); }
std::string pp(const data::fbag_container& x, bool arg0) { return data::pp< data::fbag_container >(x, arg0); }
std::string pp(const data::forall& x, bool arg0) { return data::pp< data::forall >(x, arg0); }
std::string pp(const data::forall_binder& x, bool arg0) { return data::pp< data::forall_binder >(x, arg0); }
std::string pp(const data::fset_container& x, bool arg0) { return data::pp< data::fset_container >(x, arg0); }
std::string pp(const data::function_sort& x, bool arg0) { return data::pp< data::function_sort >(x, arg0); }
std::string pp(const data::function_symbol& x, bool arg0) { return data::pp< data::function_symbol >(x, arg0); }
std::string pp(const data::lambda& x, bool arg0) { return data::pp< data::lambda >(x, arg0); }
std::string pp(const data::lambda_binder& x, bool arg0) { return data::pp< data::lambda_binder >(x, arg0); }
std::string pp(const data::list_container& x, bool arg0) { return data::pp< data::list_container >(x, arg0); }
std::string pp(const data::machine_number& x, bool arg0) { return data::pp< data::machine_number >(x, arg0); }
std::string pp(const data::set_comprehension& x, bool arg0) { return data::pp< data::set_comprehension >(x, arg0); }
std::string pp(const data::set_comprehension_binder& x, bool arg0) { return data::pp< data::set_comprehension_binder >(x, arg0); }
std::string pp(const data::set_container& x, bool arg0) { return data::pp< data::set_container >(x, arg0); }
std::string pp(const data::sort_expression& x, bool arg0) { return data::pp< data::sort_expression >(x, arg0); }
std::string pp(const data::structured_sort& x, bool arg0) { return data::pp< data::structured_sort >(x, arg0); }
std::string pp(const data::structured_sort_constructor& x, bool arg0) { return data::pp< data::structured_sort_constructor >(x, arg0); }
std::string pp(const data::structured_sort_constructor_argument& x, bool arg0) { return data::pp< data::structured_sort_constructor_argument >(x, arg0); }
std::string pp(const data::untyped_data_parameter& x, bool arg0) { return data::pp< data::untyped_data_parameter >(x, arg0); }
std::string pp(const data::untyped_identifier& x, bool arg0) { return data::pp< data::untyped_identifier >(x, arg0); }
std::string pp(const data::untyped_identifier_assignment& x, bool arg0) { return data::pp< data::untyped_identifier_assignment >(x, arg0); }
std::string pp(const data::untyped_possible_sorts& x, bool arg0) { return data::pp< data::untyped_possible_sorts >(x, arg0); }
std::string pp(const data::untyped_set_or_bag_comprehension& x, bool arg0) { return data::pp< data::untyped_set_or_bag_comprehension >(x, arg0); }
std::string pp(const data::untyped_set_or_bag_comprehension_binder& x, bool arg0) { return data::pp< data::untyped_set_or_bag_comprehension_binder >(x, arg0); }
std::string pp(const data::untyped_sort& x, bool arg0) { return data::pp< data::untyped_sort >(x, arg0); }
std::string pp(const data::untyped_sort_variable& x, bool arg0) { return data::pp< data::untyped_sort_variable >(x, arg0); }
std::string pp(const data::variable& x, bool arg0) { return data::pp< data::variable >(x, arg0); }
std::string pp(const data::where_clause& x, bool arg0) { return data::pp< data::where_clause >(x, arg0); }
data::data_equation normalize_sorts(const data::data_equation& x, const data::sort_specification& sortspec) { return data::normalize_sorts< data::data_equation >(x, sortspec); }
data::data_equation_list normalize_sorts(const data::data_equation_list& x, const data::sort_specification& sortspec) { return data::normalize_sorts< data::data_equation_list >(x, sortspec); }
void normalize_sorts(data::data_equation_vector& x, const data::sort_specification& sortspec) { data::normalize_sorts< data::data_equation_vector >(x, sortspec); }
data::data_expression normalize_sorts(const data::data_expression& x, const data::sort_specification& sortspec) { return data::normalize_sorts< data::data_expression >(x, sortspec); }
data::sort_expression normalize_sorts(const data::sort_expression& x, const data::sort_specification& sortspec) { return data::normalize_sorts< data::sort_expression >(x, sortspec); }
data::variable_list normalize_sorts(const data::variable_list& x, const data::sort_specification& sortspec) { return data::normalize_sorts< data::variable_list >(x, sortspec); }
data::data_expression translate_user_notation(const data::data_expression& x) { return data::translate_user_notation< data::data_expression >(x); }
data::data_equation translate_user_notation(const data::data_equation& x) { return data::translate_user_notation< data::data_equation >(x); }
std::set<data::sort_expression> find_sort_expressions(const data::data_equation& x) { return data::find_sort_expressions< data::data_equation >(x); }
std::set<data::sort_expression> find_sort_expressions(const data::data_expression& x) { return data::find_sort_expressions< data::data_expression >(x); }
std::set<data::sort_expression> find_sort_expressions(const data::sort_expression& x) { return data::find_sort_expressions< data::sort_expression >(x); }
std::set<data::variable> find_all_variables(const data::data_expression& x) { return data::find_all_variables< data::data_expression >(x); }
std::set<data::variable> find_all_variables(const data::data_expression_list& x) { return data::find_all_variables< data::data_expression_list >(x); }
std::set<data::variable> find_all_variables(const data::function_symbol& x) { return data::find_all_variables< data::function_symbol >(x); }
std::set<data::variable> find_all_variables(const data::variable& x) { return data::find_all_variables< data::variable >(x); }
std::set<data::variable> find_all_variables(const data::variable_list& x) { return data::find_all_variables< data::variable_list >(x); }
std::set<data::variable> find_free_variables(const data::data_expression& x) { return data::find_free_variables< data::data_expression >(x); }
std::set<data::variable> find_free_variables(const data::data_expression_list& x) { return data::find_free_variables< data::data_expression_list >(x); }
std::set<data::function_symbol> find_function_symbols(const data::data_equation& x) { return data::find_function_symbols< data::data_equation >(x); }
std::set<core::identifier_string> find_identifiers(const data::variable_list& x) { return data::find_identifiers< data::variable_list >(x); }
bool search_variable(const data::data_expression& x, const data::variable& v) { return data::search_variable< data::data_expression >(x, v); }
//--- end generated data overloads ---//

std::string pp(const std::set<variable>& x, bool precedence_aware) { return data::pp< std::set<variable> >(x, precedence_aware); }

sort_expression data_expression::sort() const
{
  using namespace atermpp;
  // This implementation is currently done in this class, because there
  // is no elegant solution of distributing the implementation of the
  // derived classes (as we need to support requesting the sort of a
  // data_expression we do need to provide an implementation here).
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  if (is_machine_number(*this))
  {
    return sort_machine_word::machine_word();
  }
#endif
  if (is_variable(*this))
  {
    const auto& v = atermpp::down_cast<variable>(*this);
    return v.sort();
  }
  else if (is_function_symbol(*this))
  {
    const auto& f = atermpp::down_cast<function_symbol>(*this);
    return f.sort();
  }
  else if (is_abstraction(*this))
  {
    if (is_forall(*this) || is_exists(*this))
    {
      return sort_bool::bool_();
    }
    else if (is_lambda(*this))
    {
      const auto& v_variables = atermpp::down_cast<atermpp::term_list<aterm> >((*this)[1]);
      sort_expression_vector s;
      for (const auto & v_variable : v_variables)
      {
        s.push_back(down_cast<sort_expression>(v_variable[1])); // Push the sort.
      }
      return function_sort(sort_expression_list(s.begin(),s.end()), atermpp::down_cast<data_expression>((*this)[2]).sort());
    }
    else
    {
      assert(is_set_comprehension(*this) || is_bag_comprehension(*this) || is_untyped_set_or_bag_comprehension(*this));
      const auto& v_variables  = atermpp::down_cast<atermpp::term_list<aterm> >((*this)[1]);
      assert(v_variables.size() == 1);

      if (is_bag_comprehension(*this))
      {
        return container_sort(bag_container(), atermpp::down_cast<const sort_expression>(v_variables.front()[1]));
      }
      else // If it is not known whether the term is a set or a bag, it returns the type of a set, as there is
           // no setbag type. This can only occur for terms that are not propertly type checked.
      {
        return container_sort(set_container(), atermpp::down_cast<sort_expression>(v_variables.front()[1]));
      }
    }
  }
  else if (is_where_clause(*this))
  {
    return atermpp::down_cast<data_expression>((*this)[0]).sort();
  }
  else if (is_untyped_identifier(*this))
  {
    return untyped_sort();
  }

  assert(is_application(*this));
  const auto& head = atermpp::down_cast<const data_expression>((*this)[0]);
  sort_expression s(head.sort());
  if (is_function_sort(s))
  {
    const auto& fs = atermpp::down_cast<function_sort>(s);
    assert(fs.domain().size()+1==this->size());
    return (fs.codomain());
  }
  return s;
}

std::set<data::variable> substitution_variables(const mutable_map_substitution<>& sigma)
{
  std::set<data::variable> result;
  for (const auto& i: sigma)
  {
    data::find_free_variables(i.second, std::inserter(result, result.end()));
  }
  return result;
}

variable_list free_variables(const data_expression& x)
{
  std::set<variable> v = find_free_variables(x);
  return variable_list(v.begin(), v.end());
}

namespace detail {

sort_expression parse_sort_expression(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("SortExpr");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  sort_expression result = data_expression_actions(p).parse_SortExpr(node);
  return result;
}

variable_list parse_variables(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("VarSpec");
  bool partial_parses = false;
  std::string var_text("var " + text);
  core::parse_node node = p.parse(var_text, start_symbol_index, partial_parses);
  variable_list result = data_specification_actions(p).parse_VarSpec(node);
  return result;
}

data_expression parse_data_expression(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("DataExpr");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  core::warn_and_or(node);
  data_expression result = data_expression_actions(p).parse_DataExpr(node);
  return result;
}

data_specification parse_data_specification_new(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("DataSpec");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  untyped_data_specification untyped_dataspec = data_specification_actions(p).parse_DataSpec(node);
  data_specification result = untyped_dataspec.construct_data_specification();
  return result;
}

variable_list parse_variable_declaration_list(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("VarsDeclList");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  variable_list result = data_specification_actions(p).parse_VarsDeclList(node);
  return result;
}

} // namespace detail

std::pair<basic_sort_vector, alias_vector> parse_sort_specification(const std::string& text)
{
  core::parser p(parser_tables_mcrl2, core::detail::ambiguity_fn, core::detail::syntax_error_fn);
  unsigned int start_symbol_index = p.start_symbol_index("SortSpec");
  bool partial_parses = false;
  core::parse_node node = p.parse(text, start_symbol_index, partial_parses);
  std::vector<atermpp::aterm> elements = detail::data_specification_actions(p).parse_SortSpec(node);
  basic_sort_vector sorts;
  alias_vector aliases;
  for (const atermpp::aterm& x: elements)
  {
    if (is_basic_sort(x))
    {
      sorts.push_back(atermpp::down_cast<basic_sort>(x));
    }
    else if (is_alias(x))
    {
      aliases.push_back(atermpp::down_cast<alias>(x));
    }
  }
  auto result = std::make_pair(sorts, aliases);
  return result;
}

} // namespace data

} // namespace mcrl2

