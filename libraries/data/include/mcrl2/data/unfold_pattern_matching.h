// Author(s): Ruud Koolen, Thomas Neele
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_DATA_UNFOLD_PATTERN_MATCHING_H
#define MCRL2_DATA_UNFOLD_PATTERN_MATCHING_H

#include "mcrl2/data/join.h"
#include "mcrl2/data/replace.h"
#include "mcrl2/data/substitutions/data_expression_assignment.h"
#include "mcrl2/data/representative_generator.h"
#include "mcrl2/data/substitutions/map_substitution.h"
#include "mcrl2/data/substitutions/variable_substitution.h"

namespace mcrl2::data
{

namespace detail
{
/**
 * \brief A rule describes a partially pattern-matched rewrite rule.
 * \details match_criteria is a set of data_expression pairs (A, B)
 * where A is a data_expression over variables in the left hand
 * side of a function definition, and B is a pattern consisting of
 * constructor applications over free variables.
 *
 * The rule match_criteria = { (A, B), (C, D) }, condition = E, rhs = R
 * describes the following rewrite proposition:
 *
 *  "if the data expression A pattern-matches to pattern B, and the
 *   data expression C pattern-matches to pattern D, and the condition
 *   E is true after substituting the proper pattern-matching variables,
 *   then the right hand side R applies (again with substitution of
 *   pattern-matched variables."
 *
 * Pattern matching can then be performed by deconstructing the patterns
 * in the right hand sides of match_criteria, and rewriting rules accordingly.
 * As an example, the following rewrite rule:
 *
 *   is_even(n) -> sign_of_list_sum(n |> l) = sign_of_list_sum(l)
 *
 * Can be represented as the following rule:
 *
 *   match_criteria = { v1 -> n |> l }, condition = is_even(n), rhs = sign_of_list_sum(l)
 *
 * Which, after one step of pattern matching, gets simplified to the following rule:
 *
 *   match_criteria = { head(v1) -> n, tail(v1) -> l }, condition = is_even(n), rhs = sign_of_list_sum(l)
 */
struct rule
{
  std::map<data_expression, data_expression> match_criteria;
  data_expression rhs;
  data_expression condition;
  variable_list variables;

  rule(const std::map<data_expression, data_expression>& mc,
       const data_expression& r,
       const data_expression& c,
       const variable_list& v)
  : match_criteria(mc)
  , rhs(r)
  , condition(c)
  , variables(v)
  {
    assert(condition.sort() == sort_bool::bool_());
  }

};

inline
std::ostream& operator<<(std::ostream& out, const rule& r)
{
  return out << core::detail::print_list(r.variables) << ". " << r.condition << " -> " << core::detail::print_map(r.match_criteria) << " = " << r.rhs;
}


/// @brief Get the top level function symbol f if expr is of the shape f or f(t1,...,tn)
/// @param expr The expression from which to extract the top function symbol
/// @return The top function symbol, or function_symbol() if it has none
inline
function_symbol get_top_fs(const data_expression& expr)
{
  if (is_function_symbol(expr))
  {
    return atermpp::down_cast<function_symbol>(expr);
  }
  if (is_application(expr))
  {
    const application& expr_appl = atermpp::down_cast<application>(expr);
    if (is_function_symbol(expr_appl.head()))
    {
      return atermpp::down_cast<function_symbol>(expr_appl.head());
    }
  }
  return function_symbol();
}

/**
 * \brief For a list of rules with equal left hand sides of match_criteria and only variables in
 * the right hand sides of match_criteria, construct a right hand side based on the
 * conditions and right hand sides of the rules.
 */
inline
data_expression construct_condition_rhs(const std::vector<rule>& rules, const data_expression& representative)
{
  // TODO: Check whether else_clause is equivalent to false. Can we use the enumerator for this?
  if (rules.size() == 0)
  {
    return representative;
  }

  data_expression result;
  for (const rule& r: rules)
  {
    std::map<variable, data_expression> substitution_map;
    for (const auto& [var_expr, pattern]: r.match_criteria)
    {
      assert(is_variable(pattern));
      substitution_map[atermpp::down_cast<variable>(pattern)] = var_expr;
    }
    map_substitution<std::map<variable, data_expression> > substitution(substitution_map);

    data_expression condition = replace_variables_capture_avoiding(r.condition, substitution);
    data_expression rhs = replace_variables_capture_avoiding(r.rhs, substitution);

    if (result == data_expression())
    {
      result = rhs;
    }
    else
    {
      result = lazy::if_(condition, rhs, result);
    }
  }

  return result;
}

/**
 * \brief For a list of rules with equal left hand sides of match_criteria, construct a right hand side.
 */
template <typename StructInfo>
data_expression construct_rhs(
  StructInfo& ssf,
  representative_generator& gen,
  const std::vector<rule>& rules,
  const sort_expression& sort
)
{
  if (rules.empty())
  {
    return construct_condition_rhs(rules, gen(sort));
  }

  /*
   * For each matching rule LHS, check whether pattern matching needs to happen on that LHS.
   *
   * We prefer to match on an expression such that *all* rules perform matching on that expression.
   * For example, for a ruleset representing the following rewrite rules:
   *
   *   optionally_remove_first_element(false, l) = l
   *   optionally_remove_first_element(true, []) = []
   *   optionally_remove_first_element(true, x |> l) = l
   *
   * We prefer to pattern match on the first parameter, not the second. Pattern matching on the
   * second parameter will (recursively) happen only on the true-branch, not the false-branch.
   * Pattern matching on the second parameter is still possible if there is no other option, but
   * it requires splitting the first rewrite rule into different cases for different
   * constructors -- suboptimal.
   */
  data_expression matching_target;
  enum class MatchType
  {
    NONE, INCOMPLETE, PARTIAL, FULL
  };
  MatchType matching_type = MatchType::NONE;
  // iterate over matches, this is the same for all rules, so just use rules[0]
  // More formally, we have rules[i].match_criteria->first is equivalent for all i.
  for (const auto& [var_expr, _]: rules[0].match_criteria)
  {
    bool variable_seen = false;
    std::set<function_symbol> constructors_seen;
    for (const rule& r: rules)
    {
      const data_expression& pattern = r.match_criteria.at(var_expr);
      if (is_variable(pattern))
      {
        variable_seen = true;
      }
      else
      {
        // either function symbol or application of function symbol
        function_symbol fs = get_top_fs(pattern);
        // this also means fs != function_symbol()
        assert(ssf.is_constructor(fs));
        constructors_seen.insert(fs);
      }
    }

    MatchType new_matching_type;
    if (constructors_seen.empty())
    {
      // No pattern matching is possible on this variable.
      new_matching_type = MatchType::NONE;
    }
    else if (variable_seen)
    {
      // There are both rules that match on this variable and rules that do not.
      // That's better than an incomplete matching but worse than a full matching.
      new_matching_type = MatchType::PARTIAL;
    }
    else if (constructors_seen.size() != ssf.get_constructors(var_expr.sort()).size())
    {
      // There are constructors for which there are no rules.
      // Thus, we have an incomplete function definition, that needs to be completed artificially.
      // A partial matching would be a better choice.
      new_matching_type = MatchType::INCOMPLETE;
    }
    else
    {
      // There is a matching rule for each constructor, and no rule with a plain variable.
      // This variable is a perfect pattern matching candidate.
      new_matching_type = MatchType::FULL;
    }

    if (new_matching_type > matching_type)
    {
      matching_target = var_expr;
      matching_type = new_matching_type;
    }
    if (matching_type == MatchType::FULL)
    {
      break;
    }
  }

  if (matching_type == MatchType::NONE)
  {
    // No constructor-based matching needs to happen.
    // All that needs to happen is incorporating the rule conditions.
    return construct_condition_rhs(rules, gen(sort));
  }

  /*
   * For each constructor, find the set of rules that apply to it, rewritten to match the constructor.
   */
  // auto type below depends on the type of ssf
  const auto& match_constructors = ssf.get_constructors(matching_target.sort());
  std::map<function_symbol, std::vector<rule> > constructor_rules;
  for (const rule& r: rules)
  {
    const data_expression& pattern = r.match_criteria.at(matching_target);
    if (is_function_symbol(pattern) || is_application(pattern))
    {
      /*
       * For a rule with a constructor pattern, strip the constructor and
       * introduce patterns for the constructor parameters.
       */
      function_symbol constructor = get_top_fs(pattern);
      assert(constructor != function_symbol());
      assert(utilities::detail::contains(match_constructors, constructor));

      data_expression_vector parameters;
      if (is_application(pattern))
      {
        const application& pattern_appl = atermpp::down_cast<application>(pattern);
        parameters.insert(parameters.end(), pattern_appl.begin(), pattern_appl.end());
      }

      rule rule = r;
      rule.match_criteria.erase(matching_target);
      // To prevent creating expressions of the shape head(l) |> tail(l), we perform a substitution here
      // This is only safe if there are no binders in the right-hand side
      std::set<data_expression> subexpr = find_data_expressions(rule.rhs);
      if (!parameters.empty() && std::none_of(subexpr.begin(), subexpr.end(), [](const data_expression& e) { return is_abstraction(e); }))
      {
        data_expression_assignment sigma(pattern,matching_target);
        rule.rhs = replace_data_expressions(rule.rhs, sigma, true);
        rule.condition = replace_data_expressions(rule.condition, sigma, true);
      }
      for (std::size_t j = 0; j < parameters.size(); j++)
      {
        function_symbol projection_function = ssf.get_projection_funcs(constructor)[j];
        data_expression lhs_expression = application(projection_function, matching_target);
        rule.match_criteria[lhs_expression] = parameters[j];
      }
      constructor_rules[constructor].push_back(rule);
    }
    else
    {
      /*
       * For a rule with a variable pattern that nonetheless needs to pattern match
       * against the possible constructors for that variable, copy the rule for each
       * possible constructor. Substitute the original un-matched term for the pattern
       * variable that disappears.
       */
      assert(is_variable(pattern));
      variable_substitution sigma(atermpp::down_cast<variable>(pattern), matching_target);
      data_expression condition = replace_variables_capture_avoiding(r.condition, sigma);
      data_expression rhs = replace_variables_capture_avoiding(r.rhs, sigma);

      for (const function_symbol& f: match_constructors)
      {
        rule rule(r.match_criteria, rhs, condition, r.variables);
        rule.match_criteria.erase(matching_target);

        set_identifier_generator generator;
        for (const variable& v: r.variables)
        {
          generator.add_identifier(v.name());
        }

        if (is_function_sort(f.sort()))
        {
          function_sort sort(f.sort());
          std::size_t index = 0;
          for (const sort_expression& s: sort.domain())
          {
            variable variable(generator("v"), s);
            function_symbol projection_function = ssf.get_projection_funcs(f)[index];
            data_expression lhs_expression = application(projection_function, matching_target);
            rule.match_criteria[lhs_expression] = variable;
            index++;
          }
        }

        constructor_rules[f].push_back(rule);
      }
    }
  }

  /*
   * Construct an rhs of the form if(is_cons1, rhs_cons1, if(is_cons2, rhs_cons2, ...)) or equivalent
   * The exact form depends on the implementation in ssf
   */
  data_expression_vector rhs;
  for (const auto& f: match_constructors)
  {
    rhs.push_back(construct_rhs(ssf, gen, constructor_rules[f], sort));
  }
  return ssf.create_cases(matching_target, rhs);
}

} // namespace detail

/**
 * \brief Check whether the given rewrite rule can be classified as a pattern matching rule.
 * \details That is, its arguments are constructed only out of unique variable occurrences and
 * constructor function symbols and constructor function applications.
 */
template <typename StructInfo>
bool is_pattern_matching_rule(StructInfo& ssf, const data_equation& rewrite_rule)
{
  // For this rewrite rule to be usable in pattern matching, its lhs must only contain
  // constructor functions and variables that occur at most once.

  std::set<data_expression> subexpressions = find_data_expressions(rewrite_rule.lhs());
  subexpressions.erase(rewrite_rule.lhs());
  if (is_application(rewrite_rule.lhs()))
  {
    subexpressions.erase(application(rewrite_rule.lhs()).head());
  }

  bool all_pattern = std::all_of(subexpressions.begin(), subexpressions.end(), [&ssf](const data_expression& expr) {
    return
      is_variable(expr) ||
      (is_function_symbol(expr) && ssf.is_constructor(atermpp::down_cast<function_symbol>(expr))) ||
      (is_application(expr) && is_function_symbol(atermpp::down_cast<application>(expr).head()) &&
                          ssf.is_constructor(function_symbol(atermpp::down_cast<application>(expr).head())));
  });
  if (!all_pattern)
  {
    return false;
  }
  if (std::all_of(subexpressions.begin(), subexpressions.end(), [](const data_expression& x){ return is_variable(x); }))
  {
    // Each argument is a variable, this is just an ordinarily defined function
    return false;
  }

  // Check whether each variable occurs at most once
  std::set<variable> variable_set;
  std::multiset<variable> variable_multiset;
  find_all_variables(rewrite_rule.lhs(), std::inserter(variable_set, variable_set.end()));
  find_all_variables(rewrite_rule.lhs(), std::inserter(variable_multiset, variable_multiset.end()));
  return variable_set.size() == variable_multiset.size();
}

/**
 * \brief This converts a collection of rewrite rules for a give function symbol into a
 * one-rule specification of the function, using recogniser and projection functions
 * to implement pattern matching.
 * \details For example, the collection of rewrite rules below:
 *
 *   sign_of_list_sum([]) = false;
 *   is_even(n) -> sign_of_list_sum(n |> l) = sign_of_list_sum(l);
 *   !is_even(n) -> sign_of_list_sum(n |> l) = !sign_of_list_sum(l);
 *
 * gets translated into the following function specification:
 *
 *   sign_of_list_sum(l) = if(is_emptylist(l), false,
 *                            if(is_even(head(l)), sign_of_list_sum(tail(l)),
 *                               !sign_of_list_sum(tail(l))))
 *
 * Two complications can arise. The rewrite rule set may contain rules that do not
 * pattern-match on the function parameters, such as 'not(not(b)) = b'; rules like
 * these are discarded.
 * More problematically, the set of rewrite rules may not be complete, or may not
 * easily be proven complete; in the example above, if the rewriter cannot rewrite
 * 'is_even(n) || !is_even(n)' to 'true', the translation cannot tell that the
 * rewrite rule set is complete.
 * In cases where a (non-constructor )function invocation can genuinely not be
 * rewritten any further, the MCRL2 semantics are unspecified (TODO check this);
 * the translation resolves this situation by returning an arbitrary value in this
 * case. Thus, in practice, the function definion above might well be translated
 * into the following:
 *
 *   sign_of_list_sum(l) = if(is_emptylist(l), false,
 *                            if(is_even(head(l)), sign_of_list_sum(tail(l)),
 *                               if(!is_even(head(l)), !sign_of_list_sum(tail(l)),
 *                                  false)))
 *
 * Where 'false' is a representative of sort_bool.
 */
template <typename StructInfo>
data_equation unfold_pattern_matching(
  const function_symbol& mapping,
  const data_equation_vector& rewrite_rules,
  StructInfo& ssf,
  representative_generator& gen,
  set_identifier_generator& id_gen
)
{
  sort_expression codomain = mapping.sort().target_sort();
  variable_vector temp_par;
  if (is_function_sort(mapping.sort()))
  {
    const function_sort& sort = atermpp::down_cast<function_sort>(mapping.sort());
    for (const sort_expression& s: sort.domain())
    {
      temp_par.emplace_back(id_gen("x"), s);
    }
  }
  variable_list new_parameters(temp_par.begin(), temp_par.end());

  // Create a rule for each data_equation
  std::vector<detail::rule> rules;
  for (const data_equation& eq: rewrite_rules)
  {
    assert(is_pattern_matching_rule(ssf, eq));

    std::map<data_expression, data_expression> match_criteria;
    if (is_application(eq.lhs()))
    {
      const application& lhs_appl = atermpp::down_cast<application>(eq.lhs());

      assert(lhs_appl.head() == mapping);
      assert(new_parameters.size() == lhs_appl.size());

      // Simultaneously iterate the parameters defined by the mapping and this
      // left-hand side to determine how matching occurs
      auto mappar_it = new_parameters.begin();
      auto lhspar_it = lhs_appl.begin();
      while(mappar_it != new_parameters.end())
      {
        match_criteria[*mappar_it] = *lhspar_it;
        ++mappar_it, ++lhspar_it;
      }

      assert(lhspar_it == lhs_appl.end());
    }

    detail::rule rule(match_criteria, eq.rhs(), eq.condition(), eq.variables());
    rules.push_back(rule);
  }
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  assert(rules.size() != 0 || sort_pos::is_most_significant_digit_function_symbol(mapping));
#else 
  assert(rules.size() != 0);
#endif

  data_expression new_lhs(application(mapping, new_parameters));
  data_expression new_rhs(construct_rhs(ssf, gen, rules, codomain));
  return data_equation(new_parameters, new_lhs, new_rhs);
}

} // namespace mcrl2::data

#endif
