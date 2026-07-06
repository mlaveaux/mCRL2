#!/usr/bin/env python3

#~ Copyright 2011 Wieger Wesselink.
#~ Distributed under the Boost Software License, Version 1.0.
#~ (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

#pylint: disable=line-too-long
#pylint: disable=missing-function-docstring
import re
import sys
from typeguard import typechecked

from mcrl2_utility import insert_text_in_file

MCRL2_ROOT: str = '../../'

REWRITE_TEXT: str = '''/// \\\\brief Rewrites all embedded expressions in an object x
/// \\\\param x an object containing expressions
/// \\\\param R a rewriter
template <typename T, typename Rewriter>
void rewrite(T& x,
             Rewriter R
            )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_rewrite_data_expressions_builder<NAMESPACE::data_expression_builder>(R).update(x);
}

/// \\\\brief Rewrites all embedded expressions in an object x
/// \\\\param x an object containing expressions
/// \\\\param R a rewriter
/// \\\\return the rewrite result
template <typename T, typename Rewriter>
T rewrite(const T& x,
          Rewriter R
         )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::make_rewrite_data_expressions_builder<NAMESPACE::data_expression_builder>(R).apply(result, x);
  return result;
}

/// \\\\brief Rewrites all embedded expressions in an object x, and applies a substitution to variables on the fly
/// \\\\param x an object containing expressions
/// \\\\param R a rewriter
/// \\\\param sigma a substitution
template <typename T, typename Rewriter, data::IsSubstitution Substitution>
void rewrite(T& x,
             Rewriter R,
             const Substitution& sigma
            )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_rewrite_data_expressions_with_substitution_builder<NAMESPACE::data_expression_builder>(R, sigma).update(x);
}

/// \\\\brief Rewrites all embedded expressions in an object x, and applies a substitution to variables on the fly
/// \\\\param x an object containing expressions
/// \\\\param R a rewriter
/// \\\\param sigma a substitution
/// \\\\return the rewrite result
template <typename T, typename Rewriter, data::IsSubstitution Substitution>
T rewrite(const T& x,
          Rewriter R,
          const Substitution& sigma
         )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::make_rewrite_data_expressions_with_substitution_builder<NAMESPACE::data_expression_builder>(R, sigma).apply(result, x);
  return result;
}
'''

SUBSTITUTE_FUNCTION_TEXT: str = '''template <typename T, typename Substitution>
void replace_sort_expressions(T& x,
                              const Substitution& sigma,
                              bool innermost
                             )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_replace_sort_expressions_builder<NAMESPACE::sort_expression_builder>(sigma, innermost).update(x);
}

template <typename T, typename Substitution>
T replace_sort_expressions(const T& x,
                           const Substitution& sigma,
                           bool innermost
                          )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::make_replace_sort_expressions_builder<NAMESPACE::sort_expression_builder>(sigma, innermost).apply(result, x);
  return result;
}

template <typename T, data::IsSubstitution Substitution>
void replace_data_expressions(T& x,
                              const Substitution& sigma,
                              bool innermost
                             )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::make_replace_data_expressions_builder<NAMESPACE::data_expression_builder>(sigma, innermost).update(x);
}

template <typename T, data::IsSubstitution Substitution>
T replace_data_expressions(const T& x,
                           const Substitution& sigma,
                           bool innermost
                          )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::make_replace_data_expressions_builder<NAMESPACE::data_expression_builder>(sigma, innermost).apply(result, x);
  return result;
}


template <typename T, data::IsSubstitution Substitution>
void replace_variables(T& x,
                       const Substitution& sigma
                      )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  core::make_update_apply_builder<NAMESPACE::data_expression_builder>(sigma).update(x);
}

template <typename T, data::IsSubstitution Substitution>
T replace_variables(const T& x,
                    const Substitution& sigma
                   )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  core::make_update_apply_builder<NAMESPACE::data_expression_builder>(sigma).apply(result, x);
  return result;
}

/* Replace all variables, including those in binders and the left hand side of assignments */
template <typename T, data::IsSubstitution Substitution>
void replace_all_variables(T& x,
                           const Substitution& sigma
                          )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  core::make_update_apply_builder<NAMESPACE::sort_expression_builder>(sigma).update(x);
}

template <typename T, data::IsSubstitution Substitution>
T replace_all_variables(const T& x,
                        const Substitution& sigma
                       )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  core::make_update_apply_builder<NAMESPACE::sort_expression_builder>(sigma).apply(result, x);
  return result;
}

/// \\\\brief Applies the substitution sigma to x.
/// \\\\pre { The substitution sigma must have the property that FV(sigma(x)) is included in {x} for all variables x. }
template <typename T, data::IsSubstitution Substitution>
void replace_free_variables(T& x,
                            const Substitution& sigma
                           )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  assert(data::is_simple_substitution(sigma));
  data::detail::make_replace_free_variables_builder<NAMESPACE::data_expression_builder, NAMESPACE::add_data_variable_builder_binding>(sigma).update(x);
}

/// \\\\brief Applies the substitution sigma to x.
/// \\\\pre { The substitution sigma must have the property that FV(sigma(x)) is included in {x} for all variables x. }
template <typename T, data::IsSubstitution Substitution>
T replace_free_variables(const T& x,
                         const Substitution& sigma
                        )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  assert(data::is_simple_substitution(sigma));
  T result;
  data::detail::make_replace_free_variables_builder<NAMESPACE::data_expression_builder, NAMESPACE::add_data_variable_builder_binding>(sigma).apply(result, x);
  return result;
}

/// \\\\brief Applies the substitution sigma to x, where the elements of bound_variables are treated as bound variables.
/// \\\\pre { The substitution sigma must have the property that FV(sigma(x)) is included in {x} for all variables x. }
template <typename T, data::IsSubstitution Substitution, typename VariableContainer>
void replace_free_variables(T& x,
                            const Substitution& sigma,
                            const VariableContainer& bound_variables
                           )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  assert(data::is_simple_substitution(sigma));
  data::detail::make_replace_free_variables_builder<NAMESPACE::data_expression_builder, NAMESPACE::add_data_variable_builder_binding>(sigma).update(x, bound_variables);
}

/// \\\\brief Applies the substitution sigma to x, where the elements of bound_variables are treated as bound variables.
/// \\\\pre { The substitution sigma must have the property that FV(sigma(x)) is included in {x} for all variables x. }
template <typename T, data::IsSubstitution Substitution, typename VariableContainer>
T replace_free_variables(const T& x,
                         const Substitution& sigma,
                         const VariableContainer& bound_variables
                        )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  assert(data::is_simple_substitution(sigma));
  T result;
  data::detail::make_replace_free_variables_builder<NAMESPACE::data_expression_builder, NAMESPACE::add_data_variable_builder_binding>(sigma).apply(result, x, bound_variables);
  return result;
}
'''

REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT: str = '''/// \\\\brief Applies sigma as a capture avoiding substitution to x.
/// \\\\param x The object to which the subsitution is applied.
/// \\\\param sigma A substitution.
/// \\\\param id_generator An identifier generator that generates names that do not appear in x and sigma
template <typename T, data::IsSubstitution Substitution>
void replace_variables_capture_avoiding(T& x,
                                        Substitution& sigma,
                                        data::set_identifier_generator& id_generator
)
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::capture_avoiding_substitution_updater<Substitution> sigma1(sigma, id_generator);
  data::detail::apply_replace_capture_avoiding_variables_builder<NAMESPACE::data_expression_builder, NAMESPACE::detail::add_capture_avoiding_replacement>(sigma1).update(x);
}

/// \\\\brief Applies sigma as a capture avoiding substitution to x.
/// \\\\param x The object to which the substiution is applied.
/// \\\\param sigma A substitution.
/// \\\\param id_generator An identifier generator that generates names that do not appear in x and sigma
template <typename T, data::IsSubstitution Substitution>
T replace_variables_capture_avoiding(const T& x,
                                     Substitution& sigma,
                                     data::set_identifier_generator& id_generator
)
  requires std::is_base_of_v<atermpp::aterm, T>
{
  data::detail::capture_avoiding_substitution_updater<Substitution> sigma1(sigma, id_generator);
  T result;
  data::detail::apply_replace_capture_avoiding_variables_builder<NAMESPACE::data_expression_builder, NAMESPACE::detail::add_capture_avoiding_replacement>(sigma1).apply(result, x);
  return result;
}

/// \\\\brief Applies sigma as a capture avoiding substitution to x.
/// \\\\param x The object to which the subsitution is applied.
/// \\\\param sigma A substitution.
template <typename T, data::IsSubstitution Substitution>
void replace_variables_capture_avoiding(T& x,
                                        Substitution& sigma
)
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::set_identifier_generator id_generator;
  id_generator.add_identifiers(NAMESPACE::find_identifiers(x));
  for (const data::variable& v: substitution_variables(sigma))
  {
    id_generator.add_identifier(v.name());
  }
  NAMESPACE::replace_variables_capture_avoiding(x, sigma, id_generator);
}

/// \\\\brief Applies sigma as a capture avoiding substitution to x.
/// \\\\param x The object to which the substiution is applied.
/// \\\\param sigma A substitution.
template <typename T, data::IsSubstitution Substitution>
T replace_variables_capture_avoiding(const T& x,
                                     Substitution& sigma
)
  requires std::is_base_of_v<atermpp::aterm, T>
{
  data::set_identifier_generator id_generator;
  id_generator.add_identifiers(NAMESPACE::find_identifiers(x));
  for (const data::variable& v: substitution_variables(sigma))
  {
    id_generator.add_identifier(v.name());
  }
  return NAMESPACE::replace_variables_capture_avoiding(x, sigma, id_generator);
}
'''

REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT: str = ''' /// \\\\brief Applies sigma as a capture avoiding substitution to x using an identifier generator.
/// \\\\details This substitution function is much faster than replace_variables_capture_avoiding, but
///          it requires an identifier generator that generates strings for fresh variables. These
///          strings must be unique in the sense that they have not been used for other variables.
/// \\\\param x The object to which the subsitution is applied.
/// \\\\param sigma A mutable substitution of which it can efficiently be checked whether a variable occurs in its
///              right hand side. The class maintain_variables_in_rhs is useful for this purpose.
/// \\\\param id_generator A generator that generates unique strings, not yet used as variable names.

template <typename T, data::IsSubstitution Substitution, typename IdentifierGenerator>
void replace_variables_capture_avoiding_with_an_identifier_generator(T& x,
                       Substitution& sigma,
                       IdentifierGenerator& id_generator
                      )
  requires (!std::is_base_of_v<atermpp::aterm, T>)
{
  data::detail::apply_replace_capture_avoiding_variables_builder_with_an_identifier_generator<NAMESPACE::data_expression_builder, NAMESPACE::detail::add_capture_avoiding_replacement_with_an_identifier_generator>(sigma, id_generator).update(x);
}

/// \\\\brief Applies sigma as a capture avoiding substitution to x using an identifier generator..
/// \\\\details This substitution function is much faster than replace_variables_capture_avoiding, but
///          it requires an identifier generator that generates strings for fresh variables. These
///          strings must be unique in the sense that they have not been used for other variables.
/// \\\\param x The object to which the substiution is applied.
/// \\\\param sigma A mutable substitution of which it can efficiently be checked whether a variable occurs in its
///              right hand side. The class maintain_variables_in_rhs is useful for this purpose.
/// \\\\param id_generator A generator that generates unique strings, not yet used as variable names.
/// \\\\return The result is the term x to which sigma has been applied.
template <typename T, data::IsSubstitution Substitution, typename IdentifierGenerator>
T replace_variables_capture_avoiding_with_an_identifier_generator(const T& x,
                    Substitution& sigma,
                    IdentifierGenerator& id_generator
                   )
  requires std::is_base_of_v<atermpp::aterm, T>
{
  T result;
  data::detail::apply_replace_capture_avoiding_variables_builder_with_an_identifier_generator<NAMESPACE::data_expression_builder, NAMESPACE::detail::add_capture_avoiding_replacement_with_an_identifier_generator>(sigma, id_generator).apply(result, x);
  return result;
}
'''

FIND_VARIABLES_FUNCTION_TEXT: str = '''/// \\brief Writes all variables that occur in an object to an output iterator.
/// \\param[in] x an object containing variables
/// \\param[in,out] o an output iterator to which all variables occurring in x are written.
template <typename T, typename OutputIterator>
void find_all_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_all_variables_traverser<NAMESPACE::variable_traverser>(o).apply(x);
}

/// \\brief Returns all variables that occur in an object.
/// \\param[in] x an object containing variables
/// \\return All variables that occur in the object x
template <typename T>
std::set<data::variable> find_all_variables(const T& x)
{
  std::set<data::variable> result;
  NAMESPACE::find_all_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \\brief Writes all free variables that occur in an object to an output iterator.
/// \\param[in] x an object containing variables
/// \\param[in,out] o an output iterator to which all variables occurring in x are added.
template <typename T, typename OutputIterator>
void find_free_variables(const T& x, OutputIterator o)
{
  data::detail::make_find_free_variables_traverser<NAMESPACE::data_expression_traverser, NAMESPACE::add_data_variable_traverser_binding>(o).apply(x);
}

/// \\brief Writes all free variables that occur in an object to an output iterator.
/// \\param[in] x an object containing variables
/// \\param[in,out] o an output iterator to which all variables occurring in x are written.
/// \\param[in] bound a container of variables
template <typename T, typename OutputIterator, typename VariableContainer>
void find_free_variables_with_bound(const T& x, OutputIterator o, const VariableContainer& bound)
{
  data::detail::make_find_free_variables_traverser<NAMESPACE::data_expression_traverser, NAMESPACE::add_data_variable_traverser_binding>(o, bound).apply(x);
}

/// \\brief Returns all free variables that occur in an object.
/// \\param[in] x an object containing variables
/// \\return All free variables that occur in the object x
template <typename T>
std::set<data::variable> find_free_variables(const T& x)
{
  std::set<data::variable> result;
  NAMESPACE::find_free_variables(x, std::inserter(result, result.end()));
  return result;
}

/// \\brief Returns all free variables that occur in an object.
/// \\param[in] x an object containing variables
/// \\param[in] bound a container of variables
/// \\return All free variables that occur in the object x
template <typename T, typename VariableContainer>
std::set<data::variable> find_free_variables_with_bound(const T& x, VariableContainer const& bound)
{
  std::set<data::variable> result;
  NAMESPACE::find_free_variables_with_bound(x, std::inserter(result, result.end()), bound);
  return result;
}

/// \\brief Writes all identifiers that occur in an object to an output iterator.
/// \\param[in] x an object containing identifiers
/// \\param[in,out] o an output iterator to which all identifiers occurring in x are written.
template <typename T, typename OutputIterator>
void find_identifiers(const T& x, OutputIterator o)
{
  data::detail::make_find_identifiers_traverser<NAMESPACE::identifier_string_traverser>(o).apply(x);
}

/// \\brief Returns all identifiers that occur in an object.
/// \\param[in] x an object containing identifiers
/// \\return All identifiers that occur in the object x
template <typename T>
std::set<core::identifier_string> find_identifiers(const T& x)
{
  std::set<core::identifier_string> result;
  NAMESPACE::find_identifiers(x, std::inserter(result, result.end()));
  return result;
}

/// \\brief Writes all sort expressions that occur in an object to an output iterator.
/// \\param[in] x an object containing sort expressions
/// \\param[in,out] o an output iterator to which all sort expressions occurring in x are written.
template <typename T, typename OutputIterator>
void find_sort_expressions(const T& x, OutputIterator o)
{
  data::detail::make_find_sort_expressions_traverser<NAMESPACE::sort_expression_traverser>(o).apply(x);
}

/// \\brief Returns all sort expressions that occur in an object.
/// \\param[in] x an object containing sort expressions
/// \\return All sort expressions that occur in the object x
template <typename T>
std::set<data::sort_expression> find_sort_expressions(const T& x)
{
  std::set<data::sort_expression> result;
  NAMESPACE::find_sort_expressions(x, std::inserter(result, result.end()));
  return result;
}

/// \\brief Writes all function symbols that occur in an object to an output iterator.
/// \\param[in] x an object containing function symbols
/// \\param[in,out] o an output iterator to which all function symbols occurring in x are written.
template <typename T, typename OutputIterator>
void find_function_symbols(const T& x, OutputIterator o)
{
  data::detail::make_find_function_symbols_traverser<NAMESPACE::data_expression_traverser>(o).apply(x);
}

/// \\brief Returns all function symbols that occur in an object.
/// \\param[in] x an object containing function symbols
/// \\return All function symbols that occur in the object x
template <typename T>
std::set<data::function_symbol> find_function_symbols(const T& x)
{
  std::set<data::function_symbol> result;
  NAMESPACE::find_function_symbols(x, std::inserter(result, result.end()));
  return result;
}
'''

@typechecked
def generate_code(filename: str, namespace: str, label: str, text: str) -> bool:
    text = re.sub('NAMESPACE', namespace, text)
    return insert_text_in_file(filename, text, f'generated {namespace} {label} code')

@typechecked
def generate_rewrite_functions() -> bool:
    result: bool = True
    result = generate_code(MCRL2_ROOT + 'libraries/data/include/mcrl2/data/rewrite.h'                  , 'data'            , 'rewrite', REWRITE_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/lps/include/mcrl2/lps/rewrite.h'                    , 'lps'             , 'rewrite', REWRITE_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/rewrite.h', 'action_formulas' , 'rewrite', REWRITE_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/rewrite.h', 'regular_formulas', 'rewrite', REWRITE_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/rewrite.h', 'state_formulas'  , 'rewrite', REWRITE_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pbes/include/mcrl2/pbes/rewrite.h'                  , 'pbes_system'     , 'rewrite', REWRITE_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pres/include/mcrl2/pres/rewrite.h'                  , 'pres_system'     , 'rewrite', REWRITE_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/process/include/mcrl2/process/rewrite.h'            , 'process'         , 'rewrite', REWRITE_TEXT) and result
    return result

@typechecked
def generate_replace_functions() -> bool:
    result: bool = True
    result = generate_code(MCRL2_ROOT + 'libraries/data/include/mcrl2/data/replace.h'                  , 'data'            , 'replace', SUBSTITUTE_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/lps/include/mcrl2/lps/replace.h'                    , 'lps'             , 'replace', SUBSTITUTE_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace.h', 'action_formulas' , 'replace', SUBSTITUTE_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace.h', 'regular_formulas', 'replace', SUBSTITUTE_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace.h', 'state_formulas'  , 'replace', SUBSTITUTE_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pbes/include/mcrl2/pbes/replace.h'                  , 'pbes_system'     , 'replace', SUBSTITUTE_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pres/include/mcrl2/pres/replace.h'                  , 'pres_system'     , 'replace', SUBSTITUTE_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/process/include/mcrl2/process/replace.h'            , 'process'         , 'replace', SUBSTITUTE_FUNCTION_TEXT) and result
    return result

@typechecked
def generate_replace_capture_avoiding_functions() -> bool:
    result: bool = True
    result = generate_code(MCRL2_ROOT + 'libraries/data/include/mcrl2/data/replace_capture_avoiding.h'                  , 'data'            , 'replace_capture_avoiding', REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/lps/include/mcrl2/lps/replace_capture_avoiding.h'                    , 'lps'             , 'replace_capture_avoiding', REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace_capture_avoiding.h', 'action_formulas' , 'replace_capture_avoiding', REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace_capture_avoiding.h', 'regular_formulas', 'replace_capture_avoiding', REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace_capture_avoiding.h', 'state_formulas'  , 'replace_capture_avoiding', REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pbes/include/mcrl2/pbes/replace_capture_avoiding.h'                  , 'pbes_system'     , 'replace_capture_avoiding', REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pres/include/mcrl2/pres/replace_capture_avoiding.h'                  , 'pres_system'     , 'replace_capture_avoiding', REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/process/include/mcrl2/process/replace_capture_avoiding.h'            , 'process'         , 'replace_capture_avoiding', REPLACE_CAPTURE_AVOIDING_FUNCTION_TEXT) and result
    return result

@typechecked
def generate_replace_capture_avoiding_with_identifier_generator_functions() -> bool:
    result: bool = True
    result = generate_code(MCRL2_ROOT + 'libraries/data/include/mcrl2/data/replace_capture_avoiding_with_an_identifier_generator.h'                  , 'data'            , 'replace_capture_avoiding_with_identifier_generator', REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/lps/include/mcrl2/lps/replace_capture_avoiding_with_an_identifier_generator.h'                    , 'lps'             , 'replace_capture_avoiding_with_identifier_generator', REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace_capture_avoiding_with_an_identifier_generator.h', 'action_formulas' , 'replace_capture_avoiding_with_identifier_generator', REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace_capture_avoiding_with_an_identifier_generator.h', 'regular_formulas', 'replace_capture_avoiding_with_identifier_generator', REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/replace_capture_avoiding_with_an_identifier_generator.h', 'state_formulas'  , 'replace_capture_avoiding_with_identifier_generator', REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pbes/include/mcrl2/pbes/replace_capture_avoiding_with_an_identifier_generator.h'                  , 'pbes_system'     , 'replace_capture_avoiding_with_identifier_generator', REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pres/include/mcrl2/pres/replace_capture_avoiding_with_an_identifier_generator.h'                  , 'pres_system'     , 'replace_capture_avoiding_with_identifier_generator', REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/process/include/mcrl2/process/replace_capture_avoiding_with_an_identifier_generator.h'            , 'process'         , 'replace_capture_avoiding_with_identifier_generator', REPLACE_CAPTURE_AVOIDING_WITH_IDENTIFIER_GENERATOR_FUNCTION_TEXT) and result
    return result

@typechecked
def generate_find_functions() -> bool:
    result: bool = True
    result = generate_code(MCRL2_ROOT + 'libraries/data/include/mcrl2/data/find.h'                  , 'data'            , 'find', FIND_VARIABLES_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/lps/include/mcrl2/lps/find.h'                    , 'lps'             , 'find', FIND_VARIABLES_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/find.h', 'action_formulas' , 'find', FIND_VARIABLES_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/find.h', 'regular_formulas', 'find', FIND_VARIABLES_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/modal_formula/include/mcrl2/modal_formula/find.h', 'state_formulas'  , 'find', FIND_VARIABLES_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pbes/include/mcrl2/pbes/find.h'                  , 'pbes_system'     , 'find', FIND_VARIABLES_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/pres/include/mcrl2/pres/find.h'                  , 'pres_system'     , 'find', FIND_VARIABLES_FUNCTION_TEXT) and result
    result = generate_code(MCRL2_ROOT + 'libraries/process/include/mcrl2/process/find.h'            , 'process'         , 'find', FIND_VARIABLES_FUNCTION_TEXT) and result
    return result

@typechecked
def main() -> int:
    result: bool = True
    result = generate_rewrite_functions() and result
    result = generate_replace_functions() and result
    result = generate_replace_capture_avoiding_functions() and result
    result = generate_replace_capture_avoiding_with_identifier_generator_functions() and result
    result = generate_find_functions() and result
    return int(not result) # 0 result indicates successful execution

if __name__ == "__main__":
    sys.exit(main())
