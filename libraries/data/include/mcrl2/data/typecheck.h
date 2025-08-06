// Author(s): Jan Friso Groote, Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/typecheck.h
/// \brief add your file description here.

#ifndef MCRL2_DATA_TYPECHECK_H
#define MCRL2_DATA_TYPECHECK_H

#include <iostream>

#include "mcrl2/atermpp/aterm_io.h"
#include "mcrl2/data/data_specification.h"
#include "mcrl2/data/detail/variable_context.h"
#include "mcrl2/data/sort_type_checker.h"

namespace mcrl2::data
{

class data_type_checker: public sort_type_checker
{
  protected:
    mutable bool was_warning_upcasting = false; // This variable is used to limit the number of upcasting warnings.

    std::map<core::identifier_string,sort_expression_list> system_constants;   //name -> Set(sort expression)
    std::map<core::identifier_string,sort_expression_list> system_functions;   //name -> Set(sort expression)
    std::map<core::identifier_string,sort_expression> user_constants;          //name -> sort expression
    std::map<core::identifier_string,sort_expression_list> user_functions;     //name -> Set(sort expression)
    data_specification type_checked_data_spec;

  public:
    /** \brief     make a data type checker.
     *             Throws a mcrl2::runtime_error exception if the data_specification is not well typed.
     *  \param[in] data_spec A data specification that does not need to have been type checked.
     *  \return    A data expression where all untyped identifiers have been replace by typed ones.
     **/
    data_type_checker(const data_specification& data_spec);

    /** \brief     Type checks a variable.
     *             Throws an mcrl2::runtime_error exception if the variable is not well typed.
     *  \details   A variable is not well typed if its name clashes with the name of a declared function, when its sort does not exist, or when
     *             the variable is used in its context with a different sort.
     *  \param[in] v A variables that is to be type checked.
     *  \param[in] context Information about the context of the variable.
     **/
    void operator()(const variable& v, const detail::variable_context& context) const;

    /** \brief     Type checks a variable list.
     *             Throws an mcrl2::runtime_error exception if the variables are not well typed.
     *  \details   A variable is not well typed if its name clashes with the name of a declared function, when its sort does not exist, or when
     *             a variable occurs in the context. Furthermore, variables cannot occur multiple times in a variable list.
     *  \param[in] l A list of variables that must be type checked.
     *  \param[in] context Information about the context of the variables in the list.
     **/
    void operator()(const variable_list& l, const detail::variable_context& context) const;

    /** \brief     Yields a type checked data specification, provided typechecking was successful.
     *             If not successful an exception is thrown.
     *  \return    a data specification where all untyped identifiers have been replace by typed ones.
     **/
    const data_specification operator()() const;

    /** \brief     Yields a type checked equation list, and sets the types in the equations right.
     *             If not successful an exception is thrown.
     *  \param[in] eqns The list of equations that is type checked and updated.
     **/
    void operator()(data_equation_vector& eqns);

    data_expression typecheck_data_expression(const data_expression& x,
                                              const sort_expression& expected_sort,
                                              const detail::variable_context& variable_context
                                             )
    {
      data_expression x1 = x;
      TraverseVarConsTypeD(variable_context, x1, expected_sort);
      x1 = data::normalize_sorts(x1, get_sort_specification());
      if (x1.sort() != expected_sort)
      {
        x1 = upcast_numeric_type(x1, expected_sort, variable_context);
      }
      return x1;
    }

    assignment typecheck_assignment(const assignment& x, const detail::variable_context& variable_context)
    {
      sort_type_checker::check_sort_is_declared(x.lhs().sort());
      data_expression rhs = typecheck_data_expression(x.rhs(), x.lhs().sort(), variable_context);
      return assignment(x.lhs(), rhs);
    }

    assignment_list typecheck_assignment_list(const assignment_list& assignments, const detail::variable_context& variable_context)
    {
      // check for name clashes
      std::set<core::identifier_string> names;
      for (const assignment& a: assignments)
      {
        const core::identifier_string& name = a.lhs().name();
        if (names.find(name) != names.end())
        {
          throw mcrl2::runtime_error("duplicate variable names in assignments: " + data::pp(assignments) + ")");
        }
        names.insert(name);
      }

      // typecheck the assignments
      assignment_vector result;
      for (const assignment& a: assignments)
      {
        result.push_back(typecheck_assignment(a, variable_context));
      }
      return assignment_list(result.begin(), result.end());
    }

    const data_specification& typechecked_data_specification() const
    {
      return type_checked_data_spec;
    }

    void print_context() const
    {
      auto const& sortspec = get_sort_specification();
      std::cout << "--- basic sorts ---" << std::endl;
      for (auto const& x: sortspec.user_defined_sorts())
      {
        std::cout << x << std::endl;
      }
      std::cout << "--- aliases ---" << std::endl;
      for (auto const& x: sortspec.user_defined_aliases())
      {
        std::cout << x << std::endl;
      }
      std::cout << "--- user constants ---" << std::endl;
      for (const auto& user_constant: user_constants)
      {
        std::cout << user_constant.first << " -> " << user_constant.second << std::endl;
      }
      std::cout << "--- user functions ---" << std::endl;
      for (const auto& user_function: user_functions)
      {
        std::cout << user_function.first << " -> " << user_function.second << std::endl;
      }
    }

  protected:
    /** \brief     Type check a data expression.
     *  Throws a mcrl2::runtime_error exception if the expression is not well typed.
     *  \param[in] data_expr A data expression that has not been type checked.
     *  \param[in] context The variable context in which this term must be typechecked.
     *  \return    a data expression where all untyped identifiers have been replace by typed ones.
     **/
    data_expression operator()(const data_expression& data_expr,
                               const detail::variable_context& context) const;

    void read_sort(const sort_expression& SortExpr);
    void read_constructors_and_mappings(const function_symbol_vector& constructors, const function_symbol_vector& mappings, const function_symbol_vector& normalized_constructors);
    void add_function(const data::function_symbol& f, const std::string& msg, bool allow_double_decls=false);
    void add_constant(const data::function_symbol& f, const std::string& msg);
    void initialise_system_defined_functions();
    void add_system_constant(const data::function_symbol& f);
    void add_system_function(const data::function_symbol& f);
    void add_system_constants_and_functions(const std::vector<data::function_symbol>& v);
    bool TypeMatchA(const sort_expression& Type_in, const sort_expression& PosType_in, sort_expression& result) const;
    bool TypeMatchL(const sort_expression_list& TypeList, const sort_expression_list& PosTypeList, sort_expression_list& result) const;
    sort_expression UnwindType(const sort_expression& Type) const;
    variable UnwindType(const variable& v) const;
    template <class T>
    atermpp::term_list<T> UnwindType(const atermpp::term_list<T>& l)
    {
      std::vector<T> result;
      for(typename atermpp::term_list<T>::const_iterator i=l.begin(); i!=l.end(); ++i)
      {
        result.push_back(UnwindType(*i));
      }
      return atermpp::term_list<T>(result.begin(),result.end());
    }

    sort_expression TraverseVarConsTypeD(
                        const detail::variable_context& DeclaredVars,
                        data_expression& DataTerm,
                        const sort_expression& PosType,
                        const bool strictly_ambiguous=true,
                        const bool warn_upcasting=false,
                        const bool print_cast_error=true) const;

    /* sort_expression TraverseVarConsTypeD(const std::map<core::identifier_string,sort_expression>& DeclaredVars,
                                         data_expression& t1,
                                         const sort_expression& t2); */

    sort_expression TraverseVarConsTypeDN(
                           const detail::variable_context& DeclaredVars,
                           data_expression& DataTerm,
                           sort_expression PosType,
                           const bool strictly_ambiguous=true,
                           const std::size_t nFactPars=std::string::npos,
                           const bool warn_upcasting=false,
                           const bool print_cast_error=true) const;

    bool InTypesA(const sort_expression& Type, sort_expression_list Types) const;
    bool EqTypesA(const sort_expression& Type1, const sort_expression& Type2) const;
    bool InTypesL(const sort_expression_list& Type, atermpp::term_list<sort_expression_list> Types) const;
    bool EqTypesL(sort_expression_list Type1, sort_expression_list Type2) const;
    bool MaximumType(const sort_expression& Type1, const sort_expression& Type2, sort_expression& result) const;
    sort_expression ExpandNumTypesUp(sort_expression Type) const;
    sort_expression_list ExpandNumTypesUpL(const sort_expression_list& type_list) const;
    sort_expression ExpandNumTypesDown(sort_expression Type) const;
    bool UnifyMinType(const sort_expression& Type1, const sort_expression& Type2, sort_expression& result) const;
    sort_expression determine_allowed_type(const data_expression& d, const sort_expression& proposed_type) const;
    bool MatchIf(const function_sort& type, sort_expression& result) const;
    bool MatchEqNeqComparison(const function_sort& type, sort_expression& result) const;
    bool MatchSqrt(const function_sort& type, sort_expression& result) const;
    bool MatchListOpCons(const function_sort& type, sort_expression& result) const;
    bool MatchListOpSnoc(const function_sort& type, sort_expression& result) const;
    bool MatchListOpConcat(const function_sort& type, sort_expression& result) const;
    bool MatchListOpEltAt(const function_sort& type, sort_expression& result) const;
    bool MatchListOpHead(const function_sort& type, sort_expression& result) const;
    bool MatchListOpTail(const function_sort& type, sort_expression& result) const;
    bool MatchSetOpSet2Bag(const function_sort& type, sort_expression& result) const;
    bool MatchFalseFunction(const function_sort& type, sort_expression& result) const;
    bool MatchListSetBagOpIn(const function_sort& type, sort_expression& result) const;
    bool match_fset_insert(const function_sort& type, sort_expression& result) const;
    bool match_fset_fbag_pick(const function_sort& type, sort_expression& result) const;
    bool match_fbag_cinsert(const function_sort& type, sort_expression& result) const;
    bool UnifyElementSort(sort_expression& Arg1, sort_expression& Arg2, sort_expression& result) const;
    bool MatchSetBagOpUnionDiffIntersect(const core::identifier_string& data_term_name, const function_sort& type, sort_expression& result) const;
    bool MatchSetOpSetCompl(const function_sort& type, sort_expression& result) const;
    bool MatchBagOpBag2Set(const function_sort& type, sort_expression& result) const;
    bool MatchBagOpBagCount(const function_sort& type, sort_expression& result) const;
    bool MatchFuncUpdate(const function_sort& type, sort_expression& result) const;
    bool MatchSetConstructor(const function_sort& type, sort_expression& result) const;
    bool MatchBagConstructor(const function_sort& type, sort_expression& result) const;
    bool UnArrowProd(const sort_expression_list& ArgTypes, sort_expression PosType, sort_expression& result) const;
    bool UnFSet(sort_expression PosType, sort_expression& result) const;
    bool UnFBag(sort_expression PosType, sort_expression& result) const;
    bool UnList(sort_expression PosType, sort_expression& result) const;
    void ErrorMsgCannotCast(sort_expression CandidateType, data_expression_list Arguments, sort_expression_list ArgumentTypes,std::string previous_reason) const;
    sort_expression UpCastNumericType(
                    sort_expression NeededType,
                    sort_expression Type,
                    data_expression& Par,
                    const detail::variable_context& DeclaredVars,
                    const bool strictly_ambiguous,
                    bool warn_upcasting=false,
                    const bool print_cast_error=false) const;
    void TransformVarConsTypeData(data_specification& data_spec);
    sort_expression_list GetNotInferredList(const atermpp::term_list<sort_expression_list>& TypeListList) const;
    sort_expression_list InsertType(const sort_expression_list& TypeList, const sort_expression& Type) const;
    std::pair<bool,sort_expression_list> AdjustNotInferredList(
            const sort_expression_list& PosTypeList,
            const atermpp::term_list<sort_expression_list>& TypeListList) const;
    bool IsTypeAllowedA(const sort_expression& Type, const sort_expression& PosType) const;
    bool IsTypeAllowedL(const sort_expression_list& TypeList, const sort_expression_list& PosTypeList) const;
    bool IsNotInferredL(sort_expression_list TypeList) const;
    bool strict_type_check(const data_expression& d) const;

    // for example Pos -> Nat, or Nat -> Int
    data_expression upcast_numeric_type(const data_expression& x,
                                        const sort_expression& expected_sort,
                                        const detail::variable_context& variable_context
                                       )
    {
      try
      {
        data_expression x1 = x;
        UpCastNumericType(expected_sort, x.sort(), x1, variable_context, false, false, false);
        return data::normalize_sorts(x1, get_sort_specification());
      }
      catch (mcrl2::runtime_error& e)
      {
        throw mcrl2::runtime_error(std::string(e.what()) + "\ncannot (up)cast " + data::pp(x) + " to type " + data::pp(expected_sort));
      }
    }

};

/** \brief     Type check a sort expression.
 *  Throws an exception if something went wrong.
 *  \param[in] sort_expr A sort expression that has not been type checked.
 *  \param[in] data_spec The data specification to use as context.
 *  \post      sort_expr is type checked.
 **/
inline
void typecheck_sort_expression(const sort_expression& sort_expr, const data_specification& data_spec)
{
  try
  {
    // sort_type_checker type_checker(data_spec.user_defined_sorts(), data_spec.user_defined_aliases());
    sort_type_checker type_checker(data_spec);
    type_checker(sort_expr);
  }
  catch (mcrl2::runtime_error& e)
  {
    throw mcrl2::runtime_error(std::string(e.what()) + "\nCould not type check sort " + pp(sort_expr));
  }
}

/** \brief     Type check a data expression.
 *  Throws an exception if something went wrong.
 *  \param[in] x A data expression that has not been type checked.
 *  \param[in] variables A container with variables that can occur in the data expression.
 *  \param[in] dataspec The data specification that is used for type checking.
 *  \post      data_expr is type checked.
 **/
template <typename VariableContainer>
data_expression typecheck_data_expression(const data_expression& x,
                                          const VariableContainer& variables,
                                          const data_specification& dataspec = data_specification()
                                         )
{
  try
  {
    data_type_checker typechecker(dataspec);
    detail::variable_context variable_context;
    variable_context.add_context_variables(variables, typechecker);
    return typechecker.typecheck_data_expression(x, untyped_sort(), variable_context);
  }
  catch (mcrl2::runtime_error& e)
  {
    throw mcrl2::runtime_error(std::string(e.what()) + "\nCould not type check data expression " + data::pp(x));
  }
}

/** \brief     Type check a data expression.
 *  Throws an exception if something went wrong.
 *  \param[in] x A data expression that has not been type checked.
 *  \param[in] dataspec Data specification to be used as context.
 *  \post      data_expr is type checked.
 **/
inline
data_expression typecheck_data_expression(const data_expression& x, const data_specification& dataspec = data_specification())
{
  return typecheck_data_expression(x, variable_list(), dataspec);
}

/** \brief     Type check a parsed mCRL2 data specification.
 *  Throws an exception if something went wrong.
 *  \param[in] data_spec A data specification that has not been type checked.
 *  \post      data_spec is type checked.
 **/
inline
void typecheck_data_specification(data_specification& data_spec)
{
  data_type_checker type_checker(data_spec);
  data_spec=type_checker();
}

inline
data_expression typecheck_untyped_data_parameter(data_type_checker& typechecker,
                                                 const core::identifier_string& name,
                                                 const data_expression_list& parameters,
                                                 const data::sort_expression& expected_sort,
                                                 const detail::variable_context& variable_context
                                                )
{
  if (parameters.empty())
  {
    return typechecker.typecheck_data_expression(untyped_identifier(name), expected_sort, variable_context);
  }
  else
  {
    return typechecker.typecheck_data_expression(application(untyped_identifier(name), parameters), expected_sort, variable_context);
  }
}

} // namespace mcrl2::data

#endif // MCRL2_DATA_TYPECHECK_H
