// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_PBES_DETAIL_PBES_REMOVE_COUNTEREXAMPLE_INFO_H
#define MCRL2_PBES_DETAIL_PBES_REMOVE_COUNTEREXAMPLE_INFO_H

#include "mcrl2/core/identifier_string.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/pbes_equation.h"
#include "mcrl2/pbes/builder.h"
#include "mcrl2/pbes/pbes_expression.h"
#include "mcrl2/pbes/rewriters/simplify_rewriter.h"

#include <regex>
#include <string>
 
namespace mcrl2::pbes_system::detail
{  
  
/// \cond INTERNAL_DOCS
static std::regex positive("Zpos_(\\d+)_.*");
static std::regex negative("Zneg_(\\d+)_.*");
static std::regex positive_or_negative("Z(neg|pos)_(\\d+)_.*");
/// \endcond

/// \brief Returns true iff the name is Zpos, name must be a counter example name.
inline
bool is_counter_example_positive(const core::identifier_string& name)
{
  std::smatch match;
  std::string name_str = name;
  return std::regex_match(name_str, match, positive);
}

inline
bool is_counter_example_name(const core::identifier_string& name)
{
  std::smatch match;
  std::string name_str = name;
  return std::regex_match(name_str, match, positive_or_negative);
}

inline
bool is_counter_example_variable(const propositional_variable& X)
{
  return is_counter_example_name(X.name());
}

/// \brief Guesses if the PBES equation is a counter example equation
inline
bool is_counter_example_equation(const pbes_equation& equation)
{
  return is_counter_example_variable(equation.variable());
}

/// \brief Guesses if the PBES variable instantiation is for counter example equation
inline
bool is_counter_example_instantiation(const propositional_variable_instantiation& inst)
{
  return is_counter_example_name(inst.name());
}

/// \brief Guesses if a pbes has counter example information
inline
bool has_counter_example_information(const pbes& pbesspec)
{
  for (const pbes_equation& eqn: pbesspec.equations())
  {
    if (is_counter_example_equation(eqn))
    {
      return true;
    }
  }
  return false;
}

struct subsitute_counterexample: public pbes_expression_builder<subsitute_counterexample>
{
  using super = pbes_expression_builder<subsitute_counterexample>;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  bool replace_Lplus;
  bool replace_Lminus;
    

  subsitute_counterexample(bool replace_Lplus, bool replace_Lminus)
    : replace_Lplus(replace_Lplus),
      replace_Lminus(replace_Lminus)
  {}

  template <class T>
  void apply(T& result, const propositional_variable_instantiation& x)
  {
    std::smatch match;
    if (replace_Lplus && std::regex_match(static_cast<const std::string&>(x.name()), match, positive))
    {
      result = true_();
    }
    else if (replace_Lminus && std::regex_match(static_cast<const std::string&>(x.name()), match, negative))
    {
      result = false_();
    }
    else
    {
      result = x;
    }
  }
};

/// Removes all equations and expressions related to counter examples from the input PBES.
inline 
mcrl2::pbes_system::pbes remove_counterexample_info(const pbes_system::pbes& pbes, bool remove_Lplus = true, bool remove_Lminus = true, bool remove_Lequations = true)
{
  // Remove all equations related to the countexamples
  std::vector<mcrl2::pbes_system::pbes_equation> equations;
  subsitute_counterexample f(remove_Lplus, remove_Lminus);

  /// Replace the corresponding PBES variables by true and false respectively.
  simplify_rewriter simplify;
  for (auto& equation : pbes.equations())
  {
    std::smatch match;
    bool Lplus = std::regex_match(static_cast<const std::string&>(equation.variable().name()), match, positive);
    bool Lminus = std::regex_match(static_cast<const std::string&>(equation.variable().name()), match, negative);

    if (!Lplus && !Lminus)
    {
      pbes_expression expression;
      f.apply(expression, equation.formula());

      // TODO: Technically simplify might result in too much rewriting taking place.
      equations.emplace_back(equation.symbol(), equation.variable(), simplify(expression));
    }

    if (!remove_Lequations && is_counter_example_equation(equation))
    {
      equations.emplace_back(equation);
    }
    else
    {
      if (Lplus && !remove_Lplus)
      {
        equations.emplace_back(equation);
      }

      if (Lminus && !remove_Lminus)
      {      
        equations.emplace_back(equation);
      }
    }
  }

  return mcrl2::pbes_system::pbes(pbes.data(), equations, pbes.initial_state());
}


} // mcrl2::pbes_system::detail

#endif // MCRL2_PBES_DETAIL_PBES_REMOVE_COUNTEREXAMPLE_INFO_H