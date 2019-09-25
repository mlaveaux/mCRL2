// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_CONSISTENCY_H
#define MCRL2_DATA_CONSISTENCY_H

#include "mcrl2/data/builder.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/replace.h"
#include "mcrl2/data/substitutions/mutable_indexed_substitution.h"
#include "mcrl2/data/variable.h"

#include <vector>
#include <set>

namespace mcrl2
{
namespace data
{
namespace detail
{

using equivalence_classes = std::vector<std::vector<variable>>;

template <template <typename> class Builder, typename Generator>
class linearize_builder : public Builder<linearize_builder<Builder, Generator>>
{
public:
  typedef Builder<linearize_builder<Builder, Generator>> super;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  linearize_builder(Generator& gen)
    : m_generator(gen)
  {}

  variable apply(const variable& var)
  {
    if (m_variables.find(var) != m_variables.end())
    {
      // The variable occurs in a different place.
      variable new_var(m_generator(), var.sort());

      // Add the new variable to the equivalence class and return it.
      m_equivalence_classes[var].emplace_back(new_var);
      return new_var;
    }
    else
    {
      m_variables.emplace(var);
    }

    return var;
  }

  equivalence_classes get_equivalence_classes()
  {
    // A set of sets (guaranteed no duplicates) of equivalence classes that must be checked for consistency.
    equivalence_classes result;

    for (auto& element : m_equivalence_classes)
    {
      result.emplace_back(element.second);
    }

    return result;
  }

private:
  std::set<variable> m_variables; ///< The set of variables that we have already seen.
  std::unordered_map<variable, std::vector<variable>> m_equivalence_classes; ///< For each original variable the set of variables that must be consistent.
  Generator& m_generator;
};

/// \brief Given an equation renames multiple occurrences of the same variable to different (unique) variables and returns the
///        equivalence class that must be checked for consistency.
template<typename Generator>
std::pair<data_equation, equivalence_classes> make_linear(const data_equation& equation, Generator& generator)
{
  // This makes the left-hand side linear, the right-hand side and condition can be the same as one instance of each variable did not change.
  linearize_builder<data_expression_builder, Generator> builder(generator);
  data_expression lhs = builder.apply(equation.lhs());

  // Obtain the new free variables (after renamings) to define the new equation.
  std::set<variable> variables = find_free_variables(lhs);

  // Generate new variables for each equation to ensure that they are unique.
  mutable_indexed_substitution<variable> sigma;
  for (auto& var : variables)
  {
    sigma[var] = variable(generator(), var.sort());
  }

  // Create the new data equation.
  data_equation linear_equation(variables, equation.condition(), lhs, equation.rhs());

  // Rename the meta variables in this equation to ensure uniqueness.
  data_equation renamed_equation = replace_variables(linear_equation, sigma);

  // Rename the variables in the same way in the equivalence classes.
  auto equivalence_classes = builder.get_equivalence_classes();
  for (auto& eq_class : equivalence_classes)
  {
    for (auto& var : variables)
    {
      std::replace(eq_class.begin(), eq_class.end(), var, static_cast<variable>(sigma(var)));
    }
  }

  return std::make_pair(renamed_equation, equivalence_classes);
}

/// \brief Check whether the given substitution sigma is consistent w.r.t. the given equivalence classes.
template<typename Substitution>
bool is_consistent(const equivalence_classes& classes, const Substitution& sigma)
{
   // We also need to check consistency of the matched rule.
   for (auto& equivalence_class : classes)
   {
     auto& subst = sigma(equivalence_class.front());
     for (auto& variable : equivalence_class)
     {
       if (sigma(variable) != subst)
       {
         return false;
       }
     }
   }

   return true;
}

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_CONSISTENCY_H
