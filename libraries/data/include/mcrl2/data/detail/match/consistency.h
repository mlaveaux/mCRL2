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

/// \brief The consistency class is a set of variables, must ensure unique entries for optimal performance..
using consistency_class = std::vector<variable>;

/// \brief The consistency partition is a set of consistency classes, the underlying vectors should be sets.
using consistency_partition = std::vector<consistency_class>;


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
      m_mapping[var].emplace_back(new_var);
      return new_var;
    }
    else
    {
      m_variables.emplace(var);
    }

    return var;
  }

  consistency_partition get_partition()
  {
    // A set of sets (guaranteed no duplicates) of equivalence classes that must be checked for consistency.
    consistency_partition result;

    for (const auto& [var, partition] : m_mapping)
    {
      result.emplace_back(partition);
    }

    return result;
  }

private:
  std::set<variable> m_variables; ///< The set of variables that we have already seen.
  std::unordered_map<variable, consistency_class> m_mapping; ///< For each original variable the set of variables that must be consistent.
  Generator& m_generator;
};

/// \brief Given an equation renames multiple occurrences of the same variable to different (unique) variables and returns the
///        equivalence class that must be checked for consistency.
template<typename Generator>
inline
std::pair<data_equation, consistency_partition> make_linear(const data_equation& equation, Generator& generator)
{
  // This makes the left-hand side linear, the right-hand side and condition can be the same as one instance of each variable did not change.
  linearize_builder<data_expression_builder, Generator> builder(generator);
  data_expression lhs = builder.apply(equation.lhs());

  // Obtain the new free variables (after renamings) to define the new equation.
  std::set<variable> variables = find_free_variables(lhs);

  // Generate new variables for each equation to ensure that they are unique.
  mutable_indexed_substitution<variable> sigma;
  for (const variable& var : variables)
  {
    sigma[var] = variable(generator(), var.sort());
  }

  // Create the new data equation.
  data_equation linear_equation(variables, equation.condition(), lhs, equation.rhs());

  // Rename the meta variables in this equation to ensure uniqueness.
  data_equation renamed_equation = replace_variables(linear_equation, sigma);

  // Rename the variables in the same way in the equivalence classes.
  consistency_partition partition = builder.get_partition();
  for (consistency_class& eq_class : partition)
  {
    for (const variable& var : variables)
    {
      std::replace(eq_class.begin(), eq_class.end(), var, static_cast<variable>(sigma(var)));
    }
  }

  return std::make_pair(renamed_equation, partition);
}

/// \brief Check whether the given substitution sigma is consistent w.r.t. the given equivalence classes.
template<typename Substitution>
inline
bool is_consistent(const consistency_partition& partition, const Substitution& sigma)
{
   // We also need to check consistency of the matched rule.
   for (const consistency_class& consistency_class : partition)
   {
     const auto& assigned = sigma(consistency_class.front());
     for (const variable& variable : consistency_class)
     {
       if (sigma(variable) != assigned)
       {
         return false;
       }
     }
   }

   return true;
}

using position = std::vector<std::size_t>;

/// \brief Print a position to a stream.
inline
std::ostream& operator<<(std::ostream& stream, const position& position)
{
  bool first = true;
  for (std::size_t index : position)
  {
    if (!first)
    {
      stream << "." << std::to_string(index);
    }
    else
    {
      stream << std::to_string(index);
    }
    first = false;
  }

  return stream;
}

/// \returns A position [1,2,3,4] as the string 1.2.3.4
inline
std::string to_variable_name(const position& position)
{
  std::stringstream result;
  result << "@" << position;
  return result.str();
}

/// \brief Create a variable named after the current position with no type.
inline
variable position_variable(const position& position)
{
  return variable(mcrl2::core::identifier_string(to_variable_name(position)), untyped_sort());
}

/// \brief Renames every variable to a unique name by using its position as identifier.
template<typename Substitution>
inline
void rename_variables_position(const atermpp::aterm_appl& appl, position current, Substitution& sigma)
{
  if (is_variable(appl))
  {
    // Rename the current variable to a corresponding position variable.
    sigma[static_cast<variable>(appl)] = position_variable(current);
  }
  else
  {
    // Extend the position to be one deeper into the subterm.
    current.emplace_back(0);
    for (const atermpp::aterm& argument : appl)
    {
      rename_variables_position(static_cast<const atermpp::aterm_appl&>(argument), current, sigma);
      ++current.back();
    }
  }
}

/// \brief Renames every variable to a unique name by using its position as identifier.
template<typename Substitution>
inline
void rename_variables_position(const atermpp::aterm_appl& appl, Substitution& sigma)
{
  rename_variables_position(appl, position(), sigma);
}

/// \brief Rename the variables in the data_equation such that at each position they are unique and adapt the equivalence classes.
inline
std::pair<data_equation, consistency_partition> rename_variables_unique(std::pair<data_equation, consistency_partition> rules)
{
  mutable_indexed_substitution<variable, variable> sigma;
  rename_variables_position(rules.first.lhs(), sigma);

  // Rename all variables in the resulting partition to the name indicated by sigma.
  consistency_partition result;
  for (const std::vector<variable>& set : rules.second)
  {
    std::vector<variable> new_vars;
    for (const variable& var : set)
    {
      new_vars.push_back(sigma(var));
    }
    result.push_back(new_vars);
  }

  return std::make_pair(replace_variables(rules.first, sigma), result);
}

} // namespace detail
} // namespace data
} // namespace mcrl2

#endif // MCRL2_DATA_CONSISTENCY_H
