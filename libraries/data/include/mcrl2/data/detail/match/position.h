// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_DATA_DETAIL_POSITION_H
#define MCRL2_DATA_DETAIL_POSITION_H

#include "mcrl2/atermpp/aterm_appl.h"
#include "mcrl2/data/replace.h"
#include "mcrl2/data/substitutions/mutable_map_substitution.h"
#include "mcrl2/data/variable.h"

#include <set>
#include <vector>
#include <sstream>

namespace mcrl2::data::detail
{

/// \brief A position in the term.
using position = std::vector<std::size_t>;

/// \brief Compute the set of positions of t such that t[p] is a variable for all p in fringe.
inline
void fringe_impl(const atermpp::aterm_appl& appl, position current, std::set<position>& fringe)
{
  if (is_variable(appl))
  {
    // The current position has a variable.
    fringe.insert(current);
  }
  else
  {
    // Extend the position to be one deeper into the subterm.
    current.emplace_back(0);
    for (const atermpp::aterm& argument : appl)
    {
      fringe_impl(static_cast<const atermpp::aterm_appl&>(argument), current, fringe);
      ++current.back();
    }
  }
}

/// \brief Compute the set of positions of t such that t[p] is a variable for all p in fringe(t).
inline
std::set<position> fringe(const atermpp::aterm_appl& appl)
{
  std::set<position> result;
  fringe_impl(appl, position(), result);
  return result;
}

/// \brief Returns t[pos], i.e., the term at the given position using a index to keep track of the pos.
inline
std::optional<data_expression> at_position_impl(const data_expression& t, const position& pos, std::size_t index)
{
  if (pos.empty())
  {
    // t[emptypos] = t
    return t;
  }
  else
  {
    std::size_t arg = pos[index];

    if (arg < t.size())
    {
      const data_expression& u = static_cast<const data_expression&>(static_cast<atermpp::aterm_appl>(t)[arg]);
      if (pos.size() == index + 1)
      {
        return u;
      }
      else
      {
        return at_position_impl(u, pos, index + 1);
      }
    }
  }

  return {};
}

/// \brief Returns t[pos], i.e., the term at the given position.
inline
std::optional<data_expression> at_position(const data_expression& t, const position& pos)
{
  return at_position_impl(t, pos, 0);
}

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

/// \brief Replace the position variable at the given position by the expression c.
inline
data_expression assign_at_position(const data_expression& term, const position& pos, const data_expression& c)
{
  mutable_map_substitution<std::map<variable, data_expression>> sigma;
  sigma[position_variable(pos)] = c;
  return replace_variables(term, sigma);
}

/// \returns True iff there exists l in L : (exists pos' <= pos : head(l[pos']) in V
bool has_variable_higher_impl(const data_expression& t, const position& pos, std::size_t index)
{
  // These two conditions check pos' <= pos.
  if (pos.empty() || index >= pos.size())
  {
    return false;
  }
  else if (is_variable(static_cast<data_expression>(t[pos[index]])))
  {
    return true;
  }
  else
  {
    assert(pos[index] < t.size());
    return has_variable_higher_impl(static_cast<data_expression>(t[pos[index]]), pos, index + 1);
  }
}

/// \returns True iff (exists pos' <= pos : head(l[pos']) in V
bool has_variable_higher(const data_expression& appl, const position& pos)
{
  return has_variable_higher_impl(appl, pos, 0);
}



} // namespace mcrl2::data::detail

#endif // MCRL2_DATA_DETAIL_POSITION_H
