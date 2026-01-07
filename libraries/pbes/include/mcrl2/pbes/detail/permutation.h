// Author(s): Maurice Laveaux and Menno Bartels
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/data/substitutions/mutable_map_substitution.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/replace.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/container/flat_map.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <ranges>

namespace mcrl2::pbes_system::detail
{

/// A representation of a permutation.
class permutation
{
public:
  permutation() = default;

  permutation(const boost::container::flat_map<std::size_t, std::size_t>& mapping)
    : m_mapping(mapping)
  {}

  /// Parse a permutation from a string of the shape x->y, y->z etc.
  permutation(const std::string& input)
  {
    boost::container::flat_map<std::size_t, std::size_t> mapping;

    // Remove the surrounding brackets if present.
    std::string trimmed_input = boost::trim_copy(input);
    std::string input_no_brackets = trimmed_input;
    if (!trimmed_input.empty() && trimmed_input.front() == '[' && trimmed_input.back() == ']')
    {
      input_no_brackets = trimmed_input.substr(1, trimmed_input.size() - 2);
    }

    // Parse all the commas.
    std::istringstream iss(input_no_brackets);
    std::string token;
    while (std::getline(iss, token, ','))
    {
      auto arrow_pos = token.find("->");
      if (arrow_pos == std::string::npos)
      {
        throw mcrl2::runtime_error("Invalid permutation format: " + token);
      }

      std::string from_str = boost::trim_copy(token.substr(0, arrow_pos));
      std::string to_str = boost::trim_copy(token.substr(arrow_pos + 2));

      std::size_t from = std::stoul(from_str);
      std::size_t to = std::stoul(to_str);

      if (mapping.contains(from))
      {
        throw mcrl2::runtime_error("Invalid permutation: multiple mappings for " + from_str);
      }

      mapping[from] = to;
    }

    m_mapping = mapping;
  }

  boost::container::flat_map<std::size_t, std::size_t> mapping() const { return m_mapping; }

  std::size_t operator[](std::size_t i) const
  {
    auto it = m_mapping.find(i);
    if (it != m_mapping.end())
    {
      return it->second;
    }
    return i;
  }

  // Returns true iff the permutation is the identity.
  bool is_identity() const
  {
    for (const auto& [key, value]: m_mapping)
    {
      if (key != value)
      {
        return false;
      }
    }
    return true;
  }

  // Applies the permutation to a set of indices.
  std::set<std::size_t> permute(const std::set<std::size_t>& s) const
  {
    std::set<std::size_t> result;
    for (const auto& i: s)
    {
      result.insert((*this)[i]);
    }
    return result;
  }

  /// Returns the concatenation of this permutation with another permutation.
  permutation concat(const permutation& other) const
  {
    boost::container::flat_map<std::size_t, std::size_t> new_mapping;

    for (const auto& [key, value]: m_mapping)
    {
      new_mapping[key] = other[value];
    }

    for (const auto& [key, value]: other.m_mapping)
    {
      assert(m_mapping.find(key) == m_mapping.end());
      new_mapping[key] = value;
    }

    return permutation(new_mapping);
  }

  bool operator==(const permutation& other) const { return m_mapping == other.m_mapping; }

private:
  boost::container::flat_map<std::size_t, std::size_t> m_mapping;
};

/// Iterator that generates all permutations of a given set of indices
class permutation_iterator
{
public:
  using value_type = permutation;
  using difference_type = std::ptrdiff_t;
  using pointer = const permutation*;
  using reference = const permutation&;
  using iterator_category = std::forward_iterator_tag;
  using iterator_concept = std::forward_iterator_tag;

  permutation_iterator(std::vector<std::size_t> indices)
    : m_indices(indices)
  {
    m_current_permutation = m_indices;
    next_permutation(); // Skip the identity permutation
  }

  permutation_iterator()
    : m_finished(true)
  {}

  permutation operator*() const
  {
    boost::container::flat_map<std::size_t, std::size_t> mapping;
    for (std::size_t i = 0; i < m_indices.size(); ++i)
    {
      mapping[m_indices[i]] = m_current_permutation[i];
    }
    return permutation(mapping);
  }

  permutation_iterator& operator++()
  {
    if (!next_permutation())
    {
      m_finished = true;
    }
    return *this;
  }

  permutation_iterator operator++(int)
  {
    permutation_iterator tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const permutation_iterator& other) const
  {
    if (m_finished && other.m_finished)
      return true;
    if (m_finished != other.m_finished)
      return false;
    return m_current_permutation == other.m_current_permutation;
  }

private:
  bool next_permutation()
  {
    // Edge case: no next permutation
    if (m_current_permutation.size() < 2)
    {
      return false;
    }

    // Find the largest index k such that a[k] < a[k + 1]
    int k = -1;
    for (int i = m_current_permutation.size() - 2; i >= 0; --i)
    {
      if (m_current_permutation[i] < m_current_permutation[i + 1])
      {
        k = static_cast<int>(i);
        break;
      }
    }

    if (k == -1)
    {
      return false; // No next permutation
    }

    // Find the largest index l greater than k such that a[k] < a[l]
    int l = -1;
    for (std::size_t i = m_current_permutation.size() - 1; i > k; --i)
    {
      if (m_current_permutation[k] < m_current_permutation[i])
      {
        l = static_cast<int>(i);
        break;
      }
    }

    // Swap a[k] and a[l]
    std::swap(m_current_permutation[k], m_current_permutation[l]);

    // Reverse the suffix starting at a[k + 1]
    std::reverse(m_current_permutation.begin() + k + 1, m_current_permutation.end());

    return true;
  }

  std::vector<std::size_t> m_indices;
  std::vector<std::size_t> m_current_permutation;
  bool m_finished = false;
};

/// Range abstraction over the iterator.
class permutation_range
{
public:
  permutation_range(const std::vector<std::size_t>& indices)
    : m_indices(indices)
  {
    std::sort(m_indices.begin(), m_indices.end());
  }

  permutation_iterator begin() const { return permutation_iterator(m_indices); }

  permutation_iterator end() const { return permutation_iterator(); }

private:
  std::vector<std::size_t> m_indices;
};

static_assert(std::forward_iterator<permutation_iterator>);
static_assert(std::ranges::range<permutation_range>);

/// Returns all the permutations for the given indices.
inline permutation_range permutation_group(const std::vector<std::size_t>& indices)
{
  return permutation_range(indices);
}

/// Apply the given detail::permutation to a pbes expression
inline pbes_expression
apply_permutation(const pbes_expression& expr, const std::vector<data::variable>& parameters, const detail::permutation& pi)
{
  data::mutable_map_substitution<> sigma;
  for (std::size_t i = 0; i < parameters.size(); ++i)
  {
    sigma[parameters.at(i)] = parameters.at(pi[i]);
  }

  auto result = pbes_system::replace_variables(expr, sigma);

  result = replace_propositional_variables(result,
    [sigma, pi, parameters](const pbes_system::propositional_variable_instantiation& x) -> pbes_system::pbes_expression
    {
      std::vector<data::data_expression> new_parameters(x.parameters().size());
      for (std::size_t i = 0; i < x.parameters().size(); ++i)
      {
        new_parameters[pi[i]] = data::data_expression(*std::next(x.parameters().begin(), i));
      }
      return propositional_variable_instantiation(x.name(), data::data_expression_list(new_parameters));
    });

  mCRL2log(log::debug) << "pi(phi): \n" << expr << "\n" << result << std::endl;
  return result;
}

/// Prints the permutation as a mapping
inline std::ostream& operator<<(std::ostream& out, const permutation& p)
{
  out << "[";
  bool first = true;
  for (const auto& [key, value]: p.mapping())
  {
    if (!first)
    {
      out << ", ";
    }

    out << key << " -> " << value;
    first = false;
  }
  out << "]";

  return out;
}

} // namespace mcrl2::pbes_system