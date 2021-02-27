// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/symbolic_reachability.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_SYMBOLIC_REACHABILITY_H
#define MCRL2_LPS_SYMBOLIC_REACHABILITY_H

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <random>
#include <regex>
#include <sylvan_ldd.hpp>
#include <boost/dynamic_bitset.hpp>
#include "mcrl2/core/detail/print_utility.h"
#include "mcrl2/data/consistency.h"
#include "mcrl2/data/data_specification.h"
#include "mcrl2/data/enumerator.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/data/undefined.h"
#include "mcrl2/utilities/parse_numbers.h"
#include "mcrl2/utilities/text_utility.h"
#include "mcrl2/utilities/indexed_set.h"

namespace mcrl2 {

namespace lps {

constexpr std::uint32_t relprod_ignore = std::numeric_limits<std::uint32_t>::max(); // used by alternative_relprod/relprev

// Return a permuted copy v' of v with v'[i] = v[permutation[i]]
// N.B. the implementation is not efficient
template <typename Container>
Container permute_copy(const Container& v, const std::vector<std::size_t>& permutation)
{
  using T = typename Container::value_type;

  std::size_t n = v.size();
  assert(permutation.size() == n);
  std::vector<T> v_(v.begin(), v.end());
  std::vector<T> result(n);
  for (std::size_t i = 0; i < n; i++)
  {
    result[i] = v_[permutation[i]];
  }
  return Container(result.begin(), result.end());
}

inline
std::vector<std::set<std::size_t>> parse_summand_groups(std::string text, std::size_t n)
{
  using utilities::regex_split;
  using utilities::trim_copy;
  using utilities::detail::as_set;

  // check if the format of the group is correct
  std::string group_format = R"(\s*\d+(\s+\d+)*\s*)";
  if (!std::regex_match(text, std::regex(group_format + "(;" + group_format + ")*")))
  {
    throw mcrl2::runtime_error("The format of the groups '" + text + "' is incorrect.");
  }

  std::vector<std::set<std::size_t>> result;
  text = trim_copy(text);
  for (const std::string& numbers: regex_split(text, "\\s*;\\s*"))
  {
    result.push_back(as_set(utilities::parse_natural_number_sequence(numbers)));
  }

  // check if result is a partition of [0 .. n)
  std::set<std::size_t> all;
  std::set<std::size_t> expected;
  for (std::size_t i = 0; i < n; i++)
  {
    expected.insert(i);
  }
  for (const auto& group: result)
  {
    std::size_t size = all.size();
    all.insert(group.begin(), group.end());
    if (all.size() != size + group.size())
    {
      throw mcrl2::runtime_error("The groups in '" + text + "' are not disjoint.");
    }
  }
  if (all != expected)
  {
    throw mcrl2::runtime_error("The groups '" + text + "' do not form a partition of [0 .. " + std::to_string(n-1) + "].");
  }

  return result;
}

// puts every summand in a separate group
inline
std::vector<std::set<std::size_t>> compute_summand_groups_default(const std::vector<boost::dynamic_bitset<>>& patterns)
{
  std::vector<std::set<std::size_t>> result;
  for (std::size_t i = 0; i < patterns.size(); i++)
  {
    result.push_back({i});
  }
  return result;
}

// joins summands with exactly the same pattern
inline
std::vector<std::set<std::size_t>> compute_summand_groups_simple(const std::vector<boost::dynamic_bitset<>>& patterns)
{
  std::map<boost::dynamic_bitset<>, std::set<std::size_t>> group_map;
  for (std::size_t i = 0; i < patterns.size(); i++)
  {
    group_map[patterns[i]].insert(i);
  }
  std::vector<std::set<std::size_t>> groups;
  for (const auto& [_, group]: group_map)
  {
    groups.push_back(group);
  }
  return groups;
}

inline
std::vector<std::set<std::size_t>> compute_summand_groups(const std::string& text, const std::vector<boost::dynamic_bitset<>>& patterns)
{
  if (text == "none")
  {
    return compute_summand_groups_default(patterns);
  }
  else if (text == "simple")
  {
    return compute_summand_groups_simple(patterns);
  }
  else
  {
    return parse_summand_groups(text, patterns.size());
  }
}

inline
std::vector<std::size_t> compute_variable_order_default(std::size_t n)
{
  std::vector<std::size_t> result;
  for (std::size_t i = 0; i < n; i++)
  {
    result.push_back(i);
  }
  return result;
}

inline
std::vector<std::size_t> compute_variable_order_random(std::size_t n, bool exclude_first_variable = false)
{
  std::vector<std::size_t> result;
  for (std::size_t i = 0; i < n; i++)
  {
    result.push_back(i);
  }
  std::random_device rd;
  std::mt19937 g(rd());
  std::size_t first = exclude_first_variable ? 1 : 0;
  std::shuffle(result.begin() + first, result.end(), g);
  return result;
}

inline
std::vector<std::size_t> parse_variable_order(std::string text, std::size_t n, bool exclude_first_variable = false)
{
  using utilities::trim_copy;
  using utilities::detail::as_set;

  // check if the format of the order is correct
  std::string format = R"(\s*\d+(\s+\d+)*\s*)";
  if (!std::regex_match(text, std::regex(format)))
  {
    throw mcrl2::runtime_error("The format of the variable order '" + text + "' is incorrect.");
  }

  text = trim_copy(text);
  std::vector<std::size_t> result = utilities::parse_natural_number_sequence(trim_copy(text));

  std::set<std::size_t> expected;
  for (std::size_t i = 0; i < n; i++)
  {
    expected.insert(i);
  }
  if (as_set(result) != expected)
  {
    throw mcrl2::runtime_error("The variable order '" + text + "' is not a permutation of [0 .. " + std::to_string(n-1) + "].");
  }

  if (exclude_first_variable && result[0] != 0)
  {
    throw mcrl2::runtime_error("The variable order '" + text + "' does not start with a zero.");
  }
  return result;
}

inline
std::vector<std::size_t> compute_variable_order(const std::string& text, const std::vector<boost::dynamic_bitset<>>& patterns, bool exclude_first_variable = false)
{
  std::size_t n = patterns.front().size() / 2;
  if (text == "none")
  {
    return compute_variable_order_default(n);
  }
  else if (text == "random")
  {
    return compute_variable_order_random(n, exclude_first_variable);
  }
  else
  {
    return parse_variable_order(text, n, exclude_first_variable);
  }
}

inline
std::vector<boost::dynamic_bitset<>> reorder_read_write_patterns(const std::vector<boost::dynamic_bitset<>>& patterns, const std::vector<std::size_t>& variable_order)
{
  std::size_t n = variable_order.size();
  std::vector<boost::dynamic_bitset<>> result;
  for (const auto& pattern: patterns)
  {
    boost::dynamic_bitset<> reordered_pattern(2*n);
    for (std::size_t i = 0; i < n; i++)
    {
      std::size_t j = variable_order[i];
      reordered_pattern[2*i] = pattern[2*j];
      reordered_pattern[2*i+1] = pattern[2*j+1];
    }
    result.push_back(reordered_pattern);
  }
  return result;
}

struct symbolic_reachability_options
{
  data::rewrite_strategy rewrite_strategy = data::jitty;
  bool cached = false;
  bool chaining = false;
  bool detect_deadlocks = false;
  bool one_point_rule_rewrite = false;
  bool replace_constants_by_variables = false;
  bool remove_unused_rewrite_rules = false;
  bool no_discard = false;
  bool no_discard_read = false;
  bool no_discard_write = false;
  bool no_relprod = false;
  bool info = false;
  std::string summand_groups;
  std::string variable_order;
  std::string dot_file;
};

inline
std::ostream& operator<<(std::ostream& out, const symbolic_reachability_options& options)
{
  out << "rewrite-strategy = " << options.rewrite_strategy << std::endl;
  out << "cached = " << std::boolalpha << options.cached << std::endl;
  out << "chaining = " << std::boolalpha << options.chaining << std::endl;
  out << "detect_deadlocks = " << std::boolalpha << options.detect_deadlocks << std::endl;
  out << "one-point-rule-rewrite = " << std::boolalpha << options.one_point_rule_rewrite << std::endl;
  out << "replace-constants-by-variables = " << std::boolalpha << options.replace_constants_by_variables << std::endl;
  out << "remove-unused-rewrite-rules = " << std::boolalpha << options.remove_unused_rewrite_rules << std::endl;
  out << "no-discard = " << std::boolalpha << options.no_discard << std::endl;
  out << "no-read = " << std::boolalpha << options.no_discard_read << std::endl;
  out << "no-write = " << std::boolalpha << options.no_discard_write << std::endl;
  out << "no-relprod = " << std::boolalpha << options.no_relprod << std::endl;
  out << "info = " << std::boolalpha << options.info << std::endl;
  out << "groups = " << options.summand_groups << std::endl;
  out << "reorder = " << options.variable_order << std::endl;
  out << "dot = " << options.dot_file << std::endl;
  return out;
}

inline
std::set<std::size_t> parameter_indices(const std::set<data::variable>& parameters, const std::map<data::variable, std::size_t>& index)
{
  std::set<std::size_t> result;
  for (const data::variable& p: parameters)
  {
    result.insert(index.at(p));
  }
  return result;
}

inline
std::string print_read_write_patterns(const std::vector<boost::dynamic_bitset<>>& patterns)
{
  std::size_t n = patterns.front().size() / 2;

  auto print_used = [n](const boost::dynamic_bitset<>& pattern)
  {
    std::ostringstream out;
    for (std::size_t i = 0; i < n; i++)
    {
      bool read = pattern[2*i];
      bool write = pattern[2*i+1];
      bool used = read || write;
      out << (used ? "1" : "0");
    }
    return out.str();
  };

  // N.B. The stream operator of boost::dynamic_bitset prints the bits in reverse order(!)
  auto print_rw = [n](const boost::dynamic_bitset<>& pattern)
  {
    std::ostringstream out;
    for (std::size_t i = 0; i < n; i++)
    {
      bool read = pattern[2*i];
      bool write = pattern[2*i+1];
      out << (read ? "1" : "0");
      out << (write ? "1" : "0");
    }
    return out.str();
  };

  auto print_compact = [n](const boost::dynamic_bitset<>& pattern)
  {
    std::ostringstream out;
    for (std::size_t i = 0; i < n; i++)
    {
      bool read = pattern[2*i];
      bool write = pattern[2*i+1];
      if (read && write)
      {
        out << '+';
      }
      else if (read)
      {
        out << 'r';
      }
      else if (write)
      {
        out << 'w';
      }
      else
      {
        out << '-';
      }
    }
    return out.str();
  };

  std::ostringstream out;
  out << "used parameters" << std::endl;
  for (std::size_t i = 0; i < patterns.size(); i++)
  {
    out << std::setw(4) << i << " " << print_used(patterns[i]) << std::endl;
  }
  out << std::endl;

  out << "read/write patterns" << std::endl;
  for (std::size_t i = 0; i < patterns.size(); i++)
  {
    out << std::setw(4) << i << " " << print_rw(patterns[i]) << std::endl;
  }
  out << std::endl;

  out << "read/write patterns compacted" << std::endl;
  for (std::size_t i = 0; i < patterns.size(); i++)
  {
    out << std::setw(4) << i << " " << print_compact(patterns[i]) << std::endl;
  }
  return out.str();
}

// Set some additional read and write bits to true in patterns. This is only used for testing purposes.
template <typename Options>
void adjust_read_write_patterns(std::vector<boost::dynamic_bitset<>>& patterns, const Options& options)
{
  for (boost::dynamic_bitset<>& pattern: patterns)
  {
    std::size_t n = pattern.size() / 2;
    for (std::size_t j = 0; j < n; j++)
    {
      bool read = pattern[2*j];
      bool write = pattern[2*j + 1];
      if (options.no_discard)
      {
        read = true;
        write = true;
      }
      else if (read || write)
      {
        if (options.no_discard_read)
        {
          read = true;
        }
        if (options.no_discard_write)
        {
          write = true;
        }
      }
      pattern[2*j] = read;
      pattern[2*j + 1] = write;
    }
  }
}

inline
std::vector<boost::dynamic_bitset<>> compute_summand_group_patterns(const std::vector<boost::dynamic_bitset<>>& patterns, const std::vector<std::set<std::size_t>> groups)
{
  std::vector<boost::dynamic_bitset<>> result;
  boost::dynamic_bitset<> empty(patterns.front().size());

  for (const std::set<std::size_t>& group: groups)
  {
    if (group.empty())
    {
      continue;
    }
    boost::dynamic_bitset<> pattern = empty;
    for (std::size_t i: group)
    {
      pattern |= patterns[i];
    }
    result.push_back(pattern);
  }
  return result;
}

/// \brief A bidirectional mapping between data expressions of a given sort and numbers
class data_expression_index
{
    friend std::ostream& operator<<(std::ostream&, const data_expression_index&);

  protected:
    data::sort_expression m_sort;
    mcrl2::utilities::indexed_set<data::data_expression,
      std::hash<data::data_expression>,
      std::equal_to<data::data_expression>,
      std::allocator<data::data_expression>,
      true> m_values;

  public:
    data_expression_index(const data::sort_expression& sort)
      : m_sort(sort)
    {}

    /// \brief Returns the index of the given value. If the value is not present yet, it will be added.
    std::uint32_t index(const data::data_expression& value)
    {
      assert(value.sort() == m_sort);
      return m_values.insert(value).first;
    }

    /// \brief Returns the value corresponding to index
    const data::data_expression& value(std::uint32_t i) const
    {
      return m_values[i];
    }

    const data::sort_expression& sort() const
    {
      return m_sort;
    }

    bool has_value(const data::data_expression& value) const
    {
      return m_values.find(value) != m_values.end();
    }

    bool has_index(std::uint32_t i) const
    {
      return i < m_values.size();
    }

    std::size_t size() const
    {
      return m_values.size();
    }

    auto begin() const
    {
      return m_values.begin();
    }

    auto end() const
    {
      return m_values.end();
    }
};

inline
std::ostream& operator<<(std::ostream& out, const data_expression_index& x)
{
  out << "data index with sort = " << x.sort() << " values = " << core::detail::print_list(x.m_values);
  return out;
}

template <typename T>
std::vector<T> project(const std::vector<T>& v, const std::vector<std::size_t>& used)
{
  std::vector<T> result;
  result.reserve(used.size());
  for (std::size_t i = 0; i < used.size(); i++)
  {
    result.push_back(v[used[i]]);
  }
  return result;
}

struct summand_group
{
  struct summand
  {
    data::data_expression condition;
    data::variable_list variables; // the summand variables
    std::vector<data::data_expression> next_state; // the projected next state vector
    std::vector<int> copy; // copy node information that is needed by sylvan::ldds::relprod

    summand(const data::data_expression& condition_, const data::variable_list& variables_, const std::vector<data::data_expression>& next_state_, const std::vector<int>& copy_)
      : condition(condition_), variables(variables_), next_state(next_state_), copy(copy_)
    {}
  };

  std::vector<summand> summands; // the summands of the group
  std::vector<data::variable> read_parameters; // the read parameters
  std::vector<std::size_t> read; // indices of the read parameters
  std::vector<std::size_t> read_pos; // indices of the read parameters in a zipped transition xy
  std::vector<data::variable> write_parameters; // the write parameters
  std::vector<std::size_t> write; // indices of the write parameters
  std::vector<std::size_t> write_pos; // indices of the write parameters in a zipped transition xy
  sylvan::ldds::ldd L; // the projected transition relation
  sylvan::ldds::ldd Ldomain; // the domain of L
  sylvan::ldds::ldd Ir; // meta data needed by sylvan::ldds::relprod
  sylvan::ldds::ldd Ip; // meta data needed by sylvan::ldds::project
  double learn_time = 0.0; // The time to learn the transitions for this group.
  mcrl2::utilities::copyable_mutex m_mutex; // A mutex per summand group for updating L.

  std::pair<std::vector<std::size_t>, std::vector<std::size_t>> compute_read_write_pos() const
  {
    std::vector<std::size_t> rpos;
    std::vector<std::size_t> wpos;

    auto ri = read.begin();
    auto wi = write.begin();
    std::size_t index = 0;
    while (ri != read.end() && wi != write.end())
    {
      if (*ri <= *wi)
      {
        rpos.push_back(index++);
        ri++;
      }
      else
      {
        wpos.push_back(index++);
        wi++;
      }
    }
    while (ri != read.end())
    {
      rpos.push_back(index++);
      ri++;
    }
    while (wi != write.end())
    {
      wpos.push_back(index++);
      wi++;
    }

    return { rpos, wpos };
  }

  summand_group(const data::variable_list& process_parameters, const boost::dynamic_bitset<>& read_write_pattern)
  {
    using namespace sylvan::ldds;
    using utilities::detail::as_vector;

    std::size_t n = process_parameters.size();

    // Indicates for every position whether the parameter should be projected on (for read dependencies).
    std::vector<std::uint32_t> Ip_values;

    for (std::size_t j = 0; j < n; j++)
    {
      bool is_read = read_write_pattern[2*j];
      bool is_write = read_write_pattern[2*j + 1];
      Ip_values.push_back(is_read ? 1 : 0);
      if (is_read)
      {
        read.push_back(j);
      }
      if (is_write)
      {
        write.push_back(j);
      }
    }

    // The index after which all values in Ip_values are the same.
    int i = Ip_values.size() - 1;
    while (i > 0 && Ip_values[i] == Ip_values.back())
    {
      --i;
    }

    Ip_values = std::vector<std::uint32_t>(Ip_values.begin(), Ip_values.begin() + i + 2);
    if (Ip_values.back() == 0)
    {
      Ip_values.back() = static_cast<std::uint32_t>(-2);
    }
    else
    {
      Ip_values.back() = static_cast<std::uint32_t>(-1);
    }

    read_parameters = project(as_vector(process_parameters), read);
    write_parameters = project(as_vector(process_parameters), write);
    L = empty_set();
    Ldomain = empty_set();
    Ir = compute_meta(read, write);
    Ip = cube(Ip_values);
    std::tie(read_pos, write_pos) = compute_read_write_pos();
  }
};

inline
std::ostream& operator<<(std::ostream& out, const summand_group& x)
{
  using namespace sylvan::ldds;
  for (const auto& smd: x.summands)
  {
    out << "condition = " << smd.condition << std::endl;
    out << "variables = " << core::detail::print_list(smd.variables) << std::endl;
    out << "next state = " << core::detail::print_list(smd.next_state) << std::endl;
    out << "copy = " << core::detail::print_list(smd.copy) << std::endl;

    std::vector<std::string> assignments;
    auto vi = x.write_parameters.begin();
    auto ni = smd.next_state.begin();
    for (; vi != x.write_parameters.end(); ++vi, ++ni)
    {
      assignments.push_back(data::pp(*vi) + " := " + data::pp(*ni));
    }
    out << "assignments = " << core::detail::print_list(assignments) << std::endl;
  }
  out << "read = " << core::detail::print_list(x.read) << std::endl;
  out << "read parameters = " << core::detail::print_list(x.read_parameters) << std::endl;
  out << "write = " << core::detail::print_list(x.write) << std::endl;
  out << "write parameters = " << core::detail::print_list(x.write_parameters) << std::endl;
  out << "L = " << print_ldd(x.L) << std::endl;
  out << "Ir = " << print_ldd(x.Ir) << std::endl;
  out << "Ip = " << print_ldd(x.Ip) << std::endl;
  return out;
}

// A very inefficient implementation of relprod, that matches the specification closely
sylvan::ldds::ldd alternative_relprod(const sylvan::ldds::ldd& todo, const summand_group& R)
{
  using namespace sylvan::ldds;

  auto split = [&](const std::vector<std::uint32_t>& xy)
  {
    std::vector<std::uint32_t> x;
    std::vector<std::uint32_t> y;
    for (std::size_t j: R.read_pos)
    {
      x.push_back(xy[j]);
    }
    for (std::size_t j: R.write_pos)
    {
      y.push_back(xy[j]);
    }
    return std::make_pair(x, y);
  };

  auto match = [&](const std::vector<std::uint32_t>& x, const std::vector<std::uint32_t>& x_)
  {
    for (std::size_t j = 0; j < x_.size(); j++)
    {
      if (x[R.read[j]] != x_[j])
      {
        return false;
      }
    }
    return true;
  };

  auto replace = [&](std::vector<std::uint32_t> x, const std::vector<std::uint32_t>& y_)
  {
    for (std::size_t j = 0; j < y_.size(); j++)
    {
      if (y_[j] != relprod_ignore)
      {
        x[R.write[j]] = y_[j];
      }
    }
    return x;
  };

  ldd result = empty_set();

  auto todo_elements = ldd_solutions(todo);
  for (const std::vector<std::uint32_t>& xy: ldd_solutions(R.L))
  {
    auto [x_, y_] = split(xy);
    for (const auto& x: todo_elements)
    {
      if (match(x, x_))
      {
        auto y = replace(x, y_);
        result = union_cube(result, y);
      }
    }
  }
  return result;
}

// A very inefficient implementation of relprev, that matches the specification closely
sylvan::ldds::ldd alternative_relprev(const sylvan::ldds::ldd& Y, const summand_group& R, const sylvan::ldds::ldd& X)
{
  using namespace sylvan::ldds;

  ldd result = empty_set();
  for (const std::vector<std::uint32_t>& x_values: ldd_solutions(X))
  {
    ldd x = cube(x_values);
    ldd y = alternative_relprod(x, R);
    if (intersect(y, Y) != empty_set())
    {
      result = union_(result, x);
    }
  }
  return result;
}

std::vector<data::data_expression> ldd2state(const std::vector<data_expression_index>& data_index, const std::vector<std::uint32_t>& x)
{
  std::vector<data::data_expression> result;
  for (std::size_t i = 0; i < x.size(); i++)
  {
    if (x[i] == relprod_ignore)
    {
      result.push_back(data::undefined_data_expression());
    }
    else
    {
      result.push_back(data_index[i].value(x[i]));
    }
  }
  return result;
}

std::vector<data::data_expression> ldd2state(const std::vector<data_expression_index>& data_index, const std::vector<std::uint32_t>& x, const std::vector<std::size_t>& used)
{
  std::vector<data::data_expression> result;
  for (std::size_t i = 0; i < used.size(); i++)
  {
    if (x[i] == relprod_ignore)
    {
      result.push_back(data::undefined_data_expression());
    }
    else
    {
      result.push_back(data_index[used[i]].value(x[i]));
    }
  }
  return result;
}

inline
std::string print_state(const std::vector<data_expression_index>& data_index, const std::vector<std::uint32_t>& x)
{
  return core::detail::print_list(ldd2state(data_index, x));
}

inline
std::string print_state(const std::vector<data_expression_index>& data_index, const std::vector<std::uint32_t>& x, const std::vector<std::size_t>& used)
{
  return core::detail::print_list(ldd2state(data_index, x, used));
}

inline
std::string print_states(const std::vector<data_expression_index>& data_index, const sylvan::ldds::ldd& x)
{
  std::ostringstream out;
  auto solutions = ldd_solutions(x);
  bool multiline = solutions.size() > 1;
  std::string sep = multiline ? ",\n" : ", ";

  out << "{" << (multiline ? "\n" : " ");
  for (std::size_t i = 0; i < solutions.size(); i++)
  {
    if (i > 0)
    {
      out << sep;
    }
    out << print_state(data_index, solutions[i]);
  }
  out << (multiline ? "\n" : " ") << "}";
  return out.str();
}

std::string print_states(const std::vector<data_expression_index>& data_index, const sylvan::ldds::ldd& x, const std::vector<std::size_t>& used)
{
  std::ostringstream out;
  auto solutions = ldd_solutions(x);
  bool multiline = solutions.size() > 1;
  std::string sep = multiline ? ",\n" : ", ";

  out << "{" << (multiline ? "\n" : " ");
  for (std::size_t i = 0; i < solutions.size(); i++)
  {
    if (i > 0)
    {
      out << sep;
    }
    out << print_state(data_index, solutions[i], used);
  }
  out << (multiline ? "\n" : " ") << "}";
  return out.str();
}

std::string print_transition(const std::vector<data_expression_index>& data_index, const std::vector<std::uint32_t>& xy)
{
  std::size_t n = xy.size() / 2;
  std::vector<std::uint32_t> x(n);
  std::vector<std::uint32_t> y(n);
  for (std::size_t i = 0; i < n; i++)
  {
    x[i] = xy[2*i];
    y[i] = xy[2*i+1];
  }
  return print_state(data_index, x) + " -> " + print_state(data_index, y);
}

std::string print_transition(const std::vector<data_expression_index>& data_index, const std::uint32_t* xy, const std::vector<std::size_t>& read, const std::vector<std::size_t>& write)
{
  std::vector<std::uint32_t> x;
  std::vector<std::uint32_t> y;
  auto ri = read.begin();
  auto wi = write.begin();
  auto xyi = xy;
  while (ri != read.end() && wi != write.end())
  {
    if (*ri <= *wi)
    {
      ri++;
      x.push_back(*xyi++);
    }
    else
    {
      wi++;
      y.push_back(*xyi++);
    }
  }
  while (ri != read.end())
  {
    ri++;
    x.push_back(*xyi++);
  }
  while (wi != write.end())
  {
    wi++;
    y.push_back(*xyi++);
  }
  return print_state(data_index, x, read) + " -> " + print_state(data_index, y, write);
}

std::string print_relation(const std::vector<data_expression_index>& data_index, const sylvan::ldds::ldd& R, const std::vector<std::size_t>& read, const std::vector<std::size_t>& write)
{
  std::ostringstream out;
  for (const std::vector<std::uint32_t>& xy: ldd_solutions(R))
  {
    out << print_transition(data_index, xy.data(), read, write) << std::endl;
  }
  return out.str();
}

// Add operations on reals that are needed for the exploration.
inline
std::set<data::function_symbol> add_real_operators(std::set<data::function_symbol> s)
{
  std::set<data::function_symbol> result = std::move(s);
  result.insert(data::less_equal(data::sort_real::real_()));
  result.insert(data::greater_equal(data::sort_real::real_()));
  result.insert(data::sort_real::plus(data::sort_real::real_(), data::sort_real::real_()));
  return result;
}

inline
data::rewriter construct_rewriter(const data::data_specification& dataspec, data::rewrite_strategy rewrite_strategy, const std::set<data::function_symbol>& function_symbols, bool remove_unused_rewrite_rules)
{
  if (remove_unused_rewrite_rules)
  {
    return data::rewriter(dataspec,
                          data::used_data_equation_selector(dataspec, add_real_operators(function_symbols), std::set<data::variable>()),
                          rewrite_strategy);
  }
  else
  {
    return data::rewriter(dataspec, rewrite_strategy);
  }
}

template <typename EnumeratorElement>
void check_enumerator_solution(const EnumeratorElement& p, const summand_group&)
{
  if (p.expression() != data::sort_bool::true_())
  {
    // TODO: print the problematic expression in the same way it is done in lps2lts(?)
    throw data::enumerator_error("Expression does not rewrite to true or false: " + data::pp(p.expression()));
  }
}

template <typename Context>
void learn_successors_callback(WorkerP* worker, Task*, std::uint32_t* x, std::size_t n, void* context)
{
  using namespace sylvan::ldds;
  using enumerator_element = data::enumerator_list_element_with_substitution<>;

  auto p = reinterpret_cast<Context*>(context);
  auto& algorithm = p->first;
  auto& group = p->second;
  auto& data_index = algorithm.m_data_index;
  auto& per_worker = algorithm.m_workers[worker->worker];
  auto& sigma = per_worker.m_sigma;

  const auto& rewr = algorithm.m_rewr;
  const auto& enumerator = algorithm.m_enumerator;

  std::size_t x_size = group.read.size();
  std::size_t y_size = group.write.size();
  std::size_t xy_size = x_size + y_size;

  MCRL2_DECLARE_STACK_ARRAY(xy, std::uint32_t, xy_size);

  // add the assignments corresponding to x to sigma
  // add x to the transition xy
  for (std::size_t j = 0; j < x_size; j++)
  {
    sigma[group.read_parameters[j]] = data_index[group.read[j]].value(x[j]);
    xy[group.read_pos[j]] = x[j];
  }

  for (const auto& smd: group.summands)
  {
    data::data_expression condition = rewr(smd.condition, sigma);
    if (!data::is_false(condition))
    {
      enumerator.enumerate(enumerator_element(smd.variables, condition),
                           sigma,
                           [&](const enumerator_element& p) {
                             check_enumerator_solution(p, group);
                             p.add_assignments(smd.variables, sigma, rewr);
                             for (std::size_t j = 0; j < y_size; j++)
                             {
                               data::data_expression value = rewr(smd.next_state[j], sigma);
                               xy[group.write_pos[j]] = data::is_variable(value) ? lps::relprod_ignore : data_index[group.write[j]].index(value);
                             }
                             mCRL2log(log::debug) << "  " << print_transition(data_index, xy.data(), group.read, group.write) << std::endl;
                             
                             group.m_mutex.lock();
                             group.Ldomain = union_cube(group.Ldomain, x, x_size);
                             group.L = algorithm.m_options.no_relprod ? union_cube(group.L, xy.data(), xy_size) : union_cube_copy(group.L, xy.data(), smd.copy.data(), xy_size);
                             group.m_mutex.unlock();                             
                             return false;
                           },
                           data::is_false
      );
    }
    data::remove_assignments(sigma, smd.variables);
  }
  data::remove_assignments(sigma, group.read_parameters);

  if (options.cached)
  {
    group.Ldomain = union_cube(group.Ldomain, x, x_size);
  }
}

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_SYMBOLIC_REACHABILITY_H
