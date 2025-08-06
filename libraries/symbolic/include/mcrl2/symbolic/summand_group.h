// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_SYMBOLIC_SUMMAND_GROUP_H
#define MCRL2_SYMBOLIC_SUMMAND_GROUP_H

#ifdef MCRL2_ENABLE_SYLVAN

#include "mcrl2/core/detail/print_utility.h"
#include "mcrl2/data/data_expression.h"
#include "mcrl2/symbolic/ordering.h"
#include "mcrl2/symbolic/utility.h"

#include <boost/dynamic_bitset.hpp>

#include <functional>
#include <iostream>
#include <vector>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm_list.cxx"
#else
  import atermpp;
#endif

namespace mcrl2::symbolic
{

/// \brief Computes for both the read and write parameters their positions in the zipped transition relation.
inline std::pair<std::vector<std::size_t>, std::vector<std::size_t>> compute_read_write_pos(const std::vector<std::size_t>& read, const std::vector<std::size_t>& write);

/// \brief Compact the vector to avoid repeated values at the end.
inline std::vector<std::uint32_t> optimise_project(const std::vector<std::uint32_t>& Ip_values);

/// \brief A transition relation over a single set of read and write parameters for a group of summands.
struct summand_group
{
  /// \brief The information of a single summand relevant for symbolic exploration.
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
  std::size_t learn_calls = 0; // Number of learn transition calls.

  summand_group(const data::variable_list& process_parameters, const boost::dynamic_bitset<>& read_write_pattern, bool has_action)
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

    Ip_values = optimise_project(Ip_values);
    read_parameters = project(as_vector(process_parameters), read);
    write_parameters = project(as_vector(process_parameters), write);
    L = empty_set();
    Ldomain = empty_set();
    Ir = compute_meta(read, write, has_action); // Note, action is always added the end.
    Ip = cube(Ip_values);

    std::tie(read_pos, write_pos) = compute_read_write_pos(read, write);
  }

  // Construct a summand group from the given parameters, where read_parameters and write_parameters have the same order as process_parameters.
  summand_group(const data::variable_list& process_parameters, std::vector<data::variable> _read_parameters, std::vector<data::variable> _write_parameters, bool has_action)
  {
    using namespace sylvan::ldds;
    read_parameters = _read_parameters;
    write_parameters = _write_parameters;
    
    std::vector<std::uint32_t> Ip_values;

    std::size_t index = 0;
    std::size_t read_index = 0;
    std::size_t write_index = 0;
    for (const data::variable& parameter : process_parameters)
    {
      if (read_index < read_parameters.size() && parameter == read_parameters[read_index])
      {
        read_index++;
        read.emplace_back(index);
        Ip_values.emplace_back(1);
      }
      else
      {
        Ip_values.emplace_back(0);
      }

      if (write_index < write_parameters.size() && parameter == write_parameters[write_index])
      {
        write_index++;
        write.emplace_back(index);
      }

      ++index;
    }

    Ip_values = optimise_project(Ip_values);
    
    L = empty_set();
    Ldomain = empty_set();
    Ir = compute_meta(read, write, has_action); // Note, action is always added the end.
    Ip = cube(Ip_values);

    std::tie(read_pos, write_pos) = compute_read_write_pos(read, write);
  }
};

/// \brief Prints the information of this summand group.
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

std::pair<std::vector<std::size_t>, std::vector<std::size_t>> compute_read_write_pos(const std::vector<std::size_t>& read, const std::vector<std::size_t>& write)
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


std::vector<std::uint32_t> optimise_project(const std::vector<std::uint32_t>& Ip_values)
{
  if (Ip_values.size() > 1)
  {
    // The index after which all values in Ip_values are the same (at least two entries).
    int i = Ip_values.size() - 1;
    while (i > 0 && Ip_values[i] == Ip_values.back())
    {
      --i;
    }

    std::vector<std::uint32_t> result(Ip_values.begin(), Ip_values.begin() + i + 2);
    if (result.back() == 0)
    {
      result.back() = static_cast<std::uint32_t>(-2);
    }
    else
    {
      result.back() = static_cast<std::uint32_t>(-1);
    }
    return result;
  }

  return Ip_values;
}

} // namespace mcrl2::symbolic

#endif // MCRL2_ENABLE_SYLVAN

#endif // MCRL2_SYMBOLIC_SUMMAND_GROUP_H