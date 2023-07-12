// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_SYMBOLIC_SYMBOLIC_REACHABILITY_H
#define MCRL2_SYMBOLIC_SYMBOLIC_REACHABILITY_H

#ifdef MCRL2_ENABLE_SYLVAN

#include "mcrl2/data/consistency.h"
#include "mcrl2/data/data_specification.h"
#include "mcrl2/data/enumerator.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/data/undefined.h"
#include "mcrl2/utilities/stopwatch.h"
#include "mcrl2/symbolic/alternative_relprod.h"
#include "mcrl2/symbolic/summand_group.h"

#include <sylvan_ldd.hpp>

namespace mcrl2::symbolic {

struct symbolic_reachability_options
{
  data::rewrite_strategy rewrite_strategy = data::jitty;
  std::size_t max_workers = 0;
  std::size_t max_iterations = 0;
  bool cached = false;
  bool chaining = false;
  bool detect_deadlocks = false;
  bool one_point_rule_rewrite = false;
  bool print_exact = false;
  bool print_nodesize = false;
  bool replace_constants_by_variables = false;
  bool remove_unused_rewrite_rules = false;
  bool saturation = false;
  bool strategy = false;
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
  out << "print_nodesize = " << std::boolalpha << options.print_nodesize << std::endl;
  out << "replace-constants-by-variables = " << std::boolalpha << options.replace_constants_by_variables << std::endl;
  out << "remove-unused-rewrite-rules = " << std::boolalpha << options.remove_unused_rewrite_rules << std::endl;
  out << "saturation = " << std::boolalpha << options.saturation << std::endl;
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

/// \brief If ActionLabel is true then the multi-action will be rewritten and added to the relation.
template <typename Context, bool ActionLabel>
void learn_successors_callback(WorkerP*, Task*, std::uint32_t* x, std::size_t, void* context)
{
  using namespace sylvan::ldds;
  using enumerator_element = data::enumerator_list_element_with_substitution<>;

  auto p = reinterpret_cast<Context*>(context);
  auto& algorithm = p->first;
  auto& group = p->second;
  auto& sigma = algorithm.m_sigma;
  auto& data_index = algorithm.data_index();
  const auto& options = algorithm.m_options;
  const auto& rewr = algorithm.m_rewr;
  const auto& enumerator = algorithm.m_enumerator;
  std::size_t x_size = group.read.size();
  std::size_t y_size = group.write.size();
  std::size_t xy_size = x_size + y_size; 
  
  if constexpr (ActionLabel)
  {
    // One additional space for the action label.
    xy_size += 1;
  }

  MCRL2_DECLARE_STACK_ARRAY(xy, std::uint32_t, xy_size);

  // add the assignments corresponding to x to sigma
  // add x to the transition xy
  stopwatch learn_start;
  for (std::size_t j = 0; j < x_size; j++)
  {
    sigma[group.read_parameters[j]] = data_index[group.read[j]][x[j]];
    xy[group.read_pos[j]] = x[j];
  }

  std::size_t i = 0;
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
                               assert(value != data::undefined_data_expression());

                               // Determine whether this is a copy parameter, insert special value if that is the case.
                               xy[group.write_pos[j]] = smd.copy[group.write_pos[j]] ? relprod_ignore : data_index[group.write[j]].insert(value).first;
                             }

                             if constexpr (ActionLabel)
                             {
                               // Action is always located on the last index of the cube.
                               xy[xy_size - 1] = algorithm.action_index().insert(algorithm.rewrite_action(group.actions[i], rewr, sigma)).first;
                             }

                             mCRL2log(log::debug1) << "  " << print_transition(data_index, xy.data(), group.read, group.write) << std::endl;
                             group.L = options.no_relprod ? union_cube(group.L, xy.data(), xy_size) : union_cube_copy(group.L, xy.data(), smd.copy.data(), xy_size);
                             return false;
                           },
                           data::is_false
      );

      ++i;
    }
    data::remove_assignments(sigma, smd.variables);
  }
  data::remove_assignments(sigma, group.read_parameters);
  group.learn_calls += 1;
  group.learn_time += learn_start.seconds();

  if (options.cached)
  {
    group.Ldomain = union_cube(group.Ldomain, x, x_size);
  }
}

} // namespace mcrl2::symbolic

#endif // MCRL2_ENABLE_SYLVAN

#endif // MCRL2_SYMBOLIC_SYMBOLIC_REACHABILITY_H
