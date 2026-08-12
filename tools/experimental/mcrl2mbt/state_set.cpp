// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "state_set.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

#include "mcrl2/lps/print.h"

namespace mcrl2
{

// ─── helpers ──────────────────────────────────────────────────────────────────

/// \brief Return true when the multi-action represents an internal (tau) step.
static bool is_tau(const lps::multi_action& action)
{
  return action.actions().empty();
}

/// \brief Extract the first action label name from a non-tau multi-action.
static std::string first_label(const lps::multi_action& action)
{
  return std::string(action.actions().front().label().name());
}
/// rief The full wire-format label: lps::pp(), matching what adapters send.
static std::string wire_label(const lps::multi_action& action)
{
  return lps::pp(action);
}
// ─── state_set ────────────────────────────────────────────────────────────────

state_set::state_set(const lps::specification& spec, const lps::explorer_options& options, const data::rewriter& rewr)
  : m_rewr(rewr),
    m_explorer(spec, options, rewr)
{}

lps::state state_set::compute_initial_state() const
{
  data::mutable_indexed_substitution<> sigma;
  lps::state s0;
  const data::data_expression_list& init = m_explorer.initial_state();
  lps::make_state(s0,
    init.begin(),
    init.size(),
    [&](data::data_expression& result, const data::data_expression& x) { m_rewr(result, x, sigma); });
  return s0;
}

void state_set::reset_to_initial()
{
  m_states = {compute_initial_state()};
}

std::vector<lps::state> state_set::tau_closure(const std::vector<lps::state>& seeds, std::size_t depth)
{
  // BFS over tau edges; track visited states to avoid cycles.
  std::vector<lps::state> result = seeds;
  std::unordered_set<lps::state> visited(seeds.begin(), seeds.end());

  std::vector<lps::state> frontier = seeds;
  for (std::size_t d = 0; d < depth && !frontier.empty(); ++d)
  {
    std::vector<lps::state> next_frontier;
    for (const lps::state& s: frontier)
    {
      for (const auto& edge: m_explorer.out_edges(s))
      {
        if (is_tau(edge.action) && visited.find(edge.state) == visited.end())
        {
          visited.insert(edge.state);
          result.push_back(edge.state);
          next_frontier.push_back(edge.state);
        }
      }
    }
    frontier = std::move(next_frontier);
  }

  return result;
}

std::vector<lps::state> state_set::tau_closed_states(std::size_t tau_depth)
{
  return tau_closure(m_states, tau_depth);
}

bool state_set::is_action_enabled(const std::string& label, std::size_t tau_depth)
{
  for (const lps::state& s: tau_closed_states(tau_depth))
  {
    for (const auto& edge: m_explorer.out_edges(s))
    {
      if (!is_tau(edge.action) && wire_label(edge.action) == label)
      {
        return true;
      }
    }
  }
  return false;
}

bool state_set::apply_action(const std::string& label, std::size_t tau_depth)
{
  const std::vector<lps::state> closed = tau_closed_states(tau_depth);

  // Collect all immediate successors via `label`.
  std::vector<lps::state> post;
  for (const lps::state& s: closed)
  {
    for (const auto& edge: m_explorer.out_edges(s))
    {
      if (!is_tau(edge.action) && wire_label(edge.action) == label)
      {
        post.push_back(edge.state);
      }
    }
  }

  if (post.empty())
  {
    return false;
  }

  // Deduplicate
  std::sort(post.begin(), post.end());
  post.erase(std::unique(post.begin(), post.end()), post.end());

  m_states = tau_closure(post, tau_depth);
  return true;
}

state_set::enabled_actions state_set::get_enabled(std::size_t tau_depth, const io_classifier& classifier)
{
  enabled_actions result;
  std::unordered_set<std::string> seen_inputs;
  std::unordered_set<std::string> seen_outputs;

  for (const lps::state& s: tau_closed_states(tau_depth))
  {
    bool quiescent = true; // innocent until a tau or output edge is found

    for (const auto& edge: m_explorer.out_edges(s))
    {
      if (is_tau(edge.action))
      {
        quiescent = false;
        continue;
      }

      const std::string name = first_label(edge.action); // for io classification
      const std::string wire = wire_label(edge.action); // label returned on wire
      const auto type = classifier.classify(name);

      if (type == io_classifier::action_type::input && seen_inputs.insert(wire).second)
      {
        result.inputs.push_back(wire);
      }
      else if (type == io_classifier::action_type::output)
      {
        quiescent = false;
        if (seen_outputs.insert(wire).second)
        {
          result.outputs.push_back(wire);
        }
      }
    }

    if (quiescent)
    {
      result.quiescence = true;
      // Keep iterating — other states may still contribute inputs/outputs.
    }
  }

  return result;
}

bool state_set::refine_to_quiescent(std::size_t tau_depth, const io_classifier& classifier)
{
  const std::vector<lps::state> closed = tau_closed_states(tau_depth);

  std::vector<lps::state> quiescent;
  for (const lps::state& s: closed)
  {
    bool is_quiescent = true;
    for (const auto& edge: m_explorer.out_edges(s))
    {
      if (is_tau(edge.action) || classifier.is_output(first_label(edge.action)))
      {
        is_quiescent = false;
        break;
      }
    }
    if (is_quiescent)
    {
      quiescent.push_back(s);
    }
  }

  if (quiescent.empty())
  {
    return false;
  }

  m_states = std::move(quiescent);
  return true;
}

} // namespace mcrl2
