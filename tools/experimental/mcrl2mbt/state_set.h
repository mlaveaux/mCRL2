// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <string>
#include <vector>

#include "mcrl2/data/rewriter.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/lps/explorer.h"
#include "mcrl2/lps/explorer_options.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/state.h"

#include "io_classifier.h"

namespace mcrl2
{

/// \brief Maintains a set of LPS states S and supports IOCO-style operations.
///
/// Operations correspond to the MBT protocol semantics:
///  - tau-closure: S := tau*_k(S)
///  - post_action:  S := tau*_k(post_a(tau*_k(S)))
///  - quiescent subset: Q = { s in tau*_k(S) | no output or tau edge }
class state_set
{
public:
  using explorer_type = lps::explorer<false, false, lps::specification>;

  /// \brief Aggregated result of get_enabled().
  struct enabled_actions
  {
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    bool quiescence = false;
  };

  state_set(const lps::specification& spec, const lps::explorer_options& options, const data::rewriter& rewr);

  /// \brief Reset S to the singleton initial state.
  void reset_to_initial();

  /// \brief Check whether action `label` is enabled from tau*_k(S).
  bool is_action_enabled(const std::string& label, std::size_t tau_depth);

  /// \brief S := tau*_k(post_label(tau*_k(S))).
  /// \returns false if the action was not enabled (S unchanged).
  bool apply_action(const std::string& label, std::size_t tau_depth);

  /// \brief Return all enabled observable actions from tau*_k(S).
  enabled_actions get_enabled(std::size_t tau_depth, const io_classifier& classifier);

  /// \brief Refine S to quiescent subset. Returns true iff S is non-empty after refinement.
  ///
  /// S is updated to { s in tau*_k(S) | no output edge and no tau edge }.
  bool refine_to_quiescent(std::size_t tau_depth, const io_classifier& classifier);

  std::size_t size() const { return m_states.size(); }

private:
  /// \brief Compute tau*_k(seeds) — BFS following tau edges up to depth k.
  std::vector<lps::state> tau_closure(const std::vector<lps::state>& seeds, std::size_t depth);

  /// \brief Compute tau*_k(m_states).
  std::vector<lps::state> tau_closed_states(std::size_t tau_depth);

  /// \brief Compute the initial lps::state from the spec's initial expressions.
  lps::state compute_initial_state() const;

  data::rewriter m_rewr;
  explorer_type m_explorer;
  std::vector<lps::state> m_states;
};

} // namespace mcrl2
