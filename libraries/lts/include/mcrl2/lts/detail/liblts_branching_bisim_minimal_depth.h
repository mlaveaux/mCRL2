// Author(s): Jan Martens and Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
/// \file mcrl2/lts/detail/liblts_branching_bisim_minimal_depth.h
///
/// \brief Partition refinement algorithm for guaruanteed minimal depth
/// counter-examples.
///
/// \details
#ifndef MCRl2_LTS_LIBLTS_BRANCHING_BISIM_MINIMAL_DEPTH
#define MCRl2_LTS_LIBLTS_BRANCHING_BISIM_MINIMAL_DEPTH


#include "mcrl2/lts/detail/liblts_merge.h"
#include "mcrl2/lts/detail/liblts_scc.h"
#include "mcrl2/lts/detail/liblts_bisim_dnj.h"
#include "mcrl2/lts/lts_aut.h"
#include "mcrl2/lts/lts_dot.h"
#include "mcrl2/lts/lts_fsm.h"
#include "mcrl2/lts/lts_utilities.h"
#include "mcrl2/modal_formula/state_formula.h"

#include <fstream>

namespace mcrl2::lts::detail
{

template <class LTS_TYPE>
class branching_bisim_partitioner_minimal_depth
{
public:
  /** \brief Creates a branching bisimulation partitioner for an LTS.
   *  \details This partitioner is specifically for creating minimal depth counter-examples for branching bisimulation.
   *           It guarantees stability w.r.t. the old partition before considering new splitter blocks. This might cause
   *           this implementation to be less efficient than other partitioners.
   *  \param l Reference to the LTS.
   *  \param init_l2 reference to the initial state of lts2.
   */
  branching_bisim_partitioner_minimal_depth(LTS_TYPE& l, const std::size_t init_l2)
      : m_lts(l),
        initial_l2(init_l2)
  {

    for (transition trans : m_lts.get_transitions())
    {
      trans_out[trans.from()].emplace_back(trans.label(), trans.to());
      if (is_tau(trans.label()))
      {
        silent_out[trans.from()].insert(trans.to());
        silent_in[trans.to()].insert(trans.from());
      }
    }

    for (state_type s = 0; s < m_lts.num_states(); s++)
    {
      state2block[s] = 0;

      if (silent_out[s].size() == 0)
      {
        bottom_states.push_back(s);
      }
    }



    // Initialize the partition with a single block.
    blocks.push_back(block());
    blocks[0].state_index = 0;
    blocks[0].block_index = 0;
    blocks[0].parent_block_index = 0;
    blocks[0].level = 0;
    level2blocksidx[0].insert(0);

    // Now we can start refining the partition.
    std::size_t num_old_blocks = 0;
    std::size_t num_blocks_created = 1;
    std::size_t level = 0;

    while (num_blocks_created > num_old_blocks && in_same_class(m_lts.initial_state(), initial_l2))
    {
      level += 1;
      num_old_blocks = num_blocks_created;
      num_blocks_created = refine_partition();
      assert(level2blocksidx[level].size() == num_blocks_created);
      state2sig = std::map<state_type, signature_type>();
      mCRL2log(mcrl2::log::verbose) << "Refined partition to " << num_blocks_created
        << " blocks on level " << level << "."
        << std::endl;
    }
  }

  /** \brief Creates a state formula that distinguishes state s from state t.
   *  \details The states s and t are non branching bisimilar states. A distinguishign state formula phi is
   *           returned, which has the property that \f$s \in \llbracket\phi\rrbracket\f$ and \f$t \notin \llbracket\phi\rrbracket\f$. Based on the preprint "Minimal Depth Distinguishing Formulas without Until for
   * Branching Bisimulation", 2024 by Jan Martens and Jan Friso Groote.
   *  \param[in] s The state number for which the resulting formula should be true
   *  \param[in] t The state number for which the resulting formula should be false
   *  \return A minimal observation depth distinguishing state formula, that is often also minimum negation-depth and
   * irreducible. */
  mcrl2::state_formulas::state_formula dist_formula_mindepth(size_t s, size_t t)
  {
    assert(s == m_lts.initial_state() && t == initial_l2);
    assert(state2block[s] != state2block[t]);
    std::pair<block_index_type, block_index_type> b1b2 = min_split_blockpair(state2block[s], state2block[t]);
    return dist_formula(b1b2.first, b1b2.second);
  }

  bool in_same_class(const std::size_t s, const std::size_t t) { return state2block[s] == state2block[t]; }


private:
  LTS_TYPE& m_lts;
  state_type initial_l2;
  state_type max_state_index = 0;
  using block_index_type = std::size_t;
  using level_type = std::size_t;
  // This tuple is meant for an observation (s', a, s'') such that s -(silent)-> s' -a-> s''.
  using branching_observation_type = std::tuple<block_index_type, label_type, block_index_type>;
  using signature_type = std::set<branching_observation_type>;
  using observation = std::pair<label_type, state_type>;
  using blockpair_type = std::pair<block_index_type, block_index_type>;

  std::map<state_type, std::set<state_type>> silent_in;
  std::map<state_type, std::set<state_type>> silent_out;
  std::map<state_type, std::vector<observation>> trans_out;
  std::map<state_type, std::size_t> state2num_touched;
  std::map<state_type, block_index_type> state2block;
  std::map<state_type, signature_type> state2sig;
  std::vector<state_type> bottom_states;
  std::map<level_type, std::set<block_index_type>> level2blocksidx;

  std::map<std::pair<block_index_type, block_index_type>, mcrl2::state_formulas::state_formula> blockpair2formula;
  std::map<std::pair<block_index_type, block_index_type>, std::set<block_index_type>> blockpair2truths;

  struct block
  {
    state_type state_index = 0UL;       // The state number that represent the states in this block
    block_index_type block_index = 0UL; // The sequence number of this block.
    block_index_type parent_block_index
        = 0UL;              // Index of the parent block. If there is no parent block, this refers to the block itself.
    level_type level = 0UL; // The level of the block in the partition.
    signature_type sig;

    void swap(block& b) noexcept
    {
      std::swap(b.state_index, state_index);
      std::swap(b.block_index, block_index);
      std::swap(b.parent_block_index, parent_block_index);
      std::swap(b.level, level);
    }

    bool operator==(const block& other)
    {
      return block_index == other.block_index;
    }

    bool operator!=(const block& other)
    {
      return !(*this == other);
    }
  };
  std::vector<block> blocks;

  /*
  * Auxiliary function that computes whether a label is a tau index.
  */
  bool is_tau(label_type l)
  {
    return m_lts.is_tau(m_lts.apply_hidden_label_map(l));
  }

  signature_type get_signature(state_type s)
  {
    signature_type sig;
    // Add the block index of the state to the signature.
    sig.insert(std::make_tuple(state2block[s], m_lts.tau_label_index(), state2block[s]));
    for (state_type target : silent_out[s])
    {
      sig.insert(state2sig[target].begin(), state2sig[target].end());
    }
    for (observation t : trans_out[s])
    {
      if (!is_tau(t.first) || state2block[s] != state2block[t.second])
      {
        sig.insert(std::make_tuple(state2block[s], t.first, state2block[t.second]));
      }
    }

    state2sig[s] = sig;
    return sig;
  }

  // Refine the partition exactly one level.
  std::size_t refine_partition()
  {
    std::queue<state_type> frontier;
    // start with bottom states.
    for (auto s : bottom_states)
    {
      frontier.push(s);
    }
    std::map<signature_type, block_index_type> sig2block;
    std::map<state_type, block_index_type> state2block_new;
    std::size_t num_blocks_created = 0;
    // Compute signatures in order of bottom states and reachability order.
    while (!frontier.empty())
    {
      state_type state = frontier.front();
      frontier.pop();
      signature_type sig = get_signature(state);
      if (sig2block.find(sig) == sig2block.end())
      {
        // Create the new block
        blocks.push_back(block());
        num_blocks_created += 1;
        block_index_type new_block_id = blocks.size() - 1;
        sig2block[sig] = new_block_id;
        blocks[new_block_id].state_index = state;
        blocks[new_block_id].block_index = new_block_id;
        blocks[new_block_id].parent_block_index = state2block[state];
        blocks[new_block_id].level = blocks[state2block[state]].level + 1;
        blocks[new_block_id].sig = sig;
        level2blocksidx[blocks[new_block_id].level].insert(new_block_id);
      }

      state2block_new[state] = sig2block[sig];
      state2num_touched[state] = 0;

      for (state_type backward : silent_in[state])
      {
        size_t max_out = silent_out[backward].size();
        state2num_touched[backward] += 1;
        if (state2num_touched[backward] == max_out)
        {
          frontier.push(backward);
        }
      }
    }
    // Now we have computed the new blocks, we can redefine the partition.
    state2block = state2block_new;
    return num_blocks_created;
  }

  std::pair<block_index_type, block_index_type> min_split_blockpair(block_index_type b1, block_index_type b2)
  {
    assert(blocks[b1] != blocks[b2] && blocks[b1].level == blocks[b2].level);
    while (blocks[b1].parent_block_index != blocks[b2].parent_block_index)
    {
      b1 = blocks[b1].parent_block_index;
      b2 = blocks[b2].parent_block_index;
    }
    return std::make_pair(b1, b2);
  }

    /**
   * \brief conjunction Creates a conjunction of state formulas
   * \param terms The terms of the conjunction
   * \return The conjunctive state formula
   */
  mcrl2::state_formulas::state_formula conjunction(std::vector<mcrl2::state_formulas::state_formula>& conjunctions)
  {
    return utilities::detail::join<mcrl2::state_formulas::state_formula>(
        conjunctions.begin(),
        conjunctions.end(),
        [](mcrl2::state_formulas::state_formula a, mcrl2::state_formulas::state_formula b)
        { return mcrl2::state_formulas::and_(a, b); },
        mcrl2::state_formulas::true_());
  }

  regular_formulas::regular_formula make_tau_hat(regular_formulas::regular_formula& f)
  {
    return regular_formulas::alt(f,
        regular_formulas::trans_or_nil(regular_formulas::regular_formula(action_formulas::false_())));
  }

  /**
   * \brief create_regular_formula Creates a regular formula that represents action a
   * \details In case the action comes from an LTS in the lts format.
   * \param a The action for which to create a regular formula
   * \return The created regular formula
   */
  regular_formulas::regular_formula create_regular_formula(const mcrl2::lts::action_label_string& a) const
  {
    // Copied from Olav no idea what it exactly does.
    return mcrl2::regular_formulas::regular_formula(mcrl2::action_formulas::multi_action(
        process::action_list({process::action(process::action_label(a, {}), {})})));
  }

  /**
   * \brief create_regular_formula Creates a regular formula that represents action a
   * \details In case the action comes from an LTS in the lts format.
   * \param a The action for which to create a regular formula
   * \return The created regular formula
   */
  regular_formulas::regular_formula create_regular_formula(const mcrl2::lps::multi_action& a) const
  {
    return regular_formulas::regular_formula(action_formulas::multi_action(a.actions()));
  }

  // Updates the truths values according the dist formula and the information in blockpair2truths.
  void split_and_intersect(std::set<block_index_type>& truths, std::pair<block_index_type, block_index_type> liftedB1B2)
  {
    // Add the observation to the formula
    level_type split_level = blocks[liftedB1B2.first].level;

    std::set<block_index_type> truths_new;
    for (block_index_type b_og : truths)
    {
      block_index_type b = b_og;
      // Check if b in phi (maybe lift the level)
      while (blocks[b].level > split_level)
      {
        b = blocks[b].parent_block_index;
      }
      if (blockpair2truths[liftedB1B2].find(b) != blockpair2truths[liftedB1B2].end())
      {
        truths_new.insert(b_og);
      }
    }
    truths = truths_new;
  }

  /**
   * \brief is_dist Checks if a given conjunction correctly exludes a set of blocks.
   * \param dist_blockpairs The blockpairs that were used to generate the conjuncts.
   * \param to_dist The set of blocks that should be excluded by the conjunction, forall b\in to_dist.
   * b\not\in\sem{phi_dist_blockpairs}. \return True if the conjunction correctly excludes the set of blocks, false
   * otherwise.
   */
  bool is_dist(std::set<blockpair_type>& dist_blockpairs,
    std::set<block_index_type>& to_dist)
  {
    if (to_dist.empty())
    {
      return true;
    }
    std::set<block_index_type> truths = level2blocksidx[blocks[*to_dist.begin()].level];
    return is_dist(dist_blockpairs, to_dist, truths);
  }


  /**
  * \brief is_dist overloaded to also maintain the truth values computed at the end.
  */
  bool is_dist(const std::set<blockpair_type>& dist_blockpairs,const std::set<block_index_type>& to_dist, std::set<block_index_type>& truths)
  {
    if (to_dist.empty())
    {
      return true;
    }

    truths = level2blocksidx[blocks[*to_dist.begin()].level];
    for (blockpair_type b1_b2 : dist_blockpairs)
    {
      split_and_intersect(truths, b1_b2);
    }
    for (block_index_type b : to_dist)
    {
      if (truths.find(b) != truths.end())
      {
        return false;
      }
    }
    return true;
  }


  std::vector<mcrl2::state_formulas::state_formula> filtered_dist_conjunction(
    std::map<blockpair_type, mcrl2::state_formulas::state_formula>& Phi,
    std::set<block_index_type>& Tdist,
    std::set<block_index_type>& Truths)
  {
    std::vector<mcrl2::state_formulas::state_formula> returnPhi;
    std::set<blockpair_type> phi_pairs_og;
    for (auto& phi_pair : Phi) {
      phi_pairs_og.insert(phi_pair.first);
    }

    for (auto& phi_pair : Phi) {
      // Remove the phi_pair from the formula.
      phi_pairs_og.erase(phi_pair.first);
      // Only add conjunct, if it wouldn't be distinguishing without.
      if (!is_dist(phi_pairs_og, Tdist))
      {
        phi_pairs_og.insert(phi_pair.first);
        returnPhi.push_back(phi_pair.second);
      }
    }
    is_dist(phi_pairs_og, Tdist, Truths);

    return returnPhi;
  }

  // This function computes the distinguishing state formula for two blocks.
  // Precondition is that the blocks are not the same, are on the same level and have the same parent block.
  mcrl2::state_formulas::state_formula dist_formula(block_index_type block_index1, block_index_type block_index2)
  {
    assert(block_index1 != block_index2);
    if (blockpair2formula.find(std::make_pair(block_index1, block_index2)) != blockpair2formula.end())
    {
      // We computed this already ( probably won't happen much in practice but guarantees polynomial runtime).
      return blockpair2formula[std::make_pair(block_index1, block_index2)];
    }

    block block1 = blocks[block_index1];
    block block2 = blocks[block_index2];

    // make blocks same level
    assert(block1.level == block2.level); // This should be true, otherwise need to make them same level.
    assert(block1.parent_block_index == block2.parent_block_index); // This should be true, otherwise need to make them same level.)

    signature_type ds = block1.sig;
    signature_type dt = block2.sig;

    // Find a distinguishing observation s- tau ->> s' -a-> s''
    std::tuple<block_index_type, label_type, block_index_type> dist_obs;
    bool found_obs = false;
    for (auto path : ds)
    {
      if (!is_tau(std::get<1>(path)) or std::get<0>(path) != std::get<2>(path))
      {
        if (dt.find(path) == dt.end())
        {
          dist_obs = path;
          found_obs = true;
          break;
        }
      }
    }
    // If no such observation exists, we flip the blocks.
    if (!found_obs)
    {
      auto phi = mcrl2::state_formulas::not_(dist_formula(block2.block_index, block1.block_index));
      blockpair2formula[std::make_pair(block1.block_index, block2.block_index)] = phi;
      // Compute the truth values \sem{!\phi} = lvl2blocksidx[lvl] - \sem{\phi}
      std::set_difference(level2blocksidx[block1.level].begin(),
          level2blocksidx[block1.level].end(),
          blockpair2truths[std::make_pair(block2.block_index, block1.block_index)].begin(),
          blockpair2truths[std::make_pair(block2.block_index, block1.block_index)].end(),
          std::inserter(blockpair2truths[std::make_pair(block1.block_index, block2.block_index)],
              blockpair2truths[std::make_pair(block1.block_index, block2.block_index)].begin()));
      return phi;
    }


    // We have a distinguishing observation, start constructing the formula.
    label_type dist_label = std::get<1>(dist_obs);
    block_index_type B1 = std::get<0>(dist_obs);
    block_index_type B2 = std::get<2>(dist_obs);
    // Log the observation for debugging purposes.
    std::vector<std::pair<block_index_type, block_index_type>> T;

    for (auto path : dt)
    {
      if (std::get<1>(path) == dist_label and (!is_tau(dist_label) or std::get<0>(path) != std::get<2>(path)))
      {
        // we might have B_1 -\tau -> B_1, I don't think its a problem.
        T.emplace_back(std::get<0>(path), std::get<2>(path));
        if (is_tau(dist_label))
        {
          // the following observation: block2 -\tau->> path2 -(tau)-> path2.
          T.emplace_back(std::get<2>(path), std::get<2>(path));
        }
      }
    }
    if (is_tau(dist_label))
    {
      // the following observation: block2 -\tau->> path2 -(tau)-> path2.
      T.push_back(std::make_pair(block2.parent_block_index, block2.parent_block_index));
    }

    //Keep a copy of T for backwards filtering in postprocessing.
    std::set<blockpair_type> T_og;
    for (auto t : T)
    {
      T_og.insert(t);
    }

    // TODO: Remove this part and make T a set.
    // Sort T, such that the block with the highest distlevel in s'' get dealt with first.
    // This is a heuristic, no idea if it improves much or can be improved
    std::sort(T.begin(), T.end(),
      [this, B2](std::pair<block_index_type, block_index_type> a, std::pair<block_index_type, block_index_type> b)
        {
        if (a.second == B2) {
          return false;
        }
        if (b.second == B2)
        {
          return true;
        }
        auto alift = min_split_blockpair(a.second, B2);
        auto blift = min_split_blockpair(b.second, B2);
        return blocks[alift.first].level < blocks[blift.first].level;
      });

    // <tau*>(<dist_label> phi1 && phi2), We remember the original keys for truth values and backwards filtering.
    std::map<blockpair_type, mcrl2::state_formulas::state_formula> Phi1;
    std::map<blockpair_type, mcrl2::state_formulas::state_formula> Phi2;

    // Track truth values for the formula s -\tau -> Truths2 [phi2] -dist_label-> Truths1 [phi1]
    // This is a bit confusing, we flip order here in the path to the formula.
    std::set<block_index_type> Truths1 = level2blocksidx[block1.level - 1];
    std::set<block_index_type> Truths2 = level2blocksidx[block1.level - 1];

    while (!T.empty())
    {
      std::pair<block_index_type, block_index_type> Bt1_Bt2 = T.back();
      // We could pop_back here on T, but the computation of the truth also handles this.
      if (Bt1_Bt2.second != B2)
      {
        std::pair<block_index_type, block_index_type> liftedPair = min_split_blockpair(B2, Bt1_Bt2.second);
        Phi1[liftedPair] = dist_formula(liftedPair.first, liftedPair.second);
        // Update truthvalues for the formula <(dist_label)> phi1.
        split_and_intersect(Truths1, liftedPair);
      }
      else
      {
        std::pair<block_index_type, block_index_type> liftedPair = min_split_blockpair(B1, Bt1_Bt2.first);
        Phi2[liftedPair] = dist_formula(liftedPair.first, liftedPair.second);
        // Update truthvalues for the formula phi2.
        split_and_intersect(Truths2, liftedPair);
      }
      // Remove observations (Bt1, Bt2) from  T of which Bt2 is not in phi1 or Bt1 is not in phi2
      std::vector<std::pair<block_index_type, block_index_type>> T_new;
      for (auto bt1_bt2 : T)
      {
        // T_new only consists stuch that both Bt1 and Bt2 are still not distinguished.
        if (Truths1.find(bt1_bt2.second) != Truths1.end() && Truths2.find(bt1_bt2.first) != Truths2.end())
        {
          T_new.push_back(bt1_bt2);
        }
      }
      T = T_new;
    }

    // Backwards filtering, remove each formula and see if it is still distinguishing, starting with Phi2.
    std::set<block_index_type> Tset;
    //We assume Phi2.
    for (auto B_Bp : T_og)
    {
      if (Truths2.find(B_Bp.first) != Truths2.end()) {
        Tset.insert(B_Bp.second);
      }
    }

    std::vector<mcrl2::state_formulas::state_formula> returnPhi1 = filtered_dist_conjunction(Phi1, Tset, Truths1);


    Tset.clear();
    //Now we assume phi1.
    for (auto B_Bp : T_og)
    {
      if (Truths1.find(B_Bp.second) != Truths1.end()) {
        Tset.insert(B_Bp.first);
      }
    }
    std::vector<mcrl2::state_formulas::state_formula> returnPhi2 = filtered_dist_conjunction(Phi2, Tset, Truths2);

    //Done with formula, set the truth values for the formula.
    // Initialize the truth values for the formula
    blockpair2truths[std::make_pair(block_index1, block_index2)] = std::set<block_index_type>();
    for (block_index_type b : level2blocksidx[block1.level])
    {
      signature_type sig = blocks[b].sig;
      for (auto path : sig)
      {
        if (std::get<1>(path) == dist_label)
        {
          block_index_type Bt1 = std::get<0>(path);
          block_index_type Bt2 = std::get<2>(path);
          if (Truths1.find(Bt2) != Truths1.end() && Truths2.find(Bt1) != Truths2.end())
          {
            blockpair2truths[std::make_pair(block1.block_index, block2.block_index)].insert(b);
            break;
          }
        }
      }
    }

    // Consruct the regular formula for the diamond operator
    mcrl2::regular_formulas::regular_formula diamond = create_regular_formula(m_lts.action_label(dist_label));
    // If the action is tau, we need to add the tau_hat operator.
    if (is_tau(dist_label))
    {
      // we mimic <\hat{tau}> phi := <tau> phi || phi , by <tau+false*> phi.
      diamond = make_tau_hat(diamond);
    }

    // diamond formula <dist_label> phi1
    mcrl2::state_formulas::state_formula returnPhi = mcrl2::state_formulas::may(diamond , conjunction(returnPhi1));
    // diamond formula: <dist_label>phi1 && phi2
    if (!returnPhi2.empty())
    {
      returnPhi = mcrl2::state_formulas::and_(returnPhi, conjunction(returnPhi2));
    }

    // <tau*> (<dist_label> phi1 && phi2)
    blockpair2formula[std::make_pair(block1.block_index, block2.block_index)]
        = mcrl2::state_formulas::may(mcrl2::regular_formulas::trans_or_nil(
            create_regular_formula(m_lts.action_label(0))), returnPhi);
    return blockpair2formula[std::make_pair(block1.block_index, block2.block_index)];
  }

};

template <class LTS_TYPE>
bool destructive_branching_bisimulation_compare_minimal_depth(LTS_TYPE& l1,
    LTS_TYPE& l2,
    const std::string& counter_example_file)
{
  std::size_t init_l2 = l2.initial_state() + l1.num_states();
  mcrl2::lts::detail::merge(l1, l2);
  l2.clear(); // No use for l2 anymore.

  // First remove tau loops in case of branching bisimulation
  bool preserve_divergences = false;
  detail::scc_partitioner<LTS_TYPE> scc_part(l1);
  scc_part.replace_transition_system(preserve_divergences);

  init_l2 = scc_part.get_eq_class(init_l2);

   // Run a faster branching bisimulation algorithm as preprocessing, no preversing of loops.
  detail::bisim_partitioner_dnj branching_bisim_part(l1, true, preserve_divergences);
  if (branching_bisim_part.in_same_class(l1.initial_state(), init_l2))
  {
    return true;
  }

  init_l2 = branching_bisim_part.get_eq_class(init_l2);
  branching_bisim_part.finalize_minimized_LTS();

  //Start the new debugging to get minimal depth splitting information.
  branching_bisim_partitioner_minimal_depth<LTS_TYPE> branching_bisim_min(l1, init_l2);

  //Min-depth partitioner should agree with dnj.
  assert(!branching_bisim_min.in_same_class(l1.initial_state(), init_l2));

  // LTSs are not bisimilar, we can create a counter example.
  std::string filename = "Counterexample.mcf";
  if (!counter_example_file.empty())
  {
    filename = counter_example_file;
  }

  mcrl2::state_formulas::state_formula counter_example_formula
      = branching_bisim_min.dist_formula_mindepth(l1.initial_state(), init_l2);

  std::ofstream counter_file(filename);
  counter_file << mcrl2::state_formulas::pp(counter_example_formula);
  mCRL2log(mcrl2::log::info) << "Saved counterexample to: \"" << filename << "\"" << std::endl;
  return false;
}

} // namespace mcrl2::lts::detail
#endif // MCRl2_LTS_LIBLTS_BRANCHING_BISIM_MINIMAL_DEPTH
