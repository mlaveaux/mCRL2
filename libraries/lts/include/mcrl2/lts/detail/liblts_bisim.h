// Author(s): Muck van Weerdenburg, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts/detail/liblts_bisim.h

#ifndef MCRL2_LTS_DETAIL_LIBLTS_BISIM_H
#define MCRL2_LTS_DETAIL_LIBLTS_BISIM_H

#include <fstream>
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/lts/lts_utilities.h"
#include "mcrl2/lts/detail/liblts_scc.h"
#include "mcrl2/lts/detail/liblts_merge.h"
#include "mcrl2/lts/lts_aut.h"
#include "mcrl2/lts/lts_fsm.h"
#include "mcrl2/lts/lts_dot.h"

namespace mcrl2::lts::detail
{

template < class LTS_TYPE>
class bisim_partitioner
{

  public:
    /** \brief Creates a bisimulation partitioner for an LTS.
     *  \details This bisimulation partitioner applies the algorithm
     *  defined in J.F. Groote and F.W. Vaandrager. An efficient algorithm for branching bisimulation and stuttering
     *  equivalence. In M.S. Paterson, editor, Proceedings 17th ICALP, Warwick, volume 443 of Lecture Notes in Computer
     *  Science, pages 626-638. Springer-Verlag, 1990. The only difference is that this algorithm uses actions labels
     *  on transitions. Therefore, each list of non_inert_transitions is grouped such that all transitions with the
     *  same label are grouped together. Tau transitions (which have as label index the number of labels) occur at the
     *  beginning of this list.
     *
     *  If branching is true, then branching bisimulation is used, otherwise strong bisimulation is applied.
     *  If preserve_divergence is true, then branching must be true. In this case states with an internal tau loop
     *  are considered to be different from states without a tau loop. In this way the divergences are preserved.
     *
     *  The input transition system is not allowed to contain tau loops, except that if preserve_divergence is true
     *  tau loops on a single state are allowed as they indicate divergences. Using the scc partitioner the tau
     *  loops must first be removed before applying this algorithm.
     *  \warning Note that when compiled with optimisations, bisimulation partitioning
     *  is much faster than compiled without any optimisation. The difference can go up to a factor 10.
     *  \param[in] l Reference to the LTS. The LTS l is only changed if \ref replace_transition_system is called. */
    bisim_partitioner(LTS_TYPE& l,
        const bool branching = false,
        const bool preserve_divergence = false,
        const bool generate_counter_examples = false)
        : aut(l),
          branching(branching),
          store_counter_info(generate_counter_examples)
    {
      assert(branching || !preserve_divergence);
      mCRL2log(log::verbose) << (preserve_divergence?"Divergence preserving b":"B") <<
                  (branching?"ranching b":"") << "isimulation partitioner created for "
                  << l.num_states() << " states and " <<
                  l.num_transitions() << " transitions\n";

      if (generate_counter_examples)
      {
        outgoing_transitions = transitions_per_outgoing_state_action_pair(aut.get_transitions(), aut.hidden_label_set());
      }

      create_initial_partition(branching,preserve_divergence);
      refine_partition_until_it_becomes_stable(branching, preserve_divergence);
    }


    /** \brief Destroys this partitioner. */
    ~bisim_partitioner()=default;

    /** \brief Replaces the transition relation of the current lts by the transitions
     *         of the bisimulation reduced transition system.
     * \details Each transition (s,l,s') is replaced by a transition (t,l,t') where
     * t and t' are the equivalence classes to which classes of the LTS. If the label l is
     * internal, which is detected using the function is_tau, then it is only returned
     * if t!=t' or preserve_divergence=true. This effectively removes all inert transitions.
     * Duplicates are removed from the transitions in the new lts.
     * Also the number of states and the initial state are adapted by this method.
     *
     * \pre The bisimulation equivalence classes have been computed.
     * \param[in] branching Causes non internal transitions to be removed.
     * \param[in] preserve_divergences Preserves tau loops on states. */
    void replace_transition_system(const bool branching, const bool preserve_divergences)
    {
      // Put all the non inert transitions in a set. Add the transitions that form a self
      // loop. Such transitions only exist in case divergence preserving branching bisimulation is
      // used. A set is used to remove double occurrences of transitions.
      std::set < transition > resulting_transitions;

      for (const transition& i: aut.get_transitions())
      {
        const bool is_transition_i_hidden=aut.is_tau(aut.apply_hidden_label_map(i.label()));
        if (!branching ||
            !is_transition_i_hidden ||
            get_eq_class(i.from())!=get_eq_class(i.to()) ||
            (preserve_divergences && i.from()==i.to()))
        {
          resulting_transitions.insert(
            transition(
              get_eq_class(i.from()),
// In the line below all hidden transitions are replaced by an explicit tau. It is possible to
// always use i.label() where the resulting transition is only implicitly hidden. This has as effect that
// the labels of non inert hidden transitions are preserved (e.g. for use in counter examples) but that
// there are possibly multiple hidden transitions between states, leading to a larger number of
// transitions than strictly necessary.
              (is_transition_i_hidden?aut.tau_label_index():i.label()),
              get_eq_class(i.to())));
        }
      }
      // Remove the old transitions
      aut.clear_transitions();

      // Copy the transitions from the set into the transition system.
      for (const transition& t: resulting_transitions)
      {
        aut.add_transition(t);
      }

      // Merge the states, by setting the state labels of each state to the concatenation of the state labels of its
      // equivalence class.
      if (aut.has_state_info())   /* If there are no state labels this step can be ignored */
      {
        /* Create a vector for the new labels */
        std::vector<typename LTS_TYPE::state_label_t> new_labels(num_eq_classes());

        for(std::size_t i=aut.num_states(); i>0; )
        {
          --i;
          const std::size_t new_index=get_eq_class(i);
          new_labels[new_index] = new_labels[new_index] + aut.state_label(i);
        }

        aut.set_num_states(num_eq_classes());
        for(std::size_t i=0; i<num_eq_classes(); ++i)
        {
          aut.set_state_label(i,new_labels[i]);
        }
      }
      else
      {
        aut.set_num_states(num_eq_classes());
      }

      aut.set_initial_state(get_eq_class(aut.initial_state()));
    }

    /** \brief Gives the number of bisimulation equivalence classes of the LTS.
     *  \return The number of bisimulation equivalence classes of the LTS.
     */
    std::size_t num_eq_classes() const
    {
      return max_state_index;
    }


    /** \brief Gives the bisimulation equivalence class number of a state.
     *  \param[in] s A state number.
     *  \return The number of the bisimulation equivalence class to which \e s belongs. */
    std::size_t get_eq_class(const std::size_t s) const
    {
      assert(s<block_index_of_a_state.size());
      return blocks[block_index_of_a_state[s]].state_index;
    }


    /** \brief Returns whether two states are in the same bisimulation equivalence class.
     *  \param[in] s A state number.
     *  \param[in] t A state number.
     *  \retval true if \e s and \e t are in the same bisimulation equivalence class;
     *  \retval false otherwise. */
    bool in_same_class(const std::size_t s, const std::size_t t) const
    {
      return get_eq_class(s)==get_eq_class(t);
    }


    /** \brief Creates a state formula that distinguishes state s from state t.
     *  \details The states s and t are non bisimilar states. If they are
     *           bisimilar an exception is raised. A counter state formula phi is
     *           returned, which has the property that s |= phi and not t |= phi.
     *           Based on "Computing Distinguishing Formulas for Branching Bisimulation", 1991 by Henri Korver
     *  \param[in] s The state number for which the resulting formula should be true
     *  \param[in] t The state number for which the resulting formula should be false
     *  \return A distinguishing state formula. */
    mcrl2::state_formulas::state_formula counter_formula(std::size_t s, std::size_t t);

  private:

    using block_index_type = std::size_t;
    using state_type = std::size_t;
    using label_type = std::size_t;

    state_type max_state_index = 0;
    LTS_TYPE& aut;
    bool branching;
    bool store_counter_info;

    struct non_bottom_state
    {
      state_type state;
      std::vector < state_type > inert_transitions; // Only the target state is interesting.

      non_bottom_state(const state_type s)
        : state(s)
      {}
      non_bottom_state(const state_type s, const std::vector < state_type > &it)
        : state(s), inert_transitions(it)
      {}
    };

    struct block
    {
      state_type state_index = 0UL;              // The state number that represent the states in this block
      block_index_type block_index = 0UL;        // The sequence number of this block.
      block_index_type parent_block_index = 0UL; // Index of the parent block.
      // If there is no parent block, this refers to the block
      // itself.
      std::vector < state_type > bottom_states; // The non bottom states must be ordered
      // on tau reachability. The deepest
      // states occur first in the vector.
      std::vector < non_bottom_state > non_bottom_states;
      // The non_inert transitions are grouped per label. The (non-inert) transitions
      // with tau labels are at the end of this vector.
      std::vector < transition > non_inert_transitions;

      void swap(block& b) noexcept
      {
        state_type state_index1=b.state_index;
        b.state_index=state_index;
        state_index=state_index1;

        block_index_type block_index1=b.block_index;
        b.block_index=block_index;
        block_index=block_index1;

        block_index_type parent_block_index1=b.parent_block_index;
        b.parent_block_index=parent_block_index;
        parent_block_index=parent_block_index1;

        bottom_states.swap(b.bottom_states);
        non_bottom_states.swap(b.non_bottom_states);
        non_inert_transitions.swap(b.non_inert_transitions);
      }
    };

    std::vector < block > blocks;

    // std::vector < bool > block_is_active;       // Indicates whether this is still a block in the partition.
    // Blocks that are split become inactive.
    std::vector < state_type > block_index_of_a_state;
    std::vector < bool > block_flags;
    std::vector < bool > block_is_in_to_be_processed;
    std::vector < bool > state_flags;

    std::vector< block_index_type > to_be_processed;
    std::vector< block_index_type > BL;

    // Counter example information based on "Computing Distinguishing Formulas for Branching Bisimulation", 1991 by Henri Korver
    // Given a block B that is split using action alpha and block B' into blocks L and R, with P as the partition just before this split, then
    // right_child[B] = R
    std::map < block_index_type, block_index_type > right_child;
    // split_by_action[B] = alpha,
    std::map < block_index_type, label_type > split_by_action;
    // split_by_block[B] = B',
    std::map < block_index_type, block_index_type > split_by_block;
    // r_alpha[R] = { C \in P | \exists_{s \in R, s' \in C} : s -a-> s'} and
    std::map < block_index_type, std::set < block_index_type > > r_alpha;
    // r_tau[R] = { C \in P | \exists_{s \in R, s' \in C} : s -tau-> s' && C != B}.
    std::map < block_index_type, std::set < block_index_type > > r_tauP;
    // map from source state and action to target state, makes generating counterexample info easier
    outgoing_transitions_per_state_action_t outgoing_transitions;

    // A counter for creating fresh variables
    int fresh_var_counter = 0;


    void create_initial_partition(const bool branching,
                                  const bool preserve_divergences)

    {
      to_be_processed.clear();

      block initial_partition;

      // First store the bottom and non bottom states.
      sort_transitions(aut.get_transitions(), aut.hidden_label_set(), mcrl2::lts::src_lbl_tgt);

      state_type last_non_stored_state_number=0;
      bool bottom_state=true;
      std::vector < state_type > current_inert_transitions;

      {
        // Reserve enough space, such that no reallocations of the vector are required when adding transitions.
        // For this purpose, first the number of inert transitions must be counted, to avoid reserving too much
        // space. This for instance leads to a waste of memory (terabytes for reducing 30M states), especially,
        // when calculating ia strong bisimulation reduction.
        std::size_t initial_partition_non_inert_counter=0;
        std::size_t current_inert_transition_counter=0;
        const std::vector<transition> & trans=aut.get_transitions();
        for (auto t: trans)
        {
          if (branching && aut.is_tau(aut.apply_hidden_label_map(t.label())))
          {
            if (preserve_divergences && t.from()==t.to())
            {
              initial_partition_non_inert_counter++;
            }
            else
            {
              current_inert_transition_counter++;
            }
          }
        }
        current_inert_transitions.reserve(initial_partition_non_inert_counter);
        initial_partition.non_inert_transitions.reserve(current_inert_transition_counter);
      }

      const std::vector<transition>& trans=aut.get_transitions();
      for (std::vector<transition>::const_iterator r=trans.begin(); r!=trans.end(); ++r)
      {
        const transition t=* r;

        if (branching && aut.is_tau(aut.apply_hidden_label_map(t.label())))
        {
          if (preserve_divergences && t.from()==t.to())
          {
            initial_partition.non_inert_transitions.push_back(transition(t.from(),aut.tau_label_index(),t.to()));
          }
          else
          {
            current_inert_transitions.push_back(t.to());
            bottom_state=false;
          }
        }
        std::vector<transition>::const_iterator next_i=r;
        ++next_i;
        if (next_i==trans.end() || t.from()!=next_i->from())
        {
          // store the current from state
          for (; last_non_stored_state_number<t.from(); ++last_non_stored_state_number)
          {
            initial_partition.bottom_states.push_back(last_non_stored_state_number);
          }

          if (bottom_state)
          {
            initial_partition.bottom_states.push_back(t.from());
          }
          else
          {
            initial_partition.non_bottom_states.push_back(non_bottom_state(t.from()));
            current_inert_transitions.swap(initial_partition.non_bottom_states.back().inert_transitions);
            bottom_state=true;
          }
          assert(last_non_stored_state_number==t.from());
          last_non_stored_state_number++;
        }
      }

      for (; last_non_stored_state_number<aut.num_states(); ++last_non_stored_state_number)
      {
        initial_partition.bottom_states.push_back(last_non_stored_state_number);
      }

      // Sort the non-bottom states such that the deepest states occur first.
      // Raise an exception if there is a non trivial tau loop.

      order_on_tau_reachability(initial_partition.non_bottom_states);

      // Store the non-inert transitions (i.e. the non tau transitions)
      sort_transitions(aut.get_transitions(), aut.hidden_label_set(), mcrl2::lts::lbl_tgt_src);
      const std::vector<transition> & trans1=aut.get_transitions();
      for (auto t: trans1)
      {
        if (!branching || !aut.is_tau(aut.apply_hidden_label_map(t.label())))
        {
          // Note that by sorting the transitions first, the non_inert_transitions are grouped per label.
          initial_partition.non_inert_transitions.push_back(t);
        }
      }

      initial_partition.block_index=0;
      initial_partition.state_index=0;
      max_state_index=1;
      initial_partition.parent_block_index=0;
      blocks.push_back(block());
      blocks.back().swap(initial_partition);
      block_index_of_a_state=std::vector < block_index_type >(aut.num_states(),0);
      block_flags.push_back(false);
      state_flags=std::vector < bool >(aut.num_states(),false);
      block_is_in_to_be_processed.push_back(false);
      to_be_processed.clear();
      BL.clear();
    } // end create_initial_partition

// Refine the partition until the partition has become stable
#ifndef NDEBUG
    void refine_partition_until_it_becomes_stable(const bool branching, const bool preserve_divergence)
#else
    void refine_partition_until_it_becomes_stable(const bool, const bool)
#endif
    {
#ifndef NDEBUG
      std::size_t consistency_check_counter=1;
      std::size_t consistency_check_barrier=1;
#endif
      bool partition_is_unstable=true; // This boolean indicates that the partition becomes unstable
      // because an inert transition becomes non inert.
      while (!to_be_processed.empty() || partition_is_unstable)
      {
#ifndef NDEBUG
        // Avoid checking too often. This is too time consuming, even in debug mode.
        consistency_check_counter++;
        if (consistency_check_counter>=consistency_check_barrier)
        {
          consistency_check_counter=0;
          consistency_check_barrier=consistency_check_barrier*2;
          check_internal_consistency_of_the_partitioning_data_structure(branching, preserve_divergence);
        }
#endif
        if (to_be_processed.empty() && partition_is_unstable)
        {
          // Put all blocks in to_be_processed;
          for (block_index_type i=0; i< blocks.size(); ++i)
          {
            to_be_processed.push_back(i);
          }
          block_is_in_to_be_processed=std::vector < bool >(blocks.size(),true);
          partition_is_unstable=false;
        }

        const block_index_type splitter_index=to_be_processed.back();
        to_be_processed.pop_back();
        block_is_in_to_be_processed[splitter_index]=false;

        // Split with the splitter block, unless it is to_be_processed as we have to reconsider it
        // completely anyhow at some later point.

        // Walk through the non_inert_transitions via a natural number, because
        // the blocks and the non_inert_transitions can move around in memory. This makes it
        // unsafe to use interators.

        for (std::size_t i=0; i<blocks[splitter_index].non_inert_transitions.size(); ++i)
        {
          // The flag of the starting state of *i is raised and its block is added to BL;

          assert(blocks[splitter_index].non_inert_transitions[i].from()<aut.num_states());
          state_flags[blocks[splitter_index].non_inert_transitions[i].from()]=true;
          const block_index_type marked_block_index=block_index_of_a_state[blocks[splitter_index].non_inert_transitions[i].from()];
          if (block_flags[marked_block_index]==false)
          {
            block_flags[marked_block_index]=true;
            BL.push_back(marked_block_index);
          }

          // If the label of the next action is different, we must carry out the splitting.
          if ((i+1)==blocks[splitter_index].non_inert_transitions.size() ||
              aut.apply_hidden_label_map(blocks[splitter_index].non_inert_transitions[i].label())!=
              aut.apply_hidden_label_map(blocks[splitter_index].non_inert_transitions[i+1].label()))
          {
            // We consider BL which contains references to all blocks from which a state from splitter
            // can be reached. If not all flags of the non bottom states in a block in BL are set, the
            // non flagged non bottom states are moved to a new block.

            split_the_blocks_in_BL(partition_is_unstable,aut.apply_hidden_label_map(blocks[splitter_index].non_inert_transitions[i].label()),splitter_index);

          }

        }
      }
#ifndef NDEBUG
      check_internal_consistency_of_the_partitioning_data_structure(branching, preserve_divergence);
#endif
      block_flags.clear();
      block_is_in_to_be_processed.clear();
      state_flags.clear();
      to_be_processed.clear();
      BL.clear();
    }

    void refine_partion_with_respect_to_divergences();

    void split_the_blocks_in_BL(
      bool& partition_is_unstable,
      const label_type splitter_label,
      const block_index_type splitter_block)
    {
      for (block_index_type i1: BL)
      {
        block_flags[i1] = false;
        std::vector < state_type > flagged_states;
        std::vector < state_type > non_flagged_states;
        std::vector < state_type > i1_bottom_states;
        i1_bottom_states.swap(blocks[i1].bottom_states);

        for (state_type bottom_state: i1_bottom_states)
        {
          if (state_flags[bottom_state])
          {
            // state is flagged.
            flagged_states.push_back(bottom_state);
          }
          else
          {
            // state is not flagged. It will be moved to a new block.
            non_flagged_states.push_back(bottom_state);
            block_index_of_a_state[bottom_state] = blocks.size();
          }
        }
        assert(!flagged_states.empty()||!blocks[i1].non_bottom_states.empty()||i1_bottom_states.empty());
        block_index_type reset_state_flags_block = i1;

        if (!non_flagged_states.empty())
        {
          // There are flagged and non flagged states. So, the block must be split.
          // Move the unflagged states to the new block.

          if (mCRL2logEnabled(log::debug))
          {
            const std::size_t m = static_cast<std::size_t>(std::pow(10.0, std::floor(std::log10(static_cast<double>((blocks.size()+1)/2)))));
            if ((blocks.size()+1)/2 % m==0)
            {
              mCRL2log(log::debug) << "Bisimulation partitioner: create block " << (blocks.size()+1)/2 << std::endl;
            }
          }

          // Create a first new block.
          blocks.push_back(block());
          block_index_type new_block1=blocks.size()-1;
          blocks.back().state_index=max_state_index;
          max_state_index++;
          blocks.back().block_index=new_block1;
          blocks.back().parent_block_index = i1;

          non_flagged_states.swap(blocks.back().bottom_states);
          // Put the indices of first split block to to_be_processed.
          to_be_processed.push_back(blocks.back().block_index);
          block_is_in_to_be_processed.push_back(true);

          // Create a second new block.
          blocks.push_back(block());
          block_index_type new_block2=blocks.size()-1;
          blocks.back().state_index = blocks[i1].state_index;
          blocks.back().block_index=new_block2;
          reset_state_flags_block=new_block2;
          blocks.back().parent_block_index = i1;

          // Move the flagged states to the second block, and let the block index of these states refer to this block.
          flagged_states.swap(blocks.back().bottom_states);
          std::vector < state_type > &reference_to_flagged_states_of_block2=blocks.back().bottom_states;
          for (state_type j: reference_to_flagged_states_of_block2)
          {
            block_index_of_a_state[j] = new_block2;
          }

          // NOTES: new_block1 = R, new_block2 = L
          // Store counter formula info
          if (store_counter_info)
          {
            right_child[i1] = new_block1;
            split_by_action[i1] = splitter_label;
            split_by_block[i1] = splitter_block;

            // compute r_alpha[new_block1]
            std::set<block_index_type> reachable_blocks = {};
            for (state_type source_state : blocks[new_block1].bottom_states)
            {
              for (outgoing_transitions_per_state_action_t::const_iterator i = outgoing_transitions.lower_bound(std::pair<state_type, label_type>(source_state, splitter_label));
                i != outgoing_transitions.upper_bound(std::pair<state_type, label_type>(source_state, splitter_label)); ++i)
              {
                state_type target_state = to(i);
                block_index_type target_block = block_index_of_a_state[target_state];
                reachable_blocks.insert(target_block == new_block1 || target_block == new_block2 ? i1 : target_block);
              }
            }
            r_alpha[new_block1] = reachable_blocks;

            if (branching) {
              // compute r_tau[new_block1]
              std::set<block_index_type> tau_reachable_blocks = {};
              for (state_type source_state : blocks[new_block1].bottom_states)
              {
                for (label_type lab = 0; lab < aut.num_action_labels(); ++lab)
                {
                  if (aut.is_tau(aut.apply_hidden_label_map(lab)))
                  {
                    for (outgoing_transitions_per_state_action_t::const_iterator i = outgoing_transitions.lower_bound(std::pair<state_type, label_type>(source_state, lab));
                      i != outgoing_transitions.upper_bound(std::pair<state_type, label_type>(source_state, lab)); ++i)
                    {
                      state_type target_state = to(i);
                      block_index_type target_block = block_index_of_a_state[target_state];
                      if (!(target_block == i1 || target_block == new_block1 || target_block == new_block2))
                      {
                        tau_reachable_blocks.insert(target_block);
                      }
                    }
                  }
                }
              }
              r_tauP[new_block1] = tau_reachable_blocks;
            }
          }

          // Put the indices of second split block to to_be_processed.
          to_be_processed.push_back(blocks.back().block_index);
          block_is_in_to_be_processed.push_back(true);

          // reset the flag of block *i1, which is being split.
          block_is_in_to_be_processed[i1] = false;

          // The flag fields of the new blocks is set to false;
          block_flags.push_back(false);
          block_flags.push_back(false);

          // Declare already some space for transitions, such that we can
          // put inert transitions that become non inert in there when investigating
          // the non bottom states. After investigating the non bottom states,
          // the transitions are split over the vectors below.

          std::vector < transition > flagged_non_inert_transitions;
          std::vector < transition > non_flagged_non_inert_transitions;
          // Reserve enough space for transitions to be copied. Otherwise, resizing may lead to
          // lot of unneccesary copying...
          flagged_non_inert_transitions.reserve(blocks[i1].non_inert_transitions.size());
          non_flagged_non_inert_transitions.reserve(blocks[i1].non_inert_transitions.size());

          // Next we scan the non-bottom states of *i1. If for some non-bottom state the flag is not raised
          // and if none of the outgoing P-inert transitions leads to a state in the old block then this
          // state becomes a non bottom state of B2.

          std::vector < non_bottom_state > flagged_non_bottom_states;
          std::vector < non_bottom_state > non_flagged_non_bottom_states;
          std::vector < non_bottom_state > i1_non_bottom_states;
          i1_non_bottom_states.swap(blocks[i1].non_bottom_states);
          for (typename std::vector < non_bottom_state >::iterator k=i1_non_bottom_states.begin();
               k!=i1_non_bottom_states.end(); ++k)
          {
            const std::vector < state_type > &inert_transitions=k->inert_transitions;
            if (!state_flags[k->state])
            {
              bool all_transitions_end_in_unflagged_block=true;
              for (std::vector < state_type > :: const_iterator l=inert_transitions.begin();
                   all_transitions_end_in_unflagged_block && l!=inert_transitions.end(); ++l)
              {
                if (block_index_of_a_state[*l]!= new_block1)
                {
                  block_index_of_a_state[*l]=new_block2;
                  all_transitions_end_in_unflagged_block=false;
                }
              }
              if (all_transitions_end_in_unflagged_block)
              {
                // Move *k to the non flagged block. Swap the inert transitions to avoid copying.
                non_bottom_state s(k->state);
                s.inert_transitions.swap(k->inert_transitions);
                non_flagged_non_bottom_states.push_back(s);
                block_index_of_a_state[k->state]=new_block1;
                continue;
              }
            }
            // Move *k to the flagged block; note that the transitions can have become
            // non-inert. So, investigate them separately.
            std::vector < state_type > remaining_inert_transitions;
            for (state_type inert_transition: inert_transitions)
            {
              if (block_index_of_a_state[inert_transition] == new_block1)
              {
                // The transition *l (*k,tau_label,*l) becomes a non inert transition in the new
                // block.
                non_flagged_non_inert_transitions.push_back(
                  transition(k->state, aut.tau_label_index(), inert_transition));
              }
              else
              {
                // The transition represented by *l remains an inert transition.
                block_index_of_a_state[inert_transition] = new_block2;
                remaining_inert_transitions.push_back(inert_transition);
              }
            }
            if (remaining_inert_transitions.empty()) // The last outgoing inert tau transition just became non inert.
              // k->state has become a bottom state. Otherwise it remains
              // a non bottom state.
            {
              blocks[new_block2].bottom_states.push_back(k->state);
              block_index_of_a_state[k->state]=new_block2;
              partition_is_unstable=true;
            }
            else
            {
              flagged_non_bottom_states.push_back(non_bottom_state(k->state,remaining_inert_transitions));
              block_index_of_a_state[k->state]=new_block2;
            }
          }
          non_flagged_non_bottom_states.swap(blocks[new_block1].non_bottom_states);
          flagged_non_bottom_states.swap(blocks[new_block2].non_bottom_states);

          // Finally the non-inert transitions are distributed over both blocks in the obvious way.
          // Note that this must be done after all states are properly put into a new block.

          assert(i1 < blocks.size());
          std::vector < transition > i1_non_inert_transitions;
          i1_non_inert_transitions.swap(blocks[i1].non_inert_transitions);
          for (auto& i1_non_inert_transition: i1_non_inert_transitions)
          {
            if (block_index_of_a_state[i1_non_inert_transition.to()] == new_block1)
            {
              non_flagged_non_inert_transitions.push_back(i1_non_inert_transition);
            }
            else
            {
              block_index_of_a_state[i1_non_inert_transition.to()] = new_block2;
              flagged_non_inert_transitions.push_back(i1_non_inert_transition);
            }
          }

          // Only put the parts that we need to store. Avoid that reserved space for
          // current_inert_transition_counter and flagged_non_inert_transitions is stored
          // in the non_inert_transitions in both blocks. Therefore, we do not use a copy
          // constructor, but put elements in there one by one.
          // The two lines below were old code, wasting memory bandwidth.
          // non_flagged_non_inert_transitions.swap(blocks[new_block1].non_inert_transitions);
          // flagged_non_inert_transitions.swap(blocks[new_block2].non_inert_transitions);

          blocks[new_block1].non_inert_transitions.reserve(non_flagged_non_inert_transitions.size());
          for (const auto& non_flagged_non_inert_transition: non_flagged_non_inert_transitions)
          {
            blocks[new_block1].non_inert_transitions.push_back(non_flagged_non_inert_transition);
          }

          blocks[new_block2].non_inert_transitions.reserve(flagged_non_inert_transitions.size());
          for (const auto& flagged_non_inert_transition: flagged_non_inert_transitions)
          {
            blocks[new_block2].non_inert_transitions.push_back(flagged_non_inert_transition);
          }
        }
        else
        {
          // Nothing changed, so put the bottom states back again.
          i1_bottom_states.swap(blocks[i1].bottom_states);
        }
        // reset the state flags
        std::vector < state_type > &flagged_states_to_be_unflagged=blocks[reset_state_flags_block].bottom_states;
        for (state_type j: flagged_states_to_be_unflagged)
        {
          state_flags[j] = false;
        }

        std::vector < non_bottom_state > &flagged_states_to_be_unflagged1=blocks[reset_state_flags_block].non_bottom_states;
        for (typename std::vector < non_bottom_state >::const_iterator j=flagged_states_to_be_unflagged1.begin();
             j!=flagged_states_to_be_unflagged1.end(); ++j)
        {
          state_flags[j->state]=false;
        }
      }
      BL.clear();
    }

    void order_recursively_on_tau_reachability(
      const state_type s,
      std::map < state_type, std::vector < state_type > >& inert_transition_map,
      std::vector < non_bottom_state >& new_non_bottom_states,
      std::set < state_type >& visited)
    {
      if (inert_transition_map.count(s)>0) // The state s is a bottom state. We need not to investigate these.
      {
        if (visited.count(s)==0)
        {
          visited.insert(s);
          std::vector < state_type> &inert_transitions=inert_transition_map[s];
          for (state_type inert_transition: inert_transitions)
          {
            order_recursively_on_tau_reachability(inert_transition,
              inert_transition_map,
              new_non_bottom_states,
              visited);
          }
          new_non_bottom_states.push_back(non_bottom_state(s));
          inert_transitions.swap(new_non_bottom_states.back().inert_transitions);
        }
      }
    }

    void order_on_tau_reachability(std::vector < non_bottom_state > &non_bottom_states)
    {
      std::set < state_type > visited;
      std::map < state_type, std::vector < state_type > > inert_transition_map;
      for (typename std::vector < non_bottom_state >::iterator i=non_bottom_states.begin();
           i!=non_bottom_states.end(); ++i)
      {
        i->inert_transitions.swap(inert_transition_map[i->state]);
      }
      std::vector < non_bottom_state > new_non_bottom_states;

      for (typename std::vector < non_bottom_state >::const_iterator i=non_bottom_states.begin();
           i!=non_bottom_states.end(); ++i)
      {
        order_recursively_on_tau_reachability(i->state, inert_transition_map, new_non_bottom_states, visited);
      }
      new_non_bottom_states.swap(non_bottom_states);
    }

    /**
     * \brief conjunction Creates a conjunction of state formulas
     * \param[in] terms The terms of the conjunction
     * \return The conjunctive state formula
     */
    mcrl2::state_formulas::state_formula conjunction(std::set<mcrl2::state_formulas::state_formula> terms) const
    {
      return utilities::detail::join<mcrl2::state_formulas::state_formula>(terms.begin(), terms.end(),
        [](mcrl2::state_formulas::state_formula a, mcrl2::state_formulas::state_formula b)
        { return mcrl2::state_formulas::and_(a, b); }, mcrl2::state_formulas::true_());
    }

    /**
     * \brief create_regular_formula Creates a regular formula that represents action a
     * \details In case the action comes from an LTS in the lts format.
     * \param[in] a The action for which to create a regular formula
     * \return The created regular formula
     */
    regular_formulas::regular_formula create_regular_formula(const mcrl2::lps::multi_action& a) const
    {
      return regular_formulas::regular_formula(action_formulas::multi_action(a.actions()));
    }

    /**
     * \brief create_regular_formula Creates a regular formula that represents action a
     * \details In case the action comes from an LTS in the aut or fsm format.
     * \param[in] a The action for which to create a regular formula
     * \return The created regular formula
     */
    regular_formulas::regular_formula create_regular_formula(const mcrl2::lts::action_label_string& a) const
    {
      return mcrl2::regular_formulas::regular_formula(mcrl2::action_formulas::multi_action(
        process::action_list({ process::action(process::action_label(a, {}), {}) })));
    }

    /**
     * \brief until_formula Creates a state formula that corresponds to the until operator phi1<a>phi2 from HMLU
     * \details This operator intuitively means: "phi1 holds while stuttering until we can do an a-step after which phi2 holds"
     *          In the operators of the mu-calculus that mCRL2 supports we can define this as:
     *              phi2 || (mu X.phi1 && (\<tau\>X || \<a\>phi2))  if a = tau
     *              mu X.phi1 && (\<tau\>X || \<a\>phi2)            else
     * \param[in] phi1 The first state formula for the until operator
     * \param[in] a The action for the until operator
     * \param[in] phi2 The second state formula for the until operator
     * \return A state formula that corresponds to the until operator phi1<a>phi2 from HMLU
     */
    mcrl2::state_formulas::state_formula until_formula(const mcrl2::state_formulas::state_formula& phi1, const label_type& a,
                                                       const mcrl2::state_formulas::state_formula& phi2)
    {
      std::string var = "X" + std::to_string(fresh_var_counter++);
      mcrl2::state_formulas::state_formula tauStep =
          mcrl2::state_formulas::may(regular_formulas::regular_formula(action_formulas::multi_action()), mcrl2::state_formulas::variable(var, {}));
      mcrl2::state_formulas::state_formula lastStep = mcrl2::state_formulas::may(create_regular_formula(aut.action_label(a)), phi2);

      mcrl2::state_formulas::state_formula until =
          mcrl2::state_formulas::mu(var, {}, mcrl2::state_formulas::and_(phi1, mcrl2::state_formulas::or_(tauStep, lastStep)));

      if (aut.is_tau(aut.apply_hidden_label_map(a)))
      {
        until = mcrl2::state_formulas::or_(phi2, until);
      }
      return until;
    }

    /* \brief Creates a state formula that distinguishes states in block B1 from state in block B2.
       \details Based on "Computing Distinguishing Formulas for Branching Bisimulation", 1991 by Henri Korver.
                Variable names used below correspond to definitions in this paper.
       \param[in] B1 The block on which the resulting formula must be true.
       \param[in] B2 The block on which the resulting formula must be false.
       \return A distinguishing state formula.
    */
    mcrl2::state_formulas::state_formula counter_formula_aux(const block_index_type B1, const block_index_type B2)
    {
      // First find the smallest block containing both B1 and B2.
      // Find all blocks containing B1.
      std::set < block_index_type > blocks_containing_B1;
      block_index_type B1p = B1;
      blocks_containing_B1.insert(B1p);
      while (blocks[B1p].parent_block_index != B1p)
      {
        B1p = blocks[B1p].parent_block_index;
        blocks_containing_B1.insert(B1p);
      }

      // Go up the block tree from B2 until we find a block that contains B1.
      block_index_type DB = B2;
      while (blocks_containing_B1.count(DB) == 0)
      {
        assert(blocks[DB].parent_block_index != DB);
        DB = blocks[DB].parent_block_index;
      }

      // Now DB is the smallest block containing both B1 and B2.

      const block_index_type R = right_child.at(DB);
      const label_type a = split_by_action.at(DB);
      const block_index_type Bp = split_by_block.at(DB);

      // We distinguish Bp with every block in r_alpha[R]
      std::set<mcrl2::state_formulas::state_formula> Gamma2;
      for (block_index_type PPB : r_alpha.at(R))
      {
        Gamma2.insert(counter_formula_aux(Bp, PPB));
      }
      mcrl2::state_formulas::state_formula Phi2 = conjunction(Gamma2);

      mcrl2::state_formulas::state_formula Phi;
      if (branching)
      {
        // We distinguish DB with every block in r_tauP[R]
        std::set<mcrl2::state_formulas::state_formula> Gamma1;
        for (block_index_type PPB : r_tauP.at(R))
        {
          Gamma1.insert(counter_formula_aux(DB, PPB));
        }
        mcrl2::state_formulas::state_formula Phi1 = conjunction(Gamma1);

        // In case a == tau, we also need to distinguish Bp with R
        if (aut.is_tau(aut.apply_hidden_label_map(a)))
        {
          Phi2 = mcrl2::state_formulas::and_(Phi2, counter_formula_aux(Bp, R));
        }

        Phi = until_formula(Phi1, a, Phi2);
      }
      else
      {
        Phi = mcrl2::state_formulas::may(create_regular_formula(aut.action_label(a)), Phi2);
      }

      if (blocks_containing_B1.count(R) == 0)
      {
        return Phi;
      }
      else
      {
        return mcrl2::state_formulas::not_(Phi);
      }
    }


#ifndef NDEBUG
    // The method below is intended to check the consistency of the internal data
    // structure. Its sole purpose is to detect programming errors. It has no
    // side effects on the data structure. If a problem occurs, execution halts with
    // an assert.

    void check_internal_consistency_of_the_partitioning_data_structure(
      const bool branching,
      const bool preserve_divergence) const
    {
      state_type total_number_of_states=0;
      std::size_t total_number_of_transitions=0;

      assert(!blocks.empty());
      std::set < block_index_type > block_indices;

      assert(block_index_of_a_state.size()==aut.num_states());
      for (typename std::vector < block >::const_iterator i=blocks.begin();
           i!=blocks.end(); ++i)
      {
        // Check the block_index.
        assert(i->block_index<blocks.size());
        assert(block_indices.count(i->block_index)==0);
        block_indices.insert(i->block_index);

        // Check the bottom states.
        const std::vector < state_type > &i_bottom_states=i->bottom_states;

        for (state_type bottom_state: i_bottom_states)
        {
          total_number_of_states++;
          assert(bottom_state < aut.num_states());
          // Check that the block number of the state is maintained properly.
          assert(block_index_of_a_state[bottom_state] == i->block_index);
        }

        // Check the non bottom states. In particular check that there is no tau loop
        // in these non bottom states.
        const std::vector < non_bottom_state > &i_non_bottom_states=i->non_bottom_states;
        std::set < state_type > visited;
        std::set < state_type > local_bottom_states;
        for (typename std::vector < non_bottom_state >::const_iterator j=i_non_bottom_states.begin();
             j!=i_non_bottom_states.end(); ++j)
        {
          local_bottom_states.insert(j->state);
        }

        for (typename std::vector < non_bottom_state >::const_iterator j=i_non_bottom_states.begin();
             j!=i_non_bottom_states.end(); ++j)
        {
          total_number_of_states++;
          assert(j->state<aut.num_states());
          // Check that the block number of the state is maintained properly.
          assert(block_index_of_a_state[j->state]==i->block_index);
          const std::vector < state_type > &j_inert_transitions=j->inert_transitions;
          for (state_type inert_transition: j_inert_transitions)
          {
            total_number_of_transitions++;
            assert(inert_transition < aut.num_states());
            // Check that the inert transitions are well ordered.
            assert(visited.count(inert_transition) > 0 || local_bottom_states.count(inert_transition) == 0);
          }
          visited.insert(j->state);
        }

        // Check the non_inert_transitions. It is required that the transitions
        // are grouped per label, and that tau transitions must be inert.

        const std::vector < transition > &i_non_inert_transitions=i->non_inert_transitions;
        std::set < label_type > observed_action_labels;
        for (std::vector < transition >::const_iterator j=i_non_inert_transitions.begin();
             j!=i_non_inert_transitions.end(); ++j)
        {
          total_number_of_transitions++;
          assert(j->to()<aut.num_states());
          assert(j->from()<aut.num_states());

          // Check proper grouping of action labels.
          std::vector < transition >::const_iterator j_next=j;
          j_next++;
          if (j_next==i_non_inert_transitions.end() || (aut.apply_hidden_label_map(j->label())!=aut.apply_hidden_label_map(j_next->label())))
          {
            assert(observed_action_labels.count(aut.apply_hidden_label_map(j->label()))==0);
            observed_action_labels.insert(aut.apply_hidden_label_map(j->label()));
          }

          // Check whether tau transition in non inert transition vector is inert.
          if (!preserve_divergence && branching && aut.is_tau(aut.apply_hidden_label_map(j->label())))
          {
            assert(j->to()!=j->from());
          }

          // Check whether the target state of the transition is in the current block.
          assert(block_index_of_a_state[j->to()]==i->block_index);
        }
      }

      // Check total number of states and transitions.
      assert(total_number_of_states==aut.num_states());
      assert(total_number_of_transitions==aut.num_transitions());

      // Check block_index_of_a_state
      assert(block_index_of_a_state.size()==aut.num_states());
      for (block_index_type i: block_index_of_a_state)
      {
        assert(blocks[i].block_index == i);
      }

      // Check block_flags that the block flags are all set to false
      for (bool block_flag: block_flags)
      {
        assert(!block_flag);
      }

      // Check that state_flags are all false.
      for (bool state_flag: state_flags)
      {
        assert(!state_flag);
      }

      // Check to_be_processed
      // Check block_is_in_to_be_processed
      std::vector < bool > temporary_block_is_in_to_be_processed(blocks.size(),false);

      for (block_index_type i: to_be_processed)
      {
        temporary_block_is_in_to_be_processed[i] = true;
      }
      for (state_type i=0; i<blocks.size(); ++i)
      {
        assert(!block_is_in_to_be_processed[i] || temporary_block_is_in_to_be_processed[i]);
      }

      // Check that BL is empty.
      assert(BL.empty());

    }

#endif // not NDEBUG

};


/** \brief Reduce transition system l with respect to strong or (divergence preserving) branching bisimulation.
 * \param[in,out] l The transition system that is reduced.
 * \param[in] branching If true branching bisimulation is applied, otherwise strong bisimulation.
 * \param[in] preserve_divergences Indicates whether loops of internal actions on states must be preserved. If false
 *            these are removed. If true these are preserved.  */
template<class LTS_TYPE>
void bisimulation_reduce(LTS_TYPE& l, bool branching = false, bool preserve_divergences = false);

/** \brief Checks whether the two initial states of two lts's are strong or branching bisimilar.
 * \details This lts and the lts l2 are not usable anymore after this call.
 *          The space consumption is O(n) and time is O(nm). It uses the branching bisimulation
 *          algorithm by Groote and Vaandrager from 1990.
 * \param[in,out] l1 A first transition system.
 * \param[in,out] l2 A second transition system.
 * \param[in] branching If true branching bisimulation is used, otherwise strong bisimulation is applied.
 * \param[in] preserve_divergences If true and branching is true, preserve tau loops on states.
 * \param[in] generate_counter_examples Whether to generate a counter example
 * \param[in] counter_example_file The file to store the counter example in
 * \param[in] structured_output
 * \retval True iff the initial states of the current transition system and l2 are (divergence preserving) (branching) bisimilar */
template<class LTS_TYPE>
bool destructive_bisimulation_compare(LTS_TYPE& l1,
  LTS_TYPE& l2,
  bool branching = false,
  bool preserve_divergences = false,
  bool generate_counter_examples = false,
  const std::string& counter_example_file = "",
  bool structured_output = false);

/** \brief Checks whether the two initial states of two lts's are strong or branching bisimilar.
 * \details The current transitions system and the lts l2 are first duplicated and subsequently
 *          reduced modulo bisimulation. If memory space is a concern, one could consider to
 *          use destructive_bisimulation_compare. This routine uses the Groote-Vaandrager
 *          branching bisimulation routine. It runs in O(mn) and uses O(n) memory where n is the
 *          number of states and m is the number of transitions.
 * \param[in,out] l1 A first transition system.
 * \param[in,out] l2 A second transistion system.
 * \param[in] branching If true branching bisimulation is used, otherwise strong bisimulation is applied.
 * \param[in] preserve_divergences If true and branching is true, preserve tau loops on states.
 * \param[in] generate_counter_examples Whether to generate a counter example
 * \param[in] counter_example_file The file to store the counter example in
 * \param[in] structured_output
 * \retval True iff the initial states of the current transition system and l2 are (divergence preserving) (branching) bisimilar */
template<class LTS_TYPE>
bool bisimulation_compare(const LTS_TYPE& l1,
  const LTS_TYPE& l2,
  bool branching = false,
  bool preserve_divergences = false,
  bool generate_counter_examples = false,
  const std::string& counter_example_file = "",
  bool structured_output = false);

template < class LTS_TYPE>
mcrl2::state_formulas::state_formula bisim_partitioner<LTS_TYPE>::counter_formula(
  const std::size_t s, const std::size_t t)
{
  if (get_eq_class(s)==get_eq_class(t))
  {
    throw mcrl2::runtime_error("Requesting a counter state formula for two bisimilar states. Such a state formula is not useful.");
  }

  return counter_formula_aux(block_index_of_a_state[s], block_index_of_a_state[t]);
}


template < class LTS_TYPE>
void bisimulation_reduce(LTS_TYPE& l,
                         const bool branching /*=false */,
                         const bool preserve_divergences /*=false */)
{
  // First, remove tau loops in case of branching bisimulation.
  if (branching)
  {
    scc_reduce(l,preserve_divergences);
  }

  // Secondly, apply the branching bisimulation reduction algorithm. If there are no tau's,
  // this will automatically yield strong bisimulation.
  detail::bisim_partitioner<LTS_TYPE> bisim_part(l, branching, preserve_divergences);

  // Assign the reduced LTS
  bisim_part.replace_transition_system(branching,preserve_divergences);
}

template < class LTS_TYPE>
bool bisimulation_compare(
  const LTS_TYPE& l1,
  const LTS_TYPE& l2,
  const bool branching /* =false*/,
  const bool preserve_divergences /*=false*/,
  const bool generate_counter_examples /*= false*/,
  const std::string& /* counter_example_file = ""*/,
  const bool structured_output /* = false */)
{
  LTS_TYPE l1_copy(l1);
  LTS_TYPE l2_copy(l2);
  return destructive_bisimulation_compare(l1_copy,l2_copy,branching,preserve_divergences,
                                          generate_counter_examples, structured_output);
}

/** \brief Checks whether the two initial states of two lts's are strong or branching bisimilar.
 * \details This lts and the lts l2 are not usable anymore after this call.
 *          The space consumption is O(n) and time is O(nm). It uses the branching bisimulation
 *          algorithm by Groote and Vaandrager from 1990.
 * \param[in/out] l1 A first transition system.
 * \param[in/out] l2 A second transition system.
 * \param[in] branching If true branching bisimulation is used, otherwise strong bisimulation is applied.
 * \param[in] preserve_divergences If true and branching is true, preserve tau loops on states.
 * \param[in] generate_counter_examples Whether to generate a counter example
 * \param[in] counter_example_file The file to store the counter example in
 * \param[in] structured_output
 * \retval True iff the initial states of the current transition system and l2 are (divergence preserving) (branching) bisimilar */
template < class LTS_TYPE>
bool destructive_bisimulation_compare(
  LTS_TYPE& l1,
  LTS_TYPE& l2,
  const bool branching /* =false*/,
  const bool preserve_divergences /*=false*/,
  const bool generate_counter_examples /* = false */,
  const std::string& counter_example_file /*= ""*/,
  const bool /*structured_output = false */)
{
  std::size_t init_l2 = l2.initial_state() + l1.num_states();
  mcrl2::lts::detail::merge(l1,l2);
  l2.clear(); // No use for l2 anymore.

  // First remove tau loops in case of branching bisimulation.
  if (branching)
  {
    detail::scc_partitioner<LTS_TYPE> scc_part(l1);
    scc_part.replace_transition_system(preserve_divergences);
    init_l2 = scc_part.get_eq_class(init_l2);
  }

  detail::bisim_partitioner<LTS_TYPE> bisim_part(l1, branching, preserve_divergences, generate_counter_examples);
  if (generate_counter_examples && !bisim_part.in_same_class(l1.initial_state(), init_l2))
  {
    mcrl2::state_formulas::state_formula counter_example_formula = bisim_part.counter_formula(l1.initial_state(), init_l2);

    std::string filename = "Counterexample.mcf";
    if (!counter_example_file.empty())
    {
      filename = counter_example_file;
    }
    std::ofstream counter_file(filename);
    counter_file << mcrl2::state_formulas::pp(counter_example_formula);
    counter_file.close();
    mCRL2log(mcrl2::log::info) << "Saved counterexample to: \"" << filename << "\"" << std::endl;

  }
  return bisim_part.in_same_class(l1.initial_state(),init_l2);
}


}

#endif // MCRL2_LTS_DETAIL_LIBLTS_BISIM_H
