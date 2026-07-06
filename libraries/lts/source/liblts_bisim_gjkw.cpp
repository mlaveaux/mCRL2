// Author(s): David N. Jansen, Radboud Universiteit, Nijmegen, The Netherlands
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/// \file liblts_bisim_gjkw.cpp
///
/// \brief O(m log n)-time stuttering equivalence algorithm
///
/// \details This file implements the efficient partition refinement algorithm
/// by Groote / Jansen / Keiren / Wijs to calculate the stuttering equivalence
/// quotient of a Kripke structure.  (Labelled transition systems are converted
/// to Kripke structures before the main algorithm).
/// The file accompanies the planned publication in the ACM Trans. Comput. Log.
/// Log. special issue for TACAS 2016, to appear in 2017.
///
/// \author David N. Jansen, Radboud Universiteit, Nijmegen, The Netherlands

#include "mcrl2/lts/detail/liblts_bisim_gjkw.h"
#include "mcrl2/lts/detail/coroutine.h"
#include "mcrl2/lts/lts_aut.h"
#include "mcrl2/lts/lts_fsm.h"
#include "mcrl2/lts/lts_utilities.h"



namespace mcrl2::lts::detail
{
namespace bisim_gjkw
{





/* ************************************************************************* */
/*                                                                           */
/*                   R E F I N A B L E   P A R T I T I O N                   */
/*                                                                           */
/* ************************************************************************* */





block_t* block_t::refinable_first = nullptr;
state_type block_t::nr_of_blocks = 0;
constln_t* constln_t::nontrivial_first = nullptr;
                                                                                #ifndef NDEBUG
                                                                                    // These variables are only accessed in debug mode.  In release mode,
                                                                                    // accessing them would lead to a linker error.
                                                                                    state_info_const_ptr state_info_entry::s_i_begin;
                                                                                    state_info_const_ptr state_info_entry::s_i_end;
                                                                                    permutation_const_iter_t block_t::perm_begin;
                                                                                #endif
/// \brief refine the block (the blue subblock is smaller)
/// \details This function is called after a refinement function has found
/// that the blue subblock is the smaller one.  It creates a new block for
/// the blue states.
/// \param blue_nonbottom_end iterator past the last blue non-bottom state
/// \returns pointer to the new (blue) block
block_t* block_t::split_off_blue(permutation_iter_t const blue_nonbottom_end)
{                                                                               assert(unmarked_nonbottom_end() >= blue_nonbottom_end);
                                                                                assert(unmarked_nonbottom_begin() <= blue_nonbottom_end);
                                                                                assert(0 != unmarked_bottom_size());
    permutation_iter_t const splitpoint = blue_nonbottom_end +
                                                        static_cast<permutation_iter_t::difference_type>(unmarked_bottom_size()); assert(splitpoint < end());  assert(begin() < splitpoint);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    unsigned const max_counter = check_complexity::log_n -
    /* It is not necessary to reset the nottoblue counters; these counters   */                   check_complexity::ilog2((state_type) (splitpoint - begin()));
    /* are anyway only valid for the maybe-blue states.                      */     assert((state_type) (splitpoint - begin()) <= size()/2);
                                                                                #endif
    if (state_type swapcount = std::min(unmarked_bottom_size(),
                   (state_type) (marked_nonbottom_end() - blue_nonbottom_end));
                                                                0 != swapcount)
    {
        // vector swap the states:
        permutation_iter_t pos1 = blue_nonbottom_end;
        permutation_iter_t pos2 = unmarked_bottom_end();
        state_info_ptr const temp = *pos1;
        for (;;)
        {
            --pos2;                                                             mCRL2complexity(*pos2, add_work(check_complexity::
                                                                                               Move_Blue_or_Red_to_a_new_block_NewB_swap_3_29, max_counter), );
            *pos1 = *pos2;
            (*pos1)->pos = pos1;
            ++pos1;
            if (0 == --swapcount)
            {
              break;
            }
            *pos2 = *pos1;
            (*pos2)-> pos = pos2;
        }
        *pos2 = temp;
        (*pos2)->pos = pos2;
    }

    // create a new block for the blue states
    block_t* const NewB = new block_t(constln(), begin(), splitpoint);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    NewB->work_counter = work_counter;
                                                                                #endif
    if (BLOCK_NO_SEQNR != seqnr())
    {
        NewB->assign_seqnr();
    }
    // NewB->set_begin(begin());
    NewB->set_bottom_begin(blue_nonbottom_end);
    NewB->set_marked_nonbottom_begin(blue_nonbottom_end);
    // NewB->set_marked_bottom_begin(splitpoint);
    /* NewB->set_end(splitpoint);                                            */ assert(NewB->to_constln.empty());
    /* NewB->set_inert_begin(?);                                             */ mCRL2complexity(NewB, add_work(check_complexity::
    /* NewB->set_inert_end(?);                                               */             Move_Blue_or_Red_to_a_new_block_NewB_pointer_3_29, max_counter), );
    for(permutation_iter_t s_iter=NewB->begin(); NewB->end()!=s_iter; ++s_iter)
    {                                                                           // mCRL2complexity(*s_iter, ...) -- optimized to the above call.
        (*s_iter)->block = NewB;
    }

    // adapt the old block: it only keeps the red states
    // set_end(end());
    set_bottom_begin(marked_bottom_begin());
    // set_marked_bottom_begin(marked_bottom_begin());
    set_marked_nonbottom_begin(marked_bottom_begin());
    set_begin(splitpoint);

    constln()->make_nontrivial();

    return NewB;
}

block_t* block_t::split_off_red(permutation_iter_t const red_nonbottom_begin)
{                                                                               assert(marked_nonbottom_begin() == red_nonbottom_begin);
                                                                                assert(marked_nonbottom_end() >= red_nonbottom_begin);
                                                                                assert(0 != marked_size());
    permutation_iter_t const splitpoint = red_nonbottom_begin +
                                                        unmarked_bottom_size(); assert(begin() < splitpoint);  assert(splitpoint < end());
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    unsigned const max_counter = check_complexity::log_n -
    /* It is not necessary to reset the nottoblue counters; these counters   */                     check_complexity::ilog2((state_type) (end() - splitpoint));
    /* are anyway only valid for the maybe-blue states.                      */     assert((state_type) (end() - splitpoint) <= size() / 2);
                                                                                #endif
    if (state_type swapcount = std::min(unmarked_bottom_size(),
                  (state_type) (marked_nonbottom_end() - red_nonbottom_begin));
                                                                0 != swapcount)
    {
        // vector swap the states:
        permutation_iter_t pos1 = red_nonbottom_begin;
        permutation_iter_t pos2 = unmarked_bottom_end();
        state_info_ptr const temp = *pos1;
        for (;;)
        {                                                                       mCRL2complexity(*pos1, add_work(check_complexity::
                                                                                               Move_Blue_or_Red_to_a_new_block_NewB_swap_3_29, max_counter), );
            *pos1 = *--pos2;
            (*pos1)->pos = pos1;
            ++pos1;
            if (0 == --swapcount)
            {
              break;
            }
            *pos2 = *pos1;
            (*pos2)->pos = pos2;
        }
        *pos2 = temp;
        (*pos2)->pos = pos2;
    }
    // create a new block for the red states
    block_t* const NewB = new block_t(constln(), splitpoint, end());
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    NewB->work_counter = work_counter;
                                                                                #endif
    if (BLOCK_NO_SEQNR != seqnr())
    {
        NewB->assign_seqnr();
    }
    // NewB->set_end(end());
    NewB->set_marked_bottom_begin(marked_bottom_begin());
    NewB->set_bottom_begin(marked_bottom_begin());
    NewB->set_marked_nonbottom_begin(marked_bottom_begin());
    /* NewB->set_begin(splitpoint);                                          */ assert(NewB->to_constln.empty());
    /* NewB->set_inert_begin(?);                                             */ mCRL2complexity(NewB, add_work(check_complexity::
    /* NewB->set_inert_end(?);                                               */             Move_Blue_or_Red_to_a_new_block_NewB_pointer_3_29, max_counter), );
    for(permutation_iter_t s_iter=NewB->begin(); NewB->end()!=s_iter; ++s_iter)
    {                                                                           // mCRL2complexity(*s_iter, ...) -- optimized to the above call.
        (*s_iter)->block = NewB;
    }

    // adapt the old block: it only keeps the blue states
    // set_begin(begin());
    set_marked_nonbottom_begin(red_nonbottom_begin);
    set_bottom_begin(red_nonbottom_begin);
    set_marked_bottom_begin(splitpoint);
    set_end(splitpoint);

    constln()->make_nontrivial();

    return NewB;
}
                                                                                #ifndef NDEBUG
                                                                                    /// \brief print a slice of the partition (typically a block)
                                                                                    /// \details If the slice indicated by the parameters is not empty, the
                                                                                    /// states in this slice will be printed.
                                                                                    /// \param message text printed as a title if the slice is not empty
                                                                                    /// \param B       block that is being printed (it is checked whether
                                                                                    ///                states belong to this block)
                                                                                    /// \param begin   iterator to the beginning of the slice
                                                                                    /// \param end     iterator past the end of the slice
                                                                                    void part_state_t::print_block(const char* const message,
                                                                                                                      const block_t* const B,
                                                                                                                      permutation_const_iter_t begin,
                                                                                                                      permutation_const_iter_t const end) const
                                                                                    {
                                                                                        if (0 != end - begin)
                                                                                        {
                                                                                            mCRL2log(log::debug) << "\t\t" << message
                                                                                                                           << (1 < end-begin ? "s:\n" : ":\n");
                                                                                            do
                                                                                            {
                                                                                                mCRL2log(log::debug) << "\t\t\t"
                                                                                                                                       << (*begin)->debug_id();
                                                                                                if (B != (*begin)->block)
                                                                                                {
                                                                                                    mCRL2log(log::debug) << ", inconsistent: "
                                                                                                                   "points to " << (*begin)->block->debug_id();
                                                                                                }
                                                                                                if (begin != (*begin)->pos)
                                                                                                {
                                                                                                    mCRL2log(log::debug)
                                                                                                               << ", inconsistent pointer to state_info_entry";
                                                                                                }
                                                                                                mCRL2log(log::debug) << '\n';
                                                                                            }
                                                                                            while (++begin != end);
                                                                                        }
                                                                                    }

                                                                                    /// \brief print the partition as a tree (per constellation and block)
                                                                                    /// \details The function prints all constellations (in order); for each
                                                                                    /// constellation it prints the blocks it consists of; and for each block,
                                                                                    /// it lists its states, separated into nonbottom and bottom states.
                                                                                    /// \param part_tr partition for the transitions
                                                                                    void part_state_t::print_part(const part_trans_t& part_tr) const
                                                                                    {
                                                                                        assert(permutation.begin() < permutation.end());
                                                                                        for (const constln_t* C = permutation.front()->constln(); ; )
                                                                                        {
                                                                                            mCRL2log(log::debug) << C->debug_id() << ":\n";
                                                                                            assert(C->begin() < C->end());
                                                                                            for (const block_t* B = (*C->begin())->block; ; )
                                                                                            {
                                                                                                mCRL2log(log::debug) << "\t" << B->debug_id();
                                                                                                if (C != B->constln())
                                                                                                {
                                                                                                    mCRL2log(log::debug) << ", inconsistent: "
                                                                                                                      "points to " << B->constln()->debug_id();
                                                                                                }
                                                                                                mCRL2log(log::debug) << ":\n";
                                                                                                print_block("Non-bottom state", B, B->nonbottom_begin(),
                                                                                                                                           B->nonbottom_end());
                                                                                                print_block("Bottom state", B, B->bottom_begin(),
                                                                                                                                              B->bottom_end());
                                                                                                mCRL2log(log::debug) << "\t\tThis block has ";
                                                                                                if (B->inert_end() == part_tr.B_to_C_begin())
                                                                                                {
                                                                                                    mCRL2log(log::debug)
                                                                                                               << "no transitions to its own constellation.\n";
                                                                                                    assert(B->inert_begin() == B->inert_end());
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                    assert(B->inert_end() ==
                                                                                                                         B->inert_end()[-1].B_to_C_slice->end);
                                                                                                    assert(B->inert_end()[-1].B_to_C_slice->from_block() == B);
                                                                                                    assert(B->inert_end()[-1].B_to_C_slice->to_constln() == C);
                                                                                                    mCRL2log(log::debug) << B->inert_end() -
                                                                                                                         B->inert_end()[-1].B_to_C_slice->begin
                                                                                                    <<" transition(s) to its own constellation,\n\t\tof which "
                                                                                                                     << B->inert_end() - B->inert_begin()
                                                                                                                     << (1 == B->inert_end() - B->inert_begin()
                                                                                                                           ? " is inert.\n" : " are inert.\n");
                                                                                                }
                                                                                                // go to next block
                                                                                                if (C->end()
                                                                                                    == B->end())
                                                                                                {
                                                                                                  break;
                                                                                                }
                                                                                                B = (*B->end())->block;
                                                                                            }
                                                                                            // go to next constellation
                                                                                            if (permutation.end()
                                                                                                == C->end())
                                                                                            {
                                                                                              break;
                                                                                            }
                                                                                            C = (*C->end())->constln();
                                                                                        }
                                                                                    }

                                                                                    /// \brief print all transitions
                                                                                    /// \details For each state (in order), its outgoing transitions are
                                                                                    /// listed, sorted by goal constellation.  The function also indicates
                                                                                    /// where the current constellation pointer of the state points at.
                                                                                    void part_state_t::print_trans() const
                                                                                    {
                                                                                        fixed_vector<state_info_entry>::const_iterator const state_info_end =
                                                                                                                                          state_info.end() - 1;
                                                                                        for (fixed_vector<state_info_entry>::const_iterator state_iter =
                                                                                                state_info.begin(); state_info_end != state_iter; ++state_iter)
                                                                                        {
                                                                                            // print transitions out of state
                                                                                            succ_const_iter_t succ_constln_iter = state_iter->succ_begin();
                                                                                            if (state_iter->succ_end() != succ_constln_iter)
                                                                                            {
                                                                                                mCRL2log(log::debug) << state_iter->debug_id()
                                                                                                                                                      << ":\n";
                                                                                                assert(state_iter->succ_begin() <=
                                                                                                                               state_iter->inert_succ_begin());
                                                                                                assert(state_iter->inert_succ_begin() <=
                                                                                                                                 state_iter->inert_succ_end());
                                                                                                assert(state_iter->succ_begin()==state_iter->inert_succ_begin()
                                                                                                    || *state_iter->inert_succ_begin()[-1].target->constln() <=
                                                                                                                                       *state_iter->constln());
                                                                                                assert(state_iter->succ_begin()==state_iter->inert_succ_begin()
                                                                                                    || state_iter->inert_succ_begin()[-1].target->block !=
                                                                                                                                            state_iter->block);
                                                                                                assert(state_iter->inert_succ_begin() ==
                                                                                                                                state_iter->inert_succ_end() ||
                                                                                                          (state_iter->inert_succ_begin()->target->block ==
                                                                                                                                           state_iter->block &&
                                                                                                           state_iter->inert_succ_end()[-1].target->block ==
                                                                                                                                           state_iter->block));
                                                                                                assert(state_iter->succ_end() == state_iter->inert_succ_end()||
                                                                                                             *state_iter->constln() <
                                                                                                             *state_iter->inert_succ_end()->target->constln());
                                                                                                do
                                                                                                {
                                                                                                    // print transitions to a constellation
                                                                                                    mCRL2log(log::debug) << "\ttransitions to "
                                                                                                            << succ_constln_iter->target->constln()->debug_id()
                                                                                                                                                      << ":\n";
                                                                                                    succ_const_iter_t s_iter = succ_constln_iter;
                                                                                                    // set succ_constln_iter to the end of the transitions to
                                                                                                    // this constellation
                                                                                                    succ_constln_iter=succ_entry::slice_end(succ_constln_iter);
                                                                                                    for ( ;s_iter != succ_constln_iter ;++s_iter)
                                                                                                    {
                                                                                                        mCRL2log(log::debug) << "\t\tto "
                                                                                                                                 << s_iter->target->debug_id();
                                                                                                        if (state_iter->inert_succ_begin() <= s_iter &&
                                                                                                                         s_iter < state_iter->inert_succ_end())
                                                                                                        {
                                                                                                            mCRL2log(log::debug) << " (inert)";
                                                                                                        }
                                                                                                        if (state_iter->current_constln() == s_iter)
                                                                                                        {
                                                                                                            mCRL2log(log::debug)
                                                                                                                                      << " <- current_constln";
                                                                                                        }
                                                                                                        mCRL2log(log::debug) << '\n';
                                                                                                        assert(s_iter->B_to_C->pred->succ == s_iter);
                                                                                                        assert(s_iter->B_to_C->pred->source == &*state_iter);
                                                                                                    }
                                                                                                }
/* ************************************************************************* */                 while (state_iter->succ_end() != succ_constln_iter);
/*                                                                           */                 if (state_iter->current_constln() == state_iter->succ_end())
/*                           T R A N S I T I O N S                           */                 {
/*                                                                           */                     mCRL2log(log::debug)
/* ************************************************************************* */                                                  << "\t\t<- current_constln\n";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                #endif // ifndef NDEBUG
/// \brief handle the transitions from the splitter to its own constellation
/// \details split_inert_to_C splits the B_to_C slice of block SpB to its own
/// constellation into two slices: one for the inert and one for the non-inert
/// transitions.  It is called with SpB just after a constellation is split, as
/// the transitions from SpB to itself (= the inert transitions) now go to a
/// different constellation than the other transitions from SpB to its old
/// constellation.  It does, however, not adapt the other transition arrays to
/// reflect that noninert and inert transitions from block SpB would go to
/// different constellations.
///
/// Its time complexity is O(1+min{|out_noninert(SpB-->C)|, |inert_out(SpB)|}).
/// \param SpB  pointer to the splitter block
void part_trans_t::split_inert_to_C(block_t* const SpB)
{
    // if there are no inert transitions
    if (SpB->inert_begin() == SpB->inert_end())
    {
        if (SpB->inert_end() != B_to_C.begin())
        {
            // There are noninert transitions from SpB to its old
            // constellation:  they all go to SpC.
            SpB->SetFromRed(SpB->inert_end()[-1].B_to_C_slice);
            // There are no more transitions from SpB to its own constellation
            SpB->set_inert_begin(B_to_C.begin());
            SpB->set_inert_end(B_to_C.begin());
        }
        return;
    }
    B_to_C_desc_iter_t const slice = SpB->inert_begin()->B_to_C_slice;
    // if all transitions are inert
    if (slice->begin == SpB->inert_begin())
    {
        return;
    }

    // now the slice actually has to be split
    B_to_C_desc_iter_t new_slice;
    // select the smaller number of swaps to decide which part should be the
    // new one:
    if(SpB->inert_begin() - slice->begin < slice->end - SpB->inert_begin())
    {
        // fewer noninert transitions
        SpB->to_constln.emplace_front(slice->begin, SpB->inert_begin());
        new_slice = SpB->to_constln.begin();
        /* SpB->SetFromRed(new_slice);                                       */ assert(new_slice->from_block() == SpB);
        slice->begin = SpB->inert_begin();
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    new_slice->work_counter = slice->work_counter;
                                                                                    // We actually change the pointers in the new slice (i. e. the pointers of
                                                                                    // the non-inert transitions) because there are fewer of them; however, we
                                                                                    // still have to assign this work to the inert transitions.
                                                                                    mCRL2complexity(slice, add_work(check_complexity::
                                                                                              Register_that_inert_transitions_from_s_go_to_NewC_B_to_C_2_17,
                                                                                              check_complexity::log_n-check_complexity::ilog2(SpB->size())), );
                                                                                #endif
    }
    else
    {
        // fewer (or equal) inert transitions
        SpB->to_constln.emplace_back(SpB->inert_begin(), slice->end);
        new_slice = std::prev(SpB->to_constln.end());
        slice->end = SpB->inert_begin();
        SpB->SetFromRed(slice);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    new_slice->work_counter = slice->work_counter;
                                                                                    mCRL2complexity(new_slice, add_work(check_complexity::
                                                                                              Register_that_inert_transitions_from_s_go_to_NewC_B_to_C_2_17,
                                                                                              check_complexity::log_n-check_complexity::ilog2(SpB->size())), );
                                                                                #endif
    }
    // set the slice pointers of the smaller part to the new slice:
    for (B_to_C_iter_t iter = new_slice->begin; new_slice->end != iter; ++iter)
    {                                                                           // mCRL2complexity(...) -- optimized to the above calls.
                                                                                assert(B_to_C.end() > iter);  assert(iter->pred->succ->B_to_C == iter);
        iter->B_to_C_slice = new_slice;
    }
}


/// \brief transition target is moved to a new constellation
/// \details part_trans_t::change_to_C has to be called after a transition
/// target has changed its constellation.  The member function will adapt the
/// transition data structure.  It assumes that the transition is non-inert and
/// that the new constellation does not (yet) have _inert_ incoming
/// transitions.  It returns the boundary between transitions to SpC and
/// transitions to NewC in the state's outgoing transition array.
///
/// In debug mode, two additional parameters `SpC` (splitter constellation) and
/// `NewC` (new constellation) are available for assertions.
/// \param pred_iter  transition that has to be changed
/// \param first_transition_of_state  This is the first transition of the
///                                   state, so a new constln slice is started.
/// \param first_transition_of_block  This is the first transition of the
///                                   block, so a new B_to_C slice has to be
///                                   allocated.
succ_iter_t part_trans_t::change_to_C(pred_iter_t const pred_iter,              ONLY_IF_DEBUG( constln_t* const SpC, constln_t* const NewC, )
    bool const first_transition_of_state, bool const first_transition_of_block)
{                                                                               assert(pred_iter<pred.end()); assert(pred_iter->succ->B_to_C->pred==pred_iter);
    // adapt the B_to_C array:
    // always move the transition to the beginning of the slice (this will make
    // it easier because inert transitions are stored at the end of a slice).
    B_to_C_iter_t const old_B_to_C_pos = pred_iter->succ->B_to_C;
    B_to_C_desc_iter_t const old_B_to_C_slice = old_B_to_C_pos->B_to_C_slice;
    B_to_C_iter_t const new_B_to_C_pos = old_B_to_C_slice->begin;               assert(new_B_to_C_pos < B_to_C.end());
    B_to_C_desc_iter_t new_B_to_C_slice;                                        assert(new_B_to_C_pos->pred->succ->B_to_C == new_B_to_C_pos);
    if (first_transition_of_block)
    {
        // create a new slice in B_to_C for the transitions from RfnB to NewC
        block_t* const RfnB = pred_iter->source->block;
        RfnB->to_constln.emplace_back(new_B_to_C_pos, new_B_to_C_pos);
        new_B_to_C_slice = std::prev(RfnB->to_constln.end());
                                                                                #ifndef NDEBUG
                                                                                    new_B_to_C_slice->work_counter = old_B_to_C_slice->work_counter;
                                                                                #endif
        RfnB->SetFromRed(old_B_to_C_slice);
    }
    else
    {
        new_B_to_C_slice = new_B_to_C_pos[-1].B_to_C_slice;
    }
    ++new_B_to_C_slice->end;
    ++old_B_to_C_slice->begin;                                                  assert(new_B_to_C_slice->end == old_B_to_C_slice->begin);
    swap_B_to_C(pred_iter->succ, new_B_to_C_pos->pred->succ);
    new_B_to_C_pos->B_to_C_slice = new_B_to_C_slice;
    if (old_B_to_C_slice->begin == old_B_to_C_slice->end)
    {
        // this was the last transition from RfnB to SpC
        block_t* const RfnB = pred_iter->source->block;
        if (RfnB->inert_end() == old_B_to_C_slice->begin)
        {                                                                       assert(RfnB->inert_begin() == RfnB->inert_end());
            // this was the last transition from RfnB to its own constellation
            RfnB->set_inert_begin(B_to_C.begin());
            RfnB->set_inert_end(B_to_C.begin());
        }                                                                       assert(RfnB->to_constln.begin() == old_B_to_C_slice);
        RfnB->to_constln.erase(old_B_to_C_slice);
    }
    // adapt the outgoing transition array:
    // move the transition to the beginning
    succ_iter_t const old_out_pos = pred_iter->succ;                            assert(succ.end() > old_out_pos);
                                                                                assert(old_out_pos->B_to_C->pred->succ == old_out_pos);
    /* move to beginning                                                     */ assert(*NewC < *SpC);
    succ_iter_t const new_out_pos = old_out_pos->slice_begin();
    swap_out(pred_iter, new_out_pos->B_to_C->pred);
    if (succ_iter_t const old_slice_end = succ_entry::slice_end(old_out_pos);
                                               old_slice_end > new_out_pos + 1)
    {
        // the following assignment might assign an illegal transition if the
        // old slice becomes empty -- but then it doesn't hurt, because it will
        // be overwritten below.  However, Visual Studio C complains.
        old_slice_end[-1].set_slice_begin_or_before_end(new_out_pos + 1);       assert(old_slice_end - 1 == new_out_pos->slice_begin_or_before_end());
    }                                                                           else
                                                                                {   assert(old_slice_end - 1 == new_out_pos);
                                                                                    assert(old_slice_end-1==old_slice_end[-1].slice_begin_or_before_end());
                                                                                }
    if (first_transition_of_state)
    {
        new_out_pos->set_slice_begin_or_before_end(new_out_pos);
    }
    else
    {                                                                           assert(pred_iter->source->current_constln() == new_out_pos);
        new_out_pos->set_slice_begin_or_before_end(
                                  new_out_pos[-1].slice_begin_or_before_end());
    }
    return new_out_pos + 1;
}


/// \brief Split outgoing transitions of a state in the splitter
/// \details split_s_inert_out splits the outgoing transitions from `s` to its
/// own constellation into two:  the inert transitions become transitions to
/// the new constellation of which `s` is now part;  the non-inert transitions
/// remain transitions to `OldC` (the old constellation, only available in debug
/// mode).
/// Its time complexity is O(1 + min { |out_notau(s)|, |out_tau(s)| }).
/// \result true iff the state also has transitions to OldC
bool part_trans_t::split_s_inert_out(state_info_ptr s                           ONLY_IF_DEBUG( , constln_t* OldC )
                                                     )
{
    constln_t* NewC = s->constln();                                             assert(*NewC < *OldC); assert(NewC->sort_key + OldC->size() == OldC->sort_key);
                                                                                assert(OldC->end() == NewC->begin() || NewC->end() == OldC->begin());
                                                                                succ_iter_t split
                                                                                    = s->inert_succ_begin();
                                                                                succ_iter_t to_C_end
                                                                                    = s->inert_succ_end();
                                                                                assert(
                                                                                    s->succ_begin() == to_C_end
                                                                                    || to_C_end[-1]
                                                                                               .slice_begin_or_before_end()
                                                                                           < to_C_end);
                                                                                succ_iter_t to_C_begin
                                                                                    = s->succ_begin() == to_C_end
                                                                                          ? to_C_end
                                                                                          : to_C_end[-1]
                                                                                                .slice_begin_or_before_end();
                                                                                //< If s has no transitions to OldC at
                                                                                //all, then to_C_begin may be the
                                                                                // beginning of the constln_slice for
                                                                                // transitions to another constellation.
                                                                                // We will check that later.
                                                                                bool result = to_C_begin < split;                                           assert(to_C_begin <= split);  assert(split <= to_C_end);
                                                                                assert(succ.end() == split || split->B_to_C->pred->succ == split);
                                                                                assert(succ.end() == to_C_end || to_C_end->B_to_C->pred->succ == to_C_end);
                                                                                assert(succ.end() == to_C_begin ||
                                                                                                             to_C_begin->B_to_C->pred->succ == to_C_begin);
                                                                                if (!result)
                                                                                {
                                                                                  ;
                                                                                }
                                                                                else if (split < to_C_end)
                                                                                {
                                                                                  // s has both inert and non-inert
                                                                                  // transitions
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
        /* the out-transitions of s also have to be swapped.                 */     unsigned const max_counter = check_complexity::log_n -
        /* Actually only B_to_C and the target need to be swapped, as the    */                                          check_complexity::ilog2(NewC->size());
        /* constln_slices are (still) identical.                             */     assert(*NewC < *OldC);
                                                                                #endif
        trans_type swapcount = std::min(to_C_end - split, split - to_C_begin);
        split = to_C_end - split + to_C_begin;                                  assert(0 != swapcount);

        succ_iter_t pos1 = to_C_begin;
        succ_iter_t pos2 = to_C_end;
        state_info_ptr temp_target = pos1->target;
        B_to_C_iter_t temp_B_to_C = pos1->B_to_C;
        for (;;)
        {
            --pos2;                                                             mCRL2complexity(pos2->B_to_C->pred, add_work(check_complexity::
                                                                                               Register_that_inert_transitions_from_s_go_to_NewC_swap_2_17,
                                                                                                                                               max_counter), );
            pos1->target = pos2->target;
            pos1->B_to_C = pos2->B_to_C;
            pos1->B_to_C->pred->succ = pos1;
            ++pos1;
            if (0 == --swapcount)
            {
              break;
            }
            pos2->target = pos1->target;
            pos2->B_to_C = pos1->B_to_C;
            pos2->B_to_C->pred->succ = pos2;
        }
        pos2->target = temp_target;
        pos2->B_to_C = temp_B_to_C;
        pos2->B_to_C->pred->succ = pos2;

        s->set_inert_succ_begin_and_end(to_C_begin, split);                     assert(s->succ_begin() == s->inert_succ_begin() ||
                                                                                                         *s->inert_succ_begin()[-1].target->constln() < *NewC);
                                                                                assert(s->inert_succ_begin()->target->block == s->block);
                                                                                assert(s->inert_succ_end()[-1].target->block == s->block);
                                                                                assert(s->inert_succ_end() < s->succ_end());
        /* set the pointer to the slice for the inert transitions.           */ assert(OldC == s->inert_succ_end()->target->constln());
        to_C_end[-1].set_slice_begin_or_before_end(split);
        split[-1].set_slice_begin_or_before_end(to_C_begin);
        for (succ_iter_t succ_iter = split - 1; to_C_begin != succ_iter; )
        {
            --succ_iter;                                                        mCRL2complexity(succ_iter->B_to_C->pred,add_work(check_complexity::
                                                                                               Register_that_inert_transitions_from_s_go_to_NewC_succ_2_17,
                                                                                                                                               max_counter), );
                                                                                assert(succ.end() > succ_iter);
                                                                                assert(succ_iter->B_to_C->pred->succ == succ_iter);
            succ_iter->set_slice_begin_or_before_end(split - 1);
        }
    }
    else if (*to_C_begin->target->constln() > *NewC)
    {
        // s has (noninert) transitions to OldC, but no (inert)
        /* transitions to NewC.                                              */ assert(to_C_begin->target->constln() == OldC);
        s->set_inert_succ_begin_and_end(to_C_begin, to_C_begin);                assert(s->succ_begin() == s->inert_succ_begin() ||
                                                                                                         *s->inert_succ_begin()[-1].target->constln() < *NewC);
                                                                                assert(s->succ_end() > s->inert_succ_end());  assert(to_C_end == split);
                                                                                assert(OldC == s->inert_succ_end()->target->constln());
    }
    else
    {
        // s actually hasn't got transitions to OldC at all, so `result`
        // should not be true.
        result = false;
    }
                                                                                #ifndef NDEBUG
                                                                                    if (succ_iter_t succ_iter = s->succ_begin(); succ_iter != s->succ_end())
                                                                                    {
                                                                                        do
                                                                                        {
                                                                                            assert(succ_iter->B_to_C->pred->source == s);
                                                                                        }
                                                                                        while (++succ_iter != s->succ_end());
                                                                                        assert(s->succ_begin() == s->inert_succ_begin() ||
                                                                                                         *s->inert_succ_begin()[-1].target->constln() < *NewC);
                                                                                        assert(s->inert_succ_begin() == s->inert_succ_end() ||
                                                                                                          (s->inert_succ_begin()->target->block == s->block &&
                                                                                                           s->inert_succ_end()[-1].target->block == s->block));
                                                                                        assert(s->inert_succ_end() == s->succ_end() ||
                                                                                                              *NewC < *s->inert_succ_end()->target->constln());
                                                                                    }
                                                                                #endif // ifndef NDEBUG
    return result;
}


/// \brief handle B_to_C slices after a new blue block has been created
/// \details part_trans_t::new_blue_block_created splits the B_to_C-slices to
/// reflect that some transitions now start in the new block NewB.  They can no
/// longer be in the same slice as the transitions that start in the old block.
/// Its time complexity is O(1 + |out(NewB)|).
/// \param RfnB  the old block that has been refined (now the red subblock)
/// \param NewB  the new block (the blue subblock)
void part_trans_t::new_blue_block_created(block_t* const RfnB,
                                                           block_t* const NewB)
{                                                                               assert(RfnB->constln()==NewB->constln());  assert(NewB->end()==RfnB->begin());
    NewB->set_inert_begin_and_end(B_to_C.begin(), B_to_C.begin());
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    unsigned const max_counter = check_complexity::log_n -
                                                                                                                         check_complexity::ilog2(NewB->size());
                                                                                    mCRL2complexity(NewB, add_work(check_complexity::
    /* for all outgoing transitions of NewB                                  */                   Move_Blue_or_Red_to_a_new_block_states_3_29, max_counter), );
                                                                                #endif
    for(permutation_iter_t s_iter=NewB->begin(); NewB->end()!=s_iter; ++s_iter)
    {                                                                           // mCRL2complexity(*s_iter, ...) -- optimized to the above call
        for (succ_iter_t succ_iter = (*s_iter)->succ_begin();
                               (*s_iter)->succ_end() != succ_iter; ++succ_iter)
        {                                                                       mCRL2complexity(succ_iter->B_to_C->pred,add_work(check_complexity::
                                                                                                    Move_Blue_or_Red_to_a_new_block_succ_3_29, max_counter), );
            /* Move the transition to a new slice:                           */ assert(succ.end()>succ_iter); assert(succ_iter->B_to_C->pred->succ==succ_iter);
            B_to_C_iter_t const old_pos = succ_iter->B_to_C;
            B_to_C_desc_iter_t const old_B_to_C_slice = old_pos->B_to_C_slice;
            B_to_C_iter_t const after_new_pos = old_B_to_C_slice->end;          assert(B_to_C.end() > old_pos);  assert(after_new_pos > old_pos);
                                                                                assert(old_pos->pred->succ->B_to_C == old_pos);
                                                                                assert(B_to_C.end() == after_new_pos ||
                                                                                                           after_new_pos->pred->succ->B_to_C == after_new_pos);
            B_to_C_desc_iter_t new_B_to_C_slice;
            if (B_to_C.end() == after_new_pos ||
                    after_new_pos->pred->source->block != NewB ||
                        after_new_pos->pred->succ->target->constln() !=
                                                  succ_iter->target->constln())
            {
                // create a new B_to_C-slice
                // this can only happen when the first transition from
                // *s_iter to a new constellation is handled.
                if (!old_B_to_C_slice->needs_postprocessing())
                {
                    // (During primary refinement [Line 2.23]: There is no need
                    // to put the FromRed-slice of the new blue block into
                    // position like in new_red_block_created(), as only the
                    // red subblock will be refined further in Line 2.24.)
                    // During postprocessing:
                    // the old_B_to_C_slice does not need postprocessing, so
                    // its corresponding new slice should be in a similar
                    // position near the beginning of the list of slices.
                    NewB->to_constln.emplace_front(after_new_pos,
                                                                after_new_pos); // new_B_to_C_slice is not yet fully initialised, therefore
                    new_B_to_C_slice = NewB->to_constln.begin();                // the assertion fails:
                                                                                // assert(new_B_to_C_slice->from_block() == NewB);
                }
                else
                {
                    // During postprocessing:
                    // the old_B_to_C_slice needs postprocessing, so also the
                    // new_B_to_C_slice will need postprocessing.
                    NewB->to_constln.emplace_back(after_new_pos,after_new_pos);
                    new_B_to_C_slice = std::prev(NewB->to_constln.end());
                }
                                                                                #ifndef NDEBUG
                                                                                    new_B_to_C_slice->work_counter = old_B_to_C_slice->work_counter;
                                                                                #endif
                if (RfnB->inert_end() == after_new_pos)
                {
                    /* this is the first transition from NewB to its own     */ assert(NewB->inert_begin() == B_to_C.begin());
                    /* constellation.  Adapt the pointers accordingly.       */ assert(NewB->inert_end() == B_to_C.begin());
                    NewB->set_inert_end(after_new_pos);
                    NewB->set_inert_begin(after_new_pos);
                }
            }
            else
            {
                // the slice at after_new_pos is already the correct one
                new_B_to_C_slice = after_new_pos->B_to_C_slice;
            }
            --new_B_to_C_slice->begin;
            --old_B_to_C_slice->end;                                            assert(new_B_to_C_slice->begin == old_B_to_C_slice->end);
            B_to_C_iter_t new_pos = after_new_pos - 1;                          assert(B_to_C.end() > new_pos);
                                                                                assert(new_pos->pred->succ->B_to_C == new_pos);
            if (RfnB->inert_end() == after_new_pos)
            {
                // The transition goes from NewB to the constellation of
                // RfnB and NewB.
                if (RfnB->inert_begin() <= old_pos)
                {
                    // The transition is inert and has to be moved over the
                    // non-inert transitions of NewB.
                    NewB->set_inert_begin(NewB->inert_begin() - 1);
                    new_pos = NewB->inert_begin();
                    // old_pos --> new_pos --> new_B_to_C_slice->begin -->
                    // old_pos
                    swap3_B_to_C(succ_iter, new_pos->pred->succ,
                                    new_B_to_C_slice->begin->pred->succ);
                }
                else
                {
                    // The transition is non-inert, but it has to be moved
                    // over the inert transitions of RfnB.
                    RfnB->set_inert_begin(RfnB->inert_begin() - 1);
                    // old_pos --> new_pos --> RfnB->inert_begin() -> old_pos
                    swap3_B_to_C(succ_iter, new_pos->pred->succ,
                                    RfnB->inert_begin()->pred->succ);
                }
                RfnB->set_inert_end(RfnB->inert_end() - 1);
                if (old_B_to_C_slice->begin == old_B_to_C_slice->end)
                {
                    // This was the last transition from RfnB to its own
                    // constellation.
                    RfnB->set_inert_begin(B_to_C.begin());
                    RfnB->set_inert_end(B_to_C.begin());

                    RfnB->to_constln.erase(old_B_to_C_slice);
                }
            }
            else
            {
                // The transition goes from NewB to a constellation that
                // does not contain RfnB or NewB.  No special treatment is
                // required.
                swap_B_to_C(succ_iter, new_pos->pred->succ);
                if (old_B_to_C_slice->begin == old_B_to_C_slice->end)
                {
                    RfnB->to_constln.erase(old_B_to_C_slice);
                }
            }
            new_pos->B_to_C_slice = new_B_to_C_slice;
        }
    }
                                                                                #ifndef NDEBUG
                                                                                    if (RfnB->inert_begin() == RfnB->inert_end() &&
                                                                                            RfnB->inert_end() != B_to_C.begin() &&
                                                                                                (RfnB->inert_end()[-1].pred->source->block != RfnB ||
                                                                                                 RfnB->inert_end()[-1].pred->succ->target->constln()
                                                                                                                                           != RfnB->constln()))
                                                                                    {
                                                                                        assert(0 && "The old block has no transitions to its own "
                                                                                                           "constellation, but its inert_begin and "
                                                                                                           "inert_end pointers are not set to B_to_C.begin()");
                                                                                    }
                                                                                    if (RfnB->inert_end() != B_to_C.begin())
                                                                                    {
                                                                                        assert(RfnB->inert_end()[-1].pred->source->block == RfnB);
                                                                                        assert(RfnB->inert_end()[-1].pred->succ->target->constln() ==
                                                                                                                                              RfnB->constln());
                                                                                    }
                                                                                    if (NewB->inert_begin() == NewB->inert_end() &&
                                                                                            NewB->inert_end() != B_to_C.begin() &&
                                                                                                (NewB->inert_end()[-1].pred->source->block != NewB ||
                                                                                                 NewB->inert_end()[-1].pred->succ->target->constln()
                                                                                                                                           != NewB->constln()))
                                                                                    {
                                                                                        assert(0 && "The new block has no transitions to its own "
                                                                                                           "constellation, but its inert_begin and "
                                                                                                           "inert_end pointers are not set to B_to_C.begin()");
                                                                                    }
                                                                                    if (NewB->inert_end() != B_to_C.begin())
                                                                                    {
                                                                                        assert(NewB->inert_end()[-1].pred->source->block == NewB);
                                                                                        assert(NewB->inert_end()[-1].pred->succ->target->constln() ==
                                                                                                                                              NewB->constln());
                                                                                    }
                                                                                #endif
}


/// \brief handle B_to_C slices after a new red block has been created
/// \details part_trans_t::new_red_block_created splits the B_to_C-slices to
/// reflect that some transitions now start in the new block NewB.  They can no
/// longer be in the same slice as the transitions that start in the old block.
/// Its time complexity is O(1 + |out(NewB)|).
/// \param RfnB  the old block that has been refined (now the blue subblock)
/// \param NewB  the new block (the red subblock)
/// \param postprocessing  true iff the refinement happened during
///                        postprocessing.  (Otherwise, the refinement should
///                        preserve the information about `FromRed()`).
void part_trans_t::new_red_block_created(block_t* const RfnB,
                                block_t* const NewB, bool const postprocessing)
{                                                                               assert(RfnB->constln()==NewB->constln());  assert(NewB->begin()==RfnB->end());
    NewB->set_inert_begin_and_end(B_to_C.begin(), B_to_C.begin());
    bool old_fromred_invalid = false;
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    unsigned const max_counter = check_complexity::log_n -
                                                                                                                         check_complexity::ilog2(NewB->size());
                                                                                    mCRL2complexity(NewB, add_work(check_complexity::
    /* for all outgoing transitions of NewB                                  */                   Move_Blue_or_Red_to_a_new_block_states_3_29, max_counter), );
                                                                                #endif
    for(permutation_iter_t s_iter=NewB->begin(); NewB->end()!=s_iter; ++s_iter)
    {                                                                           // mCRL2complexity(*s_iter, ...) -- optimized to the above call
        for (succ_iter_t succ_iter = (*s_iter)->succ_begin();
                               (*s_iter)->succ_end() != succ_iter; ++succ_iter)
        {                                                                       mCRL2complexity(succ_iter->B_to_C->pred,add_work(check_complexity::
                                                                                                    Move_Blue_or_Red_to_a_new_block_succ_3_29, max_counter), );
            /* Move the transition to a new slice:                           */ assert(succ.end()>succ_iter); assert(succ_iter->B_to_C->pred->succ==succ_iter);
            B_to_C_iter_t const old_pos = succ_iter->B_to_C;
            B_to_C_desc_iter_t const old_B_to_C_slice = old_pos->B_to_C_slice;
            B_to_C_iter_t new_pos = old_B_to_C_slice->begin;                    assert(B_to_C.end() > old_pos); assert(old_pos->pred->succ->B_to_C == old_pos);
                                                                                assert(B_to_C.end() > new_pos); assert(new_pos->pred->succ->B_to_C == new_pos);
            B_to_C_desc_iter_t new_B_to_C_slice;                                assert(new_pos <= old_pos);
            if (B_to_C.begin() == new_pos ||
                    new_pos[-1].pred->source->block != NewB ||
                        new_pos[-1].pred->succ->target->constln() !=
                                            succ_iter->target->constln())
            {
                // create a new B_to_C-slice
                // this can only happen when the first transition from
                // *s_iter to a new constellation is handled.
                if (!postprocessing ? !old_fromred_invalid &&
                                  RfnB->to_constln.begin() == old_B_to_C_slice
                                   : !old_B_to_C_slice->needs_postprocessing())
                {
                    // During primary refinement (Line 2.23):
                    // The old B_to_C_slice is in the FromRed-position, i. e.
                    // it contains the transitions to SpC.  So the new slice
                    // also contains the transitions to SpC.
                    // During postprocessing:
                    // The old_B_to_C_slice does not need postprocessing, so
                    // its corresponding new slice should be in a similar
                    // position near the beginning of the list of slices.
                    NewB->to_constln.emplace_front(new_pos, new_pos);           // new_B_to_C_slice is not yet fully initialised, therefore
                    new_B_to_C_slice = NewB->to_constln.begin();                // the assertion fails:
                                                                                // assert(new_B_to_C_slice->from_block() == NewB);
                }
                else
                {
                    // During primary refinement (Line 2.23):
                    // The old_B_to_C_slice is not in the FromRed-position. The
                    // corresponding slice of NewB should be moved into a
                    // position that does not change a potential FromRed-slice
                    // there.
                    // During postprocessing:
                    // The old_B_to_C_slice needs postprocessing, so also the
                    // new_B_to_C_slice will need postprocessing.
                    NewB->to_constln.emplace_back(new_pos, new_pos);
                    new_B_to_C_slice = std::prev(NewB->to_constln.end());
                }
                                                                                #ifndef NDEBUG
                                                                                    new_B_to_C_slice->work_counter = old_B_to_C_slice->work_counter;
                                                                                #endif
                if (RfnB->inert_end() == old_B_to_C_slice->end)
                {
                    /* this is the first transition from NewB to its own     */ assert(NewB->inert_begin() == B_to_C.begin());
                    /* constellation.  Adapt the pointers accordingly.       */ assert(NewB->inert_end() == B_to_C.begin());
                    NewB->set_inert_end(new_pos);
                    NewB->set_inert_begin(new_pos);
                }
            }
            else
            {
                // the slice before new_pos is already the correct one
                new_B_to_C_slice = new_pos[-1].B_to_C_slice;
            }
            ++new_B_to_C_slice->end;
            ++old_B_to_C_slice->begin;                                          assert(new_B_to_C_slice->end == old_B_to_C_slice->begin);
            if (RfnB->inert_end() == old_B_to_C_slice->end)
            {
                // The transition goes from NewB to the constellation of
                // RfnB and NewB.
                NewB->set_inert_end(NewB->inert_end() + 1);
                if (RfnB->inert_begin() <= old_pos)
                {
                    // The transition is inert and has to be moved over the
                    // non-inert transitions of RfnB.
                    // old_pos --> new_pos --> RfnB->inert_begin() --> old_pos
                    swap3_B_to_C(succ_iter, new_pos->pred->succ,
                                              RfnB->inert_begin()->pred->succ);
                    RfnB->set_inert_begin(RfnB->inert_begin() + 1);
                }
                else
                {
                    // The transition is non-inert, but it has to be moved
                    // over the inert transitions of NewB.
                    new_pos = NewB->inert_begin();
                    NewB->set_inert_begin(NewB->inert_begin() + 1);
                    // old_pos --> new_pos --> NewB->inert_end() - 1 -> old_pos
                    swap3_B_to_C(succ_iter, new_pos->pred->succ,
                                             NewB->inert_end()[-1].pred->succ);
                }
                if (old_B_to_C_slice->begin == old_B_to_C_slice->end)
                {
                    // This was the last transition from RfnB to its own
                    // constellation.
                    RfnB->set_inert_begin(B_to_C.begin());
                    RfnB->set_inert_end(B_to_C.begin());

                    if (!old_fromred_invalid &&
                                  RfnB->to_constln.begin() == old_B_to_C_slice)
                    {
                        old_fromred_invalid = true;
                    }
                    RfnB->to_constln.erase(old_B_to_C_slice);
                }
            }
            else
            {
                // The transition goes from NewB to a constellation that
                // does not contain RfnB or NewB.  No special treatment is
                // required.
                swap_B_to_C(succ_iter, new_pos->pred->succ);
                if (old_B_to_C_slice->begin == old_B_to_C_slice->end)
                {
                    if (!old_fromred_invalid &&
                                  RfnB->to_constln.begin() == old_B_to_C_slice)
                    {
                        old_fromred_invalid = true;
                    }
                    RfnB->to_constln.erase(old_B_to_C_slice);
                }
            }
            new_pos->B_to_C_slice = new_B_to_C_slice;
        }
    }
                                                                                #ifndef NDEBUG
                                                                                    if (RfnB->inert_begin() == RfnB->inert_end() &&
                                                                                            RfnB->inert_end() != B_to_C.begin() &&
                                                                                                (RfnB->inert_end()[-1].pred->source->block != RfnB ||
                                                                                                 RfnB->inert_end()[-1].pred->succ->target->constln()
                                                                                                                                           != RfnB->constln()))
                                                                                    {
                                                                                        assert(0 && "The old block has no transitions to its own "
                                                                                                           "constellation, but its inert_begin and "
                                                                                                           "inert_end pointers are not set to B_to_C.begin()");
                                                                                    }
                                                                                    if (RfnB->inert_end() != B_to_C.begin())
                                                                                    {
                                                                                        assert(RfnB->inert_end()[-1].pred->source->block == RfnB);
                                                                                        assert(RfnB->inert_end()[-1].pred->succ->target->constln() ==
                                                                                                                                              RfnB->constln());
                                                                                    }
                                                                                    if (NewB->inert_begin() == NewB->inert_end() &&
                                                                                            NewB->inert_end() != B_to_C.begin() &&
                                                                                                (NewB->inert_end()[-1].pred->source->block != NewB ||
                                                                                                 NewB->inert_end()[-1].pred->succ->target->constln()
                                                                                                                                           != NewB->constln()))
                                                                                    {
                                                                                        assert(0 && "The new block has no transitions to its own "
                                                                                                           "constellation, but its inert_begin and "
                                                                                                           "inert_end pointers are not set to B_to_C.begin()");
                                                                                    }
                                                                                    if (NewB->inert_end() != B_to_C.begin())
                                                                                    {
                                                                                        assert(NewB->inert_end()[-1].pred->source->block == NewB);
                                                                                        assert(NewB->inert_end()[-1].pred->succ->target->constln() ==
                                                                                                                                              NewB->constln());
                                                                                    }
                                                                                #endif
}
                                                                                #ifndef NDEBUG
                                                                                    /// \brief assert that the data structure is consistent and stable
                                                                                    /// \details The data structure is tested against a large number of
                                                                                    /// assertions to ensure that everything is consistent, e. g. pointers that
                                                                                    /// should point to successors of state s actually point to a transition
                                                                                    /// that starts in s.
                                                                                    ///
                                                                                    /// Additionally, it is asserted that the partition is stable. i. e. every
                                                                                    /// bottom state in every block can reach exactly the constellations in the
                                                                                    /// list of constellations that should be reachable from it, and every
                                                                                    /// nonbottom state can reach a subset of them.
                                                                                    void part_trans_t::assert_stability(const part_state_t& part_st) const
                                                                                    {
                                                                                        #ifdef MCRL2_GJKW_PARANOID_CHECK
                                                                                            std::vector<state_info_const_ptr> part_st_predecessors;
                                                                                        #endif
                                                                                        // count the nontrivial constellations (to check later whether every
                                                                                        // nontrivial constellation is reachable from the first nontrivial
                                                                                        // constellation)
                                                                                        state_type nr_of_nontrivial_constellations = 0;
                                                                                        const constln_t* C = constln_t::get_some_nontrivial();
                                                                                        if (nullptr != C)
                                                                                        {
                                                                                            for (;;)
                                                                                            {
                                                                                                ++nr_of_nontrivial_constellations;
                                                                                                if (C->get_nontrivial_next()
                                                                                                    == C)
                                                                                                {
                                                                                                  break;
                                                                                                }
                                                                                                C = C->get_nontrivial_next();
                                                                                                assert(nullptr != C);
                                                                                            }
                                                                                        }

                                                                                        // for all constellations C do
                                                                                        C = part_st.permutation.front()->constln();
                                                                                        assert(C->begin() == part_st.permutation.begin());
                                                                                        for (;;)
                                                                                        {
                                                                                            // assert some properties of constellation C
                                                                                            assert(C->begin() < C->end());
                                                                                            assert(C->postprocess_begin == C->postprocess_end);
                                                                                            unsigned const max_C = check_complexity::log_n -
                                                                                                                            check_complexity::ilog2(C->size());
                                                                                            const block_t* B = (*C->begin())->block;
                                                                                            assert(C->begin() == B->begin());
                                                                                            if (C->is_trivial())
                                                                                            {
                                                                                                // a trivial constellation contains exactly one block
                                                                                                assert(B->end() == C->end());
                                                                                            }
                                                                                            else
                                                                                            {
                                                                                                // a nontrivial constellation contains at least two blocks
                                                                                                assert(B->end() < C->end());
                                                                                                --nr_of_nontrivial_constellations;
                                                                                            }
                                                                                            // for all blocks B in C do
                                                                                            for (;;)
                                                                                            {
                                                                                                // assert some properties of block B
                                                                                                assert(B->constln() == C);
                                                                                                assert(B->nonbottom_begin() <= B->nonbottom_end());
                                                                                                assert(B->marked_nonbottom_begin()==B->marked_nonbottom_end());
                                                                                                assert(B->bottom_begin() < B->bottom_end());
                                                                                                assert(B->marked_bottom_begin() == B->marked_bottom_end());
                                                                                                assert(!B->is_refinable());
                                                                                                assert(B->inert_begin() <= B->inert_end());
                                                                                                unsigned const max_B = check_complexity::log_n -
                                                                                                                            check_complexity::ilog2(B->size());
                                                                                                mCRL2complexity(B, no_temporary_work(max_C, max_B), );
                                                                                                // count inert transitions in block
                                                                                                state_type nr_of_inert_successors =
                                                                                                                               B->inert_end()-B->inert_begin();
                                                                                                permutation_const_iter_t s_iter;
                                                                                                #ifdef MCRL2_GJKW_PARANOID_CHECK
                                                                                                    state_type nr_of_inert_predecessors=nr_of_inert_successors;
                                                                                                    // make sure that every non-bottom state can reach some
                                                                                                    // bottom  state  in  the  block.   This  is  done  using
                                                                                                    // a  simple  graph  algorithm  for  reachability.   The
                                                                                                    // algorithm  should  mark  every  source  of  an  inert
                                                                                                    // transition exactly once.  For this, we misuse the field
                                                                                                    // state_info_entry::notblue:  It is == STATE_TYPE_MAX if
                                                                                                    // and only if the state has been marked.

                                                                                                    // Because we run through the bottom states and their
                                                                                                    // predecessor transitions anyway, we also verify a few
                                                                                                    // other properties, in particular everything we want to
                                                                                                    // verify about inert predecessor transitions.

                                                                                                    // for all bottom states s in B do
                                                                                                    assert(part_st_predecessors.empty());
                                                                                                    part_st_predecessors.reserve(B->nonbottom_end() -
                                                                                                                                         B->nonbottom_begin());
                                                                                                    s_iter = B->bottom_begin();
                                                                                                    do
                                                                                                    {
                                                                                                        state_info_const_ptr const s = *s_iter;
                                                                                                        nr_of_inert_predecessors -=
                                                                                                                   s->inert_pred_end() - s->inert_pred_begin();
                                                                                                        // for all inert predecessors pred of s do
                                                                                                        for (pred_const_iter_t pred_iter=s->inert_pred_begin();
                                                                                                                  pred_iter < s->inert_pred_end(); ++pred_iter)
                                                                                                        {
                                                                                                            // assert some properties of the predecessor
                                                                                                            // transition
                                                                                                            assert(pred_iter->succ->B_to_C->B_to_C_slice->
                                                                                                                                            from_block() == B);
                                                                                                            assert(pred_iter->succ->B_to_C->B_to_C_slice->
                                                                                                                                            to_constln() == C);
                                                                                                            assert(pred_iter->succ->target == s);

                                                                                                            state_info_const_ptr const pred=pred_iter->source;
                                                                                                            // if pred is not yet marked as predecessor then
                                                                                                            if (STATE_TYPE_MAX != pred->notblue)
                                                                                                            {
                                                                                                                // assert some properties of the predecessor
                                                                                                                // state
                                                                                                                assert(pred->block == B);
                                                                                                                assert(pred->pos < B->nonbottom_end());
                                                                                                                // mark pred as predecessor
                                                                                                                const_cast<state_type&>(pred->notblue) =
                                                                                                                                                STATE_TYPE_MAX;
                                                                                                                // add pred to the list of predecessors
                                                                                                                part_st_predecessors.push_back(pred);
                                                                                                            // end if
                                                                                                            }
                                                                                                        // end for
                                                                                                        }
                                                                                                    // end for
                                                                                                    }
                                                                                                    while(++s_iter < B->bottom_end());

                                                                                                    // Now that we have collected the predecessors of the
                                                                                                    // bottom states, we have to extend this set to their
                                                                                                    // indirect predecessors.

                                                                                                    // for all states s in the list of predecessors do
                                                                                                    for (std::vector<state_info_const_ptr>::iterator s_iter =
                                                                                                                  part_st_predecessors.begin();
                                                                                                                 s_iter < part_st_predecessors.end(); ++s_iter)
                                                                                                    {
                                                                                                        state_info_const_ptr const s = *s_iter;
                                                                                                        nr_of_inert_predecessors -=
                                                                                                                   s->inert_pred_end() - s->inert_pred_begin();
                                                                                                        // for all inert predecessors pred of s do
                                                                                                        for (pred_const_iter_t pred_iter=s->inert_pred_begin();
                                                                                                                  pred_iter < s->inert_pred_end(); ++pred_iter)
                                                                                                        {
                                                                                                            // assert some properties of the predecessor
                                                                                                            // transition
                                                                                                            assert(pred_iter->succ->B_to_C->B_to_C_slice->
                                                                                                                                            from_block() == B);
                                                                                                            assert(pred_iter->succ->B_to_C->B_to_C_slice->
                                                                                                                                            to_constln() == C);
                                                                                                            assert(pred_iter->succ->target == s);

                                                                                                            state_info_const_ptr const pred=pred_iter->source;
                                                                                                            // if pred is not yet marked as predecessor then
                                                                                                            if (STATE_TYPE_MAX != pred->notblue)
                                                                                                            {
                                                                                                                // assert some properties of the predecessor
                                                                                                                // state
                                                                                                                assert(pred->block == B);
                                                                                                                assert(pred->pos < B->nonbottom_end());
                                                                                                                // mark pred as predecessor
                                                                                                                const_cast<state_type&>(pred->notblue) =
                                                                                                                                                STATE_TYPE_MAX;
                                                                                                                // add pred to the list of predecessors
                                                                                                                part_st_predecessors.push_back(pred);
                                                                                                            // end if
                                                                                                            }
                                                                                                        // end for
                                                                                                        }
                                                                                                    // end for
                                                                                                    }
                                                                                                    assert(part_st_predecessors.size() == (std::size_t)
                                                                                                                  (B->nonbottom_end() - B->nonbottom_begin()));
                                                                                                    part_st_predecessors.clear();
                                                                                                    assert(0 == nr_of_inert_predecessors);
                                                                                                    // now all nonbottom states should have
                                                                                                    // s->notblue == STATE_TYPE_MAX.
                                                                                                #endif
                                                                                                // verify to_constln list:
                                                                                                bool to_own_constln = false;
                                                                                                for (B_to_C_desc_const_iter_t iter = B->to_constln.begin();
                                                                                                                           B->to_constln.end() != iter; ++iter)
                                                                                                {
                                                                                                    assert(iter->from_block() == B);
                                                                                                    mCRL2complexity(iter, no_temporary_work(
                                                                                                       check_complexity::log_n -
                                                                                                       check_complexity::ilog2(iter->to_constln()->size())), );
                                                                                                    if (iter->to_constln() == C)
                                                                                                    {
                                                                                                        assert(!to_own_constln);
                                                                                                        to_own_constln = true;
                                                                                                        assert(iter->begin <= B->inert_begin());
                                                                                                        assert(iter->end == B->inert_end());
                                                                                                    }
                                                                                                }
                                                                                                if (!to_own_constln)
                                                                                                {
                                                                                                    assert(B->nonbottom_begin() == B->nonbottom_end());
                                                                                                    assert(B->inert_begin() == B_to_C_begin());
                                                                                                    assert(B->inert_end() == B_to_C_begin());
                                                                                                }

                                                                                                // for all states s in B do
                                                                                                s_iter = B->begin();
                                                                                                do
                                                                                                {
                                                                                                    state_info_const_ptr const s = *s_iter;
                                                                                                    // assert some properties of state s
                                                                                                    assert(s->pos == s_iter);
                                                                                                    assert(s->block == B);
                                                                                                    assert(s->pred_begin() <= s->inert_pred_begin());
                                                                                                    assert(s->inert_pred_begin() <= s->inert_pred_end());
                                                                                                    assert(s->inert_pred_end() == s->pred_end());
                                                                                                    assert(s->succ_begin() <= s->inert_succ_begin());
                                                                                                    assert(s->inert_succ_begin() <= s->inert_succ_end());
                                                                                                    assert(s->inert_succ_end() <= s->succ_end());
                                                                                                    assert(s->succ_begin() == s->current_constln() ||
                                                                                                             s->succ_end() == s->current_constln() ||
                                                                                                                  *s->current_constln()[-1].target->constln() <
                                                                                                                     *s->current_constln()->target->constln());
                                                                                                    mCRL2complexity(s, no_temporary_work(max_B,
                                                                                                                               s_iter >= B->bottom_begin()), );
                                                                                                    // count reachable constellations
                                                                                                    state_type nr_of_reachable_constlns = B->to_constln.size();

                                                                                                    // for all constln-slices of successors of s do
                                                                                                    if (succ_const_iter_t succ_iter = s->succ_begin();
                                                                                                                                     succ_iter < s->succ_end())
                                                                                                    {
                                                                                                        for (;;)
                                                                                                        {
                                                                                                            succ_const_iter_t before_end =
                                                                                                                        succ_iter->slice_begin_or_before_end();
                                                                                                            succ_const_iter_t slice_end = before_end + 1;
                                                                                                            assert(succ_iter < slice_end);
                                                                                                            assert(before_end->slice_begin_or_before_end() ==
                                                                                                                                                    succ_iter);
                                                                                                            const constln_t* const targetC =
                                                                                                                                  succ_iter->target->constln();
                                                                                                            // for all noninert transitions in the
                                                                                                            // constln-slice do
                                                                                                            if (targetC == C)
                                                                                                            {
                                                                                                                assert(s->inert_succ_end() == slice_end);
                                                                                                                slice_end = s->inert_succ_begin();
                                                                                                            }
                                                                                                            for (; succ_iter < slice_end; ++succ_iter)
                                                                                                            {
                                                                                                                // assert some properties of the successor
                                                                                                                // transition
                                                                                                                assert(succ_iter->target->block != B);
                                                                                                                assert(succ_iter->target->constln()==targetC);
                                                                                                                if (succ_iter != before_end)
                                                                                                                {
                                                                                                                    assert(succ_iter->
                                                                                                                      slice_begin_or_before_end()==before_end);
                                                                                                                }
                                                                                                                assert(succ_iter->B_to_C < B->inert_begin() ||
                                                                                                                          B->inert_end() <= succ_iter->B_to_C);
                                                                                                                assert(succ_iter->B_to_C->B_to_C_slice->
                                                                                                                                            from_block() == B);
                                                                                                                assert(succ_iter->B_to_C->B_to_C_slice->
                                                                                                                                      to_constln() == targetC);
                                                                                                                assert(succ_iter->B_to_C->pred->succ ==
                                                                                                                                                    succ_iter);
                                                                                                                assert(succ_iter->B_to_C->pred->source == s);
                                                                                                                mCRL2complexity(succ_iter->B_to_C->pred,
                                                                                                                    no_temporary_work(max_B, check_complexity::
                                                                                                                       log_n - check_complexity::ilog2(
                                                                                                                         succ_iter->target->constln()->size()),
                                                                                                                   check_complexity::log_n - check_complexity::
                                                                                                                       ilog2(succ_iter->target->block->size()),
                                                                                                                               s_iter >= B->bottom_begin()), );
                                                                                                            // end for
                                                                                                            }
                                                                                                            // if we have reached the inert transitions then
                                                                                                            if (targetC == C)
                                                                                                            {
                                                                                                                // for all inert transitions in the
                                                                                                                // constln-slice do
                                                                                                                for (slice_end = s->inert_succ_end();
                                                                                                                            succ_iter < slice_end; ++succ_iter)
                                                                                                                {
                                                                                                                    // assert some properties of inert
                                                                                                                    // transitions
                                                                                                                    assert(succ_iter->target->block == B);
                                                                                                                    if (succ_iter != before_end)
                                                                                                                    {
                                                                                                                        assert(before_end == succ_iter->
                                                                                                                                  slice_begin_or_before_end());
                                                                                                                    }
                                                                                                                    assert(B->inert_begin() <=
                                                                                                                                            succ_iter->B_to_C);
                                                                                                                    assert(succ_iter->B_to_C < B->inert_end());
                                                                                                                    assert(succ_iter->B_to_C->B_to_C_slice->
                                                                                                                                            from_block() == B);
                                                                                                                    assert(succ_iter->B_to_C->B_to_C_slice->
                                                                                                                                      to_constln() == targetC);
                                                                                                                    assert(succ_iter->B_to_C->pred->succ ==
                                                                                                                                                    succ_iter);
                                                                                                                    assert(succ_iter->B_to_C->pred->source==s);
                                                                                                                    mCRL2complexity(succ_iter->B_to_C->pred,
                                                                                                                                no_temporary_work(max_B,
                                                                                                                                       max_C, max_B, false), );
                                                                                                                // end for
                                                                                                                }
                                                                                                            // end if
                                                                                                            }
                                                                                                            else
                                                                                                            {
                                                                                                                --nr_of_reachable_constlns;
                                                                                                            }
                                                                                                    // end for
                                                                                                            if (s->succ_end()
                                                                                                                <= succ_iter)
                                                                                                            {
                                                                                                              break;
                                                                                                            }
                                                                                                            assert(0 != nr_of_reachable_constlns);
                                                                                                            assert(*targetC < *succ_iter->target->constln());
                                                                                                        }
                                                                                                    }
                                                                                                    // if s is a nonbottom state then
                                                                                                    if (s_iter < B->bottom_begin())
                                                                                                    {
                                                                                                        assert(s->inert_succ_begin() < s->inert_succ_end());
                                                                                                        nr_of_inert_successors -=
                                                                                                                   s->inert_succ_end() - s->inert_succ_begin();
                                                                                                        // the following assertion is necessary because s must
                                                                                                        // have a transition to its own constellation (namely
                                                                                                        // an inert one), but this constln_slice is not
                                                                                                        // counted.
                                                                                                        assert(0 != nr_of_reachable_constlns);

                                                                                                        #ifdef MCRL2_GJKW_PARANOID_CHECK
                                                                                                            // assert that s can reach a bottom state
                                                                                                            assert(STATE_TYPE_MAX == s->notblue);
                                                                                                            const_cast<state_type&>(s->notblue) = 1;
                                                                                                        #endif
                                                                                                    }
                                                                                                    else
                                                                                                    {
                                                                                                        // (s is a bottom state.)
                                                                                                        // assert that not too few constellations are
                                                                                                        // reachable.
                                                                                                        assert((state_type) to_own_constln ==
                                                                                                                                     nr_of_reachable_constlns);
                                                                                                        assert(s->inert_succ_begin() == s->inert_succ_end());
                                                                                                        // the following assertions are necessary because it
                                                                                                        // could be that the state has no transition to its own
                                                                                                        // constellation.
                                                                                                        assert(s->succ_begin() == s->inert_succ_begin() ||
                                                                                                           *s->inert_succ_begin()[-1].target->constln() <= *C);
                                                                                                        assert(s->inert_succ_end() == s->succ_end() ||
                                                                                                                 *s->inert_succ_end()->target->constln() > *C);
                                                                                                    }

                                                                                                    // for all noninert predecessors of s do
                                                                                                    for (pred_const_iter_t pred_iter=s->noninert_pred_begin();
                                                                                                               pred_iter < s->noninert_pred_end(); ++pred_iter)
                                                                                                    {
                                                                                                        // assert some properties of the predecessor
                                                                                                        assert(pred_iter->succ->B_to_C->B_to_C_slice->
                                                                                                                                            from_block() != B);
                                                                                                        assert(pred_iter->succ->B_to_C->B_to_C_slice->
                                                                                                                                            to_constln() == C);
                                                                                                        assert(pred_iter->succ->target == s);
                                                                                                        assert(pred_iter->source->block != B);
                                                                                                    // end for
                                                                                                    }
                                                                                                    #ifndef MCRL2_GJKW_PARANOID_CHECK
                                                                                                        // for all inert predecessors of s do
                                                                                                        for (pred_const_iter_t pred_iter=s->inert_pred_begin();
                                                                                                                  pred_iter < s->inert_pred_end(); ++pred_iter)
                                                                                                        {
                                                                                                            // assert some properties of the predecessor
                                                                                                            assert(pred_iter->succ->B_to_C->B_to_C_slice->
                                                                                                                                            from_block() == B);
                                                                                                            assert(pred_iter->succ->B_to_C->B_to_C_slice->
                                                                                                                                            to_constln() == C);
                                                                                                            assert(pred_iter->succ->target == s);
                                                                                                            assert(pred_iter->source->block == B);
                                                                                                            assert(pred_iter->source->pos<B->nonbottom_end());
                                                                                                        // end for
                                                                                                        }
                                                                                                    #endif
                                                                                                // end for (all states s in B)
                                                                                                }
                                                                                                while (++s_iter < B->end());
                                                                                                assert(0 == nr_of_inert_successors);
                                                                                            // end for (all blocks B in C)
                                                                                                /* *************************************************************************
                                                                                                 */
                                                                                                if (B->end()
                                                                                                    == C->end())
                                                                                                {
                                                                                                  break;
                                                                                                }
/*                                                                           */                 assert(B->end() < C->end());
/*                            A L G O R I T H M S                            */                 assert(B->end() == (*B->end())->block->begin());
/*                                                                           */                 B = (*B->end())->block;
/* ************************************************************************* */             }
                                                                                        // end for (all constellations C)
                                                                                            if (C->end()
                                                                                                == part_st.permutation
                                                                                                    .end())
                                                                                            {
                                                                                              break;
                                                                                            }
                                                                                            assert(C->end() < part_st.permutation.end());
                                                                                            assert(C->end() == (*C->end())->constln()->begin());
                                                                                            C = (*C->end())->constln();
/*===========================================================================*/         }
/*                           initialisation helper                           */         assert(0 == nr_of_nontrivial_constellations);
/*===========================================================================*/         return;
                                                                                    }
                                                                                #endif
/// \brief constructor of the helper class
template<class LTS_TYPE>
bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>::
bisim_partitioner_gjkw_initialise_helper(LTS_TYPE& l, bool const branching,
                                                bool const preserve_divergence)
  : aut(l),
    nr_of_states(l.num_states()),
    orig_nr_of_states(l.num_states()),
    nr_of_transitions(l.num_transitions()),
    noninert_out_per_state(l.num_states(), 0),
    inert_out_per_state(l.num_states(), 0),
    noninert_in_per_state(l.num_states(), 0),
    inert_in_per_state(l.num_states(), 0),
    noninert_out_per_block(1, 0),
    inert_out_per_block(1, 0),
    states_per_block(1, l.num_states())
{

    mCRL2log(log::verbose) << "O(m log n) "
                    << (preserve_divergence ? "Divergence preserving b" : "B")
                    << (branching ? "ranching b" : "")
                    << "isimulation partitioner created for " << l.num_states()
                    << " states and " << l.num_transitions()
                    << " transitions [GJKW 2017]\n";
    // Iterate over the transitions and collect new states
    for (const transition& t: aut.get_transitions())
    {
        if (!branching || !aut.is_tau(aut.apply_hidden_label_map(t.label())) ||
                                   (preserve_divergence && t.from() == t.to()))
        {
            // (possibly) create new state
            Key const k(aut.apply_hidden_label_map(t.label()), t.to());
            std::pair<typename std::unordered_map<Key, state_type,
                KeyHasher>::iterator,bool> extra_state = extra_kripke_states.
                                    insert(std::make_pair(k, nr_of_states));
            if (extra_state.second)
            {
                noninert_in_per_state.push_back(0);
                inert_in_per_state.push_back(0);
                noninert_out_per_state.push_back(0);
                inert_out_per_state.push_back(0);

                // (possibly) create new block
                std::pair<std::unordered_map<label_type, state_type>::iterator,
                    bool> const action_block = action_block_map.insert(
                          std::make_pair(aut.apply_hidden_label_map(t.label()),
                                                     states_per_block.size()));
                if (action_block.second)
                {
                    noninert_out_per_block.push_back(0);
                    inert_out_per_block.push_back(0);
                    states_per_block.push_back(0);
                }

                ++noninert_in_per_state[t.to()];
                ++noninert_out_per_state[nr_of_states];
                ++noninert_out_per_block[action_block.first->second];
                ++states_per_block[action_block.first->second];
                ++nr_of_states;
                ++nr_of_transitions;
            }
            ++noninert_in_per_state[extra_state.first->second];
            ++noninert_out_per_state[t.from()];
            ++noninert_out_per_block[0];
        }
        else
        {
            ++inert_in_per_state[t.to()];
            if (1 == ++inert_out_per_state[t.from()])
            {
                // this is the first inert outgoing transition of t.from()
                ++nr_of_nonbottom_states;
            }
            ++inert_out_per_block[0];
        }
    }
    mCRL2log(log::verbose) << "Number of extra states: "
                                         << extra_kripke_states.size() << "\n";
                                                                                #ifndef NDEBUG
                                                                                    check_complexity::init(nr_of_states);
                                                                                #endif
}


/// \brief initialise the state in part_st and the transitions in part_tr
template<class LTS_TYPE>
inline void bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>::
init_transitions(part_state_t& part_st, part_trans_t& part_tr,
                          bool const branching, bool const preserve_divergence)
{                                                                               assert(part_st.state_size() == get_nr_of_states());
                                                                                assert(part_tr.trans_size() == get_nr_of_transitions());
    // initialise blocks and B_to_C slices
    permutation_iter_t begin = part_st.permutation.begin();
    constln_t* const constln = new constln_t(get_nr_of_states(), begin,
                              part_st.permutation.end(), part_tr.B_to_C_end());
    if (1 < states_per_block.size())
    {
        constln->make_nontrivial();
    }
    std::vector<block_t*> blocks(states_per_block.size());
    B_to_C_iter_t B_to_C_begin = part_tr.B_to_C.begin();
    for (state_type B = 0; B < states_per_block.size(); ++B)
    {
        permutation_iter_t const end = begin + states_per_block[B];
        blocks[B] = new block_t(constln, begin, end);
        if (0 == noninert_out_per_block[B] && 0 == inert_out_per_block[B])
        {
            blocks[B]->set_inert_begin_and_end(part_tr.B_to_C.begin(),
                                                       part_tr.B_to_C.begin()); assert(blocks[B]->to_constln.empty());
        }
        else
        {
            blocks[B]->set_inert_begin_and_end(B_to_C_begin +
                                                     noninert_out_per_block[B],
                B_to_C_begin + noninert_out_per_block[B] +
                                                       inert_out_per_block[B]);
            blocks[B]->to_constln.emplace_back(B_to_C_begin,
                                                       blocks[B]->inert_end());
            B_to_C_desc_iter_t const slice =
                                        std::prev(blocks[B]->to_constln.end()); assert(B_to_C_begin < slice->end);
            for (; slice->end != B_to_C_begin; ++B_to_C_begin)
            {
                B_to_C_begin->B_to_C_slice = slice;
            }
        }
        begin = end;
    }                                                                           assert(part_st.permutation.end() == begin);
    /* only block 0 has a sequence number and non-bottom states:             */ assert(part_tr.B_to_C.end() == B_to_C_begin);
    blocks[0]->assign_seqnr();
    blocks[0]->set_bottom_begin(blocks[0]->begin() + nr_of_nonbottom_states);
    blocks[0]->set_marked_nonbottom_begin(blocks[0]->bottom_begin());

    // initialise states and succ slices
    part_st.state_info.front().set_pred_begin(part_tr.pred.begin());
    part_st.state_info.front().set_succ_begin(part_tr.succ.begin());
    for (state_type s = 0; get_nr_of_states() != s; ++s)
    {
        part_st.state_info[s].set_pred_end(part_st.state_info[s].pred_begin() +
                             noninert_in_per_state[s] + inert_in_per_state[s]);
        part_st.state_info[s].set_inert_pred_begin(part_st.state_info[s].
                                      pred_begin() + noninert_in_per_state[s]);
        // part_st.state_info[s+1].set_pred_begin(part_st.state_info[s].
        //                                                         pred_end());

        succ_iter_t succ_iter = part_st.state_info[s].succ_begin();
        succ_iter_t succ_end = succ_iter +
                            noninert_out_per_state[s] + inert_out_per_state[s];
        part_st.state_info[s].set_succ_end(succ_end);
        part_st.state_info[s].set_current_constln(succ_end);
        part_st.state_info[s].set_inert_succ_begin_and_end(part_st.
            state_info[s].succ_begin() + noninert_out_per_state[s], succ_end);
        if (succ_iter < succ_end)
        {
            --succ_end;
            for (; succ_iter < succ_end; ++succ_iter)
            {
                succ_iter->set_slice_begin_or_before_end(succ_end);
            }                                                                   assert(succ_iter == succ_end);
            succ_end->set_slice_begin_or_before_end(
                                           part_st.state_info[s].succ_begin());
        }
        else
        {
          assert(succ_end == succ_iter);
        }
        if (s < aut.num_states())
        {
            // s is not an extra Kripke state.  It is in block 0.
            part_st.state_info[s].block = blocks[0];
            if (0 != inert_out_per_state[s])
            {
                /* non-bottom state:                                         */ assert(0 != nr_of_nonbottom_states);
                --nr_of_nonbottom_states;
                part_st.state_info[s].pos = blocks[0]->begin() +
                                                        nr_of_nonbottom_states;
            }
            else
            {                                                                   // The following assertion is incomplete; only the second
                // bottom state:                                                // assertion (after the assignment) makes sure that not too
                                                                                // many states become part of this slice.
                                                                                assert(0 != states_per_block[0]);
                --states_per_block[0];
                part_st.state_info[s].pos = blocks[0]->begin() +
                                                           states_per_block[0]; assert(part_st.state_info[s].pos >= blocks[0]->bottom_begin());
            }
            *part_st.state_info[s].pos = &part_st.state_info[s];
            // part_st.state_info[s].notblue = 0;
        }
    }

    // initialise transitions (and finalise extra Kripke states)
    for (const transition& t: aut.get_transitions())
    {
        if (!branching || !aut.is_tau(aut.apply_hidden_label_map(t.label())) ||
                                   (preserve_divergence && t.from() == t.to()))
        {
            // take transition through an extra intermediary state
            Key const k(aut.apply_hidden_label_map(t.label()), t.to());
            state_type const extra_state = extra_kripke_states[k];
            if (0 != noninert_out_per_state[extra_state])
            {
                state_type const extra_block =
                       action_block_map[aut.apply_hidden_label_map(t.label())];
                // now initialise extra_state correctly
                part_st.state_info[extra_state].block = blocks[extra_block];    assert(0 != states_per_block[extra_block]);
                --states_per_block[extra_block];
                part_st.state_info[extra_state].pos = blocks[extra_block]->
                                       begin() + states_per_block[extra_block];
                *part_st.state_info[extra_state].pos =
                                              &part_st.state_info[extra_state];
                // part_st.state_info[extra_state].notblue = 0;

                // state extra_state has exactly one outgoing transition,
                // namely a noninert transition to to t.to().  It has to be
                /* initialised now.                                          */ assert(0 != noninert_in_per_state[t.to()]);
                --noninert_in_per_state[t.to()];
                pred_iter_t const t_pred =
                             part_st.state_info[t.to()].noninert_pred_begin() +
                                                 noninert_in_per_state[t.to()];
                --noninert_out_per_state[extra_state];                          assert(0 == noninert_out_per_state[extra_state]);
                succ_iter_t const t_succ =
                                  part_st.state_info[extra_state].succ_begin(); assert(0 != noninert_out_per_block[extra_block]);
                B_to_C_iter_t const t_B_to_C =
                                         blocks[extra_block]->inert_begin() -
                                         noninert_out_per_block[extra_block]--;
                t_pred->source = &part_st.state_info[extra_state];
                t_pred->succ = t_succ;
                t_succ->target = &part_st.state_info[t.to()];
                t_succ->B_to_C = t_B_to_C;
                // t_B_to_C->B_to_C_slice = (already initialised);
                t_B_to_C->pred = t_pred;
            }
            /* noninert transition from t.from() to extra_state              */ assert(0 != noninert_in_per_state[extra_state]);
            --noninert_in_per_state[extra_state];
            pred_iter_t const t_pred =
                        part_st.state_info[extra_state].noninert_pred_begin() +
                                            noninert_in_per_state[extra_state]; assert(0 != noninert_out_per_state[t.from()]);
            --noninert_out_per_state[t.from()];
            succ_iter_t const t_succ=part_st.state_info[t.from()].succ_begin()+
                                              noninert_out_per_state[t.from()]; assert(0 != noninert_out_per_block[0]);
            B_to_C_iter_t const t_B_to_C = blocks[0]->inert_begin() -
                                                   noninert_out_per_block[0]--;

            t_pred->source = &part_st.state_info[t.from()];
            t_pred->succ = t_succ;
            t_succ->target = &part_st.state_info[extra_state];
            t_succ->B_to_C = t_B_to_C;
            // t_B_to_C->B_to_C_slice = (already initialised);
            t_B_to_C->pred = t_pred;
        }
        else
        {
            /* inert transition from t.from() to t.to()                      */ assert(0 != inert_in_per_state[t.to()]);
            --inert_in_per_state[t.to()];
            pred_iter_t const t_pred =
                                part_st.state_info[t.to()].inert_pred_begin() +
                                                    inert_in_per_state[t.to()]; assert(0 != inert_out_per_state[t.from()]);
            --inert_out_per_state[t.from()];
            succ_iter_t const t_succ =
                              part_st.state_info[t.from()].inert_succ_begin() +
                                                 inert_out_per_state[t.from()]; assert(0 != inert_out_per_block[0]);
            --inert_out_per_block[0];
            B_to_C_iter_t const t_B_to_C = blocks[0]->inert_begin() +
                                                        inert_out_per_block[0];

            t_pred->source = &part_st.state_info[t.from()];
            t_pred->succ = t_succ;
            t_succ->target = &part_st.state_info[t.to()];
            t_succ->B_to_C = t_B_to_C;
            // t_B_to_C->B_to_C_slice = (already initialised);
            t_B_to_C->pred = t_pred;
        }
    }
    noninert_out_per_state.clear(); inert_out_per_state.clear();
    noninert_in_per_state.clear();  inert_in_per_state.clear();
    noninert_out_per_block.clear(); inert_out_per_block.clear();
    states_per_block.clear();

    aut.clear_transitions();

    mCRL2log(log::verbose) << "Size of the resulting Kripke structure: "
                               << get_nr_of_states() << " states and "
                               << get_nr_of_transitions() << " transitions.\n";
}

/// \brief Replaces the transition relation of the current LTS by the
/// transitions of the bisimulation-reduced transition system.
/// \details Each transition (s, l, s') is replaced by a transition (t, l, t'),
/// where t and t' are the equivalence classes of s and s', respectively.  If
/// the label l is internal, then the transition is only added if t != t' or
/// preserve_divergence == true.  This effectively removes all inert
/// transitions.  Duplicates are removed from the transitions in the new LTS.
///
/// Note that the number of states nor the initial state are not adapted by
/// this method.  These must be set separately.
///
/// The code is very much inspired by liblts_bisim_gw.h, which was written by
/// Anton Wijs.
///
/// \pre The bisimulation equivalence classes have been computed.
/// \note In debug mode, an additional parameter `branching` is available to
///       indicate that non-internal transitions should be removed.
template <class LTS_TYPE>
void bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>::
         replace_transition_system(const part_state_t& part_st,                 ONLY_IF_DEBUG( const bool branching, )
                                   const bool preserve_divergence)
{
    std::unordered_map <state_type, Key> to_lts_map;
    // obtain a map from state to <action, state> pair from extra_kripke_states
    for (typename std::unordered_map<Key, state_type, KeyHasher>::iterator it =
            extra_kripke_states.begin(); it != extra_kripke_states.end(); ++it)
    {
        to_lts_map.insert(std::make_pair(it->second, it->first));
    }
    extra_kripke_states.clear();

    label_type const tau_label = aut.tau_label_index();
    // In the following loop, we visit a bottom state of each block and take
    // its transitions.  As the partition is (assumed to be) stable, in this
    // way we visit each transition of a lumped state exactly once.
    for (permutation_const_iter_t s_iter = part_st.permutation.begin();
                                 part_st.permutation.end() != s_iter; ++s_iter)
    {
        block_t *B = (*s_iter)->block;
        // forward to last state of block, i. e. to a bottom state:
        s_iter = B->end() - 1;                                                  assert(B->bottom_end() > s_iter);  assert(B->bottom_begin() <= s_iter);
                                                                                assert(B->end() == B->constln()->end());
        state_type const s_eq = B->seqnr();
        if (BLOCK_NO_SEQNR == s_eq)
        {
            break;
        }

        for (succ_const_iter_t succ_iter = (*s_iter)->succ_begin();
                                          (*s_iter)->succ_end() != succ_iter; )
        {
            if (state_type t_eq = succ_iter->target->block->seqnr();
                                                        BLOCK_NO_SEQNR != t_eq)
            {                                                                   assert(branching);
                // We have a transition that originally was inert.
                if (s_eq == t_eq)
                {
                    // The transition is still inert.
                    if (!preserve_divergence)
                    {
                        // As we do not preserve divergence, we do not add it.
                        // Nor will we add other transitions to the same
                        // constellation.
                        succ_iter = succ_entry::slice_end(succ_iter);
                        continue;
                    }
                    if (*s_iter != succ_iter->target)
                    {
                        // The transition was not a self-loop to start with.
                        // So we do not add it either.
                        ++succ_iter;
                        continue;
                    }
                }                                                               assert(std::cmp_not_equal(-1, tau_label));
                aut.add_transition(transition(s_eq, tau_label, t_eq));
            }
            else
            {
                state_type const tgt_id = succ_iter->target -
                                                   &part_st.state_info.front();
                // We have a non-inert transition to an intermediary state.
                // Look up the label and where the transition from the
                // intermediary state goes.
                Key const k = to_lts_map.find(tgt_id)->second;
                t_eq = part_st.state_info[k.second].block->seqnr();             assert(BLOCK_NO_SEQNR != t_eq);
                aut.add_transition(transition(s_eq, k.first, t_eq));            // The target state could also be found through the pointer
                                                                                // structure (but we also need the labels, which are not stored
                                                                                // in the refinable partition):
                                                                                assert(&part_st.state_info[k.second]==succ_iter->target->succ_begin()->target);
                                                                                assert(1 == succ_iter->target->succ_end() - succ_iter->target->succ_begin());
            }
            // Skip over other transitions from the same state to the same
            // constellation -- they would be mapped to the same resulting
            // transition.
            succ_iter = succ_entry::slice_end(succ_iter);
        }
    }

    // Merge the states, by setting the state labels of each state to the concatenation of the state labels of its
    // equivalence class.

    if (aut.has_state_info())   /* If there are no state labels this step can be ignored */
    {
      /* Create a vector for the new labels */
      std::vector<typename LTS_TYPE::state_label_t> new_labels(block_t::nr_of_blocks);

      for(std::size_t i=aut.num_states(); i>0; )
      {
        --i;
        const std::size_t new_index=part_st.block(i)->seqnr();  /* get_eq_class(i) */
        new_labels[new_index] = new_labels[new_index] + aut.state_label(i);
      }

      aut.set_num_states(block_t::nr_of_blocks);
      for(std::size_t i=0; i<block_t::nr_of_blocks; ++i)
      {
        aut.set_state_label(i,new_labels[i]);
      }
    }
    else
    {
      aut.set_num_states(block_t::nr_of_blocks);
    }

    aut.set_initial_state(part_st.block(aut.initial_state())->seqnr());
}

} // end namespace bisim_gjkw



/*=============================================================================
=            dbStutteringEquivalence -- Algorithm 2 in [GJKW 2017]            =
=============================================================================*/



template <class LTS_TYPE>
void bisim_partitioner_gjkw<LTS_TYPE>::create_initial_partition_gjkw(
                          bool const branching, bool const preserve_divergence)
{
    // 2.2: P := P_0, i. e. the initial, cycle-free partition; C = {S}
    // and
    // 2.3: Initialise all temporary data
    init_helper.init_transitions(part_st, part_tr, branching,
                                                          preserve_divergence);
}


template <class LTS_TYPE>
void bisim_partitioner_gjkw<LTS_TYPE>::
                                refine_partition_until_it_becomes_stable_gjkw()
{
                                                                                #ifndef NDEBUG
                                                                                    if (mCRL2logEnabled(log::debug))
                                                                                    {
                                                                                        part_st.print_part(part_tr);
                                                                                        part_st.print_trans();
                                                                                    }

                                                                                    part_tr.assert_stability(part_st);
                                                                                #endif
    // 2.4: while C contains a nontrivial constellation SpC do
    while (nullptr != bisim_gjkw::constln_t::get_some_nontrivial())
    {                                                                           // check_complexity::add_work is called below, after SpB has been found
        bisim_gjkw::constln_t* const SpC =
                                  bisim_gjkw::constln_t::get_some_nontrivial();
        // 2.5: Choose a small splitter block SpB subset of SpC from P,
        //      i.e. |SpB| <= 1/2*|SpC|
        // and
        // 2.6: Create a new constellation NewC and move SpB from SpC to NewC
        // and
        // 2.7: C := partition C where SpB is removed from SpC and NewC is
        //      added
        bisim_gjkw::block_t* const SpB = SpC->split_off_small_block();
        bisim_gjkw::constln_t* const NewC = SpB->constln();
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    unsigned const max_counter = check_complexity::log_n -
        /*-------------------- find predecessors of SpB ---------------------*/                                           check_complexity::ilog2(SpB->size());
                                                                                #endif
        /* 2.8: Mark block SpB as refinable                                  */ mCRL2complexity(SpB, add_work(check_complexity::
        /* 2.9: Mark all states of SpB as predecessors                       */               while_C_contains_a_nontrivial_constellation_2_4, max_counter), );
            /* we deviate from the published algorithm: only states with a   */ assert(nullptr == bisim_gjkw::block_t::get_some_refinable());
            // transition to SpC will be marked.  SpB will be split separately.
        // 2.17: Register that the transitions from s to inert_out(s) go to
        //       NewC (instead of SpC)
            // (Before we enter the loop, we already adapt the ``B_to_C''
            // transition array.)
        part_tr.split_inert_to_C(SpB);
        // 2.10: for all s in SpB do
            /* We have to walk through the states from end to beginning so   */ mCRL2complexity(SpB,
            /* we can mark state s if necessary.  Marking will move s to a   */              add_work(check_complexity::for_all_s_in_SpB_2_10, max_counter), );
            // position that has already been visited.
        for (bisim_gjkw::permutation_iter_t s_iter = SpB->end();
                                                      SpB->begin() != s_iter; )
        {
            --s_iter;
            bisim_gjkw::state_info_ptr const s = *s_iter;                       // mCRL2complexity(s, ...) -- optimized to the above call.
            // 2.11: for all s_prime in noninert_in(s) do
            for (bisim_gjkw::pred_iter_t pred_iter = s->noninert_pred_begin();
                              s->noninert_pred_end() != pred_iter; ++pred_iter)
            {                                                                   mCRL2complexity(pred_iter, add_work(
                                                                                             check_complexity::for_all_s_prime_in_pred_s_2_11, max_counter), );
                                                                                assert(part_tr.pred_end() > pred_iter);
                                                                                assert(pred_iter->succ->B_to_C->pred == pred_iter);
                bisim_gjkw::state_info_ptr const s_prime = pred_iter->source;
                // 2.12: Mark the block of s_prime as refinable
                bool const first_transition_of_block =
                                              s_prime->block->make_refinable();
                // 2.13: Mark s_prime as predecessor of SpB
                bool const first_transition_of_state =
                                                 s_prime->block->mark(s_prime);
                // 2.14: Register that s_prime->s goes to NewC (instead of SpC)
                // and
                // 2.15: Store whether s' still has some transition to SpC\SpB
                s_prime->set_current_constln(part_tr.change_to_C(pred_iter,     ONLY_IF_DEBUG( SpC, NewC, )
                        first_transition_of_state, first_transition_of_block));
            // 2.16: end for
            }
            // 2.17: Register that the transitions from s to inert_out(s) go to
            //       NewC (instead of SpC)
                // (Here, we only adapt the ``succ'' transition array.)
            if (part_tr.split_s_inert_out(s                                     ONLY_IF_DEBUG( , SpC )
                                           ))
            {
                // 2.18: Store whether s still has some transition to SpC\SpB
                    // we deviate from the explanation of the published
                    // algorithm: we store this information by marking states
                    // with transitions to SpC.
                SpB->mark(s);
            }
        // 2.19: end for
        }

        // second pass through the predecessors to correct the pointers to the
        // constln_slices
        for (auto s : *SpB)
        {
            for (bisim_gjkw::pred_iter_t pred_iter = s->noninert_pred_begin();
                              s->noninert_pred_end() != pred_iter; ++pred_iter)
            {
                bisim_gjkw::succ_iter_t succ = pred_iter->succ;
                bisim_gjkw::succ_iter_t before_end =
                                      pred_iter->source->current_constln() - 1;
                if (succ != before_end)
                {                                                               assert(succ<before_end); assert(succ>=before_end->slice_begin_or_before_end());
                                                                                assert(succ->slice_begin_or_before_end() ==
                                                                                                                      before_end->slice_begin_or_before_end());
                    succ->set_slice_begin_or_before_end(before_end);
                }
            }
        }
                                                                                #ifndef NDEBUG
                                                                                    // The following tests cannot be executed during the above loops
                                                                                    // because a state s_prime may have multiple transitions to SpB.
                                                                                    for (auto s : *SpB)
                                                                                    {
                                                                                        for (bisim_gjkw::pred_const_iter_t pred_iter=s->noninert_pred_begin();
                                                                                                              s->noninert_pred_end() != pred_iter; ++pred_iter)
                                                                                        {
                                                                                            bisim_gjkw::state_info_const_ptr const s_prime = pred_iter->source;
                                                                                            // check consistency of s_prime->current_constln()
                                                                                            assert(s_prime->succ_begin() == s_prime->current_constln() ||
                                                                                                    *s_prime->current_constln()[-1].target->constln() <= *SpC);
                                                                                            assert(s_prime->succ_end() == s_prime->current_constln() ||
                                                                                                       *SpC <= *s_prime->current_constln()->target->constln());
                                                                                                // s_prime must have a transition to the new constellation
                                                                                            assert(s_prime->succ_begin() < s_prime->current_constln());
                                                                                            assert(s_prime->current_constln()[-1].target->constln() == NewC);
                                                                                            // check consistency of s_prime->inert_succ_begin() and
                                                                                            // s_prime->inert_succ_end()
                                                                                            assert(s_prime->succ_begin() == s_prime->inert_succ_begin() ||
                                                                                                      *s_prime->inert_succ_begin()[-1].target->constln() <=
                                                                                                                                          *s_prime->constln());
                                                                                            assert(s_prime->succ_begin() == s_prime->inert_succ_begin() ||
                                                                                                           s_prime->inert_succ_begin()[-1].target->block !=
                                                                                                                                               s_prime->block);
                                                                                            assert(s_prime->inert_succ_begin()==s_prime->inert_succ_end()||
                                                                                                      (s_prime->inert_succ_begin()->target->block ==
                                                                                                                                            s_prime->block &&
                                                                                                       s_prime->inert_succ_end()[-1].target->block ==
                                                                                                                                              s_prime->block));
                                                                                            assert(s_prime->succ_end() == s_prime->inert_succ_end() ||
                                                                                                      *s_prime->constln() <
                                                                                                                *s_prime->inert_succ_end()->target->constln());
                                                                                        }
                                                                                        // check consistency of s->inert_succ_begin() and
                                                                                        // s->inert_succ_end()
                                                                                        assert(s->succ_begin() == s->inert_succ_begin() ||
                                                                                                 *s->inert_succ_begin()[-1].target->constln() < *s->constln());
                                                                                        assert(s->inert_succ_begin() == s->inert_succ_end() ||
        /*------------------ stabilise the partition again ------------------*/                       (s->inert_succ_begin()->target->block == s->block &&
                                                                                                           s->inert_succ_end()[-1].target->block == s->block));
        /* deviation from the published algorithm: we first refine the       */         assert(s->succ_end() == s->inert_succ_end() ||
        /* splitter according to the marking of states (marked states have a */                       *s->constln() < *s->inert_succ_end()->target->constln());
        /* transition to SpC).                                               */     }
                                                                                #endif // ifndef NDEBUG
        if (0 != SpB->marked_size())
        {
            if (1 == SpB->size())
            {
                // if the block only contains a single state, the refinement
                // would be trivial anyway (and it cannot find new bottom
                // states).
                SpB->set_marked_bottom_begin(SpB->bottom_end());
            }
            else
            {
                bisim_gjkw::block_t* RedB = refine(SpB, SpC, nullptr, false);
                if (0 != RedB->unmarked_bottom_size())
                {
                    postprocess_new_bottom(RedB);
                }
                else
                {
                    RedB->set_marked_bottom_begin(RedB->bottom_end());
                }
            }
        }

        // 2.20: for all refinable blocks RfnB do
        while (nullptr != bisim_gjkw::block_t::get_some_refinable())
        {
            bisim_gjkw::block_t* const RfnB =
                                     bisim_gjkw::block_t::get_some_refinable();
            // 2.21: Mark block RfnB as non-refinable
            RfnB->make_nonrefinable();
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    // The work in this loop has to be assigned to the B_to_C-
                                                                                    // slice that contains the transitions from RfnB to NewC.
                                                                                    // Note that this is not FromRed, so we have to find it in a
                                                                                    // different way.  First find a marked state:
                                                                                    bisim_gjkw::permutation_iter_t s_mark_iter = RfnB->marked_bottom_begin();
                                                                                    if (s_mark_iter == RfnB->marked_bottom_end())
                                                                                    {
                                                                                        s_mark_iter = RfnB->marked_nonbottom_begin();
                                                                                        assert(s_mark_iter < RfnB->marked_nonbottom_end());
                                                                                    }
                                                                                    // Now, using the current_constln pointer of s, find a
                                                                                    // transition to NewC.
                                                                                    bisim_gjkw::succ_iter_t to_NewC =
                                                                                                                  std::prev((*s_mark_iter)->current_constln());
                                                                                    assert((*s_mark_iter)->succ_begin() <= to_NewC);
                                                                                    assert(to_NewC->target->constln() == NewC);
                                                                                    mCRL2complexity(to_NewC->B_to_C->B_to_C_slice, add_work(
                                                                                         check_complexity::for_all_refinable_blocks_RfnB_2_20, max_counter), );
                                                                                #endif
            if (1 == RfnB->size())
            {
                // if the block only contains a single state, the refinement
                // would be trivial anyway (and it cannot find new bottom
                // states).
                RfnB->set_marked_bottom_begin(RfnB->bottom_end());
                continue;
            }
            // 2.22: <RedB, BlueB> := Refine(RfnB, NewC, {marked states in
            //                                                       RfnB}, {})
            bisim_gjkw::block_t* RedB = refine(RfnB, NewC, nullptr, false);
            // 2.23: if RedB contains new bottom states then
            if (0 != RedB->unmarked_bottom_size())
            {
                // 2.24: RedB := PostprocessNewBottom(RedB, BlueB)
                RedB = postprocess_new_bottom(RedB);
                if (nullptr == RedB)
                {
                  continue;
                }
            // 2.25: end if
            }
            // 2.30: Unmark all states of the original RfnB as predecessors
                // (first part) This needs to be done to make sure the call to
                // Refine in line 2.26 does not regard some states as marked.
            else
            {
                RedB->set_marked_bottom_begin(RedB->bottom_end());
            }
                                                                                assert(0 == RedB->marked_size());
            if (1 == RedB->size() || nullptr == RedB->FromRed(SpC))
            {
                // If the block only contains a single state or has no
                // transitions to SpC, the refinement would be trivial anyway
                // (and it cannot find new bottom states).
                continue;
            }

            // 2.26: <RedB, BlueB> := Refine(RedB, SpC\SpB, {}, {transitions
            //                                           from RedB to SpC\SpB})
            RedB = refine(RedB, SpC, RedB->FromRed(SpC), false
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    , NewC
                                                                                #endif
                                                              );
            // 2.27: if RedB contains new bottom states then
            if (0 != RedB->unmarked_bottom_size())
            {
                // 2.28: PostprocessNewBottom(RedB, BlueB)
                postprocess_new_bottom(RedB);
            // 2.29: endif
            }
            // 2.30: Unmark all states of the original RfnB as predecessors
                // if postprocess_new_bottom is called, that procedure already
                // unmarks the new bottom states.
            else
            {
                RedB->set_marked_bottom_begin(RedB->bottom_end());
            }                                                                   assert(0 == RedB->marked_size());
        // 2.31: end for
        }
                                                                                #ifndef NDEBUG
                                                                                    if (mCRL2logEnabled(log::debug))
                                                                                    {
                                                                                        part_st.print_part(part_tr);
                                                                                        part_st.print_trans();
                                                                                    }

    /* 2.32: end while                                                       */     part_tr.assert_stability(part_st);
                                                                                #endif
    }
    // 2.33: return P
        // (this happens implicitly, through the bisim_partitioner_gjkw object
        // data)
    mCRL2log(log::verbose) << "number of blocks in the quotient: "
                                  << bisim_gjkw::block_t::nr_of_blocks << '\n';
}
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    namespace bisim_gjkw
                                                                                    {
/*=============================================================================
=                     Refine -- Algorithm 3 in [GJKW 2017]                    =
=============================================================================*/
                                                                                    /// \brief moves temporary counters to normal ones if the blue block is
                                                                                    /// smaller
                                                                                    /// \details When a refinement has finished and the blue block turns out to
                                                                                    /// be smaller, this function moves the corresponding temporary work to the
                                                                                    /// normal counters and cancels the work on the red state counters.
                                                                                    /// \param BlueB  pointer to the blue block
                                                                                    /// \param RedB   pointer to the red block
                                                                                    /// \param NewC   (only needed if called after the refinement in line 2.26)
                                                                                    ///               pointer to the constellation NewC, i. e. the
                                                                                    ///               constellation that was used as the basis of marking
                                                                                    ///               states.  If blue_is_smaller() is called from another
                                                                                    ///               place, this parameter should be nullptr.
                                                                                    static void blue_is_smaller(block_t* BlueB, block_t* RedB,
                                                                                                                                         const constln_t* NewC)
                                                                                    {
                                                                                        assert(nullptr != RedB);
                                                                                        if (nullptr != BlueB)
                                                                                        {
                                                                                            assert(BlueB->size() <= RedB->size());
                                                                                            unsigned const max_NewB = nullptr == BlueB ? 0
                                                                                              : check_complexity::log_n-check_complexity::ilog2(BlueB->size());
                                                                                            for (auto s : *BlueB)
                                                                                            {
                                                                                                mCRL2complexity(s, finalise_work(check_complexity::
                                                                                                                  while_Test_is_not_empty_3_6l_s_is_blue_3_11l,
                                                                                                                  check_complexity::refine_bottom_state_3_6l,
                                                                                                                  max_NewB), );
                                                                                                mCRL2complexity(s, finalise_work(check_complexity::
                                                                                                                   while_Blue_contains_unvisited_states_3_15l,
                                                                                                                   check_complexity::refine_visited_state_3_15,
                                                                                                                   max_NewB), );
                                                                                                for (succ_iter_t succ = s->succ_begin(); s->succ_end() != succ;
                                                                                                                                                        ++succ)
                                                                                                {
                                                                                                    mCRL2complexity(succ->B_to_C->pred,finalise_work(
                                                                                                           check_complexity::
                                                                                                           if_s_prime_has_transition_to_SpC_3_23l,
                                                                                                           check_complexity::
                                                                                                           refine_outgoing_transition_3_6_or_23l, max_NewB), );
                                                                                                }
                                                                                                for (pred_iter_t pred = s->pred_begin(); s->pred_end() != pred;
                                                                                                                                                        ++pred)
                                                                                                {
                                                                                                    mCRL2complexity(pred, finalise_work(check_complexity::
                                                                                                             for_all_s_prime_in_pred_s_setminus_Red_3_18l,
                                                                                                             check_complexity::refine_incoming_transition_3_18,
                                                                                                             max_NewB), );
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                        // cancel work counters for the red states, and also measure the work
                                                                                        // done in the blue coroutine on states that turned out to be red.
                                                                                        unsigned const max_NewC = nullptr == NewC ? 0
                                                                                             : check_complexity::log_n - check_complexity::ilog2(NewC->size());
                                                                                        for (auto s : *RedB)
                                                                                        {
                                                                                            mCRL2complexity(s, cancel_work(check_complexity::
                                                                                                                 while_Red_contains_unvisited_states_3_15r), );
                                                                                            for (succ_iter_t succ=s->succ_begin(); s->succ_end()!=succ; ++succ)
                                                                                            {
                                                                                                mCRL2complexity(succ->B_to_C->pred, cancel_work(
                                                                                                         check_complexity::while_FromRed_is_not_empty_3_6r), );
                                                                                                // the following counter should only be set for transitions to
                                                                                                // the splitter constellation.  For others transitions, no
                                                                                                // (temporary) work should have been done.
                                                                                                mCRL2complexity(succ->B_to_C->pred, finalise_work(
                                                                                                            check_complexity::
                                                                                                            while_Test_is_not_empty_3_6l_s_is_red_3_9l,
                                                                                                            check_complexity::
                                                                                                            refine_outgoing_transition_to_marked_state_3_6l,
                                                                                                            succ->target->constln() == NewC ? max_NewC : 0), );
                                                                                                // the following counter only gets 1 during postprocessing, at
                                                                                                // a time when state s is already stored as bottom state.
                                                                                                mCRL2complexity(succ->B_to_C->pred, finalise_work(
                                                                                                     check_complexity::
                                                                                                     while_Test_is_not_empty_3_6l_s_is_red_3_9l_postprocessing,
                                                                                                     check_complexity::
                                                                                                     refine_outgoing_transition_postprocess_new_bottom_3_6l,
                                                                                                     1), );
                                                                                                // the following counter only gets 1 before postprocessing, at
                                                                                                // the same time as when it is discovered to be a new bottom
                                                                                                // state
                                                                                                mCRL2complexity(succ->B_to_C->pred, finalise_work(
                                                                                                    check_complexity::if_s_prime_has_transition_to_SpC_3_23l,
                                                                                                    check_complexity::
                                                                                                    refine_outgoing_transition_from_new_bottom_3_23l, 1), );
                                                                                            }
                                                                                            for (pred_iter_t pred=s->pred_begin(); s->pred_end()!=pred; ++pred)
                                                                                            {
                                                                                                mCRL2complexity(pred, cancel_work(check_complexity::
                                                                                                                           for_all_s_prime_in_pred_s_3_18r), );
                                                                                            }
                                                                                        }
                                                                                        // check the balance between useful and cancelled work:
                                                                                        check_complexity::check_temporary_work();
                                                                                    }

                                                                                    /// \brief moves temporary counters to normal ones if the red block is
                                                                                    /// smaller
                                                                                    /// \details When a refinement has finished and the red block turns out to
                                                                                    /// be smaller, this function moves the corresponding temporary work to the
                                                                                    /// normal counters and cancels the work on the blue state counters.
                                                                                    /// \param BlueB  pointer to the blue block
                                                                                    /// \param RedB   pointer to the red block
                                                                                    static void red_is_smaller(block_t* BlueB, block_t* RedB)
                                                                                    {
                                                                                        assert(nullptr != BlueB);
                                                                                        assert(nullptr != RedB);
                                                                                        assert(BlueB->size() >= RedB->size());
                                                                                        for (auto s : *BlueB)
                                                                                        {
                                                                                            mCRL2complexity(s, cancel_work(check_complexity::
                                                                                                              while_Test_is_not_empty_3_6l_s_is_blue_3_11l), );
                                                                                            mCRL2complexity(s, cancel_work(check_complexity::
                                                                                                                while_Blue_contains_unvisited_states_3_15l), );
                                                                                            for (succ_iter_t succ=s->succ_begin(); s->succ_end()!=succ; ++succ)
                                                                                            {
                                                                                                mCRL2complexity(succ->B_to_C->pred, cancel_work(
                                                                                                                  check_complexity::
                                                                                                                  if_s_prime_has_transition_to_SpC_3_23l), );
                                                                                            }
                                                                                            for (pred_iter_t pred=s->pred_begin(); s->pred_end()!=pred; ++pred)
                                                                                            {
                                                                                                mCRL2complexity(pred, cancel_work(check_complexity::
                                                                                                              for_all_s_prime_in_pred_s_setminus_Red_3_18l), );
                                                                                            }
                                                                                        }
                                                                                        unsigned const max_NewB = check_complexity::log_n -
                                                                                                                         check_complexity::ilog2(RedB->size());
                                                                                        for (auto s : *RedB)
                                                                                        {
                                                                                            mCRL2complexity(s, finalise_work(check_complexity::
                                                                                                     while_Red_contains_unvisited_states_3_15r,
                                                                                                     check_complexity::refine_visited_state_3_15, max_NewB), );
                                                                                            for (succ_iter_t succ=s->succ_begin(); s->succ_end()!=succ; ++succ)
                                                                                            {
                                                                                                mCRL2complexity(succ->B_to_C->pred, finalise_work(
                                                                                                           check_complexity::while_FromRed_is_not_empty_3_6r,
                                                                                                           check_complexity::
                                                                                                           refine_outgoing_transition_3_6_or_23l, max_NewB), );
                                                                                                mCRL2complexity(succ->B_to_C->pred, cancel_work(
                                                                                                                check_complexity::
                                                                                                                while_Test_is_not_empty_3_6l_s_is_red_3_9l), );
                                                                                                mCRL2complexity(succ->B_to_C->pred, cancel_work(
                                                                                                 check_complexity::
                                                                                                 while_Test_is_not_empty_3_6l_s_is_red_3_9l_postprocessing), );
                                                                                                mCRL2complexity(succ->B_to_C->pred, cancel_work(
                                                                                                                  check_complexity::
                                                                                                                  if_s_prime_has_transition_to_SpC_3_23l), );
                                                                                            }
                                                                                            for (pred_iter_t pred=s->pred_begin(); s->pred_end()!=pred; ++pred)
                                                                                            {
                                                                                                mCRL2complexity(pred, finalise_work(check_complexity::
                                                                                                             for_all_s_prime_in_pred_s_3_18r,
                                                                                                             check_complexity::refine_incoming_transition_3_18,
                                                                                                             max_NewB), );
                                                                                            }
                                                                                        }
                                                                                        check_complexity::check_temporary_work();
                                                                                    }

                                                                                    } // end namespace bisim_gjkw
                                                                                #endif
/// \brief refine a block into the part that can reach `SpC` and the part that
/// cannot
///
/// \details At the beginning of a call to `refine()`, the block RfnB is sliced
/// into these parts:
///
///     |          unmarked         |marked |         unmarked         |marked|
///     |         non-bottom        |non-bot|          bottom          |bottom|
///
/// The procedure will partition the block into two parts, called ``blue'' and
/// ``red''.  The red part can reach the splitter constellation `SpC`, the blue
/// part cannot.  Upon calling `refine()`, it is known that all marked states
/// are red.
///
/// If states with a strong transition to `SpC` have not yet been found, the
/// parameter `FromRed` indicates all such transitions, so `refine()` can
/// easily find their sources.  It is required that either some states are
/// marked or a non-nullptr `FromRed` is supplied.  In the latter case, also
/// the _current constellation_ pointer of all unmarked bottom states is set in
/// a way that it makes easy to check whether the state has a transition to
/// `SpC`.  (Some non-bottom states may have their _current constellation_
/// pointer set similarly; while this is not required, it still speeds up
/// `refine()`.)
///
/// `refine()` works by starting two coroutines to find states in the parts;
/// as soon as the first coroutine finishes, the refinement can be completed.
///
/// Refine slices these states further as follows:
///
///                                          refine_blue()    refine_red()
///              shared variable            local variable   local variable
///          notblue_initialised_end          visited_end     visited_begin
///                     v                          v                v
///     | blue  |notblue|notblue|    red    | blue |unknown|not blue|   red   |
///     |non-bot|  > 0  |undef'd| non-bottom|bottom|bottom | bottom | bottom  |
///             ^               ^                          ^
///       refine_blue()    block.unmarked_          block.unmarked_
///      local variable    nonbottom_end()           bottom_end()
///     blue_nonbottom_end
///
/// - *blue non-bottom:* states of which it is known that they are blue
/// - *notblue > 0:* states that have an inert transition to some blue state,
///   but for some transition it is not yet known whether they go to a red or a
///   blue state
/// - *notblue undefined:* for no inert transition it is known whether it goes
///   to a red or a blue state
/// - *red non-bottom:* states of which it is known that they are red
/// - *blue bottom:* states of which it is known that they are blue
/// - *unknown bottom:* states that have not yet been checked
/// - *not blue bottom:* states of which the blue coroutine has found that they
///   are red, but the red coroutine has not yet handled them
/// - *red bottom:* states of which everybody knows that they are red
///
/// Note that the slices of red non-bottom and the red bottom states may become
/// larger than the marked (non-)bottom states have been at the start of
/// `refine()`.  Upon termination of one of the two coroutines, states whose
/// colour has not yet been determined become known to be of the colour of the
/// unfinished coroutine.
///
/// \param RfnB             the block that has to be refined
/// \param SpC              the splitter constellation
/// \param FromRed          the set of transitions from `RfnB` to `SpC`
/// \param postprocessing   true iff `refine()` is called during postprocessing
///                         new bottom states
/// \note In debug mode, an additional parameter `NewC` is available:
///       (only used for measuring the time complexity, and only if called from
///       line 2.26) the constellation that was used as the basis of marking
///       states.  If not called from line 2.26, pass `nullptr`.
/// \result a pointer to the block that contains the red part of `RfnB`.

template <class LTS_TYPE>
bisim_gjkw::block_t* bisim_partitioner_gjkw<LTS_TYPE>::refine(
       bisim_gjkw::block_t* const RfnB, const bisim_gjkw::constln_t* const SpC,
       const bisim_gjkw::B_to_C_descriptor* const FromRed,
       bool const postprocessing
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                  , const bisim_gjkw::constln_t* NewC
                                                                                #endif
       )
{
                                                                                #ifndef NDEBUG
                                                                                    if (nullptr != FromRed)
                                                                                    {
                                                                                        assert(FromRed->from_block() == RfnB);
                                                                                        assert(FromRed->to_constln() == SpC);
                                                                                        assert(nullptr != SpC);
                                                                                        assert(0 == RfnB->marked_size());
                                                                                    }
                                                                                    else
                                                                                    {
                                                                                        assert(0 != RfnB->marked_size());
                                                                                    }
                                                                                    assert(1 < RfnB->size());
                                                                                #endif
    // 3.2: if RfnB subseteq SpC then return RfnB
    if (RfnB->constln() == SpC)
    {
        // Mark all bottom states to indicate there are no new bottom states.
        RfnB->set_marked_nonbottom_begin(RfnB->nonbottom_end());
        RfnB->set_marked_bottom_begin(RfnB->bottom_begin());
        return RfnB;
    }
    // 3.3: Test := {bottom states}\Red, ...
    // 3.4: Spend the same amount of work on either coroutine:
    // and
    // 3.39: RedB := RfnB  or  RedB := NewB , respectively
                                                                                #ifndef NDEBUG
                                                                                    bisim_gjkw::permutation_iter_t const red_end = RfnB->end();
                                                                                #endif
    bisim_gjkw::block_t* RedB;

    // NOLINTBEGIN(cppcoreguidelines-avoid-goto,misc-static-assert) -- generated by the coroutine framework macros
    COROUTINES_SECTION
        // common variables
            bisim_gjkw::permutation_iter_t notblue_initialised_end =
                                                       RfnB->nonbottom_begin();
            bool FromRed_is_empty = nullptr == FromRed;

        // variables for the blue coroutine
            bisim_gjkw::permutation_iter_t blue_visited_end;
            bisim_gjkw::state_info_ptr blue_s;
            bisim_gjkw::pred_iter_t blue_pred_iter;
            bisim_gjkw::state_info_ptr blue_s_prime;
            bisim_gjkw::permutation_iter_t blue_blue_nonbottom_end;
            bisim_gjkw::succ_const_iter_t blue_begin;
            bisim_gjkw::succ_const_iter_t blue_end;

        // variables for the red coroutine
            bisim_gjkw::B_to_C_iter_t red_fromred_visited_begin;
            bisim_gjkw::permutation_iter_t red_visited_begin;
            bisim_gjkw::state_info_ptr red_s;
            bisim_gjkw::pred_iter_t red_pred_iter;

        // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        COROUTINE_LABELS(   (REFINE_BLUE_PREDECESSOR_HANDLED)
                            (REFINE_BLUE_TESTING)
                            (REFINE_BLUE_STATE_HANDLED)
                            (REFINE_BLUE_COLLECT_BOTTOM)
                            (REFINE_RED_COLLECT_FROMRED)
                            (REFINE_RED_PREDECESSOR_HANDLED)
                            (REFINE_RED_STATE_HANDLED))

/*--------------------------- handle blue states ----------------------------*/

/// The blue coroutine finds blue states in a block, i. e. states that cannot
/// reach the splitter.  It is one of the two coroutines that together
/// implement the fast refinement of block `RfnB` into states that can reach
/// the splitter `SpC` (red states) and those that cannot (blue states).
///
/// This coroutine assumes that `RfnB` is nontrivial, i. e. contains at least
/// two states.  `refine()` may be called in two ways:  either the initially
/// red states are marked, or a set of transitions from `RfnB` to `SpC` is
/// given.  In the first case, all unmarked bottom states are blue.  In the
/// second case, however, the first thing this coroutine has to do is to decide
/// which bottom states are blue and which ones are not.  In that case, it
/// assumes that the `current_constln` pointer of bottom states is pointing to
/// a place where a transition to `SpC` could be inserted.  If there is already
/// such a transition just before or after that memory location, the bottom
/// state is actually red; otherwise, it is blue.
///
/// After having found the blue bottom states, the coroutine proceeds to find
/// blue non-bottom states by walking through the predecessors of the blue
/// states.  As soon as can be determined that all inert outgoing transitions
/// of a state lead to (other) blue states, this state is also blue -- with one
/// exception:  if this state has a (non-inert) transition to the splitter.  If
/// `refine()` was called in the second way, one has to check the outgoing
/// non-inert transitions to find out whether such a transition exists.
///
/// As soon as it becomes clear that the blue subblock will be larger, this
/// coroutine is aborted.  Otherwise, at the end `RfnB` is actually split into
/// its two parts, and the coroutine closes with updating the inert-ness of
/// transitions and finding new bottom states.  New bottom states will be
/// unmarked bottom states.  (It may also happen that the blue subblock is
/// empty;  in that case, `RfnB` is not split.)
///
/// The coroutine implements the left-hand side of Algorithm 3 in [GJKW 2017].

        COROUTINE                                                               assert(RfnB->nonbottom_begin() == notblue_initialised_end);
            // we have to decide which unmarked bottom states are blue.  So we
            // walk over all of them and check whether they have a transition
            // to SpC or not.
            // 3.3: ..., Blue := {}
            blue_visited_end = RfnB->unmarked_bottom_begin();
            // 3.5l: whenever |Blue| > |RfnB|/2 do  Abort this coroutine
                // nothing needs to be done now, as |Blue| = 0 here.

    /*  -  -  -  -  -  -  - collect blue bottom states -  -  -  -  -  -  -  */

            // 3.6l: while Test is not empty and FromRed is not empty do
                // We use the variable blue_visited_end in this loop to
                // indicate the boundary between blue states (namely those
                // in the interval [RfnB->unmarked_bottom_begin(),
                // blue_visited_end) ) and Test states (namely those in
                // [blue_visited_end, RfnB->unmarked_bottom_end()) ).
            COROUTINE_WHILE (REFINE_BLUE_COLLECT_BOTTOM,
                                RfnB->unmarked_bottom_end() > blue_visited_end)
            {
                if (FromRed_is_empty)
                {
                    // 3.14l: Blue := Blue union Test
                        // done implicitly: we now regard all unmarked
                        // bottom states as blue, i. e. the whole
                        // interval [RfnB->unmarked_bottom_begin(),
                        // RfnB->unmarked_bottom_end()).
                    // 3.5l: whenever |Blue|>|RfnB|/2 do  Abort this coroutine
                        // In this case, Test may not yet be empty, so here we
                        // have to check the condition of Line 3.5l.
                    if (RfnB->unmarked_bottom_size() > RfnB->size() / 2)
                    {
                        ABORT_THIS_COROUTINE();
                    }
                    break;
                }
                // 3.7l: Choose s in Test
                blue_s = *blue_visited_end;
                // 3.8l: if s --> SpC then
                if (blue_s->surely_has_transition_to(SpC))
                {
                    // 3.9l: Move s from Test to Red
                    // The state s is not blue.  Move it to the slice of
                    // non-blue bottom states.
                    RfnB->set_marked_bottom_begin(
                                              RfnB->marked_bottom_begin() - 1);
                    bisim_gjkw::swap_permutation(blue_visited_end,
                                                  RfnB->marked_bottom_begin());
                    // 3.5r: whenever |Red|>|RfnB|/2 do Abort the red coroutine
                    if (RfnB->marked_size() > RfnB->size() / 2)
                    {
                        ABORT_OTHER_COROUTINE();
                    }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    // The state is red.  Depending on the context of refinement, the work done
                                                                                    // on this state is attributed to one or another element of the Kripke
                                                                                    // structure.
                                                                                    if (postprocessing)
                                                                                    {
                                                                                        // refine() has been called from line 4.14.  During postprocessing, the
                                                                                        // work is ascribed to the transitions from s to SpC.
                                                                                        // Using the current_constln pointer of s, find a transition to SpC.
                                                                                        // The work has to be ascribed to the transitions from s to SpC.
                                                                                        bisim_gjkw::succ_iter_t to_SpC = blue_s->current_constln();
                                                                                        if (to_SpC == blue_s->succ_end() || to_SpC->target->constln() != SpC)
                                                                                        {
                                                                                            assert(blue_s->succ_begin() < to_SpC);
                                                                                            --to_SpC;
                                                                                            assert(to_SpC->target->constln() == SpC);
                                                                                        }
                                                                                        bisim_gjkw::succ_entry::slice_add_work_to_transns(to_SpC,
                                                                                                 check_complexity::
                                                                                                 while_Test_is_not_empty_3_6l_s_is_red_3_9l_postprocessing, 1);
                                                                                    }
                                                                                    else
                                                                                    {
                                                                                        // refine() has been called from line 2.26, and RfnB is a block with a
                                                                                        // transition to SpB.  The state has been marked because of these
                                                                                        // transitions (which are not transitions to SpC!).
                                                                                        // code similar to the case if (postprocessing).
                                                                                        bisim_gjkw::succ_iter_t to_NewC = blue_s->current_constln();
                                                                                        if (to_NewC == blue_s->succ_end() || to_NewC->target->constln() == SpC)
                                                                                        {
                                                                                            assert(blue_s->succ_begin() < to_NewC);
                                                                                            --to_NewC;
                                                                                            assert(to_NewC->target->constln() != SpC);
                                                                                        }
                                                                                        bisim_gjkw::succ_entry::slice_add_work_to_transns(to_NewC,
                                                                                              check_complexity::while_Test_is_not_empty_3_6l_s_is_red_3_9l, 1);
                                                                                    }
                                                                                #endif
                // 3.10l: else
                }
                else
                {                                                               assert(blue_s->surely_has_no_transition_to(SpC));
                    // 3.11l: Move s from Test to Blue
                    ++blue_visited_end;
                    // 3.5l: whenever |Blue|>|RfnB|/2 do  Abort this coroutine
                    if ((state_type) (blue_visited_end -
                             RfnB->unmarked_bottom_begin()) > RfnB->size() / 2)
                    {
                        ABORT_THIS_COROUTINE();
                    }                                                           mCRL2complexity(blue_s, add_work(
                                                                                         check_complexity::while_Test_is_not_empty_3_6l_s_is_blue_3_11l, 1), );
                // 3.12l: end if
                }
            // 3.13l: end while
            }
            END_COROUTINE_WHILE;                                                assert(RfnB->constln() != SpC);

            if (0 == RfnB->unmarked_bottom_size())
            {
                // all bottom states are red, so there cannot be any blue
                // states.
                RfnB->set_marked_nonbottom_begin(RfnB->marked_nonbottom_end());
                // RfnB->set_marked_bottom_begin(RfnB->bottom_begin());
                RedB = RfnB;
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    blue_is_smaller(nullptr, RedB, NewC);
                                                                                #endif
                TERMINATE_COROUTINE_SUCCESSFULLY();
            }

    /*  -  -  -  -  -  -  -  - visit blue states -  -  -  -  -  -  -  -  */

            // 3.15l: while Blue contains unvisited states do
            blue_visited_end = RfnB->unmarked_bottom_begin();
            blue_blue_nonbottom_end = RfnB->unmarked_nonbottom_begin();
            COROUTINE_DO_WHILE (REFINE_BLUE_STATE_HANDLED,
                                   blue_visited_end != blue_blue_nonbottom_end)
            {
                // 3.16l: Choose an unvisited s in Blue
                blue_s = *blue_visited_end;                                     assert(blue_visited_end < blue_blue_nonbottom_end ||
                                                                                                       (RfnB->unmarked_bottom_begin() <= blue_visited_end &&
                                                                                                              RfnB->unmarked_bottom_end() > blue_visited_end));
                                                                                assert(RfnB->unmarked_nonbottom_begin() <= blue_blue_nonbottom_end);
                                                                                assert(blue_blue_nonbottom_end <= notblue_initialised_end);
                /* 3.17l: Mark s as visited                                  */ assert(RfnB->unmarked_nonbottom_end() >= notblue_initialised_end);
                ++blue_visited_end;
                // 3.18l: for all s_prime in inert_in(s) \ Red do
                COROUTINE_FOR (REFINE_BLUE_PREDECESSOR_HANDLED,
                    blue_pred_iter = blue_s->inert_pred_begin(),
                    blue_s->inert_pred_end()!=blue_pred_iter, ++blue_pred_iter)
                {
                    blue_s_prime = blue_pred_iter->source;
                    if (blue_s_prime->pos >= RfnB->marked_nonbottom_begin())
                    {                                                           mCRL2complexity(blue_pred_iter, add_work(
                                                                                         check_complexity::for_all_s_prime_in_pred_s_setminus_Red_3_18l, 1), );
                        continue;
                    }
                    // 3.19l: if notblue(s_prime) undefined then
                    if (blue_s_prime->pos >= notblue_initialised_end)
                    {
                        // 3.20l: notblue(s_prime) := |inert_out(s_prime)|
                        blue_s_prime->notblue = blue_s_prime->inert_succ_end()-
                                              blue_s_prime->inert_succ_begin();
                        bisim_gjkw::swap_permutation(blue_s_prime->pos,
                                                      notblue_initialised_end);
                        ++notblue_initialised_end;
                    // 3.21l: end if
                    }
                    // 3.22l: notblue(s_prime) := notblue(s_prime) - 1
                    --blue_s_prime->notblue;
                    // 3.23l: if notblue(s_prime) == 0 && ...
                    if (0 != blue_s_prime->notblue)
                    {                                                           mCRL2complexity(blue_pred_iter, add_work(
                                                                                         check_complexity::for_all_s_prime_in_pred_s_setminus_Red_3_18l, 1), );
                        continue;
                    }
                    // 3.23l: ... && (FromRed == {} ||
                    //          out_noninert(s_prime) intersect SpC == {}) then
                    if (!FromRed_is_empty)
                    {
                        if (blue_s_prime->surely_has_transition_to(SpC))
                        {                                                       mCRL2complexity(blue_pred_iter, add_work(
                                                                                         check_complexity::for_all_s_prime_in_pred_s_setminus_Red_3_18l, 1), );
                            continue;
                        }
                        if (!blue_s_prime->surely_has_no_transition_to(SpC))
                        {
                            // It is not yet known whether s_prime has a
                            // transition to SpC or not.  Execute the slow test
                            // now.
                            COROUTINE_FOR (REFINE_BLUE_TESTING,
                                (blue_begin = blue_s_prime->succ_begin(),
                                 blue_end   = blue_s_prime->succ_end()   ),
                                               blue_begin < blue_end, (void) 0)
                            {
                                // binary search for transitions from s_prime
                                // to constellation SpC.
                                bisim_gjkw::succ_const_iter_t const mid =
                                      blue_begin + (blue_end - blue_begin) / 2;
                                if (*SpC <= *mid->target->constln())
                                {
                                    blue_end = mid->slice_begin();
                                }
                                if (*mid->target->constln() <= *SpC)
                                {
                                    blue_begin =
                                        bisim_gjkw::succ_entry::slice_end(mid);
                                }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    bisim_gjkw::succ_entry::slice_add_work_to_transns(mid,
                                                                                                check_complexity::if_s_prime_has_transition_to_SpC_3_23l, 1);
                                                                                #endif
                            }
                            END_COROUTINE_FOR;
                            if (blue_begin != blue_end)
                            {                                                   mCRL2complexity(blue_pred_iter, add_work(
                                                                                         check_complexity::for_all_s_prime_in_pred_s_setminus_Red_3_18l, 1), );
                                continue;
                            }
                        }
                    }
                    // 3.24l: Blue := Blue union {s_prime}
                    bisim_gjkw::swap_permutation(blue_s_prime->pos,
                                                      blue_blue_nonbottom_end);
                    ++blue_blue_nonbottom_end;
                    // 3.5l: whenever |Blue|>|RfnB| / 2 do Abort this coroutine
                    if (blue_blue_nonbottom_end -
                               RfnB->unmarked_nonbottom_begin() +
                               RfnB->unmarked_bottom_size() > RfnB->size() / 2)
                    {
                        ABORT_THIS_COROUTINE();
                    }
                    // 3.25l: end if
                        // this is implicit in the `continue` statements above.
                /* 3.26l: end for                                            */ mCRL2complexity(blue_pred_iter, add_work(
                                                                                         check_complexity::for_all_s_prime_in_pred_s_setminus_Red_3_18l, 1), );
                }
                END_COROUTINE_FOR;
            /* 3.27l: end while                                              */ mCRL2complexity(blue_s, add_work(
                                                                                           check_complexity::while_Blue_contains_unvisited_states_3_15l, 1), );
                if (RfnB->unmarked_bottom_end() == blue_visited_end)
                {
                    blue_visited_end = RfnB->unmarked_nonbottom_begin();
                }
            }
            END_COROUTINE_DO_WHILE;

    /*  -  -  -  -  -  -  -  - split off blue block -  -  -  -  -  -  -  -  */

            // 3.28l: Abort the other coroutine
            ABORT_OTHER_COROUTINE();
            // All non-blue states are red.
            // 3.29l: Move Blue to a new block NewB
            // and
            // 3.30l: Destroy all temporary data
            bisim_gjkw::block_t* const NewB =
                                 RfnB->split_off_blue(blue_blue_nonbottom_end);
            part_tr.new_blue_block_created(RfnB, NewB);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    unsigned const max_counter = check_complexity::log_n -
                                                                                                                         check_complexity::ilog2(NewB->size());
                                                                                    mCRL2complexity(NewB,
                                                                                            add_work(check_complexity::for_all_s_in_NewB_3_31, max_counter), );
                                                                                #endif
            // 3.31l: for all s in NewB do
            for (bisim_gjkw::permutation_iter_t s_iter = NewB->begin();
                                               NewB->end() != s_iter; ++s_iter)
            {
                blue_s = *s_iter;                                               // mCRL2complexity(s, ...) -- optimized to the above call.
                // 3.32l: for all s_prime in inert_in(s) \ NewB do
                for (blue_pred_iter = blue_s->inert_pred_begin();
                    blue_s->inert_pred_end()!=blue_pred_iter; ++blue_pred_iter)
                {                                                               mCRL2complexity(blue_pred_iter, add_work(
                                                                                            check_complexity::for_all_s_prime_in_pred_s_3_32l, max_counter), );
                                                                                assert(part_tr.pred_end() > blue_pred_iter);
                                                                                assert(blue_pred_iter->succ->B_to_C->pred == blue_pred_iter);
                    blue_s_prime = blue_pred_iter->source;
                    if (blue_s_prime->block == NewB) { continue; }              assert(blue_s_prime->block == RfnB);
                    // 3.33l: s_prime --> s is no longer inert
                    part_tr.make_noninert(blue_pred_iter->succ);
                    // 3.34l: if |inert_out(s_prime)| == 0 then
                    if (blue_s_prime->inert_succ_begin() ==
                                                blue_s_prime->inert_succ_end())
                    {
                        // 3.35l: s_prime is a new bottom state
                        RfnB->set_marked_nonbottom_begin(
                                                     RfnB->bottom_begin() - 1);
                        RfnB->set_bottom_begin(RfnB->marked_nonbottom_begin());
                        swap_permutation(blue_s_prime->pos,
                                                         RfnB->bottom_begin()); // assert("blue_s_prime has a transition to SpC");
                    // 3.36l: end if
                    }
                // 3.37l: end for
                }
            // 3.38l: end for
            }

            RedB = RfnB;
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    blue_is_smaller(NewB, RedB, NewC);
                                                                                #endif
        END_COROUTINE

/*---------------------------- handle red states ----------------------------*/

/// The red coroutine finds red states in a block, i. e. states that can reach
/// the splitter.  It is one of the two coroutines that together implement
/// the fast refinement of block `RfnB` into states that can reach the splitter
/// (red states) and those that cannot (blue states).
///
/// This coroutine assumes that `RfnB` is nontrivial, i. e. contains at least
/// two states, and that there is at least one red state in `RfnB`.  The red
/// states can be indicated in one of two ways: either some states of `RfnB`
/// are marked, or the parameter `FromRed` indicates the set of transitions
/// from `RfnB` to the splitter.  Exactly one of the two must be chosen:  If
/// and only if some states are marked, `FromRed` should be `nullptr`.
///
/// If `FromRed` is not `nullptr`, the coroutine first finds initially red
/// states, i. e. states with a transition to the splitter, and marks them.
/// Note that it may happen that some non-bottom states (or even only
/// non-bottom states) may be initially red as well.  When all initially red
/// states are marked, it proceeds to find other red states by walking through
/// the inert predecessors of red states.  Every such predecessor is also a red
/// state.
///
/// As soon as it becomes clear that the red subblock will be larger, this
/// coroutine is aborted.  Otherwise, at the end `RfnB` is actually split into
/// its two parts, and the coroutine closes with updating the inert-ness of
/// transitions and finding new bottom states.  New bottom states will be
/// unmarked bottom states.
///
/// The coroutine implements the right-hand side of Algorithm 3 in [GJKW 2017].

        COROUTINE
            // 3.5r: whenever |Red| > |RfnB|/2 then  Abort this coroutine
        if (RfnB->marked_size() > RfnB->size() / 2)
        {
          ABORT_THIS_COROUTINE();
        }
            // The red block contains at most RfnB->size() -
            // RfnB->unmarked_bottom_size() + FromRed->size() states.  If that
            // is <= RfnB->size() / 2, we could abort the other coroutine
            // immediately.  We don't do it here because we want to investigate
            // the effect of this and similar heuristics more systematically.

    /*  -  -  -  -  -  -  - collect states from FromRed -  -  -  -  -  -  -  */

            // 3.6r: while FromRed != {} do
            if (nullptr != FromRed)
            {                                                                   assert(FromRed->from_block() == RfnB);
                red_fromred_visited_begin = FromRed->end;                       assert(RfnB->inert_end() <= FromRed->begin ||
                                                                                                        (RfnB->inert_begin() >= red_fromred_visited_begin &&
                                                                                                               RfnB->inert_end() > red_fromred_visited_begin));
                COROUTINE_WHILE (REFINE_RED_COLLECT_FROMRED,
                                   FromRed->begin != red_fromred_visited_begin)
                {
                    // 3.10r (order of lines changed): FromRed := FromRed \ {s --> t}
                        // We can change the order of lines because the
                        // coroutine is not interrupted in between.
                    --red_fromred_visited_begin;
                    // 3.7r: Choose s --> t in FromRed
                    bisim_gjkw::state_info_ptr const red_s =
                                       red_fromred_visited_begin->pred->source; assert(RfnB->begin() <= red_s->pos);
                    /* 3.8r: Test := Test \ {s}                              */ assert(red_s->pos < RfnB->end());
                    // and
                    // 3.9r: Red := Red union {s}
                    if (red_s->pos < notblue_initialised_end)
                    {
                        // The non-bottom state has a transition to a blue
                        // state, so notblue is initialised; however, now it is
                        // revealed to be red anyway.
                        bisim_gjkw::swap_permutation(red_s->pos,
                                                    --notblue_initialised_end);
                    }
                    if (RfnB->mark(red_s) &&
                    // 3.5r: whenever |Red| > |RfnB|/2 do  Abort this coroutine
                                        RfnB->marked_size() > RfnB->size() / 2)
                    {
                        ABORT_THIS_COROUTINE();
                    }
                /* 3.13r: end while                                          */ mCRL2complexity(red_fromred_visited_begin->pred,
                                                                                             add_work(check_complexity::while_FromRed_is_not_empty_3_6r, 1), );
                }
                END_COROUTINE_WHILE;

                // The shared variable FromRed_is_empty is set to true as soon
                // as FromRed should be considered empty.  (From that moment
                // on, no slow tests are needed any more.)
                FromRed_is_empty = true;
            }                                                                   assert(FromRed_is_empty);  assert(0 != RfnB->marked_size());

    /*  -  -  -  -  -  -  -  - visit red states -  -  -  -  -  -  -  -  */

            red_visited_begin = RfnB->marked_bottom_end();
            if (RfnB->marked_bottom_begin() == red_visited_begin)
            {
                // It may happen that all found states are non-bottom states.
                // (In that case, some of these states will become new bottom
                // states.)
                red_visited_begin = RfnB->marked_nonbottom_end();
            }
            // 3.15r: while Red contains unvisited states do
            COROUTINE_DO_WHILE(REFINE_RED_STATE_HANDLED,
                           RfnB->marked_nonbottom_begin() != red_visited_begin)
            {
                // 3.17r (order of lines changed): Mark s as visited
                --red_visited_begin;                                            assert(RfnB->marked_bottom_begin() <= red_visited_begin ||
                                                                                                     (RfnB->marked_nonbottom_begin() <= red_visited_begin &&
                                                                                                            RfnB->marked_nonbottom_end() > red_visited_begin));
                // 3.16r: Choose an unvisited s in Red
                red_s = *red_visited_begin;
                // 3.18r: for all s_prime in inert_in(s) do
                COROUTINE_FOR (REFINE_RED_PREDECESSOR_HANDLED,
                     red_pred_iter = red_s->inert_pred_begin(),
                     red_s->inert_pred_end() != red_pred_iter, ++red_pred_iter)
                {
                    bisim_gjkw::state_info_ptr const s_prime =
                                                         red_pred_iter->source;
                    // 3.24r: Red := Red union {s_prime}
                    if (s_prime->pos < notblue_initialised_end)
                    {
                        // The state has a transition to a blue state, so
                        // notblue is initialised; however, now it is revealed
                        // to be red anyway.
                        bisim_gjkw::swap_permutation(s_prime->pos,
                                                    --notblue_initialised_end);
                    }
                    if (RfnB->mark_nonbottom(s_prime) &&
                    // 3.5r: whenever |Red| > |RfnB|/2 do  Abort this coroutine
                                        RfnB->marked_size() > RfnB->size() / 2)
                    {
                        ABORT_THIS_COROUTINE();
                    }
                /* 3.26r: end for                                            */ mCRL2complexity(red_pred_iter,
                                                                                             add_work(check_complexity::for_all_s_prime_in_pred_s_3_18r, 1), );
                }
                END_COROUTINE_FOR;
            /* 3.27r: end while                                              */ mCRL2complexity(red_s, add_work(
                                                                                            check_complexity::while_Red_contains_unvisited_states_3_15r, 1), );
                if (RfnB->marked_bottom_begin() == red_visited_begin)
                {
                    red_visited_begin = RfnB->marked_nonbottom_end();
                }
            }
            END_COROUTINE_DO_WHILE;

    /*  -  -  -  -  -  -  -  - split off red block -  -  -  -  -  -  -  -  */

            // 3.28r: Abort the other coroutine
            ABORT_OTHER_COROUTINE();
            // All non-red states are blue.
            // 3.29r: Move Red to a new block RedB
            // and
            // 3.30r: Destroy all temporary data
            RedB = RfnB->split_off_red(RfnB->marked_nonbottom_begin());
            part_tr.new_red_block_created(RfnB, RedB, postprocessing);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    unsigned const max_counter = check_complexity::log_n -
                                                                                                                         check_complexity::ilog2(RedB->size());
                                                                                    mCRL2complexity(RedB,
                                                                                            add_work(check_complexity::for_all_s_in_NewB_3_31, max_counter), );
                                                                                #endif
            // 3.31r: for all non-bottom s in RedB do
                // we have to run through the states backwards because
                // otherwise, we might miss out some states.
            for (bisim_gjkw::permutation_iter_t s_iter = RedB->nonbottom_end();
                                           RedB->nonbottom_begin() != s_iter; )
            {
                --s_iter;
                red_s = *s_iter;                                                // mCRL2complexity(s, ...) -- optimized to the above call.
                // 3.32r: for all s_prime in inert_out(s) \ RedB do
                for (bisim_gjkw::succ_iter_t
                               succ_iter = red_s->inert_succ_begin();
                               red_s->inert_succ_end()!=succ_iter; ++succ_iter)
                {                                                               mCRL2complexity(succ_iter->B_to_C->pred, add_work(
                                                                                            check_complexity::for_all_s_prime_in_succ_s_3_32r, max_counter), );
                                                                                assert(part_tr.succ_end() > succ_iter);
                                                                                assert(succ_iter->B_to_C->pred->succ == succ_iter);
                    if (bisim_gjkw::state_info_ptr const s_prime =
                                                             succ_iter->target;
                                                        s_prime->block == RedB)
                    {
                        continue;
                    }
                    else
                    {
                      assert(s_prime->block == RfnB);
                    }
                    // 3.33r: s --> s_prime is no longer inert
                    part_tr.make_noninert(succ_iter);
                // 3.34r: end for
                }
                // 3.35r: if |inert_out(s)| == 0 then
                if (red_s->inert_succ_begin() == red_s->inert_succ_end())
                {
                    // 3.36r: s is a new bottom state
                    RedB->set_marked_nonbottom_begin(RedB->bottom_begin() - 1);
                    RedB->set_bottom_begin(RedB->marked_nonbottom_begin());
                    swap_permutation(red_s->pos, RedB->bottom_begin());         // assert("red_s has a transition to SpC");
                // 3.37r: end if
                }
            // 3.38r: end for
            }                                                                   assert(RfnB->end() < red_end);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    red_is_smaller(RfnB, RedB);
                                                                                #endif
        END_COROUTINE
    END_COROUTINES_SECTION
    // NOLINTEND(cppcoreguidelines-avoid-goto,misc-static-assert)


    /* 3.41 P := partition P where NewB is added and the states in NewB are  */ assert(RedB->end() == red_end);
    //                             removed from RfnB
        // done through adapting the data structures in the coroutines.
    // 3.42: return <RedB, BlueB> (with old and new bottom states separated)
        // The separation is indicated by marking the old bottom states.
        // Only RedB is returned explicitly, as most calls to Refine only use
        // RedB.
    return RedB;
}



/*=============================================================================
=              PostprocessNewBottom -- Algorithm 4 in [GJKW 2017]             =
=============================================================================*/



/// \brief function object to compare two constln_t pointers based on their
/// contents
/// \details This function object is used to create a sorted set of
/// constellations, namely those that are reachable from the block of new
/// bottom states in Algorithm 4.
class constln_ptr_less
{
  public:
    bool operator() (const bisim_gjkw::constln_t* a,
                                          const bisim_gjkw::constln_t* b) const
    {
        return *a < *b;
    }
};

using R_map_t = std::set<bisim_gjkw::constln_t*, constln_ptr_less>;

/// \brief Split a block with new bottom states as needed
/// \details The function splits RedB by checking whether all new bottom
/// states can reach the constellations that RedB can reach.
///
/// When this function starts, it assumes that the old bottom states of RedB
/// are marked and the new ones are not.  It is an error if RedB does not
/// contain any new bottom states.
/// The function implements Algorithm 4 of [GJKW 2017].  It first separates the
/// old from the new bottom states; then it walks through all constellations
/// that can be reached from the new bottom states to separate them into
/// smaller, stable subblocks.  The return value is the block containing the
/// old bottom states, resulting from the first separation.
/// \param RedB  block containing new bottom states that need to be stabilised
/// \returns the block containing the old bottom states (and every state in
///          RedB that can reach some old bottom state through inert
///          transitions)
template <class LTS_TYPE>
bisim_gjkw::block_t* bisim_partitioner_gjkw<LTS_TYPE>::postprocess_new_bottom(
                                            bisim_gjkw::block_t* RedB
                                            /* , bisim_gjkw::block_t* BlueB */)
{                                                                               assert(0 != RedB->unmarked_bottom_size());
                                                                                assert(RedB->marked_nonbottom_begin() == RedB->marked_nonbottom_end());
    /*------- collect constellations reachable from new bottom states -------*/

    // 4.3: <ResultB, RfnB> := Refine(RedB, cosplit(RedB, BlueB), {old bottom
    //                                                     states in RedB}, {})
    bisim_gjkw::block_t* ResultB;
    bisim_gjkw::block_t* RfnB;
                                                                                #ifndef NDEBUG
                                                                                    bisim_gjkw::permutation_const_iter_t blue_begin = RedB->begin();
                                                                                #endif
    if (0 != RedB->marked_bottom_size())
    {
        // postprocessing == false is required here because we need to keep
        // ResultB->FromRed().
        ResultB = refine(RedB, nullptr, nullptr, false);                        assert(nullptr != ResultB);  assert(0 == ResultB->unmarked_bottom_size());
        // 4.26: Destroy all temporary data
            // As part of this line, we unmark all states.
        ResultB->set_marked_bottom_begin(ResultB->bottom_end());                assert(ResultB->begin() > blue_begin);
        RfnB = ResultB->begin()[-1]->block;
    }
    else
    {
        RfnB = RedB;
        ResultB = nullptr;
    }                                                                           assert(ResultB != RfnB);  assert(RfnB->begin() == blue_begin);
    // do not refine a trivial block
    if (1 == RfnB->size())
    {                                                                           assert(0 == RfnB->marked_size());
        return ResultB;
    }
    // 4.2: Create an empty search tree R of constellations
    R_map_t R;
Line_4_4:
    /* 4.4: for all constellations C not in R reachable from RfnB do         */ assert(0 == RfnB->marked_bottom_size());
    while (!RfnB->to_constln.front().needs_postprocessing())
    {
        bisim_gjkw::B_to_C_desc_iter_t const new_slice =
                                                      RfnB->to_constln.begin(); // try to assign the work to the transitions from bottom states in RfnB to C.
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    if (!new_slice->add_work_to_bottom_transns(check_complexity::
                                                                                                       for_all_transitions_from_bottom_states_a_priori_4_4, 1))
                                                                                    {
                                                                                        // It didn't work out -- now assign it temporarily to the
                                                                                        // B_to_C slice itself.
                                                                                        mCRL2complexity(new_slice, add_work(check_complexity::
                                                                                                        for_all_constellations_C_not_in_R_from_RfnB_4_4, 1), );
                                                                                    }
                                                                                #endif
        bisim_gjkw::constln_t* const C = new_slice->to_constln();
        // 4.5: Add C to R
        if (!R.insert(C).second) //< complexity log(n)
        {                                                                       assert(0 && "The constellation already was in R");
        }
        // 4.6: Register that the transitions from RfnB to C need
        /*      postprocessing                                               */ assert(C->postprocess_begin == C->postprocess_end);
        C->postprocess_begin = new_slice->begin;
        C->postprocess_end = new_slice->end;                                    assert(C->postprocess_begin < C->postprocess_end);
        // move the set of transitions to the end of the list
        RfnB->to_constln.splice(RfnB->to_constln.end(), RfnB->to_constln,
                                                                    new_slice);
    // 4.7: end for
    }
    // 4.8: for all bottom states s in RfnB do
    bisim_gjkw::permutation_iter_t s_iter = RfnB->unmarked_bottom_begin();
    do
    {
        bisim_gjkw::state_info_ptr const s = *s_iter;                           mCRL2complexity(s,
                                                                                         add_work(check_complexity::for_all_bottom_states_s_in_RfnB_4_8, 1), );
        // 4.9: Set the current constellation pointer of s to the first
        //      constellation it can reach
        s->set_current_constln(s->succ_begin());
    // 4.10: end for
    }
    while (RfnB->unmarked_bottom_end() != ++s_iter);

    /*--------------- stabilise w. r. t. found constellations ---------------*/

    // 4.11: for all constellations SpC in R (in order) do
    while (!R.empty())
    {                                                                           // the work has to be assigned to the transitions from (new) bottom
        bisim_gjkw::constln_t* const SpC = *R.begin();                          // states to constellation SpC.  We do that assignment in the inner
        // 4.12: for all blocks B with transitions to SpC that need             // loop.  Ensure that the inner loop is executed:
        /*               postprocessing do                                   */ assert(SpC->postprocess_begin != SpC->postprocess_end);
        do
        {
            // add_work is called indirectly: after refining, work is assigned
            // to every (new) bottom state in the refined block.
            bisim_gjkw::B_to_C_iter_t const B_iter = SpC->postprocess_begin;    assert(part_tr.B_to_C_end() > B_iter);
            bisim_gjkw::block_t* const B = B_iter->pred->source->block;         assert(B_iter->pred->succ->B_to_C == B_iter);
            bisim_gjkw::B_to_C_desc_iter_t FromRed = B_iter->B_to_C_slice;      assert(FromRed->begin == B_iter);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    bool postproc_a_posteriori = !FromRed->
                                                                                        add_work_to_bottom_transns(check_complexity::
                                                                                                      for_all_transitions_that_need_postproc_a_priori_4_12, 1);
                                                                                #endif
            // 4.13: Delete the transitions from B to SpC from those that need
            //       postprocessing
            SpC->postprocess_begin = FromRed->end;                              assert(SpC->postprocess_begin <= SpC->postprocess_end);
                                                                                assert(!B->to_constln.empty());
            if (B->to_constln.begin() != FromRed)
            {                                                                   assert(FromRed->from_block() == B);
                B->to_constln.splice(B->to_constln.begin(), B->to_constln,
                                                                      FromRed);
            }
            /* do not refine a trivial block                                 */ assert(0 == B->marked_size());
            if (1 == B->size())
            {                                                                   assert(!postproc_a_posteriori);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    if (0 != FromRed->work_counter.get_work_counter_4_4())
                                                                                    {
                                                                                        if (!FromRed->add_work_to_bottom_transns(check_complexity::
                                                                                                   for_all_transitions_from_bottom_states_a_posteriori_4_4, 1))
                                                                                        {  assert(0);  }
                                                                                        FromRed->work_counter.reset_work_counter_4_4();
                                                                                    }
                                                                                #endif
                continue;
            }
            // refine() may destroy FromRed:  if the red subblock is smaller,
            // the transitions in FromRed are moved to a new B_to_C slice.
            // Therefore we select some transition from it, to restore the
            // slice afterwards.  Any transition should do.
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    bisim_gjkw::pred_iter_t FromRed_save = FromRed->begin->pred;
                                                                                #endif
            // 4.14: <RedB, BlueB> := Refine(B, SpC, {}, {transitions to SpC})
            RedB = refine(B, SpC, &*FromRed, true);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    FromRed = FromRed_save->succ->B_to_C->B_to_C_slice;
                                                                                    if (postproc_a_posteriori && !FromRed->add_work_to_bottom_transns(
                                                                                                  check_complexity::
                                                                                                  for_all_transitions_that_need_postproc_a_posteriori_4_12, 1))
                                                                                    {
                                                                                        assert(0 && "no new bottom states to assign work to");
                                                                                    }
                                                                                    // try to move the work that was not yet moved to the
                                                                                    // transitions from (new) bottom states now.
                                                                                    if (0 != FromRed->work_counter.get_work_counter_4_4())
                                                                                    {
                                                                                        if (!FromRed->add_work_to_bottom_transns(check_complexity::
                                                                                                   for_all_transitions_from_bottom_states_a_posteriori_4_4, 1))
                                                                                        {
                                                                                            assert(0 && "no new bottom states to assign work to");
                                                                                        }
                                                                                        FromRed->work_counter.reset_work_counter_4_4();
                                                                                    }
                                                                                #endif
            // 4.15: for all old bottom states s in RedB do
            for (bisim_gjkw::permutation_iter_t
                                 s_iter = RedB->marked_bottom_begin();
                                 RedB->marked_bottom_end() != s_iter; ++s_iter)
            {
                bisim_gjkw::state_info_ptr const s = *s_iter;
                // 4.16: Advance the current constellation pointer of s to the
                //       next constellation it can reach
                if (s->current_constln() < s->succ_end() &&
                                s->current_constln()->target->constln() == SpC)
                {
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    bisim_gjkw::succ_entry::slice_add_work_to_transns(s->current_constln(),
                                                                                                 check_complexity::for_all_old_bottom_states_s_in_RedB_4_15,1);
                                                                                #endif
                    s->set_current_constln(bisim_gjkw::succ_entry::slice_end(
                                                        s->current_constln()));
                }                                                               else
                                                                                {   mCRL2complexity(s, add_work(check_complexity::
                                                                                                      for_all_old_bottom_states_s_in_RedB_selfloop_4_15, 1), );
                                                                                    assert(B == RedB /* we should test: BlueB is empty */);
                                                                                    assert(B->constln() == SpC);
                                                                                }
                                                                                assert(s->succ_end() == s->current_constln() ||
            /* 4.17: end for                                                 */                               *SpC < *s->current_constln()->target->constln());
            }
            /* 4.18: if RedB contains new bottom states then                 */ assert(RedB->marked_nonbottom_begin() == RedB->marked_nonbottom_end());
            if (0 != RedB->unmarked_bottom_size())
            {
                // 4.19: <_, RfnB> := Refine(RedB, cosplit(RedB, BlueB), {old
                //                           bottom states in RedB}, {})
                                                                                #ifndef NDEBUG
                                                                                    blue_begin = RedB->begin();
                                                                                #endif
                if (0 != RedB->marked_bottom_size())
                {
                    RedB = refine(RedB, nullptr, nullptr, true);                assert(nullptr != RedB);  assert(0 == RedB->unmarked_bottom_size());
                    // 4.26: Destroy all temporary data
                        // As part of this line, we unmark all states.
                    RedB->set_marked_bottom_begin(RedB->bottom_end());          assert(RedB->begin() > blue_begin);
                    RfnB = RedB->begin()[-1]->block;
                }
                else
                {
                    RfnB = RedB;
                                                                                #ifndef NDEBUG
                                                                                    RedB = nullptr;
                                                                                #endif
                }                                                               assert(RedB != RfnB);
                                                                                assert(RfnB->begin() == blue_begin);
                // do not refine a trivial block
                if (1 == RfnB->size())
                {                                                               assert(0 == RfnB->marked_size());
                    continue;
                }
                // 4.20: Register that the transitions from RfnB to SpC need
                //       postprocessing
                bisim_gjkw::B_to_C_desc_iter_t new_slice =
                                       SpC->postprocess_begin[-1].B_to_C_slice; assert(new_slice->from_block() == RfnB);
                                                                                assert(new_slice->to_constln() == SpC);
                SpC->postprocess_begin = new_slice->begin;                      assert(new_slice->end <= SpC->postprocess_end);
                // move the set of transitions to the end of the list
                RfnB->to_constln.splice(RfnB->to_constln.end(),
                                                  RfnB->to_constln, new_slice);
                // 4.21: Restart the procedure (but keep R),
                //       i. e. go to Line 4.4
                goto Line_4_4;
            // 4.22: end if
            }
            // 4.26: Destroy all temporary data
                // As part of this line, we unmark all states.
            RedB->set_marked_bottom_begin(RedB->bottom_end());
        // 4.23: end for
        }
        while (SpC->postprocess_begin != SpC->postprocess_end);
        // 4.24: Delete SpC from R
        R.erase(R.begin());
    // 4.25: end for
    }
    /* 4.26: Destroy all temporary data                                      */ assert(R.empty());

    return ResultB;
}



/*=============================================================================
=                       explicit instantiation requests                       =
=============================================================================*/



namespace bisim_gjkw
{

template class bisim_partitioner_gjkw_initialise_helper<lts_lts_t>;
template class bisim_partitioner_gjkw_initialise_helper<lts_aut_t>;
template class bisim_partitioner_gjkw_initialise_helper<lts_fsm_t>;

} // end namespace bisim_gjkw

template class bisim_partitioner_gjkw<lts_lts_t>;
template class bisim_partitioner_gjkw<lts_aut_t>;
template class bisim_partitioner_gjkw<lts_fsm_t>;

} // end namespace detail
// end namespace lts
// end namespace mcrl2
