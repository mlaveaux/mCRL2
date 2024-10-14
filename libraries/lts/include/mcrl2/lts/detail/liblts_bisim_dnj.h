// Author(s): David N. Jansen, Institute of Software, Chinese Academy of
// Sciences, Beijing, PR China
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/// \file lts/detail/liblts_bisim_dnj.h
///
/// \brief O(m log n)-time branching bisimulation algorithm
///
/// \details This file implements an efficient partition refinement algorithm
/// for labelled transition systems inspired by Groote / Jansen / Keiren / Wijs
/// (2017) and Valmari (2009) to calculate the branching bisimulation classes
/// of a labelled transition system.  Different from the 2017 article, it does
/// not translate the LTS to a Kripke structure, but works on the LTS directly.
/// In this way the memory use can be reduced.  The algorithm is described in
/// the conference publication
///
/// David N. Jansen, Jan Friso Groote, Jeroen J.A. Keiren, Anton Wijs:
/// An O(m log n) algorithm for branching bisimilarity on labelled transition
/// systems.  In: Armin Biere, David Parker (Eds.): Tools and Algorithms for the
/// Construction and Analysis of Systems: TACAS.  (LNCS 12070).  Springer: Cham,
/// 2020.  pp 3-20.  https://doi.org/10.1007/978-3-030-45237-7_1
///
/// Partition refinement means that the algorithm maintains a partition of the
/// state space of the LTS into ``blocks''.  A block contains all states in one
/// or several branching bisimulation equivalence classes.  Blocks are being
/// refined until every block contains exactly one branching bisimulation
/// equivalence class.
///
/// The algorithm divides the non-inert transitions into *action_block-slices.*
/// An action_block-slice contains all transitions with a given action label
/// into one block.  One or more action_block-slices are joined into a *bunch.*
/// Bunches register which blocks have already been stabilised:
///
/// > Invariant:  The blocks are stable under the bunches, i. e. if a block
/// > has a transition in a bunch, then every bottom state in this block has
/// > a transition in this bunch.
///
/// However, if a bunch is non-trivial (i. e., it contains more than one
/// action_block-slice), its slices may split some blocks into finer parts:
/// not all states may have a transition in the same action_block-slice.  So, a
/// refinement step consists of moving one small action_block-slice from a
/// non-trivial bunch to a new bunch and then splitting blocks to restore the
/// Invariant.  Splitting is always done by finding a small subblock and moving
/// the states in that subblock to a new block.  Note that we always handle a
/// small part of some larger unit; that ensures that each state and transition
/// is only touched O(log n) times.
///
/// After splitting blocks, some inert transitions may have become non-inert.
/// These transitions mostly need to be inserted into their own bunch.  Also,
/// some states may have lost all their inert transitions and become new bottom
/// states.  These states need to be handled as well to re-establish the
/// Invariant.
///
/// Overall, we spend time as follows:
/// - Every  transition  is  moved  to  a  new  bunch  at  most
///   log<SUB>2</SUB>(2 * n<SUP>2</SUP>) = 2*log<SUB>2</SUB>(n) + 1 times.
///   Every move leads to O(1) work.
/// - Every state is moved to a new block at most log<SUB>2</SUB>(n) times.
///   Every move leads to work proportional to the number of its incoming and
///   outgoing transitions.
/// - Every state becomes a new bottom state at most once.  When this happens,
///   this  leads  to  work  proportional  to  the  number  of  its  outgoing
///   transitions.
/// When summing this up over all states and transitions, we get that the
/// algorithm spends time in O(m log n), where m is the number of transitions
/// and n ist the number of states.
///
/// \author David N. Jansen, Institute of Software, Chinese Academy of
/// Sciences, Beijing, PR China

// The file is best viewed on a screen or in a window that is 160 characters
// wide.  The left 80 columns contain the main text of the program.  The right
// 80 columns contain assertions and other code used for debugging.

#ifndef LIBLTS_BISIM_DNJ_H
#define LIBLTS_BISIM_DNJ_H

#include <iomanip> // for std::fixed, std::setprecision(), std::setw()
#include <ctime> // for std::clock_t, std::clock()
#include "mcrl2/lts/detail/liblts_scc.h"
#include "mcrl2/lts/detail/liblts_merge.h"
#include "mcrl2/lts/detail/coroutine.h"
#include "mcrl2/lts/detail/check_complexity.h"
#include "mcrl2/lts/detail/fixed_vector.h"
#include "mcrl2/lts/detail/simple_list.h"

#include <cstddef>   // for std::size_t

namespace mcrl2
{
namespace lts
{
namespace detail
{
                                                                                #ifndef NDEBUG
                                                                                    /// \brief include something in Debug mode
                                                                                    /// \details In a few places, we have to include an additional parameter to
                                                                                    /// a function in Debug mode.  While it is in principle possible to use
                                                                                    /// #ifndef NDEBUG ... #endif, that would lead to distributing the code
                                                                                    /// over many code lines.  This macro expands to its arguments in Debug
// state_type and trans_type are defined in check_complexity.h.                     /// mode and to nothing otherwise.
                                                                                    #define ONLY_IF_DEBUG(...) __VA_ARGS__
                                                                                #else
                                                                                    #define ONLY_IF_DEBUG(...)
                                                                                #endif
/// \brief type used to store label numbers and counts
/// \details It would be better to define it as LTS_TYPE::labels_size_type, but
/// that would require most classes to become templates.
typedef std::size_t label_type;

template <class LTS_TYPE> class bisim_partitioner_dnj;

namespace bisim_dnj
{

/// \class iterator_or_counter
/// \brief union: either iterator or counter (for initialisation)
/// \details During initialisation, we need some counters to count the number
/// of transitions for each state.  To avoid allocating separate memory for
/// these counters, we store their value in the same place where we later
/// store an iterator.
///
/// We assume that each such variable starts out as a counter and at
/// some point becomes an iterator.  That point is marked by calling
/// convert_to_iterator().  The structure may only be destroyed after that
/// call, as the destructor assumes it's an iterator.
template <class Iterator>
union iterator_or_counter
{
    /// \brief counter (used during initialisation)
    trans_type count;

    /// \brief iterator (used during main part of the algorithm)
    Iterator begin;

    /// \brief Construct the object as a counter
    iterator_or_counter()  {  count = 0;  }


    /// \brief Convert the object from counter to iterator
    void convert_to_iterator(const Iterator other)
    {
        new (&begin) Iterator(other);
    }


    /// \brief Destruct the object as an iterator
    ~iterator_or_counter()  {  begin.~Iterator();  }
};


class block_bunch_entry;
class action_block_entry;





/* ************************************************************************* */
/*                                                                           */
/*                   R E F I N A B L E   P A R T I T I O N                   */
/*                                                                           */
/* ************************************************************************* */





/// \defgroup part_state
/// \brief data structures for states
/// \details States are stored in a refinable partition data structure.  The
/// actual state information will not be moved around, but only entries in
/// a separate permutation array.  Entries is this array are grouped per
/// block, so that the states in a block can be described as a slice in the
/// permutation array.
///@{
class state_info_entry;
class permutation_entry;

/// \class permutation_t
/// \brief stores a permutation of the states, ordered by block
/// \details This is the central concept of the _refinable partition_: the
/// permutation of the states, such that states belonging to the same block are
/// adjacent.
///
/// Iterating over the states of a block will
/// therefore be done using the permutation_t array.
typedef fixed_vector<permutation_entry> permutation_t;

class block_t;
class bunch_t;

class pred_entry;
class succ_entry;

class block_bunch_slice_t;
typedef simple_list<block_bunch_slice_t>::iterator block_bunch_slice_iter_t;
typedef simple_list<block_bunch_slice_t>::const_iterator
                                                block_bunch_slice_const_iter_t;
typedef iterator_or_null_t<block_bunch_slice_t>
                                              block_bunch_slice_iter_or_null_t;

enum new_block_mode_t { new_block_is_U, new_block_is_R };


/// \class state_info_entry
/// \brief stores information about a single state
/// \details This class stores all other information about a state that the
/// partition needs.  In particular:  the block where the state belongs and the
/// position in the permutation array (i. e. the inverse of the permutation).
class state_info_entry
{
  public:
    /// \brief iterator to first inert incoming transition
    /// \details Non-inert incoming transitions of the state are stored just
    /// before the element where this iterator points to.
    ///
    /// During initialisation, this field also doubles up as a counter for
    /// the number of incoming transitions, and as the pointer to the first
    /// incoming inert transition that already has been initialised.
    iterator_or_counter<pred_entry*> pred_inert;

    /// \brief iterator to first inert outgoing transition
    /// \details Non-inert outgoing transitions of the state are stored just
    /// before the element where this iterator points to.
    ///
    /// During initialisation, this field also doubles up as a counter for the
    /// number of *inert* outgoing transitions, and as the pointer to the first
    /// outgoing inert transition that already has been initialised.
    iterator_or_counter<succ_entry*> succ_inert;

    /// \brief block where the state belongs
    /// \details During initialisation, this field is used to point at the
    /// first unused slot of the (non-inert) bledecessors, ahem, predecessors.
    /// Sorry for the mock-chinese ``typo''.  So we always assume that it
    /// starts as a pred_entry*, at some moment is converted to a block_t*, and
    /// then stays that way until it is destroyed.
    union bl_t {
        pred_entry* ed_noninert_end;
        block_t* ock;
    } bl;

    /// \brief position of the state in the permutation array
    permutation_entry* pos;

    /// \brief number of inert transitions to non-U-states
    /// \details Actually, as we cannot easily count how many inert outgoing
    /// transitions this state has, we initialize this pointer to
    /// succ_inert.begin.  Every time we find an outgoing transition to a
    /// U-state, we increase this iterator; as soon as it no longer points to
    /// an outgoing transition of this state, we have found all inert outgoing
    /// transitions.  This requires that after the inert outgoing transitions
    /// there is a transition that starts in another state, or there is a dummy
    /// transition.
    ///
    /// During initialisation, this field also doubles up as a counter for the
    /// number of *non-inert* outgoing transitions, and as the pointer to the
    /// first outgoing transition that already has been initialised.  Therefore
    /// it cannot be a `const succ_entry*`.
    iterator_or_counter<succ_entry*> untested_to_U_eqv;
                                                                                #ifndef NDEBUG
                                                                                    /// \brief print a short state identification for debugging
                                                                                    template<class LTS_TYPE>
                                                                                    std::string debug_id_short(const bisim_partitioner_dnj<LTS_TYPE>&
                                                                                                                                             partitioner) const
                                                                                    {   assert(partitioner.part_st.state_info.data() <= this);
                                                                                        assert(this < partitioner.part_st.state_info.data_end());
                                                                                        return std::to_string(this - partitioner.part_st.state_info.data());
                                                                                    }

                                                                                    /// \brief print a state identification for debugging
                                                                                    template<class LTS_TYPE>
                                                                                    std::string debug_id(const bisim_partitioner_dnj<LTS_TYPE>&
                                                                                                                                             partitioner) const
                                                                                    {
                                                                                        return "state " + debug_id_short(partitioner);
                                                                                    }
                                                                                #endif
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    mutable check_complexity::state_dnj_counter_t work_counter;
                                                                                #endif
};


/// \brief entry in the permutation array
class permutation_entry {
  public:
    /// \brief pointer to the state information data structure
    state_info_entry* st;


    /// \brief default constructor (should not be deleted)
    permutation_entry() = default;


    /// \brief move constructor
    /// \details The move constructor is called when a temporary object is
    /// created; in that case, it is not necessary to set the pos pointer.
    permutation_entry(const permutation_entry&& other) noexcept
    {
        st = other.st;
    }


    /// \brief move assignment operator
    /// \details The move assignment operator is called when an object is moved
    /// to its final place. Therefore, we have to adapt the pos pointer.  Note
    /// that std::swap also uses move assignment, so we automatically get the
    /// correct behaviour there.
    void operator=(const permutation_entry&& other) noexcept
    {
        st = other.st;
        st->pos = this;
    }
};


/// \class block_t
/// \brief stores information about a block
/// \details A block corresponds to a slice of the permutation array.  As the
/// number of blocks is initially unknown, we will allocate instances of this
/// class dynamically.
///
/// The slice in the permutation array containing the states of the block is
/// subdivided into the following subslices (in this order):
/// 1. unmarked bottom states
/// 2. marked bottom states (initially empty)
/// 3. unmarked non-bottom states
/// 4. marked non-bottom states (initially empty)
///
/// A state should be marked iff it is a predecessor of the current splitter
/// (through a strong transition).  The marking is later extended to the
/// R-states;  that are the states with a weak transition to the splitter.
///
/// (During the execution of some functions, blocks are subdivided further;
/// however, as these subdivisions are local to a single function, they are not
/// stored here.)
///
/// Note that block_t uses the default destructor; therefore, if simple_list is
/// trivially destructible, so is block_t.
class block_t
{
  public:
    /// \brief iterator to the first state of the block
    permutation_entry* begin;

    /// \brief iterator to the first marked bottom state of the block
    permutation_entry* marked_bottom_begin;

    /// \brief iterator to the first non-bottom state of the block
    permutation_entry* nonbottom_begin;

    /// \brief iterator to the first marked non-bottom state of the block
    permutation_entry* marked_nonbottom_begin;

    /// \brief iterator past the last state of the block
    permutation_entry* end;

    /// \brief list of stable block_bunch-slices with transitions from this
    /// block
    simple_list<block_bunch_slice_t> stable_block_bunch;

    /// \brief unique sequence number of this block
    /// \details After the stuttering equivalence algorithm has terminated,
    /// this number is used as a state number in the quotient LTS.
    const state_type seqnr;


    /// \brief constructor
    /// \details The constructor initialises the block to: all states are
    /// bottom states, no state is marked.
    /// \param new_begin   initial iterator to the first state of the block
    /// \param new_end     initial iterator past the last state of the block
    /// \param new_seqnr   is the sequence number of the new block
    block_t(permutation_entry* const new_begin,
                  permutation_entry* const new_end, state_type const new_seqnr)
      : begin(new_begin),
        marked_bottom_begin(new_end),
        nonbottom_begin(new_end),
        marked_nonbottom_begin(new_end),
        end(new_end),
        stable_block_bunch(),
        seqnr(new_seqnr)
    {                                                                           assert(new_begin < new_end);
    }


    /// \brief provides the number of states in the block
    state_type size() const
    {                                                                           assert(begin <= marked_bottom_begin);  assert(marked_nonbottom_begin <= end);
                                                                                assert(marked_bottom_begin <= nonbottom_begin);
                                                                                assert(nonbottom_begin <= marked_nonbottom_begin);
        return end - begin;
    }


    /// \brief provides the number of bottom states in the block
    state_type bottom_size() const
    {                                                                           assert(begin <= marked_bottom_begin);  assert(marked_nonbottom_begin <= end);
                                                                                assert(marked_bottom_begin <= nonbottom_begin);
                                                                                assert(nonbottom_begin <= marked_nonbottom_begin);
        return nonbottom_begin - begin;
    }


    /// \brief provides the number of marked bottom states in the block
    state_type marked_bottom_size() const
    {                                                                           assert(begin <= marked_bottom_begin);  assert(marked_nonbottom_begin <= end);
                                                                                assert(marked_bottom_begin <= nonbottom_begin);
                                                                                assert(nonbottom_begin <= marked_nonbottom_begin);
        return nonbottom_begin - marked_bottom_begin;
    }


    /// \brief provides the number of marked states in the block
    state_type marked_size() const
    {
        return end - marked_nonbottom_begin + marked_bottom_size();
    }


    /// \brief provides the number of unmarked bottom states in the block
    state_type unmarked_bottom_size() const
    {                                                                           assert(begin <= marked_bottom_begin);  assert(marked_nonbottom_begin <= end);
                                                                                assert(marked_bottom_begin <= nonbottom_begin);
                                                                                assert(nonbottom_begin <= marked_nonbottom_begin);
        return marked_bottom_begin - begin;
    }


    /// \brief provides the number of unmarked nonbottom states in the block
    state_type unmarked_nonbottom_size() const
    {                                                                           assert(begin <= marked_bottom_begin);  assert(marked_nonbottom_begin <= end);
                                                                                assert(marked_bottom_begin <= nonbottom_begin);
                                                                                assert(nonbottom_begin <= marked_nonbottom_begin);
        return marked_nonbottom_begin - nonbottom_begin;
    }


    /// \brief mark a non-bottom state
    /// \details Marking is done by moving the state to the slice of the marked
    /// non-bottom states of the block.
    /// \param s the non-bottom state that has to be marked
    /// \returns true iff the state was not marked before
    bool mark_nonbottom(permutation_entry* const s)
    {                                                                           assert(nonbottom_begin <= s);  assert(s < end);
                                                                                // assert(this == s->st->bl.ock); -- does not hold during initialisation
                                                                                assert(begin <= marked_bottom_begin);
                                                                                assert(marked_bottom_begin <= nonbottom_begin);
                                                                                assert(nonbottom_begin <= marked_nonbottom_begin);
        if (marked_nonbottom_begin <= s) { return false; }                      assert(marked_nonbottom_begin <= end);
        std::swap(*s, *--marked_nonbottom_begin);                               assert(nonbottom_begin <= marked_nonbottom_begin);
        return true;
    }


    /// \brief mark a state
    /// \details Marking is done by moving the state to the slice of the marked
    /// bottom or non-bottom states of the block.
    /// \param s the state that has to be marked
    /// \returns true iff the state was not marked before
    bool mark(permutation_entry* const s)
    {                                                                           assert(begin <= s);
        if (s < nonbottom_begin)                                                // assert(this == s->st->bl.ock); -- does not hold during initialisation
        {                                                                       assert(begin <= marked_bottom_begin);  assert(marked_nonbottom_begin <= end);
                                                                                assert(nonbottom_begin <= marked_nonbottom_begin);
            if (marked_bottom_begin <= s) { return false; }                     assert(marked_bottom_begin <= nonbottom_begin);
            std::swap(*s, *--marked_bottom_begin);                              assert(begin <= marked_bottom_begin);
            return true;
        }
        return mark_nonbottom(s);
    }


    /// \brief refine a block
    /// \details This function is called after a refinement function has found
    /// where to split the block into unmarked (U) and marked (R) states.
    /// It creates a new block for the smaller subblock.
    /// \param   new_block_mode indicates whether the U- or the R-block
    ///          should be the new one.  (This parameter is necessary in case
    ///          the two halves have exactly the same size.)
    /// \param   new_seqnr is the sequence number of the new block
    /// \param   new_block (if the pool allocator is used) a pointer to an
    ///          uninitialized block, where the new block will be stored.
    /// \returns pointer to the new block
                                                                                ONLY_IF_DEBUG( template<class LTS_TYPE> )
    block_t* split_off_block(enum new_block_mode_t new_block_mode,              ONLY_IF_DEBUG( const bisim_partitioner_dnj<LTS_TYPE>& partitioner, )
                                                         state_type new_seqnr);
                                                                                #ifndef NDEBUG
                                                                                    /// \brief print a block identification for debugging
                                                                                    template<class LTS_TYPE>
                                                                                    inline std::string debug_id(const bisim_partitioner_dnj<LTS_TYPE>&
                                                                                                                                             partitioner) const
                                                                                    {   assert(partitioner.part_st.permutation.data() <= begin);
                                                                                        assert(begin < end);  assert(begin <= marked_bottom_begin);
                                                                                        assert(marked_bottom_begin <= nonbottom_begin);
                                                                                        assert(nonbottom_begin <= marked_nonbottom_begin);
                                                                                        assert(marked_nonbottom_begin <= end);
                                                                                        assert(end <= partitioner.part_st.permutation.data_end());
                                                                                        return "block [" +
                                                                                           std::to_string(begin - partitioner.part_st.permutation.data()) +
                                                                                           "," + std::to_string(end - partitioner.part_st.permutation.data()) +
                                                                                                                          ") (#" + std::to_string(seqnr) + ")";
                                                                                    }
                                                                                #endif
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    mutable check_complexity::block_dnj_counter_t work_counter;
                                                                                #endif
};


/// \class part_state_t
/// \brief refinable partition data structure
/// \details This class collects all information about a partition of the
/// states.
class part_state_t
{
  public:
    /// \brief permutation array
    /// \details This is the central element of the data structure:  In this
    /// array, states that belong to the same block are stored in adjacent
    /// elements.
    permutation_t permutation;

    /// \brief array with all other information about states
    fixed_vector<state_info_entry> state_info;

    /// \brief total number of blocks with unique sequence number allocated
    /// \details Upon starting the stuttering equivalence algorithm, the number
    /// of blocks must be zero.
    state_type nr_of_blocks;

    /// \brief constructor
    /// \details The constructor allocates memory and makes the permutation and
    /// state_info arrays consistent with each other, but does not actually
    /// initialise the partition.  Immediately afterwards, the initialisation
    /// will be done in `bisim_partitioner_dnj::create_initial_partition()`.
    /// \param num_states number of states in the LTS
    part_state_t(state_type const num_states)
      : permutation(num_states),
        state_info(num_states),
        nr_of_blocks(0)
    {                                                                           assert(0 < num_states);
        permutation_entry* perm_iter(permutation.data());
                                                                                #ifdef USE_POOL_ALLOCATOR
                                                                                    static_assert(std::is_trivially_destructible<block_t>::value);
                                                                                #endif
        state_info_entry* state_iter(state_info.data());                        assert(perm_iter < permutation.data_end());
        do
        {
            state_iter->pos = perm_iter;
            perm_iter->st = state_iter;
        }
        while (++state_iter, ++perm_iter < permutation.data_end());             assert(state_iter == state_info.data_end());
    }


    #ifndef USE_POOL_ALLOCATOR
        /// \brief destructor
        /// \details The destructor also deallocates the blocks, as they are
        /// not directly referenced from anywhere.  This is only necessary if
        /// we do not use the pool allocator, as the latter will destroy the
        /// blocks wholesale.
        ~part_state_t()
        {                                                                       ONLY_IF_DEBUG( state_type deleted_blocks(0); )
            permutation_entry* perm_iter(permutation.data_end());               assert(permutation.data() < perm_iter);
            do
            {
                block_t* const B(perm_iter[-1].st->bl.ock);                     assert(B->end == perm_iter);
                perm_iter = B->begin;                                           ONLY_IF_DEBUG( ++deleted_blocks; )
                delete B;
            }
            while (permutation.data() < perm_iter);                             assert(deleted_blocks == nr_of_blocks);
        }
    #endif


    /// \brief calculate the size of the state space
    /// \returns the number of states in the LTS
    state_type state_size() const  {  return permutation.size();  }


    /// \brief find the block of a state
    /// \param s number of the state
    /// \returns a pointer to the block where state s resides in
    const block_t* block(state_type const s) const
    {
        return state_info[s].bl.ock;
    }
                                                                                #ifndef NDEBUG
                                                                                  private:
                                                                                    /// \brief print a slice of the partition (typically a block)
                                                                                    /// \details If the slice indicated by the parameters is not empty, the
                                                                                    /// message and the states in this slice will be printed.
                                                                                    /// \param message      text printed as a title if the slice is not empty
                                                                                    /// \param begin_print  iterator to the beginning of the slice
                                                                                    /// \param end_print    iterator past the end of the slice
                                                                                    /// \param partitioner  LTS partitioner (used to print more details)
                                                                                    template<class LTS_TYPE>
                                                                                    void print_block(const block_t* const B,
                                                                                                     const char* const message,
                                                                                                     const permutation_entry* begin_print,
                                                                                                     const permutation_entry* const end_print,
                                                                                                     const bisim_partitioner_dnj<LTS_TYPE>& partitioner) const
                                                                                    {   assert(B->begin <= begin_print);  assert(end_print <= B->end);
                                                                                        if (end_print == begin_print)  return;

                                                                                        mCRL2log(log::debug) << '\t' << message
                                                                                                             << (1 < end_print - begin_print ? "s:\n" : ":\n");
                                                                                        assert(begin_print < end_print);
                                                                                        do
                                                                                        {
                                                                                            mCRL2log(log::debug) << "\t\t"
                                                                                                                     << begin_print->st->debug_id(partitioner);
                                                                                            if (B != begin_print->st->bl.ock)
                                                                                            {
                                                                                                mCRL2log(log::debug) << ", inconsistent: points "
                                                                                                       "to " << begin_print->st->bl.ock->debug_id(partitioner);
                                                                                            }
                                                                                            if (begin_print != begin_print->st->pos)
                                                                                            {
                                                                                                mCRL2log(log::debug)
                                                                                                               << ", inconsistent pointer to state_info_entry";
                                                                                            }
                                                                                            mCRL2log(log::debug) << '\n';
                                                                                        }
                                                                                        while (++begin_print < end_print);
                                                                                    }
                                                                                  public:
                                                                                    /// \brief print the partition per block
                                                                                    /// \details The function prints all blocks in order.  For each block, it
                                                                                    /// lists its states, separated into nonbottom and bottom states.
                                                                                    /// \param partitioner  LTS partitioner (used to print more details)
                                                                                    template<class LTS_TYPE>
                                                                                    void print_part(const bisim_partitioner_dnj<LTS_TYPE>& partitioner) const
                                                                                    {
                                                                                        if (!mCRL2logEnabled(log::debug))  return;
                                                                                        const block_t* B(permutation.front().st->bl.ock);
                                                                                        do
                                                                                        {
                                                                                            mCRL2log(log::debug)<<B->debug_id(partitioner)<<":\n";
                                                                                            print_block(B, "Bottom state",
                                                                                                                B->begin, B->marked_bottom_begin, partitioner);
                                                                                            print_block(B, "Marked bottom state",
                                                                                                      B->marked_bottom_begin, B->nonbottom_begin, partitioner);
                                                                                            print_block(B, "Non-bottom state",
                                                                                                   B->nonbottom_begin, B->marked_nonbottom_begin, partitioner);
                                                                                            print_block(B, "Marked non-bottom state",
                                                                                                               B->marked_nonbottom_begin, B->end, partitioner);
                                                                                        // go to next block
                                                                                        }
                                                                                        while(B->end<permutation.data_end() && (B = B->end->st->bl.ock, true));
                                                                                    }

                                                                                    /// \brief asserts that the partition of states is consistent
                                                                                    /// \details It also requires that no states are marked.
                                                                                    /// \param partitioner  LTS partitioner (used to print more details)
                                                                                    template<class LTS_TYPE>
                                                                                    void assert_consistency(
                                                                                                      const bisim_partitioner_dnj<LTS_TYPE>& partitioner) const
                                                                                    {
                                                                                        const permutation_entry* perm_iter(permutation.data());
                                                                                        state_type true_nr_of_blocks(0);
                                                                                        assert(perm_iter < permutation.data_end());
                                                                                        do
                                                                                        {
                                                                                            const block_t* const block(perm_iter->st->bl.ock);
                                                                                            // block is consistent:
                                                                                            assert(block->begin == perm_iter);
                                                                                            assert(block->begin < block->marked_bottom_begin);
                                                                                            assert(block->marked_bottom_begin == block->nonbottom_begin);
                                                                                            assert(block->nonbottom_begin <= block->marked_nonbottom_begin);
                                                                                            assert(block->marked_nonbottom_begin == block->end);
                                                                                            assert(partitioner.branching||block->nonbottom_begin==block->end);
                                                                                            assert(0 <= block->seqnr);
                                                                                            assert(block->seqnr < nr_of_blocks);
                                                                                            unsigned const max_block(check_complexity::log_n -
                                                                                                                       check_complexity::ilog2(block->size()));
                                                                                            mCRL2complexity(block, no_temporary_work(max_block), partitioner);

                                                                                            // states in the block are consistent:
                                                                                            do
                                                                                            {
                                                                                                const state_info_entry* const state(perm_iter->st);
                                                                                                // assert(part_tr.pred.data() < state->pred_inert.begin);
                                                                                                assert(&state_info.back() == state ||
                                                                                                         state->pred_inert.begin <= state[1].pred_inert.begin);
                                                                                                // assert(state->pred_inert.begin < part_tr.pred.data_end());
                                                                                                // assert(state->succ_inert.begin < part_tr.succ.data_end());
                                                                                                if (perm_iter < block->nonbottom_begin)
                                                                                                {
                                                                                                    assert(&state_info.back() == state || state->
                                                                                                                succ_inert.begin <= state[1].succ_inert.begin);
                                                                                                    // assert(state->succ_inert.begin==&part_tr.succ.back() ||
                                                                                                    //    state <
                                                                                                    //     state->succ_inert.begin->block_bunch->pred->source);
                                                                                                    mCRL2complexity(state, no_temporary_work(max_block, true),
                                                                                                                                                  partitioner);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                    // assert(state->succ_inert.begin < &part_tr.succ.back());
                                                                                                    assert(&state_info.back() == state || state->
                                                                                                                 succ_inert.begin < state[1].succ_inert.begin);
                                                                                                    //assert(state ==
                                                                                                    //     state->succ_inert.begin->block_bunch->pred->source);
                                                                                                    mCRL2complexity(state, no_temporary_work(max_block, false),
                                                                                                                                                  partitioner);
                                                                                                }
                                                                                                assert(block == state->bl.ock);
                                                                                                assert(perm_iter == state->pos);
                                                                                            }
                                                                                            while (++perm_iter < block->end);
                                                                                            assert(perm_iter == block->end);
                                                                                            ++true_nr_of_blocks;
                                                                                        }
                                                                                        while (perm_iter < permutation.data_end());
                                                                                        assert(nr_of_blocks == true_nr_of_blocks);
                                                                                    }
                                                                                #endif
};

///@} (end of group part_state)




                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
/* ************************************************************************* */     static struct {
/*                                                                           */         bool operator()(const iterator_or_counter<action_block_entry*> p1,
/*                           T R A N S I T I O N S                           */                             const action_block_entry* const action_block) const
/*                                                                           */         {
/* ************************************************************************* */             return p1.begin > action_block;
                                                                                        }
                                                                                    } const action_label_greater;
                                                                                #endif


/// \defgroup part_trans
/// \brief data structures for transitions used during partition refinement
/// \details These definitions provide a partition for transition data
/// structure that can be used for the partition refinement algorithm.
///
/// Basically, transitions are stored in four arrays:
/// - `pred`: transitions grouped by goal state, to allow finding all
///   predecessors of a goal state.
///   (At the beginning and the end of the pred array, there are dummy
///   entries.)
/// - `succ`: transitions grouped by source state and bunch, to allow finding
///   all successors of a source state.  Given a transition in this array, it
///   is easy to find all transitions from the same source state in the same
///   bunch.
///   (At the beginning and the end of the succ array, there are dummy
///   entries.)
/// - `action_block`: a permutation of the transitions such that transitions
///   in the same bunch are adjacent, and within each bunch transitions with
///   the same action label and target block.
///   (Between two action_block-slices with different actions, there is a dummy
///   entry.)
/// - `block_bunch`: a permutation of the transitions such that transitions
///   from the same block in the same bunch are adjacent.
///   (At the beginning of the block_bunch array, there is a dummy entry.)
/// Entries in these four arrays are linked with each other with circular
/// iterators, so that one can find the corresponding entry in another array.
///
/// Within this sort order, inert transitions are always placed after non-inert
/// transitions.
///
/// state_info_entry and block_t (defined above) contain pointers to the slices
/// of these arrays.  For bunches and block_bunch-slices, we additionally
/// create _descriptors_ that hold some information about the slice.

///@{

/// \brief information about a transition sorted per source state
class succ_entry
{
  public:
    /// \brief circular iterator to link the four transition arrays
    block_bunch_entry* block_bunch;

    /// \brief pointer to delimit the slice of transitions in the same bunch
    /// \details For most transitions, this pointer points to the first
    /// transition that starts in the same state and belongs to the same
    /// bunch.  But if this transition is the first such transition, the
    /// pointer points to the last such transition (not one past the last, like
    /// otherwise in C and C++).
    ///
    /// For inert transitions, the value is nullptr.
    succ_entry* begin_or_before_end;


    /// \brief find the beginning of the out-slice
    succ_entry* out_slice_begin(                                                ONLY_IF_DEBUG( const fixed_vector<succ_entry>& succ )
                                );


    /// \brief find the bunch of the transition
    bunch_t* bunch() const;
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    /// \brief assign work to the transitions in an out-slice (i.e. the
                                                                                    /// transitions from one state in a specific bunch)
                                                                                    /// \details This debugging function is called to account for work that
                                                                                    /// could be assigned to any transition in the out-slice.  Just to make
                                                                                    /// sure, we therefore set the corresponding counter `ctr` for every
                                                                                    /// transition in the out-slice to `max_value`.
                                                                                    /// \param partitioner      LTS partitioner
                                                                                    /// \param out_slice_begin  pointer to the first transition in the
                                                                                    ///                         out-slice
                                                                                    /// \param ctr              type of the counter that work is assigned to
                                                                                    /// \param max_value        new value that the counter should get
                                                                                    template <class LTS_TYPE>
                                                                                    static inline void add_work_to_out_slice(
                                                                                                  const bisim_partitioner_dnj<LTS_TYPE>& partitioner,
                                                                                                  const succ_entry* out_slice_begin,
                                                                                                  enum check_complexity::counter_type ctr, unsigned max_value);
                                                                                #endif
};


/// \brief information about a transition grouped per (source block, bunch)
/// pair
class block_bunch_entry
{
  public:
    /// \brief circular iterator to link the four transition arrays
    pred_entry* pred;

    /// \brief block_bunch-slice of which this transition is part
    /// \details The slice is null iff the transition is inert.
    block_bunch_slice_iter_or_null_t slice;
};


/// \brief information about a transition sorted per target state
/// \details As I expect the transitions in this array to be moved least often,
/// I store the information on source and target state in this entry.  (It
/// could be stored in any of the four arrays describing the transition;
/// through the circular iterators, the information would be available anyway.)
class pred_entry
{
  public:
    /// \brief circular iterator to link the four transition arrays
    action_block_entry* action_block;

    /// \brief source state of the transition
    state_info_entry* source;

    /// \brief target state of the transition
    state_info_entry* target;
                                                                                #ifndef NDEBUG
                                                                                    /// \brief print a short transition identification for debugging
                                                                                    template <class LTS_TYPE>
                                                                                    std::string debug_id_short(const bisim_partitioner_dnj<LTS_TYPE>&
                                                                                                                                             partitioner) const
                                                                                    {
                                                                                        return "from " + source->debug_id_short(partitioner) +
                                                                                                                  " to " + target->debug_id_short(partitioner);
                                                                                    }

                                                                                    /// \brief print a transition identification for debugging
                                                                                    template <class LTS_TYPE>
                                                                                    std::string debug_id(const bisim_partitioner_dnj<LTS_TYPE>& partitioner)
                                                                                                                                                          const
                                                                                    {
                                                                                        // Search for the action label in partitioner.action_label
                                                                                        label_type const label(std::lower_bound(
                                                                                            partitioner.action_label.cbegin(), partitioner.action_label.cend(),
                                                                                                                      action_block, action_label_greater) -
                                                                                                                            partitioner.action_label.cbegin());
                                                                                        assert(label < partitioner.action_label.size());
                                                                                        assert(partitioner.action_label[label].begin <= action_block);
                                                                                        assert(0==label||action_block<partitioner.action_label[label-1].begin);
                                                                                        // class lts_lts_t uses a function pp() to transform the action label
                                                                                        // to a string.
                                                                                        return pp(partitioner.aut.action_label(label)) + "-transition " +
                                                                                                                                   debug_id_short(partitioner);
                                                                                    }
                                                                                #endif
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    mutable check_complexity::trans_dnj_counter_t work_counter;
                                                                                #endif
};


/// \brief information about a transition sorted per (action, target block)
/// pair
class action_block_entry
{
  public:
    /// \brief circular iterator to link the four transition arrays
    /// \details This iterator can be nullptr because we need to insert dummy
    /// elements between two action_block-slices during initialisation, to make
    /// it easier for first_move_transition_to_new_action_block() to detect
    /// whether two action_block-slices belong to the same action or not.
    succ_entry* succ;

    /// \brief pointer to delimit the slice of transitions in the same (action,
    /// block) pair
    /// \details For most transitions, this pointer points to the first
    /// transition that has the same action and goes to the same block.  But if
    /// this transition is the first such transition, the pointer points to the
    /// last such transition (not one past the last, like otherwise in C and
    /// C++).
    ///
    /// For inert transitions and dummy entries, the value is nullptr.
    action_block_entry* begin_or_before_end;


    /// \brief find the beginning of the action_block-slice
    action_block_entry* action_block_slice_begin(                               ONLY_IF_DEBUG( const action_block_entry* const action_block_begin,
                                                                                                const action_block_entry* const action_block_orig_inert_begin )
                                                 )
    {
        action_block_entry* result(begin_or_before_end);                        assert(nullptr != result);
        if (this < result)
        {                                                                       assert(this == result->begin_or_before_end);
            result = this;                                                      // The following assertion does not always hold: the function is called
        }                                                                       // immediately after a block is refined, so it may be the case that the
                                                                                // transitions are still to be moved to different slices.
                                                                                // assert(succ->block_bunch->pred->target->bl.ock ==
                                                                                //                            result->succ->block_bunch->pred->target->bl.ock);
                                                                                assert(nullptr != succ);  assert(nullptr != result->succ);
                                                                                assert(succ->bunch() == result->succ->bunch());
                                                                                assert(result == action_block_begin || nullptr == result[-1].succ ||
                                                                                            action_block_orig_inert_begin <= result ||
                                                                                            result[-1].succ->block_bunch->pred->target->bl.ock !=
                                                                                                              result->succ->block_bunch->pred->target->bl.ock);
                                                                                // assert(this has the same action as result);
        return result;
    }
};


class part_trans_t;

/// \brief bunch of transitions
/// \details Like a slice, at the end of the algorithm there will be a bunch
/// for every transition in the bisimulation quotient.  Therefore, we should
/// try to minimize the size of a bunch as much as possible.
class bunch_t
{
  public:
    /// \brief first transition in the bunch
    action_block_entry* begin;

    /// \brief pointer past the last transition in the bunch
    action_block_entry* end;

    /// \brief pointer to next non-trivial bunch (in the single-linked list) or
    /// label
    /// \details During refinement, this field stores a pointer to the next
    /// nontrivial bunch.  After refinement, it is set to the label.
    union next_nontrivial_and_label_t
    {
        /// \brief pointer to the next non-trivial bunch in the single-linked
        /// list
        /// \details This pointer is == nullptr if the bunch is trivial.  The
        /// last entry in the list points to itself so its pointer is still
        /// != nullptr.
        bunch_t* next_nontrivial;

        /// \brief action label of the transition
        label_type label;


        /// \brief constructor
        next_nontrivial_and_label_t()
        {
            next_nontrivial = nullptr;
        }
    } next_nontrivial_and_label;


    /// \brief constructor
    bunch_t(action_block_entry* const new_begin,
                                             action_block_entry* const new_end)
      : begin(new_begin),
        end(new_end),
        next_nontrivial_and_label()
    {  }


    /// \brief returns true iff the bunch is trivial
    /// \details If this bunch is the last in the list of non-trivial bunches,
    /// the convention is that the next pointer points to this bunch itself (to
    /// distinguish it from nullptr).
    bool is_trivial() const
    {
        return nullptr == next_nontrivial_and_label.next_nontrivial;
    }


    /// \brief split off a single action_block-slice from the bunch
    /// \details The function splits the current bunch after its first
    /// action_block-slice or before its last action_block-slice, whichever
    /// is smaller.  It creates a new bunch for the split-off slice and
    /// returns a pointer to the new bunch.  The caller has to adapt the
    /// block_bunch-slices.
    bunch_t* split_off_small_action_block_slice(part_trans_t& part_tr);
                                                                                #ifndef NDEBUG
                                                                                    /// \brief print a short bunch identification for debugging
                                                                                    template <class LTS_TYPE>
                                                                                    std::string debug_id_short(const bisim_partitioner_dnj<LTS_TYPE>&
                                                                                                                                             partitioner) const
                                                                                    {
                                                                                        assert(partitioner.part_tr.action_block.data() <= begin);
                                                                                        assert(end <= partitioner.part_tr.action_block_inert_begin);
                                                                                        return "bunch [" + std::to_string(begin -
                                                                                                               partitioner.part_tr.action_block.data()) + "," +
                                                                                           std::to_string(end - partitioner.part_tr.action_block.data()) + ")";
                                                                                    }

                                                                                    /// \brief print a long bunch identification for debugging
                                                                                    template <class LTS_TYPE>
                                                                                    std::string debug_id(const bisim_partitioner_dnj<LTS_TYPE>& partitioner)
                                                                                                                                                          const
                                                                                    {   assert(nullptr != end[-1].succ);
                                                                                        const action_block_entry* iter(begin);  assert(iter < end);
                                                                                        assert(nullptr != iter->succ);
                                                                                        assert(iter == iter->succ->block_bunch->pred->action_block);
                                                                                        std::string result(debug_id_short(partitioner));
                                                                                        result += " containing transition";
                                                                                        result += iter < end - 1 ? "s " : " ";
                                                                                        result += iter->succ->block_bunch->pred->debug_id_short(partitioner);
                                                                                        ++iter;
                                                                                        if (end <= iter)  return result;
                                                                                        while (nullptr == iter->succ)  ++iter;
                                                                                        assert(iter < end);
                                                                                        assert(iter == iter->succ->block_bunch->pred->action_block);
                                                                                        result += ", ";
                                                                                        result += iter->succ->block_bunch->pred->debug_id_short(partitioner);
                                                                                        if (iter < end - 3)
                                                                                        {
                                                                                            result += ", ...";
                                                                                            iter = end - 3;
                                                                                        }
                                                                                        while (++iter < end)
                                                                                        {
                                                                                            if (nullptr != iter->succ)
                                                                                            {   assert(iter == iter->succ->block_bunch->pred->action_block);
                                                                                                result += ", ";
                                                                                                result += iter->succ->block_bunch->pred->debug_id_short(
                                                                                                                                                  partitioner);
                                                                                            }
                                                                                        }
                                                                                        return result;
                                                                                    }
                                                                                #endif
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    /// \brief calculates the maximal allowed value for work counters
                                                                                    /// associated with this bunch
                                                                                    /// \details Work counters may only be nonzero if this bunch is a
                                                                                    /// single-action bunch, i. e. all its transitions have the same action
                                                                                    /// label.  Also, only then the size can be calculated as end - begin.
                                                                                    template <class LTS_TYPE>
                                                                                    int max_work_counter(const bisim_partitioner_dnj<LTS_TYPE>& partitioner)
                                                                                                                                                          const
                                                                                    {
                                                                                        // verify that the bunch only has a single action label.
                                                                                        // Search for the action label in partitioner.action_label
                                                                                        label_type const label(std::lower_bound(
                                                                                            partitioner.action_label.cbegin(), partitioner.action_label.cend(),
                                                                                                                                begin, action_label_greater) -
                                                                                                                            partitioner.action_label.cbegin());
                                                                                        assert(label < partitioner.action_label.size());
                                                                                        assert(partitioner.action_label[label].begin <= begin);
                                                                                        assert(0 == label || begin < partitioner.action_label[label-1].begin);
                                                                                        if (0 == label || end < partitioner.action_label[label - 1].begin)
                                                                                        {
                                                                                            assert(check_complexity::ilog2(end - begin) <=
                                                                                                                                      check_complexity::log_n);
                                                                                            return check_complexity::log_n-check_complexity::ilog2(end-begin);
                                                                                        }
                                                                                        return 0;
                                                                                    }

                                                                                    mutable check_complexity::bunch_dnj_counter_t work_counter;
                                                                                #endif
};


/// \brief Information about a set of transitions with the same source block,
/// in the same bunch
/// \details A block_bunch-slices contains the transitions in a single bunch
/// that start in the same block.  In the end, we want each block to be stable
/// under its block_bunch-slices, i.e. every bottom state has a transition in
/// every block_bunch-slice of the block.  Then, there will be one slice for
/// each transition in the minimised LTS, so we should try to minimize this
/// data structure as much as possible.
///
/// Also note that these slices are part of a doubly-linked list.  We cannot
/// change this to a singly-linked list because we occasionally delete an
/// element from this list.  This would be possible with a single-linked list
/// if we could infer the order of the list somehow, e.g. from the transitions
/// in a bottom state -- however, this does not work when a block loses all its
/// bottom states, i.e. when a block splits into a small U and large R but R
/// does not contain any bottom states.
class block_bunch_slice_t
{
  public:
    /// \brief pointer past the end of the transitions in the block_bunch array
    /// \details We do not need a begin pointer because we can always walk
    /// through the transitions in the slice from end to beginning.
    block_bunch_entry* end;

    /// bunch to which this slice belongs
    bunch_t* bunch;

    /// \brief pointer to the first marked transition in the block_bunch array
    /// \details If this pointer is nullptr, then the block_bunch_slice is
    /// stable.
    block_bunch_entry* marked_begin;


    /// \brief returns true iff the block_bunch-slice is registered as stable
    bool is_stable() const  {  return nullptr == marked_begin;  }


    /// \brief register that the block_bunch-slice is stable
    void make_stable()
    {                                                                           assert(!is_stable());  assert(!empty());
        marked_begin = nullptr;
    }


    /// \brief register that the block_bunch-slice is not stable
    void make_unstable()
    {                                                                           assert(is_stable());
        marked_begin = end;                                                     assert(!is_stable());
    }


    /// \brief constructor
    block_bunch_slice_t(block_bunch_entry* const new_end,
                            bunch_t* const new_bunch, bool const new_is_stable)
      : end(new_end),
        bunch(new_bunch),
        marked_begin(nullptr)
    {
        if (!new_is_stable)  make_unstable();
    }


    /// \brief returns true iff the block_bunch-slice is empty
    /// \details A block_bunch-slice should only become empty if it is
    /// unstable.
    bool empty() const
    {                                                                           // assert(std::less(part_tr.block_bunch.data(), end));
                                                                                // assert(!std::less(part_tr.block_bunch_inert_begin, end));
                                                                                // assert(part_tr.block_bunch.front().slice != this);
        return end[-1].slice != this;
    }


    /// compute the source block of the transitions in this slice
    block_t* source_block() const
    {                                                                           assert(!empty());
        return end[-1].pred->source->bl.ock;
    }
                                                                                #ifndef NDEBUG
                                                                                    /// \brief print a block_bunch-slice identification for debugging
                                                                                    template <class LTS_TYPE>
                                                                                    std::string debug_id(const bisim_partitioner_dnj<LTS_TYPE>& partitioner)
                                                                                                                                                          const
                                                                                    {
                                                                                        static struct {
                                                                                            bool operator()(const block_bunch_entry& p1,
                                                                                                                     const block_bunch_slice_t* const p2) const
                                                                                            {
                                                                                                return p1.slice != p2;
                                                                                            }
                                                                                        } const block_bunch_not_equal;

                                                                                        assert(partitioner.part_tr.block_bunch.data() < end);
                                                                                        assert(end <= partitioner.part_tr.block_bunch_inert_begin);
                                                                                        std::string const index_string(std::to_string(end -
                                                                                                                &partitioner.part_tr.block_bunch.cbegin()[1]));
                                                                                        if (empty())
                                                                                        {   //assert(!is_stable());
                                                                                            return "empty block_bunch_slice [" + index_string + "," +
                                                                                                                                            index_string + ")";
                                                                                        }
                                                                                        const block_bunch_entry* begin(
                                                                                                                 &partitioner.part_tr.block_bunch.cbegin()[1]);
                                                                                        if (trans_type bunch_size(bunch->end - bunch->begin);
                                                                                                                       (trans_type) (end - begin) > bunch_size)
                                                                                        {
                                                                                            begin = end - bunch_size;
                                                                                        }
                                                                                        begin = std::lower_bound(begin, const_cast<const block_bunch_entry*>
                                                                                                     (is_stable() || marked_begin==end ? end-1 : marked_begin),
                                                                                                                                  this, block_bunch_not_equal);
                                                                                        assert(begin->slice == this);
                                                                                        assert(begin[-1].slice != this);
                                                                                        return (is_stable() ? "stable block_bunch-slice ["
                                                                                                            : "unstable block_bunch_slice [") +
                                                                                           std::to_string(begin-&partitioner.part_tr.block_bunch.cbegin()[1]) +
                                                                                                        "," + index_string + ") containing transitions from " +
                                                                                                                   source_block()->debug_id(partitioner) +
                                                                                                                   " in " + bunch->debug_id_short(partitioner);
                                                                                    }
                                                                                #endif
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    /// \brief add work to transitions starting in bottom states
                                                                                    /// \details Sometimes an action is done whose time could be accounted for
                                                                                    /// by any transition starting in a bottom state of the block.
                                                                                    /// \param ctr          counter type to which work is assigned
                                                                                    /// \param max_value    new value of the counter
                                                                                    /// \param partitioner  LTS partitioner (to print error messages if
                                                                                    ///                     necessary)
                                                                                    template <class LTS_TYPE>
                                                                                    bool add_work_to_bottom_transns(enum check_complexity::
                                                                                                      counter_type const ctr, unsigned const max_value,
                                                                                                      const bisim_partitioner_dnj<LTS_TYPE>& partitioner) const
                                                                                    {   assert(!empty());
                                                                                        assert(1U == max_value);
                                                                                        const block_t* const block(source_block());
                                                                                        bool result(false);
                                                                                        const block_bunch_entry* block_bunch(end);
                                                                                        assert(partitioner.part_tr.block_bunch.front().slice != this);
                                                                                        assert(block_bunch[-1].slice == this);
                                                                                        do
                                                                                        {
                                                                                            --block_bunch;
                                                                                            const state_info_entry* const source(block_bunch->pred->source);
                                                                                            assert(source->bl.ock == block);
                                                                                            if (source->pos < block->nonbottom_begin /*&&
                                                                                                // the transition starts in a (new) bottom state
                                                                                                block_bunch->pred->work_counter.counters[ctr -
                                                                                                               check_complexity::TRANS_dnj_MIN] != max_value*/)
                                                                                            {
                                                                                                mCRL2complexity(block_bunch->pred, add_work(ctr, max_value),
                                                                                                                                                  partitioner);
                                                                                                (void) partitioner; // avoid unused variable warning
                                                                                                result = true;
                                                                                            }
                                                                                        }
                                                                                        while (block_bunch[-1].slice == this);
                                                                                        return result;
                                                                                    }

                                                                                    mutable check_complexity::block_bunch_dnj_counter_t work_counter;
                                                                                #endif
};


/// \brief find the beginning of the out-slice
inline succ_entry* succ_entry::out_slice_begin(                                 ONLY_IF_DEBUG( const fixed_vector<succ_entry>& succ )
                                               )
{                                                                               assert(nullptr != begin_or_before_end);
    succ_entry* result(begin_or_before_end);                                    assert(result->block_bunch->pred->action_block->succ == result);
    if (this < result)
    {                                                                           assert(nullptr != result->begin_or_before_end);
                                                                                assert(this == result->begin_or_before_end);
        result = this;                                                          assert(result->block_bunch->pred->action_block->succ == result);
    }                                                                           assert(block_bunch->pred->source == result->block_bunch->pred->source);
                                                                                // assert(this <= result); //< holds always, based on the if() above
                                                                                assert(nullptr != result->begin_or_before_end);
                                                                                assert(this <= result->begin_or_before_end);
                                                                                assert(block_bunch->slice == result->block_bunch->slice);
                                                                                assert(&succ.cbegin()[1] == result ||
                                                                                        result[-1].block_bunch->pred->source < block_bunch->pred->source ||
                                                                                                              result[-1].bunch() != block_bunch->slice->bunch);
    return result;
}


/// \brief find the bunch of a transition
inline bunch_t* succ_entry::bunch() const
{
    return block_bunch->slice->bunch;
}
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    /// \brief register that work has been done for the out-slice containing
                                                                                    /// `out_slice_begin`
                                                                                    /// \details This function should be used if work
                                                                                    /// \param partitioner      the partitioner data structure, used to write
                                                                                    ///                         diagnostic messages
                                                                                    /// \param out_slice_begin  the first transition in the out-slice
                                                                                    /// \param ctr              counter type to which work has to be assigned
                                                                                    /// \param max_value        new value of the counter
                                                                                    template <class LTS_TYPE>
                                                                                    /* static */ inline void succ_entry::add_work_to_out_slice(
                                                                                        const bisim_partitioner_dnj<LTS_TYPE>& partitioner,
                                                                                        const succ_entry* out_slice_begin,
                                                                                        enum check_complexity::counter_type const ctr,unsigned const max_value)
                                                                                    {
                                                                                        const succ_entry* const out_slice_before_end(
                                                                                                                         out_slice_begin->begin_or_before_end);
                                                                                        assert(nullptr != out_slice_before_end);
                                                                                        assert(out_slice_begin <= out_slice_before_end);
                                                                                        mCRL2complexity(out_slice_begin->block_bunch->pred,
                                                                                                                        add_work(ctr, max_value), partitioner);
                                                                                        #ifndef NDEBUG
                                                                                            while (++out_slice_begin <= out_slice_before_end)
                                                                                            {
                                                                                                // treat temporary counters specially
                                                                                                mCRL2complexity(out_slice_begin->block_bunch->pred,
                                                                                                            add_work_notemporary(ctr, max_value), partitioner);
                                                                                            }
                                                                                        #else
                                                                                            (void) partitioner; (void) out_slice_before_end;
                                                                                            // avoid unused variable warning
                                                                                        #endif
                                                                                    }
                                                                                #endif
class part_trans_t
{
  public:
    /// \brief array containing all successor entries
    /// \details The first and last entry are dummy entries, pointing to a
    /// transition from nullptr to nullptr, to make it easier to check whether
    /// there is another transition from the current state.
    fixed_vector<succ_entry> succ;

    /// \brief array containing all block_bunch entries
    /// \details The first entry is a dummy entry, pointing to a transition not
    /// contained in any slice, to make it easier to check whether there is
    /// another transition in the current block_bunch.
    fixed_vector<block_bunch_entry> block_bunch;

    /// \brief array containing all predecessor entries
    /// \details The first and last entry are dummy entries, pointing to a
    /// transition to nullptr, to make it easier to check whether there is
    /// another transition to the current state.
    fixed_vector<pred_entry> pred;

    /// \brief array containing all action_block entries
    /// \details During initialisation, the transitions are sorted according to
    /// their label. Between transitions with different labels there is a dummy
    /// entry, to make it easier to check whether there is another transition
    /// in the current action_block slice.
    ///
    /// The array may be empty if there are no transitions.
    fixed_vector<action_block_entry> action_block;

    /// \brief pointer to the first inert transition in block_bunch
    block_bunch_entry* block_bunch_inert_begin;

    /// \brief pointer to the first inert transition in action_block
    action_block_entry* action_block_inert_begin;
                                                                                #ifndef NDEBUG
                                                                                    /// \brief pointer to the first inert transition in the initial partition
                                                                                    const action_block_entry* action_block_orig_inert_begin;
                                                                                #endif
    /// \brief list of unstable block_bunch-slices
    simple_list<block_bunch_slice_t> splitter_list;
   private:
    /// \brief pointer to first non-trivial bunch
    bunch_t* first_nontrivial;

  public:
    #ifdef USE_POOL_ALLOCATOR
                                                                                static_assert(std::is_trivially_destructible<bunch_t>::value);
    #endif

    /// \brief number of new bottom states found until now.
    state_type nr_of_new_bottom_states;

    /// \brief counters to measure progress
    trans_type nr_of_bunches;
    trans_type nr_of_nontrivial_bunches;
    trans_type nr_of_action_block_slices;
    trans_type nr_of_block_bunch_slices;

    /// \brief constructor
    /// \details The constructor sets up the dummy transitions at the beginning
    /// and end of the succ, block_bunch and pred arrays.  (Dummy transitions
    /// in action_block depend on the number of transitions per action label,
    /// so they cannot be set up without knowing details about how many
    /// transitions have which label.)
    /// \param num_transitions  number of transitions of the LTS
    /// \param num_actions      number of action labels of the LTS
    part_trans_t(trans_type num_transitions,
                 trans_type num_actions)
      : succ(num_transitions + 2),
        block_bunch(num_transitions + 1),
        pred(num_transitions + 2),
        action_block(num_transitions + num_actions < 1
                                      ? 0 : num_transitions + num_actions - 1),
        block_bunch_inert_begin(block_bunch.data_end()),
        action_block_inert_begin(action_block.data_end()),
        splitter_list(),
        first_nontrivial(nullptr),
        nr_of_new_bottom_states(0),
        nr_of_bunches(0),
        nr_of_nontrivial_bunches(0),
        nr_of_action_block_slices(0),
        nr_of_block_bunch_slices(0)
    {
        succ.front().block_bunch = block_bunch.data();
        succ.back() .block_bunch = block_bunch.data();
        block_bunch.front().pred = pred.data();
        block_bunch.front().slice = nullptr;
        pred.front().source = nullptr;
        pred.front().target = nullptr;
        pred.back() .source = nullptr;
        pred.back() .target = nullptr;
    }


    /// \brief destructor
    ~part_trans_t()
    {
        #ifndef USE_POOL_ALLOCATOR
            // The destructor also deallocates the bunches, as they are not
            // directly referenced from anywhere.  This is only necessary if we
            // do not use the pool allocator, as the latter will destroy the
            // bunches wholesale.
            action_block_entry* action_block_iter(action_block.data());
            for (;;)
            {
                do
                {
                    if (action_block_inert_begin <= action_block_iter)
                    {                                                           assert(0 == nr_of_bunches);
                        return;
                    }
                }
                while (nullptr == action_block_iter->succ && (                  assert(nullptr == action_block_iter->begin_or_before_end),
                                                   ++action_block_iter, true)); assert(nullptr != action_block_iter->begin_or_before_end);
                bunch_t* const bunch(action_block_iter->succ->bunch());         assert(bunch->begin == action_block_iter);
                action_block_iter = bunch->end;
                delete bunch;                                                   ONLY_IF_DEBUG( --nr_of_bunches; )
            }
            /* unreachable                                                   */ assert(0);
        #endif
    }


    /// \brief provide some bunch from the list of non-trivial bunches
    /// \returns pointer to a bunch that is in the list of non-trivial bunches
    bunch_t* get_some_nontrivial()
    {
        return first_nontrivial;
    }


    /// \brief insert a bunch into the list of nontrivial bunches
    /// \param bunch  the bunch that has become non-trivial
    void make_nontrivial(bunch_t* const bunch)
    {                                                                           assert(1 < bunch->end - bunch->begin);  assert(bunch->is_trivial());
                                                                                // The following assertions do not necessarily hold during initialisation:
                                                                                //assert(bunch->begin <= bunch->begin->begin_or_before_end);
        bunch->next_nontrivial_and_label.next_nontrivial =
                        nullptr == first_nontrivial ? bunch : first_nontrivial; //assert(nullptr != bunch->begin->begin_or_before_end);
                                                                                //assert(nullptr != bunch->end[-1].begin_or_before_end);
                                                                                //assert(bunch->begin->begin_or_before_end <
                                                                                //                                         bunch->end[-1].begin_or_before_end);
                                                                                //assert(nullptr != end[-1].begin_or_before_end);
        first_nontrivial = bunch;                                               assert(nr_of_nontrivial_bunches < nr_of_bunches);
        ++nr_of_nontrivial_bunches;                                             //assert(end[-1].begin_or_before_end <= end);
    }


    /// \brief remove a bunch from the list of nontrivial bunches
    /// \param bunch  the bunch that has become trivial
    void make_trivial(bunch_t* const bunch)
    {                                                                           assert(!bunch->is_trivial());  assert(first_nontrivial == bunch);
        first_nontrivial =
                bunch == bunch->next_nontrivial_and_label.next_nontrivial
                  ? nullptr : bunch->next_nontrivial_and_label.next_nontrivial; assert(bunch->end - 1 == bunch->begin->begin_or_before_end);
        bunch->next_nontrivial_and_label.next_nontrivial = nullptr;             assert(0 < nr_of_nontrivial_bunches);
        --nr_of_nontrivial_bunches;                                             assert(bunch->begin == bunch->end[-1].begin_or_before_end);
    }


    /// \brief transition is moved to a new bunch
    /// \details This (and the next function) have to be called after a
    /// transition has changed its bunch.  The member function will adapt the
    /// transition data structure.  It assumes that the transition is
    /// non-inert.
    ///
    /// The work has to be done in two steps: We call the first step
    /// first_move_transition_to_new_bunch() for each transition in the new
    /// bunch, and then call second_move_transition_to_new_bunch() again for
    /// all these transitions.  The reason is that some data structures need to
    /// be finalised in the second phase.
    ///
    /// The first phase moves all transitions to their correct position in
    /// the out-slices and block_bunch-slices, but it doesn't yet create
    /// a fully correct new out-slice and block_bunch-slice.  It adapts
    /// current_out_slice of all states with a transition in the new bunch.
    /// \param action_block_iter_iter     transition that has to be changed
    /// \param bunch_T_a_Bprime           the new bunch in which the transition
    ///                                   lies
    /// \param first_transition_of_state  true iff this is the first transition
    ///                                   of the state, so a new out-slice has
    ///                                   to be allocated.
    void first_move_transition_to_new_bunch(
                                   action_block_entry* const action_block_iter,
                                   bunch_t* const bunch_T_a_Bprime,
                                   bool const first_transition_of_state)
    {

        /* -  -  -  -  -  -  -  - adapt part_tr.succ -  -  -  -  -  -  -  - */

        succ_entry* const old_succ_pos(action_block_iter->succ);                assert(nullptr != old_succ_pos);
                                                                                assert(old_succ_pos->block_bunch->pred->action_block == action_block_iter);
        succ_entry* const out_slice_begin(old_succ_pos->out_slice_begin(        ONLY_IF_DEBUG( succ )
                                                                  ));           assert(out_slice_begin->block_bunch->pred->action_block->succ ==
                                                                                                                                              out_slice_begin);
        succ_entry* const new_succ_pos(out_slice_begin->begin_or_before_end);   assert(nullptr != new_succ_pos);
                                                                                assert(out_slice_begin == new_succ_pos->begin_or_before_end);
                                                                                assert(new_succ_pos<old_succ_pos->block_bunch->pred->source->succ_inert.begin);
        /* move the transition to the end of its out-slice                   */ assert(new_succ_pos->block_bunch->pred->action_block->succ == new_succ_pos);
        if (old_succ_pos < new_succ_pos)
        {
            std::swap(old_succ_pos->block_bunch, new_succ_pos->block_bunch);
            old_succ_pos->block_bunch->pred->action_block->succ = old_succ_pos; assert(action_block_iter == new_succ_pos->block_bunch->pred->action_block);
            action_block_iter->succ = new_succ_pos;
        }                                                                       else  assert(old_succ_pos == new_succ_pos);

        // adapt the old out-slice immediately
            // If the old out-slice becomes empty, then out_slice_begin ==
            // new_succ_pos, so the two following assignments will assign the
            // same variable.  The second assignment is the relevant one.
        out_slice_begin->begin_or_before_end = new_succ_pos - 1;

        // adapt the new out-slice, as far as is possible now:
            // make the begin_or_before_end pointers of the first and last
            // transition in the slice correct immediately.  The other
            // begin_or_before_end pointers need to be corrected after all
            // transitions in the new bunch have been positioned correctly.
        if (first_transition_of_state)
        {
            new_succ_pos->begin_or_before_end = new_succ_pos;
        }
        else
        {
            succ_entry* const out_slice_before_end(
                                          new_succ_pos[1].begin_or_before_end); assert(nullptr != out_slice_before_end);
                                                                                assert(new_succ_pos < out_slice_before_end);
                                                                                assert(out_slice_before_end->block_bunch->pred->action_block->succ ==
                                                                                                                                         out_slice_before_end);
                                                                                assert(new_succ_pos + 1 == out_slice_before_end->begin_or_before_end);
            out_slice_before_end->begin_or_before_end = new_succ_pos;           assert(out_slice_before_end <
                                                                                                    new_succ_pos->block_bunch->pred->source->succ_inert.begin);
            new_succ_pos->begin_or_before_end = out_slice_before_end;           assert(bunch_T_a_Bprime == out_slice_before_end->bunch());
        }

        /* -  -  -  -  -  -  - adapt part_tr.block_bunch -  -  -  -  -  -  - */ assert(new_succ_pos == action_block_iter->succ);

        block_bunch_entry* const old_block_bunch_pos(
                                                    new_succ_pos->block_bunch); assert(old_block_bunch_pos->pred->action_block == action_block_iter);
        block_t*const source_block = old_block_bunch_pos->pred->source->bl.ock; assert(!old_block_bunch_pos->slice.is_null());
        block_bunch_slice_iter_t const old_block_bunch_slice(
                                                   old_block_bunch_pos->slice);
        block_bunch_entry* const new_block_bunch_pos(
                                               old_block_bunch_slice->end - 1); assert(nullptr != new_block_bunch_pos->pred->action_block->succ);
                                                                                assert(new_block_bunch_pos->pred->action_block->succ->block_bunch ==
                                                                                                                                          new_block_bunch_pos);
        // create or adapt the new block_bunch-slice
        block_bunch_slice_iter_t new_block_bunch_slice;
        if (new_block_bunch_pos + 1 >= block_bunch_inert_begin ||
                       (new_block_bunch_slice = (block_bunch_slice_iter_t)
                                                new_block_bunch_pos[1].slice,   assert(!new_block_bunch_pos[1].slice.is_null()),
                        bunch_T_a_Bprime != new_block_bunch_slice->bunch ||
                        source_block != new_block_bunch_slice->source_block()))
        {                                                                       assert(first_transition_of_state);
            // This is the first transition in the block_bunch-slice.
            // The old block_bunch-slice becomes unstable, and the new
            // block_bunch-slice is created unstable.

            // Note that the new block_bunch-slice should precede the old one.

            #ifdef USE_SIMPLE_LIST
                new_block_bunch_slice = splitter_list.emplace_back(
                             new_block_bunch_pos + 1, bunch_T_a_Bprime, false);
            #else
                splitter_list.emplace_back(new_block_bunch_pos + 1,
                                                      bunch_T_a_Bprime, false);
                new_block_bunch_slice = std::prev(splitter_list.end());
            #endif
            ++nr_of_block_bunch_slices;
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    new_block_bunch_slice->work_counter = old_block_bunch_slice->work_counter;
                                                                                #endif
            splitter_list.splice(splitter_list.end(),
                      source_block->stable_block_bunch, old_block_bunch_slice);
            old_block_bunch_slice->make_unstable();
        }                                                                       assert(!new_block_bunch_slice->is_stable());

        // move the transition to the end of its block_bunch-slice
        if (old_block_bunch_pos < new_block_bunch_pos)
        {
            std::swap(old_block_bunch_pos->pred, new_block_bunch_pos->pred);    assert(nullptr != old_block_bunch_pos->pred->action_block->succ);
            old_block_bunch_pos->pred->action_block->succ->block_bunch =
                                                           old_block_bunch_pos; assert(new_succ_pos == new_block_bunch_pos->pred->action_block->succ);
            new_succ_pos->block_bunch = new_block_bunch_pos;
        }                                                                       else  assert(new_block_bunch_pos == old_block_bunch_pos);
                                                                                assert(new_block_bunch_pos->slice == old_block_bunch_slice);
        new_block_bunch_pos->slice = new_block_bunch_slice;

        /* adapt the old block_bunch-slice                                   */ assert(new_block_bunch_pos + 1 == old_block_bunch_slice->marked_begin);
        old_block_bunch_slice->end = new_block_bunch_pos;
        old_block_bunch_slice->marked_begin = new_block_bunch_pos;              assert(nullptr != new_block_bunch_pos);
        if (old_block_bunch_slice->empty())
        {                                                                       assert(!old_block_bunch_slice->is_stable());
            splitter_list.erase(old_block_bunch_slice);                         assert(!new_block_bunch_slice->is_stable());
            --nr_of_block_bunch_slices;

            // Because now every bottom state has a transition in the new
            // bunch, and no state has a transition in the old bunch, there
            // is no need to refine this block.  So we make this
            // block_bunch-slice stable again.
            source_block->stable_block_bunch.splice(
                                        source_block->stable_block_bunch.end(),
                                        splitter_list, new_block_bunch_slice);
            new_block_bunch_slice->make_stable();

            // unmark the states
                // (This transition must be the last transition from
                // source_block in the new bunch, so unmarking the states now
                // will not be undone by later markings of other states.)
            source_block->marked_nonbottom_begin = source_block->end;           assert(source_block->marked_bottom_begin == source_block->begin);
            source_block->marked_bottom_begin = source_block->nonbottom_begin;
        }
    }


    /// \brief transition is moved to a new bunch, phase 2
    /// \details This (and the previous function) have to be called after a
    /// transition has changed its bunch.  The member function will adapt the
    /// transition data structure.  It assumes that the transition is
    /// non-inert.
    ///
    /// The work has to be done in two steps: We call the first step
    /// first_move_transition_to_new_bunch() for each transition in the new
    /// bunch, and then call second_move_transition_to_new_bunch() again for
    /// all these transitions.  The reason is that some data structures need to
    /// be finalised/ in the second phase.
    ///
    /// The second phase finalizes the new out-slices and block_bunch-slices
    /// that were left half-finished by the first phase.  It assumes that all
    /// block_bunch-slices are registered as stable.
    /// \param action_block_iter_iter     transition that has to be changed
    /// \param large_splitter_bunch       the large splitter_bunch that has
    ///                                   been split; the transition has moved
    ///                                   from `large_splitter_bunch` to a new,
    ///                                   small bunch.
                                                                                ONLY_IF_DEBUG( template <class LTS_TYPE> )
    void second_move_transition_to_new_bunch(
                                   action_block_entry* const action_block_iter, ONLY_IF_DEBUG( const bisim_partitioner_dnj<LTS_TYPE>& partitioner,
                                                                                                                             bunch_t* const bunch_T_a_Bprime, )
                                   bunch_t* const large_splitter_bunch)
    {                                                                           assert(nullptr != bunch_T_a_Bprime);

        /* -  -  -  -  -  -  -  - adapt part_tr.succ -  -  -  -  -  -  -  - */

            // We already moved the transition in part_tr.succ to the correct
            // place in first_move_transition_to_new_bunch(); now we have to
            // set begin_or_before_end.
        succ_entry* const new_succ_pos(action_block_iter->succ);                assert(nullptr != new_succ_pos);
                                                                                assert(new_succ_pos->block_bunch->pred->action_block == action_block_iter);
        state_info_entry* const source(
                                      new_succ_pos->block_bunch->pred->source); assert(source->pos->st == source);
                                                                                assert(new_succ_pos < source->succ_inert.begin);
                                                                                assert(source == partitioner.part_st.state_info.data() ||
                                                                                                                  source[-1].succ_inert.begin <= new_succ_pos);
                                                                                assert(nullptr != new_succ_pos->begin_or_before_end);
        succ_entry* const new_begin_or_before_end(
                       new_succ_pos->begin_or_before_end->begin_or_before_end); assert(nullptr != new_begin_or_before_end);
                                                                                assert(new_begin_or_before_end->block_bunch->pred->action_block->succ ==
                                                                                                                                      new_begin_or_before_end);
        if (new_begin_or_before_end < new_succ_pos)
        {                                                                       assert(source == partitioner.part_st.state_info.data() ||
                                                                                                       source[-1].succ_inert.begin <= new_begin_or_before_end);
            new_succ_pos->begin_or_before_end = new_begin_or_before_end;
        }
        else
        {                                                                       assert(new_begin_or_before_end == new_succ_pos);
            // This is the first or the last transition in the out-slice.
            const succ_entry* const new_before_end(
                                 new_begin_or_before_end->begin_or_before_end); assert(nullptr != new_before_end);
            if (new_begin_or_before_end <= new_before_end)
            {                                                                   assert(&partitioner.part_tr.succ.cbegin()[1] == new_begin_or_before_end ||
                /* This is the first transition in the new out-slice.        */             new_begin_or_before_end[-1].block_bunch->pred->source < source ||
                /* If there is still a transition in the old out-slice,      */                       new_begin_or_before_end[-1].bunch() != bunch_T_a_Bprime);
                /* we prepay for it.                                         */ assert(new_before_end + 1 == source->succ_inert.begin ||
                                                                                                                bunch_T_a_Bprime != new_before_end[1].bunch());
                if (source == new_succ_pos[-1].block_bunch->pred->source &&
                              new_succ_pos[-1].bunch() == large_splitter_bunch)
                {
                    // Mark one transition in the large slice
                    block_bunch_entry* const old_block_bunch_pos(
                                                 new_succ_pos[-1].block_bunch); assert(!old_block_bunch_pos->slice.is_null());
                                                                                assert(old_block_bunch_pos->pred->action_block->succ == new_succ_pos - 1);
                    block_bunch_slice_iter_t const large_splitter_slice(
                                                   old_block_bunch_pos->slice);
                    if (!large_splitter_slice->is_stable())
                    {
                        block_bunch_entry* const new_block_bunch_pos(
                                       large_splitter_slice->marked_begin - 1); assert(nullptr != new_block_bunch_pos->pred->action_block->succ);
                                                                                assert(new_block_bunch_pos->pred->action_block->succ->block_bunch ==
                                                                                                                                          new_block_bunch_pos);
                        if (old_block_bunch_pos < new_block_bunch_pos)
                        {
                            std::swap(old_block_bunch_pos->pred,
                                                    new_block_bunch_pos->pred); // assert(nullptr != old_block_bunch_pos->pred->action_block->succ);
                            old_block_bunch_pos->pred->action_block->
                                       succ->block_bunch = old_block_bunch_pos; assert(new_block_bunch_pos->pred->action_block->succ == new_succ_pos - 1);
                            new_succ_pos[-1].block_bunch = new_block_bunch_pos;
                        }
                        large_splitter_slice->marked_begin=new_block_bunch_pos; assert(nullptr != new_block_bunch_pos);
                    }                                                           else assert(1 >= source->bl.ock->size());
                }
            }                                                                   else assert(source == partitioner.part_st.state_info.data() ||
                                                                                                                source[-1].succ_inert.begin <= new_before_end);
        }
                                                                                #ifndef NDEBUG
        /* -  -  -  -  -  -  - adapt part_tr.block_bunch -  -  -  -  -  -  - */     const block_bunch_entry* new_block_bunch_pos(new_succ_pos->block_bunch);
                                                                                    assert(new_block_bunch_pos->pred->action_block->succ == new_succ_pos);
                                                                                    assert(!new_block_bunch_pos->slice.is_null());
        /* Nothing needs to be done.                                         */     block_bunch_slice_const_iter_t const new_block_bunch_slice(
                                                                                                                                   new_block_bunch_pos->slice);
                                                                                    assert(new_block_bunch_pos < new_block_bunch_slice->end);
                                                                                    assert(bunch_T_a_Bprime == new_block_bunch_slice->bunch);
                                                                                    if (new_block_bunch_pos + 1 < new_block_bunch_slice->end)  return;

                                                                                    // This transition is the last in the block_bunch-slice.  If there
                                                                                    // were some task that would need to be done exactly once per
                                                                                    // block_bunch-slice, this would be the moment.
                                                                                    do  assert(source->bl.ock == new_block_bunch_pos->pred->source->bl.ock);
                                                                                    while ((--new_block_bunch_pos)->slice == new_block_bunch_slice);
                                                                                    assert(new_block_bunch_pos <= partitioner.part_tr.block_bunch.data() ||
                                                                                                source->bl.ock != new_block_bunch_pos->pred->source->bl.ock ||
                                                                                                        bunch_T_a_Bprime != new_block_bunch_pos->slice->bunch);
                                                                                #endif
    }


  private:
    /// \brief Adapt the non-inert transitions in an out-slice to a new block
    /// \details After a block has been split, the outgoing transitions of the
    /// new block need to move to the respective block_bunch-slice of the new
    /// block.
    ///
    /// This function handles all transitions in the out-slice just before
    /// `out_slice_end`, as they all belong to the same block_bunch-slice and
    /// can be moved together.  However, transitions in `splitter_T` are
    /// excepted:  all transitions in `splitter_T` from all states are
    /// transitions of the R-subblock, so if the latter is the new block, then
    /// `splitter_T` can be moved as a whole instead of per-state.  In this
    /// case, the caller should move `splitter_T` to the list of stable
    /// block_bunch-slices of the R-subblock.
    ///
    /// The function returns the beginning of this out-slice (which can become
    /// the next out_slice_end).  It is meant to be called from the last
    /// out-slice back to the first because it inserts stable
    /// block_bunch-slices at the beginning of the list for the new block, so
    /// it would normally become ordered according to the bunch.
    /// \param out_slice_end  The transition just after the out-slice that is
    ///                       adapted
    /// \param old_block      The block in which the source state of the
    ///                       out-slice was before it was split (only needed if
    ///                       we do not use simple lists)
    /// \param splitter_T     The splitter that made this block split
    /// \returns  the beginning of this out-slice (which can become the next
    ///           out_slice_end)
                                                                                ONLY_IF_DEBUG( template <class LTS_TYPE> )
    succ_entry* move_out_slice_to_new_block(
                               succ_entry* out_slice_end,                       ONLY_IF_DEBUG( const bisim_partitioner_dnj<LTS_TYPE>& partitioner, )
                               block_t* const old_block,
                               block_bunch_slice_const_iter_t const splitter_T)
    {                                                                           assert(&succ.cbegin()[1] < out_slice_end);
        succ_entry* const out_slice_begin(
                                        out_slice_end[-1].begin_or_before_end); assert(nullptr != out_slice_begin);
                                                                                assert(out_slice_begin < out_slice_end);
                                                                                assert(out_slice_begin->block_bunch->pred->action_block->succ ==
                                                                                                                                              out_slice_begin);
        block_bunch_entry* old_block_bunch_pos(out_slice_end[-1].block_bunch);  assert(nullptr != old_block_bunch_pos->pred->action_block->succ);
                                                                                assert(!old_block_bunch_pos->slice.is_null());
        block_bunch_slice_iter_t const old_block_bunch_slice(
                                                   old_block_bunch_pos->slice); assert(old_block_bunch_pos->pred->action_block->succ->block_bunch ==
                                                                                                                                          old_block_bunch_pos);
        if (&*splitter_T == &*old_block_bunch_slice)  return out_slice_begin;

        block_bunch_entry* old_block_bunch_slice_end(
                                                   old_block_bunch_slice->end);
        state_info_entry* const source(old_block_bunch_pos->pred->source);      assert(out_slice_end <= source->succ_inert.begin);
                                                                                assert(partitioner.part_st.state_info.data() == source ||
                                                                                                                  source[-1].succ_inert.begin < out_slice_end);
        block_t* const new_block(source->bl.ock);                               assert(source == out_slice_begin->block_bunch->pred->source);
        block_bunch_slice_iter_t new_block_bunch_slice;                         assert(source->pos->st == source);
        if (old_block_bunch_slice_end >= block_bunch_inert_begin ||
            new_block != old_block_bunch_slice_end->pred->source->bl.ock ||
                (new_block_bunch_slice = (block_bunch_slice_iter_t)
                                         old_block_bunch_slice_end->slice,      assert(!old_block_bunch_slice_end->slice.is_null()),
                 old_block_bunch_slice->bunch != new_block_bunch_slice->bunch))
        {
            // the new block_bunch-slice is not suitable; create a new one and
            // insert it into the correct list.
            if (old_block_bunch_slice->is_stable())
            {
                // In most cases, but not always, the source is a bottom state.
                #ifdef USE_SIMPLE_LIST
                    new_block_bunch_slice =
                        new_block->stable_block_bunch.emplace_front(
                                           old_block_bunch_slice->end,
                                           old_block_bunch_slice->bunch, true);
                #else
                    new_block->stable_block_bunch.emplace_front(
                                           old_block_bunch_slice->end,
                                           old_block_bunch_slice->bunch, true);
                    new_block_bunch_slice =
                                         new_block->stable_block_bunch.begin();
                #endif
            }
            else
            {
                #ifdef USE_SIMPLE_LIST
                    new_block_bunch_slice = splitter_list.emplace_after(
                                          old_block_bunch_slice,
                                          old_block_bunch_slice->end,
                                          old_block_bunch_slice->bunch, false);
                #else
                    new_block_bunch_slice = splitter_list.emplace(
                                          std::next(old_block_bunch_slice),
                                          old_block_bunch_slice->end,
                                          old_block_bunch_slice->bunch, false);
                #endif
            }
            ++nr_of_block_bunch_slices;
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    new_block_bunch_slice->work_counter = old_block_bunch_slice->work_counter;
                                                                                #endif
        }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    unsigned max_counter = check_complexity::log_n -
                                                                                                                    check_complexity::ilog2(new_block->size());
                                                                                #endif
        /* move all transitions in this out-slice to the new block_bunch     */ assert(out_slice_begin < out_slice_end);
        do
        {                                                                       assert(old_block_bunch_pos == out_slice_end[-1].block_bunch);
            --out_slice_end;                                                    assert(old_block_bunch_pos->slice == old_block_bunch_slice);
                                                                                assert(source == out_slice_end->block_bunch->pred->source);
            --old_block_bunch_slice_end;                                        // assign work already now because the transition may be moved to several
                                                                                // places:
            old_block_bunch_slice_end->slice = new_block_bunch_slice;           mCRL2complexity(old_block_bunch_pos->pred, add_work(
                                                                                     check_complexity::move_out_slice_to_new_block, max_counter), partitioner);
            if (old_block_bunch_slice->is_stable() || (                         assert(!new_block_bunch_slice->is_stable()),
                 old_block_bunch_slice->marked_begin >
                                                 old_block_bunch_slice_end &&
                 (/* As the old block_bunch-slice has no marked              */ assert(nullptr != old_block_bunch_slice_end),
                  // transitions, it is enough to adapt its marked_begin
                  // and then do a simple (two-way) swap.
                  old_block_bunch_slice->marked_begin =
                                             old_block_bunch_slice_end, true)))
            {
                // The old block_bunch-slice is stable, or it has no
                // marked transitions.
                std::swap(old_block_bunch_pos->pred,
                                              old_block_bunch_slice_end->pred);
            }
            else
            {
                // The old block_bunch-slice is unstable and has marked
                // transitions.
                pred_entry* const old_pred = old_block_bunch_pos->pred;
                if (old_block_bunch_pos < old_block_bunch_slice->marked_begin)
                {
                    // The transition is not marked, but there are other
                    // marked transitions in the old block_bunch-slice.
                    // Move the transition to the non-marked part of the
                    // new block_bunch-slice.
                    block_bunch_entry* const old_marked_begin =
                                       old_block_bunch_slice->marked_begin - 1; assert(old_block_bunch_pos < old_block_bunch_slice_end);
                    old_block_bunch_slice->marked_begin = old_marked_begin;

                    old_block_bunch_pos->pred = old_marked_begin->pred;
                    old_marked_begin->pred = old_block_bunch_slice_end->pred;
                    old_block_bunch_slice_end->pred = old_pred;                 assert(nullptr != old_marked_begin->pred->action_block->succ);

                    old_marked_begin->pred->action_block->succ->
                                                block_bunch = old_marked_begin;
                }
                else
                {
                    // The transition is marked. Move to the marked part
                    // of the new block_bunch-slice.
                    block_bunch_entry* const new_marked_begin =
                                       new_block_bunch_slice->marked_begin - 1;
                    new_block_bunch_slice->marked_begin = new_marked_begin;     assert(old_block_bunch_pos < new_marked_begin ||
                                                                                                             old_block_bunch_pos == old_block_bunch_slice_end);
                    old_block_bunch_pos->pred=old_block_bunch_slice_end->pred;  assert(old_block_bunch_slice_end <= new_marked_begin);
                    old_block_bunch_slice_end->pred = new_marked_begin->pred;
                    new_marked_begin->pred = old_pred;                          assert(out_slice_end == new_marked_begin->pred->action_block->succ);

                    out_slice_end->block_bunch = new_marked_begin;
                }
            }                                                                   assert(nullptr != old_block_bunch_slice_end->pred->action_block->succ);
            old_block_bunch_slice_end->pred->action_block->succ->block_bunch =
                                                     old_block_bunch_slice_end; assert(nullptr != old_block_bunch_pos->pred->action_block->succ);
            old_block_bunch_pos->pred->action_block->succ->block_bunch =
                                                           old_block_bunch_pos;
        }
        while (out_slice_begin < out_slice_end &&
                  (old_block_bunch_pos = out_slice_end[-1].block_bunch, true));
        old_block_bunch_slice->end = old_block_bunch_slice_end;

        if (old_block_bunch_slice->empty())
        {
            if (old_block_bunch_slice->is_stable())
            {
                // If the new block is R, then the old (U) block loses
                // exactly one stable block_bunch-slice, namely the one we
                // just stabilised for (`splitter_T`).  We could perhaps
                // optimize this by moving that slice as a whole to the new
                // block -- perhaps later.
                //
                // If the new block is U, then the old (R) block loses
                // no stable block_bunch-slices if it contains any bottom
                // state.  If it doesn't contain any bottom state, it will
                // definitely keep `splitter_T`, but nothing else can be
                // guaranteed.
                //
                // So old_block_bunch_slice may be deleted, in particular
                // if the new block is U, but not exclusively.
                old_block->stable_block_bunch.erase(old_block_bunch_slice);
            }
            else
            {
                splitter_list.erase(old_block_bunch_slice);
            }
            --nr_of_block_bunch_slices;
        }
        return out_slice_begin;
    }


    /// \brief handle one transition after a block has been split
    /// \details The main task of this method is to move the transition to the
    /// correct place in the action_block slice.
    ///
    /// This function handles phase 1.  Because the new action_block-slice
    /// cannot be adapted completely until all transitions into the new block
    /// have been handled through phase 1, the next function handles them again
    /// in phase 2.
    /// \param pred_iter  transition that has to be moved
    void first_move_transition_to_new_action_block(pred_entry* const pred_iter)
    {
        action_block_entry* const old_action_block_pos(
                                                      pred_iter->action_block); assert(nullptr != old_action_block_pos->succ);
                                                                                assert(old_action_block_pos->succ->block_bunch->pred == pred_iter);
        action_block_entry* const action_block_slice_begin(
                             old_action_block_pos->action_block_slice_begin(    ONLY_IF_DEBUG( action_block.data(), action_block_orig_inert_begin )
                                                                            )); assert(nullptr != action_block_slice_begin->succ);
                                                                                assert(action_block_slice_begin->succ->block_bunch->pred->action_block ==
                                                                                                                                     action_block_slice_begin);
        action_block_entry* const new_action_block_pos(
                                action_block_slice_begin->begin_or_before_end); assert(nullptr != new_action_block_pos);
                                                                                assert(action_block_slice_begin == new_action_block_pos->begin_or_before_end);
                                                                                assert(nullptr != new_action_block_pos->succ);
                                                                                assert(new_action_block_pos->succ->block_bunch->pred->action_block ==
        /* move the transition to the end of the action_block-slice          */                                                          new_action_block_pos);
        if (old_action_block_pos < new_action_block_pos)
        {
            succ_entry* const temp(new_action_block_pos->succ);                 assert(nullptr != temp);  assert(nullptr != old_action_block_pos->succ);
            new_action_block_pos->succ = old_action_block_pos->succ;
            old_action_block_pos->succ = temp;
            temp->block_bunch->pred->action_block = old_action_block_pos;       assert(pred_iter == new_action_block_pos->succ->block_bunch->pred);
            pred_iter->action_block = new_action_block_pos;

        // adapt the old action_block-slice immediately
            action_block_slice_begin->begin_or_before_end =
                                                      new_action_block_pos - 1;
        }
        else
        {                                                                       assert(old_action_block_pos == new_action_block_pos);
            if (action_block_slice_begin < new_action_block_pos)
            {
                // The old action_block-slice is not empty, so we have to adapt
                // the pointer at the beginning.  (If it is empty, it may
                // happen that `new_action_block_pos - 1` is an illegal value.)
                action_block_slice_begin->begin_or_before_end =
                                                      new_action_block_pos - 1;
            }
            else  --nr_of_action_block_slices;
        }                                                                       assert(nullptr != new_action_block_pos->succ);
                                                                                assert(pred_iter == new_action_block_pos->succ->block_bunch->pred);
        // adapt the new action_block-slice, as far as is possible now
            // make the begin_or_before_end pointers of the first and last
            // transition in the slice correct immediately.  The other
            // begin_or_before_end pointers need to be corrected after all
            // transitions in the new bunch have been positioned correctly.
        if (new_action_block_pos + 1 >= action_block_inert_begin ||
            nullptr == new_action_block_pos[1].succ ||
            new_action_block_pos[1].succ->bunch() !=
                                         new_action_block_pos->succ->bunch() ||
            new_action_block_pos[1].succ->block_bunch->pred->target->bl.ock !=
                                                     pred_iter->target->bl.ock)
        {
            // This is the first transition that moves to this new
            // action_block-slice.
            new_action_block_pos->begin_or_before_end = new_action_block_pos;
            ++nr_of_action_block_slices;
        }
        else
        {
            action_block_entry* const action_block_slice_before_end(
                                  new_action_block_pos[1].begin_or_before_end); assert(nullptr != action_block_slice_before_end);
                                                                                assert(new_action_block_pos < action_block_slice_before_end);
                                                                                assert(nullptr != action_block_slice_before_end->succ);
                                                                                assert(action_block_slice_before_end->succ->block_bunch->pred->action_block ==
                                                                                                                                action_block_slice_before_end);
                                                                                assert(new_action_block_pos + 1 ==
                                                                                                           action_block_slice_before_end->begin_or_before_end);
            action_block_slice_before_end->begin_or_before_end =
                                                          new_action_block_pos; assert(action_block_slice_before_end->succ->block_bunch->
                                                                                                            pred->target->bl.ock == pred_iter->target->bl.ock);
            new_action_block_pos->begin_or_before_end =
                                                 action_block_slice_before_end; assert(action_block_slice_before_end < action_block_inert_begin);
        }
    }


    /// \brief handle one transition after a block has been split, phase 2
    /// \details Because the new action_block-slice cannot be adapted
    /// completely until all transitions into the new block have been handled
    /// through phase 1 see the previous function), this function handles them
    /// again in phase 2.
    /// \param pred_iter  transition that has to be moved
    void second_move_transition_to_new_action_block(
                                                   pred_entry* const pred_iter)
    {
        action_block_entry* const new_action_block_pos(
                                                      pred_iter->action_block); assert(nullptr != new_action_block_pos->succ);
                                                                                assert(new_action_block_pos->succ->block_bunch->pred == pred_iter);
        action_block_entry* const old_begin_or_before_end(
                                    new_action_block_pos->begin_or_before_end); assert(nullptr != old_begin_or_before_end);
                                                                                assert(nullptr != old_begin_or_before_end->succ);
                                                                                assert(old_begin_or_before_end->succ->block_bunch->pred->action_block ==
                                                                                                                                      old_begin_or_before_end);
        if (action_block_entry* const new_begin_or_before_end(
                                 old_begin_or_before_end->begin_or_before_end); assert(nullptr != new_begin_or_before_end),
                                                                                assert(nullptr != new_begin_or_before_end->succ),
                                                                                assert(new_begin_or_before_end->succ->block_bunch->pred->action_block ==
                                                                                                                                      new_begin_or_before_end),
                                new_begin_or_before_end < new_action_block_pos)
        {                                                                       assert(old_begin_or_before_end==new_begin_or_before_end->begin_or_before_end);
            new_action_block_pos->begin_or_before_end =
                                                       new_begin_or_before_end; assert(new_action_block_pos <= old_begin_or_before_end);
            return;
        }                                                                       else  assert(new_begin_or_before_end == new_action_block_pos);
        if (old_begin_or_before_end < new_action_block_pos)  return;

        // this is the first transition in the new action_block-slice.
        // Check whether the bunch it belongs to has become nontrivial.
        bunch_t* const bunch(new_action_block_pos->succ->bunch());
        if (!bunch->is_trivial()) { return; }                                    assert(old_begin_or_before_end + 1 == bunch->end);
        if (bunch->begin < new_action_block_pos)
        {
            make_nontrivial(bunch);
        }
    }


    /// \brief adapt data structures for a transition that has become non-inert
    /// \details If the action_block-slice and the block_bunch-slice that
    /// precede the inert transitions in the respective arrays fit, the
    /// transition is added to these arrays instead of creating a new one.
    /// This only works if:
    /// - the action_block-slice has the same target block and the same action
    ///   as old_pred_pos
    /// - the block_bunch-slice has the same source block as old_pred_pos
    /// - the bunch must contain the action_block-slice.
    /// If only the last two conditions are fulfilled, we can start a new
    /// action_block-slice in the same bunch.  (It would be best for this if
    /// the R-subblock's block_bunch-slice would be the new one, because that
    /// would generally allow to add the new non-inert transitions to the last
    /// splitter.)
    ///
    /// The state is only marked if is becomes a new bottom state.  Otherwise,
    /// the marking/unmarking of the state is unchanged.
    /// \param         old_pred_pos the transition that needs to be adapted.
    /// \param[in,out] new_noninert_block_bunch_ptr the bunch where new
    ///                             non-inert transitions have to be stored.
    ///                             If no such bunch has yet been created, it
    ///                             is nullptr; in that case, make_noninert()
    ///                             creates a new bunch.
    /// \returns true iff the state became a new bottom state
    bool make_noninert(pred_entry* const old_pred_pos,
          block_bunch_slice_iter_or_null_t* const new_noninert_block_bunch_ptr)
    {
        state_info_entry* const source(old_pred_pos->source);                   assert(source->pos->st == source);
        state_info_entry* const target(old_pred_pos->target);                   assert(target->pos->st == target);

        action_block_entry* const new_action_block_pos(
                                                   action_block_inert_begin++); assert(nullptr != new_action_block_pos->succ);
                                                                                assert(new_action_block_pos->succ->block_bunch->pred->action_block ==
                                                                                                                                         new_action_block_pos);
        succ_entry* const new_succ_pos(source->succ_inert.begin++);             assert(new_succ_pos->block_bunch->pred->action_block->succ == new_succ_pos);
        block_bunch_entry* const new_block_bunch_pos(
                                                    block_bunch_inert_begin++); assert(nullptr != new_block_bunch_pos->pred->action_block->succ);
                                                                                assert(new_block_bunch_pos->pred->action_block->succ->block_bunch ==
                                                                                                                                          new_block_bunch_pos);
        action_block_entry* const old_action_block_pos(
                                                   old_pred_pos->action_block); assert(new_action_block_pos <= old_action_block_pos);

        succ_entry* const old_succ_pos(old_action_block_pos->succ);             assert(nullptr != old_succ_pos);
        block_bunch_entry* const old_block_bunch_pos(
                                                    old_succ_pos->block_bunch); assert(old_pred_pos == old_block_bunch_pos->pred);
        pred_entry* const new_pred_pos(target->pred_inert.begin++);             assert(nullptr != new_pred_pos->action_block->succ);
                                                                                assert(new_pred_pos->action_block->succ->block_bunch->pred == new_pred_pos);

        /* adapt action_block                                                */ assert(nullptr == new_action_block_pos->begin_or_before_end);
        if (new_action_block_pos < old_action_block_pos)
        {                                                                       assert(nullptr == old_action_block_pos->begin_or_before_end);
            old_action_block_pos->succ = new_action_block_pos->succ;            assert(nullptr != old_action_block_pos->succ);
                                                                                assert(old_action_block_pos->succ->block_bunch->pred->action_block ==
                                                                                                                                         new_action_block_pos);
            old_action_block_pos->succ->block_bunch->pred->action_block =
                                                          old_action_block_pos;
        }                                                                       else  assert(new_action_block_pos == old_action_block_pos);
        new_action_block_pos->succ = new_succ_pos;                              assert(nullptr != new_succ_pos);
        // new_action_block_pos->begin_or_before_end = ...; -- see below

        /* adapt succ                                                        */ assert(nullptr == new_succ_pos->begin_or_before_end);
        if (new_succ_pos < old_succ_pos)
        {                                                                       assert(nullptr == old_succ_pos->begin_or_before_end);
            old_succ_pos->block_bunch = new_succ_pos->block_bunch;              assert(old_succ_pos->block_bunch->pred->action_block->succ == new_succ_pos);
            old_succ_pos->block_bunch->pred->action_block->succ = old_succ_pos; assert(nullptr != old_succ_pos);
        }                                                                       else  assert(new_succ_pos == old_succ_pos);
        new_succ_pos->block_bunch = new_block_bunch_pos;
        // new_succ_pos->begin_or_before_end = ...; -- see below

        /* adapt block_bunch                                                 */ assert(new_block_bunch_pos->slice.is_null());
        if (new_block_bunch_pos < old_block_bunch_pos)
        {                                                                       assert(old_block_bunch_pos->slice.is_null());
                                                                                assert(nullptr != old_block_bunch_pos->pred->action_block->succ);
            old_block_bunch_pos->pred = new_block_bunch_pos->pred;              assert(old_block_bunch_pos->pred->action_block->succ->block_bunch ==
                                                                                                                                          new_block_bunch_pos);
                                                                                assert(nullptr != old_block_bunch_pos->pred->action_block->succ);
            old_block_bunch_pos->pred->action_block->succ->block_bunch =
                                                           old_block_bunch_pos;
        }                                                                       else  assert(new_block_bunch_pos == old_block_bunch_pos);
        new_block_bunch_pos->pred = new_pred_pos;
        // new_block_bunch_pos->slice = ...; -- see below

        // adapt pred
        if (new_pred_pos < old_pred_pos)
        {
            // We need std::swap here to swap the whole content, including
            // work counters in case we measure the work.
            std::swap(*old_pred_pos, *new_pred_pos);                            assert(nullptr != old_pred_pos->action_block->succ);
                                                                                assert(old_pred_pos->action_block->succ->block_bunch->pred == new_pred_pos);
            old_pred_pos->action_block->succ->block_bunch->pred = old_pred_pos;
        }                                                                       else  assert(new_pred_pos == old_pred_pos);
        new_pred_pos->action_block = new_action_block_pos;

        /* make the state a bottom state if necessary                        */ assert(source->bl.ock->nonbottom_begin <= source->pos);
        bool became_bottom(false);                                              assert(succ.back().block_bunch->pred->source != source);
        if (source != source->succ_inert.begin->block_bunch->pred->source)
        {
            block_t* const source_block(source->bl.ock);
            // make the state a marked bottom state
            if (source->pos >= source_block->marked_nonbottom_begin)
            {
                std::swap(*source->pos,
                                      *source_block->marked_nonbottom_begin++);
            }                                                                   assert(source->pos < source_block->marked_nonbottom_begin);
            std::swap(*source->pos, *source_block->nonbottom_begin++);
            ++nr_of_new_bottom_states;
            became_bottom = true;
        }

        bunch_t* new_noninert_bunch;                                            assert(nullptr != new_action_block_pos);
        if (!new_noninert_block_bunch_ptr->is_null())
        {
            // There is already some new non-inert transition from this block.
            // So we can reuse this block_bunch and its bunch.
            // (However, it may be the case that the current transition goes to
            // another block; in the latter case, we have to create a new
            // action_block-slice.)

            // extend the bunch
            new_noninert_bunch = (*new_noninert_block_bunch_ptr)->bunch;        assert(new_action_block_pos >= new_noninert_bunch->end);
                                                                                ONLY_IF_DEBUG(
                                                                                    for (const action_block_entry*temp_action_block_pos=new_action_block_pos ;
                                                                                                            temp_action_block_pos > new_noninert_bunch->end ; )
                                                                                    {
                                                                                        assert(nullptr == (--temp_action_block_pos)->succ);
                                                                                    }
                                                                                )
            new_noninert_bunch->end = action_block_inert_begin;
            /* extend the block_bunch-slice                                  */ assert((*new_noninert_block_bunch_ptr)->end == new_block_bunch_pos);
            (*new_noninert_block_bunch_ptr)->end = block_bunch_inert_begin;
            if (!(*new_noninert_block_bunch_ptr)->is_stable())
            {                                                                   assert((*new_noninert_block_bunch_ptr)->marked_begin == new_block_bunch_pos);
                (*new_noninert_block_bunch_ptr)->marked_begin =
                                                       block_bunch_inert_begin;
            }
            new_block_bunch_pos->slice = *new_noninert_block_bunch_ptr;

            /* adapt the action_block-slice                                  */ assert(new_noninert_bunch->begin < new_action_block_pos);
            if (nullptr != new_action_block_pos[-1].succ &&
                target->bl.ock == new_action_block_pos[-1].
                                       succ->block_bunch->pred->target->bl.ock)
            {
                // the action_block-slice is suitable: extend it
                action_block_entry* const action_block_slice_begin(
                                 new_action_block_pos[-1].begin_or_before_end); assert(nullptr != action_block_slice_begin);
                                                                                assert(new_action_block_pos-1==action_block_slice_begin->begin_or_before_end);
                                                                                assert(nullptr != action_block_slice_begin->succ);
                                                                                assert(action_block_slice_begin->succ->block_bunch->pred->action_block ==
                                                                                                                                     action_block_slice_begin);
                action_block_slice_begin->begin_or_before_end =
                                                          new_action_block_pos;
                new_action_block_pos->begin_or_before_end =
                                                      action_block_slice_begin;
            }
            else
            {
                // create a new action_block-slice
                new_action_block_pos->begin_or_before_end=new_action_block_pos;
                if (new_noninert_bunch->is_trivial())
                {                                                               // Only during initialisation, it may happen that we add new non-inert
                    make_nontrivial(new_noninert_bunch);                        // transitions to a nontrivial bunch:
                }
                                                                                #ifndef NDEBUG
                                                                                    else
                                                                                    {
                                                                                        // We make sure that new_noninert_bunch is the first bunch in
                                                                                        // action_block (and because it's always the last one, it will be the
                                                                                        // only one, so there is only one bunch, as ).
                                                                                        for (const action_block_entry* iter = action_block.data();
                                                                                                                      iter < new_noninert_bunch->begin; ++iter)
                                                                                        {
                                                                                            assert(nullptr == iter->succ);
                                                                                            assert(nullptr == iter->begin_or_before_end);
                                                                                        }
                                                                                    }
                                                                                #endif
                ++nr_of_action_block_slices;
            }

            /* adapt the out-slice                                           */ assert(source != succ.front().block_bunch->pred->source);
            if (source == new_succ_pos[-1].block_bunch->pred->source &&
                                new_succ_pos[-1].bunch() == new_noninert_bunch)
            {
                // the out-slice is suitable: extend it.
                succ_entry* const out_slice_begin(
                                         new_succ_pos[-1].begin_or_before_end); assert(nullptr != out_slice_begin);
                                                                                assert(new_succ_pos - 1 == out_slice_begin->begin_or_before_end);
                out_slice_begin->begin_or_before_end = new_succ_pos;            assert(out_slice_begin->block_bunch->pred->action_block->succ ==
                                                                                                                                              out_slice_begin);
                new_succ_pos->begin_or_before_end = out_slice_begin;
                return became_bottom;
            }
        }
        else
        {
            // create a new bunch for noninert transitions
            new_noninert_bunch =
                #ifdef USE_POOL_ALLOCATOR
                    simple_list<block_bunch_slice_t>::get_pool().
                    template construct<bunch_t>
                #else
                    new bunch_t
                #endif
                              (new_action_block_pos, action_block_inert_begin);
            ++nr_of_bunches;

            // create a new block_bunch-slice
            #ifdef USE_SIMPLE_LIST
                block_bunch_slice_iter_t new_noninert_block_bunch(
                    splitter_list.emplace_back(
                          block_bunch_inert_begin, new_noninert_bunch, false));
            #else
                splitter_list.emplace_back(block_bunch_inert_begin,
                                                    new_noninert_bunch, false);
                block_bunch_slice_iter_t new_noninert_block_bunch(
                                               std::prev(splitter_list.end()));
            #endif
            ++nr_of_block_bunch_slices;
            new_block_bunch_pos->slice = new_noninert_block_bunch;
            *new_noninert_block_bunch_ptr = new_noninert_block_bunch;

            // create a new action_block-slice
            new_action_block_pos->begin_or_before_end = new_action_block_pos;
            ++nr_of_action_block_slices;
        }                                                                       assert(&succ.cbegin()[1] == new_succ_pos ||
                                                                                                new_succ_pos[-1].block_bunch->pred->source < source ||
        /* create a new out-slice                                            */                                new_succ_pos[-1].bunch() != new_noninert_bunch);
        new_succ_pos->begin_or_before_end = new_succ_pos;
        return became_bottom;
    }


  public:
    /// \brief Split all data structures after a new block has been created
    /// \details This function splits the block_bunch- and action_block-slices
    /// to reflect that some transitions now start or end in the new block.
    /// They can no longer be in the same slice as the transitions that start
    /// or end in the old block, respectively.  It also marks the transitions
    /// that have become non-inert as such and finds new bottom states.
    ///
    /// Its time complexity is O(1 + |in(new_block)| + |out(new_block)|).
    /// \param new_block       the new block
    /// \param old_block       the old block (from which new_block was split
    ///                        off)
    /// \param add_new_noninert_to_splitter  indicates to which
    ///                        block_bunch-slice new non-inert transitions
    ///                        should be added:  if this parameter is `false`,
    ///                        a new slice is created;  if it is `true`,
    ///                        new non-inert transitions are added to
    ///                        `splitter_T`.  The latter can be done iff
    ///                        `splitter_T` is the last block_bunch-slice
    /// \param splitter_T      the splitter that caused new_block and old_block
    ///                        to separate from each other
    /// \param new_block_mode  indicates whether the new block is U or R
                                                                                ONLY_IF_DEBUG( template<class LTS_TYPE> )
    void adapt_transitions_for_new_block(
                  block_t* const new_block,
                  block_t* const old_block,                                     ONLY_IF_DEBUG( const bisim_partitioner_dnj<LTS_TYPE>& partitioner, )
                  bool const add_new_noninert_to_splitter,
                  const block_bunch_slice_iter_t splitter_T,
                  enum new_block_mode_t const new_block_mode)
    {                                                                           assert(splitter_T->is_stable());
        // We begin with a bottom state so the new block gets a sorted list of
        // stable block_bunch-slices.
        permutation_entry* s_iter(new_block->begin);                            assert(s_iter < new_block->end);
        do
        {
            state_info_entry* const s(s_iter->st);                              assert(new_block == s->bl.ock);
                                                                                assert(s->pos == s_iter);
            /* -  -  -  -  -  -  adapt part_tr.block_bunch  -  -  -  -  -  - */
                                                                                assert(s != succ.front().block_bunch->pred->source);
            for (succ_entry* succ_iter(s->succ_inert.begin);
                                s == succ_iter[-1].block_bunch->pred->source; )
            {
                succ_iter = move_out_slice_to_new_block(succ_iter,              ONLY_IF_DEBUG( partitioner, )
                                                        old_block, splitter_T); assert(succ_iter->block_bunch->pred->action_block->succ == succ_iter);
                                                                                assert(s == succ_iter->block_bunch->pred->source);
                                                                                // add_work_to_out_slice(succ_iter, ...) -- subsumed in the call below
            }

            /*-  -  -  -  -  -  adapt part_tr.action_block  -  -  -  -  -  -*/
                                                                                assert(s != pred.front().target);
            for (pred_entry* pred_iter(s->pred_inert.begin);
                                                  s == (--pred_iter)->target; )
            {                                                                   assert(pred.data() < pred_iter);
                                                                                assert(nullptr != pred_iter->action_block->succ);
                                                                                assert(pred_iter->action_block->succ->block_bunch->pred == pred_iter);
                first_move_transition_to_new_action_block(pred_iter);           // mCRL2complexity(pred_iter, ...) -- subsumed in the call below
            }                                                                   // mCRL2complexity(s, ...) -- subsumed in the call at the end
        }
        while (++s_iter < new_block->end);

        if (new_block_is_R == new_block_mode)
        {                                                                       assert(splitter_T->source_block() == new_block);
            // The `splitter_T` slice moves completely from the old to the new
            // block.  We move it as a whole to the new block_bunch list.
            new_block->stable_block_bunch.splice(
                    new_block->stable_block_bunch.begin(),
                                    old_block->stable_block_bunch, splitter_T);
        }                                                                       else  assert(splitter_T->source_block() == old_block);

        // We cannot join the loop above with the one below because transitions
        // in the action_block-slices need to be handled in two phases.

        s_iter = new_block->begin;                                              assert(s_iter < new_block->end);
        do
        {
            state_info_entry* const s(s_iter->st);                              assert(s->pos == s_iter);  assert(s != pred.front().target);
            for (pred_entry* pred_iter(s->pred_inert.begin);
                                                  s == (--pred_iter)->target; )
            {                                                                   assert(pred.data() < pred_iter);
                                                                                assert(nullptr != pred_iter->action_block->succ);
                                                                                assert(pred_iter->action_block->succ->block_bunch->pred == pred_iter);
                second_move_transition_to_new_action_block(pred_iter);          // mCRL2complexity(pred_iter, ...) -- subsumed in the call below
            }                                                                   // mCRL2complexity(s, ...) -- subsumed in the call at the end
        }
        while (++s_iter < new_block->end);
                                                                                assert(0 == new_block->marked_size());  assert(0 == old_block->marked_size());
        /* -  -  -  -  -  - find new non-inert transitions -  -  -  -  -  - */  assert(block_bunch.data_end() - block_bunch_inert_begin ==
                                                                                                           action_block.data_end() - action_block_inert_begin);
        if (block_bunch_inert_begin < block_bunch.data_end())
        {
            block_bunch_slice_iter_or_null_t new_noninert_block_bunch;
            if (add_new_noninert_to_splitter)
            {
                new_noninert_block_bunch = splitter_T;
            }
            else
            {
                new_noninert_block_bunch = nullptr;
            }
            if (new_block_is_U == new_block_mode)
            {                                                                   assert(old_block == new_block->end->st->bl.ock);
                                                                                assert(new_block->end < partitioner.part_st.permutation.data_end());
                permutation_entry* target_iter(new_block->begin);               assert(target_iter < new_block->end);
                do
                {
                    state_info_entry* const s(target_iter->st);                 assert(s->pos == target_iter);
                    // check all incoming inert transitions of s, whether they
                    /* still start in new_block                              */ assert(s != pred.back().target);
                    for (pred_entry* pred_iter(s->pred_inert.begin);
                                           s == pred_iter->target; ++pred_iter)
                    {                                                           assert(pred_iter < &pred.back());
                                                                                assert(nullptr != pred_iter->action_block->succ);
                        state_info_entry* const t(pred_iter->source);           assert(pred_iter->action_block->succ->block_bunch->pred == pred_iter);
                                                                                assert(t->pos->st == t);
                        if (new_block != t->bl.ock)
                        {                                                       assert(old_block == t->bl.ock);
                            if (!make_noninert(pred_iter,
                                                    &new_noninert_block_bunch))
                                       // make_noninert() may modify *pred_iter
                            {
                                old_block->mark_nonbottom(t->pos);
                            }
                        }                                                       // mCRL2complexity(old value of *pred_iter, ...) -- overapproximated by the
                                                                                // call below
                    }                                                           // mCRL2complexity(s, ...) -- subsumed in the call at the end
                }
                while (++target_iter < new_block->end);                         assert(0 < old_block->bottom_size());
            }
            else
            {                                                                   assert(new_block_is_R == new_block_mode);
                    /* We have to be careful because make_noninert may move  */ assert(partitioner.part_st.permutation.data() < new_block->begin);
                    /* a state either forward (to the marked states) or      */ assert(old_block == new_block->begin[-1].st->bl.ock);
                    /* back (to the bottom states).                          */ assert(0 < old_block->bottom_size());
                for(permutation_entry* source_iter(new_block->nonbottom_begin);
                             source_iter < new_block->marked_nonbottom_begin; )
                {
                    state_info_entry* const s(source_iter->st);                 assert(s->pos == source_iter);
                    // check all outgoing inert transitions of s, whether they
                    /* still end in new_block.                               */ assert(succ.back().block_bunch->pred->source != s);
                    succ_entry* succ_iter(s->succ_inert.begin);                 assert(succ_iter < &succ.back());
                    bool dont_mark(true);                                       assert(s == succ_iter->block_bunch->pred->source);
                    do
                    {                                                           assert(succ_iter->block_bunch->pred->action_block->succ == succ_iter);
                        if (new_block !=
                                  succ_iter->block_bunch->pred->target->bl.ock)
                        {                                                       assert(old_block == succ_iter->block_bunch->pred->target->bl.ock);
                            dont_mark = make_noninert(
                                                  succ_iter->block_bunch->pred,
                                                    &new_noninert_block_bunch);
                        }                                                       // mCRL2complexity(succ_iter->block_bunch->pred, ...) -- overapproximated by
                                                                                // the call below
                    }
                    while (s == (++succ_iter)->block_bunch->pred->source);
                    if (dont_mark)  ++source_iter;
                    else
                    {                                                           assert(s->pos == source_iter);
                        new_block->mark_nonbottom(source_iter);
                    }                                                           assert(new_block->nonbottom_begin <= source_iter);
                                                                                // mCRL2complexity(s, ...) -- overapproximated by the call at the end
                }
            }
        }                                                                       else  assert(block_bunch_inert_begin == block_bunch.data_end());
                                                                                mCRL2complexity(new_block, add_work(check_complexity::
                                                                                    adapt_transitions_for_new_block, check_complexity::log_n -
                                                                                                     check_complexity::ilog2(new_block->size())), partitioner);
    }
                                                                                #ifndef NDEBUG
                                                                                    /// \brief print all transitions
                                                                                    /// \details Transitions are printed organised into bunches.
                                                                                    template <class LTS_TYPE>
                                                                                    void print_trans(const bisim_partitioner_dnj<LTS_TYPE>& partitioner) const
                                                                                    {
                                                                                        if (!mCRL2logEnabled(log::debug))  return;
                                                                                        // print all outgoing transitions grouped per successor and out-slice
                                                                                        const succ_entry* succ_iter(&succ.cbegin()[1]);
                                                                                        if (succ_iter >= &succ.back())
                                                                                        {
                                                                                            mCRL2log(log::debug) << "No transitions.\n";
                                                                                            return;
                                                                                        }
                                                                                        const state_info_entry* source(succ_iter->block_bunch->pred->source);
                                                                                        mCRL2log(log::debug) << source->debug_id(partitioner) << ":\n";
                                                                                        if (succ_iter->block_bunch->slice.is_null())
                                                                                        {
                                                                                            mCRL2log(log::debug) << "\tInert successors:\n";
                                                                                        }
                                                                                        block_bunch_slice_iter_or_null_t current_out_bunch(nullptr);
                                                                                        do
                                                                                        {
                                                                                            bool always_print=false;
                                                                                            if (source != succ_iter->block_bunch->pred->source)
                                                                                            {   assert(source < succ_iter->block_bunch->pred->source);
                                                                                                source = succ_iter->block_bunch->pred->source;
                                                                                                mCRL2log(log::debug)
                                                                                                    << source->debug_id(partitioner) << ":\n";
                                                                                                always_print=true;
                                                                                            }
                                                                                            if (always_print ||
                                                                                                succ_iter->block_bunch->slice != current_out_bunch)
                                                                                            {   //assert(!current_out_bunch.is_null());
                                                                                                if (succ_iter->block_bunch->slice.is_null())
                                                                                                {   assert(succ_iter == source->succ_inert.begin);
                                                                                                    mCRL2log(log::debug)<<"\tInert successors:\n";
                                                                                                    current_out_bunch = nullptr;
                                                                                                }
                                                                                                else
                                                                                                {   assert(succ_iter < source->succ_inert.begin);
                                                                                                    //assert(!current_out_bunch.is_null());
                                                                                                    //assert(current_out_bunch == splitter_list.end() ||
                                                                                                    //           current_out_bunch->bunch != succ_iter->bunch());
                                                                                                    mCRL2log(log::debug) << "\tSuccessors in "
                                                                                                      <<succ_iter->bunch()->debug_id_short(partitioner)<<":\n";
                                                                                                    current_out_bunch = succ_iter->block_bunch->slice;
                                                                                                }
                                                                                            }
                                                                                            mCRL2log(log::debug) << "\t\t"
                                                                                                << succ_iter->block_bunch->pred->debug_id(partitioner) << '\n';
                                                                                        }
                                                                                        while (++succ_iter < &succ.back());

                                                                                        // print all transitions grouped per bunch and action_block-slice
                                                                                        const action_block_entry* action_block_iter(action_block.data());
                                                                                        do  assert(action_block_iter < action_block.data_end());
                                                                                        while (nullptr == action_block_iter->succ &&
                                                                                                   (assert(nullptr == action_block_iter->begin_or_before_end),
                                                                                                    ++action_block_iter, true));
                                                                                        do
                                                                                        {
                                                                                            const action_block_entry* bunch_end;
                                                                                            const action_block_entry* action_block_slice_end;
                                                                                            assert(nullptr != action_block_iter->succ);
                                                                                            if (action_block_iter->succ->block_bunch->slice.is_null())
                                                                                            {   assert(action_block_iter == action_block_inert_begin);
                                                                                                mCRL2log(log::debug) <<"Inert transition slice [";
                                                                                                action_block_slice_end = bunch_end = action_block.data_end();
                                                                                            }
                                                                                            else
                                                                                            {
                                                                                                const bunch_t* const bunch(action_block_iter->succ->bunch());
                                                                                                assert(nullptr != bunch);
                                                                                                mCRL2log(log::debug) << bunch->debug_id_short(
                                                                                                                   partitioner) << ":\n\taction_block-slice [";
                                                                                                assert(bunch->begin == action_block_iter);
                                                                                                bunch_end = bunch->end;
                                                                                                assert(bunch_end <= action_block_inert_begin);
                                                                                                assert(nullptr != action_block_iter->begin_or_before_end);
                                                                                                action_block_slice_end =
                                                                                                                  action_block_iter->begin_or_before_end + 1;
                                                                                            }
                                                                                            assert(action_block_slice_end <= bunch_end);
                                                                                            // for all action_block-slices in bunch
                                                                                            for (;;)
                                                                                            {
                                                                                                mCRL2log(log::debug) << (action_block_iter -
                                                                                                                                    action_block.data()) << ","
                                                                                                   << (action_block_slice_end - action_block.data()) << "):\n";
                                                                                                // for all transitions in the action_block-slice
                                                                                                assert(action_block_iter < action_block_slice_end);
                                                                                                do
                                                                                                {
                                                                                                    assert(nullptr != action_block_iter->succ);
                                                                                                    mCRL2log(log::debug) << "\t\t"
                                                                                                            << action_block_iter->succ->block_bunch->
                                                                                                                           pred->debug_id(partitioner) << '\n';
                                                                                                }
                                                                                                while (++action_block_iter < action_block_slice_end);
                                                                                                // go to next action_block-slice in the same bunch
                                                                                                while (action_block_iter < bunch_end &&
                                                                                                                            nullptr == action_block_iter->succ)
                                                                                                {
                                                                                                    assert(nullptr == action_block_iter->begin_or_before_end);
                                                                                                    ++action_block_iter;
                                                                                                    assert(action_block_iter < bunch_end);
                                                                                                }
                                                                                                if (action_block_iter >= bunch_end)  break;
                                                                                                assert(nullptr != action_block_iter->begin_or_before_end);
                                                                                                action_block_slice_end =
                                                                                                                  action_block_iter->begin_or_before_end + 1;
                                                                                                mCRL2log(log::debug) << "\taction_block-slice [";
                                                                                            }
                                                                                            // go to next bunch
                                                                                            assert(action_block_iter == bunch_end);
                                                                                            while (action_block_iter < action_block.data_end() &&
                                                                                                                            nullptr == action_block_iter->succ)
                                                                                            {
                                                                                                assert(nullptr == action_block_iter->begin_or_before_end);
                                                                                                ++action_block_iter;
                                                                                            }
                                                                                        }
                                                                                        while (action_block_iter < action_block.data_end());
                                                                                    }
                                                                                #endif
};


/// \brief refine a block
/// \details This function is called after a refinement function has found
/// where to split the block into unmarked (U) and marked (R) states.
/// It creates a new block for the smaller subblock.
/// \param  new_block_mode  indicates whether the U- or the R-block should be
///                         the new one.  (This parameter is necessary in case
///                         the two halves have exactly the same size.)
/// \param  new_seqnr       is the sequence number of the new block
/// \returns pointer to the new block
                                                                                ONLY_IF_DEBUG( template<class LTS_TYPE> )
inline block_t* block_t::split_off_block(
        enum new_block_mode_t const new_block_mode,                             ONLY_IF_DEBUG( const bisim_partitioner_dnj<LTS_TYPE>& partitioner, )
                                                    state_type const new_seqnr)
{                                                                               assert(0 < marked_size());  assert(0 < unmarked_bottom_size());
    // create a new block
    block_t* new_block;
    state_type swapcount(std::min(marked_bottom_size(),
                                                   unmarked_nonbottom_size()));
    if (permutation_entry* const splitpoint(marked_bottom_begin +
                                                    unmarked_nonbottom_size()); assert(begin < splitpoint),  assert(splitpoint < end),
                                                                                assert(splitpoint->st->pos == splitpoint),
                                              new_block_is_U == new_block_mode)
    {                                                                           assert((state_type) (splitpoint - begin) <= size()/2);
        new_block =
                    #ifdef USE_POOL_ALLOCATOR
                        simple_list<block_bunch_slice_t>::get_pool().
                        template construct<block_t>
                    #else
                        new block_t
                    #endif
                                                (begin, splitpoint, new_seqnr);
        new_block->nonbottom_begin = marked_bottom_begin;

        // adapt the old block: it only keeps the R-states
        begin = splitpoint;
        nonbottom_begin = marked_nonbottom_begin;
    }
    else
    {                                                                           assert(new_block_is_R == new_block_mode);
        new_block =
                    #ifdef USE_POOL_ALLOCATOR
                        simple_list<block_bunch_slice_t>::get_pool().
                        template construct<block_t>
                    #else
                        new block_t
                    #endif
                                                  (splitpoint, end, new_seqnr);
        new_block->nonbottom_begin = marked_nonbottom_begin;                    assert((state_type) (end - splitpoint) <= size()/2);

        // adapt the old block: it only keeps the U-states
        end = splitpoint;
        nonbottom_begin = marked_bottom_begin;
    }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
    /* swap contents                                                         */     new_block->work_counter = work_counter;
                                                                                #endif
    // The structure of a block is
    // |  unmarked  |   marked   |  unmarked  |   marked   |
    // |   bottom   |   bottom   | non-bottom | non-bottom |
    // We have to swap the marked bottom with the unmarked non-bottom
    // states.
    //
    // It is not necessary to reset the untested_to_U counters; these
    // counters are anyway only valid for the states in the respective
    // slice.

    if (0 < swapcount)
    {
        // vector swap the states:
        permutation_entry* pos1(marked_bottom_begin);
        permutation_entry* pos2(marked_nonbottom_begin);                        assert(pos1 < pos2);
        permutation_entry const temp(std::move(*pos1));
        for (;;)
        {
            --pos2;                                                             assert(pos1 < pos2);
            *pos1 = std::move(*pos2);
            ++pos1;
            if (0 >= --swapcount) { break; }                                    assert(pos1 < pos2);
            *pos2 = std::move(*pos1);                                           // mCRL2complexity(new_block_is_U == new_block_mode ? pos1[-1] : *pos2, ...)
        }                                                                       // -- overapproximated by the call at the end
        *pos2 = std::move(temp);                                                // mCRL2complexity(new_block_is_U == new_block_mode ? pos1[-1] : *pos2, ...)
    }                                                                           // -- overapproximated by the call at the end
                                                                                #ifndef NDEBUG
                                                                                    { const permutation_entry* s_iter(begin);  assert(s_iter < end);
                                                                                    do  assert(s_iter->st->pos == s_iter);
                                                                                    while (++s_iter < end); }
                                                                                #endif
    // unmark all states in both blocks
    marked_nonbottom_begin = end;
    marked_bottom_begin = nonbottom_begin;
    new_block->marked_bottom_begin = new_block->nonbottom_begin;                assert(new_block->size() <= size());

    /* set the block pointer of states in the new block                      */ assert(new_block->marked_nonbottom_begin == new_block->end);
    permutation_entry* s_iter(new_block->begin);                                assert(s_iter < new_block->end);
    do
    {                                                                           assert(s_iter->st->pos == s_iter);
        s_iter->st->bl.ock = new_block;                                         // mCRL2complexity (*s_iter, ...) -- subsumed in the call below
    }
    while (++s_iter < new_block->end);                                          mCRL2complexity(new_block, add_work(check_complexity::split_off_block,
                                                                                       check_complexity::log_n - check_complexity::ilog2(new_block->size())),
                                                                                                                                                  partitioner);
    return new_block;
}


/// \brief split off a single action_block-slice from the bunch
/// \details The function splits the current bunch after its first
/// action_block-slice or before its last action_block-slice, whichever is
/// smaller.  It creates a new bunch for the split-off slice and returns a
/// pointer to the new bunch.  The caller has to adapt the block_bunch-slices.
/// \param part_tr  the data structure containing information about the
///                 partition of transitions (needed to find the list of
///                 non-trivial bunches)
/// \returns pointer to a new bunch containing one small action_block-slice
///          that was originally in this bunch
inline bunch_t* bunch_t::split_off_small_action_block_slice(
                                                         part_trans_t& part_tr)
{                                                                               assert(begin < end);  assert(nullptr != begin->succ);
                                                                                assert(nullptr != begin->begin_or_before_end);
    action_block_entry* const first_slice_end(begin->begin_or_before_end + 1);  assert(nullptr != end[-1].succ);  assert(nullptr!=end[-1].begin_or_before_end);
    action_block_entry* const last_slice_begin(end[-1].begin_or_before_end);    assert(begin < first_slice_end);  assert(first_slice_end <= last_slice_begin);
    bunch_t* bunch_T_a_Bprime;
        /* Line 2.6: Select some a in Act and B' in Pi_s such that           */ assert(last_slice_begin < end);  assert(nullptr != first_slice_end[-1].succ);
        /*            |T--a-->B'| < 1/2 |T|                                  */ assert(nullptr != last_slice_begin->succ);
    if (first_slice_end - begin > end - last_slice_begin)
    {
        // Line 2.7: Pi_t := Pi_t \ {T} union { T--a-->B', T \ T--a-->B' }
        bunch_T_a_Bprime =
                            #ifdef USE_POOL_ALLOCATOR
                                 simple_list<block_bunch_slice_t>::get_pool().
                                 template construct<bunch_t>
                            #else
                                new bunch_t
                            #endif
                                                       (last_slice_begin, end); assert(nullptr != bunch_T_a_Bprime);
        end = last_slice_begin;
        while (nullptr == end[-1].succ)
        {
            --end;                                                              assert(first_slice_end <= end);  assert(nullptr == end->begin_or_before_end);
        }                                                                       assert(nullptr != end[-1].begin_or_before_end);
        if (first_slice_end == end)  part_tr.make_trivial(this);
    }
    else
    {
        // Line 2.7: Pi_t := Pi_t \ {T} union { T--a-->B', T \ T--a-->B' }
        bunch_T_a_Bprime =
                            #ifdef USE_POOL_ALLOCATOR
                                 simple_list<block_bunch_slice_t>::get_pool().
                                 template construct<bunch_t>
                            #else
                                new bunch_t
                            #endif
                                                      (begin, first_slice_end); assert(nullptr != bunch_T_a_Bprime);
        begin = first_slice_end;
        while (nullptr == begin->succ)
        {                                                                       assert(nullptr == begin->begin_or_before_end);
            ++begin;                                                            assert(begin <= last_slice_begin);
        }                                                                       assert(nullptr != begin->begin_or_before_end);
        if (begin == last_slice_begin)  part_tr.make_trivial(this);
    }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    bunch_T_a_Bprime->work_counter = work_counter;
                                                                                #endif
    ++part_tr.nr_of_bunches;
    return bunch_T_a_Bprime;
}
///@} (end of group part_trans)

} // end namespace bisim_dnj





/* ************************************************************************* */
/*                                                                           */
/*                            A L G O R I T H M S                            */
/*                                                                           */
/* ************************************************************************* */





/// \defgroup part_refine
/// \brief classes to calculate the stutter equivalence quotient of a LTS
///@{



/*=============================================================================
=                                 main class                                  =
=============================================================================*/





/// \class bisim_partitioner_dnj
/// \brief implements the main algorithm for the branching bisimulation
/// quotient
template <class LTS_TYPE>
class bisim_partitioner_dnj
{
  private:
    /// \brief modes that determine details of how split() should work
    enum refine_mode_t{extend_from_marked_states,
                       extend_from_marked_states__add_new_noninert_to_splitter,
                       extend_from_splitter };

    /// \brief automaton that is being reduced
    LTS_TYPE& aut;
                                                                                ONLY_IF_DEBUG( public: )
    /// \brief partition of the state space into blocks
    bisim_dnj::part_state_t part_st;

    /// \brief partitions of the transitions (with bunches and
    /// action_block-slices)
    bisim_dnj::part_trans_t part_tr;
  private:
    /// \brief action label slices
    /// \details In part_tr.action_block, no information about the action label
    /// is actually stored with the transitions, to save memory.  Entry l of
    /// this array contains a pointer to the first entry in
    /// part_tr.action_block with label l.
    ///
    /// During initialisation, entry l of this array contains a counter to
    /// indicate how many non-inert transitions with action label l have been
    /// found.
    fixed_vector<bisim_dnj::iterator_or_counter<
                                bisim_dnj::action_block_entry*> > action_label;
                                                                                ONLY_IF_DEBUG( public: )
    /// \brief true iff branching (not strong) bisimulation has been requested
    bool const branching;
  private:
    /// \brief true iff divergence-preserving branching bisimulation has been
    /// requested
    /// \details Note that this field must be false if strong bisimulation has
    /// been requested.  There is no such thing as divergence-preserving strong
    /// bisimulation.
    bool const preserve_divergence;
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    friend class bisim_dnj::pred_entry;
                                                                                    friend class bisim_dnj::bunch_t;
                                                                                #endif
  public:
    /// time measurement after the end of create_initial_partition()
    std::clock_t end_initial_part;

    /// \brief constructor
    /// \details The constructor constructs the data structures and immediately
    /// calculates the partition corresponding with the bisimulation quotient.
    /// It destroys the transitions on the LTS (to save memory) but does not
    /// adapt the LTS to represent the quotient's transitions.
    /// \param new_aut                 LTS that needs to be reduced
    /// \param new_branching           If true branching bisimulation is used,
    ///                                otherwise strong bisimulation is
    ///                                applied.
    /// \param new_preserve_divergence If true and branching is true, preserve
    ///                                tau loops on states.
    bisim_partitioner_dnj(LTS_TYPE& new_aut, bool const new_branching = false,
                                    bool const new_preserve_divergence = false)
      : aut(new_aut),
        part_st(new_aut.num_states()),
        part_tr(new_aut.num_transitions(), new_aut.num_action_labels()),
        action_label(new_aut.num_action_labels()),
        branching(new_branching),
        preserve_divergence(new_preserve_divergence)
    {                                                                           assert(branching || !preserve_divergence);

        mCRL2log(log::verbose) << "Start initialisation.\n";
        create_initial_partition();                                             ONLY_IF_DEBUG( part_tr.action_block_orig_inert_begin =
                                                                                                                       part_tr.action_block_inert_begin; )
        end_initial_part = std::clock();
        mCRL2log(log::verbose) << "Start refining\n";
        refine_partition_until_it_becomes_stable();
    }


    /// \brief Calculate the number of equivalence classes
    /// \details The number of equivalence classes (which is valid after the
    /// partition has been constructed) is equal to the number of states in the
    /// bisimulation quotient.
    state_type num_eq_classes() const
    {
        return part_st.nr_of_blocks;
    }


    /// \brief Get the equivalence class of a state
    /// \details After running the minimisation algorithm, this function
    /// produces the number of the equivalence class of a state.  This number
    /// is the same as the number of the state in the minimised LTS to which
    /// the original state is mapped.
    /// \param s state whose equivalence class needs to be found
    /// \returns sequence number of the equivalence class of state s
    state_type get_eq_class(state_type const s) const
    {
        return part_st.block(s)->seqnr;
    }


    /// \brief Adapt the LTS after minimisation
    /// \details After the efficient branching bisimulation minimisation, the
    /// information about the quotient LTS is only stored in the partition data
    /// structure of the partitioner object.  This function exports the
    /// information back to the LTS by adapting its states and transitions:  it
    /// updates the number of states and adds those transitions that are
    /// mandated by the partition data structure.  If desired, it also creates
    /// a vector containing an arbritrary (example) original state per
    /// equivalence class.
    ///
    /// The main parameter and return value are implicit with this function: a
    /// reference to the LTS was stored in the object by the constructor.
    void finalize_minimized_LTS()
    {
        // The labels have already been stored in
        // next_nontrivial_and_label.label by
        // refine_partition_until_it_becomes_stable().

        // for all blocks
        const bisim_dnj::permutation_entry* s_iter(part_st.permutation.data()); assert(s_iter < part_st.permutation.data_end());
        do
        {
            const bisim_dnj::block_t* const B(s_iter->st->bl.ock);
            // for all block_bunch-slices of the block
            for (const bisim_dnj::block_bunch_slice_t& block_bunch :
                                                         B->stable_block_bunch)
            {                                                                   assert(block_bunch.is_stable());  assert(!block_bunch.empty());
                const bisim_dnj::pred_entry* const
                                                pred(block_bunch.end[-1].pred); assert(pred->source->bl.ock == B);
                                                                                assert(nullptr != pred->action_block->succ);
                /* add a transition from the source block to the goal block  */ assert(pred->action_block->succ->block_bunch->pred == pred);
                /* with the indicated label.                                 */ assert(pred->action_block->succ->block_bunch->slice == &block_bunch);
                label_type const
                     label(block_bunch.bunch->next_nontrivial_and_label.label); assert(0 <= label);  assert(label < action_label.size());
                aut.add_transition(transition(B->seqnr, label,
                                                 pred->target->bl.ock->seqnr));
            }
            s_iter = B->end;
        }
        while (s_iter < part_st.permutation.data_end());

        // Merge the states, by setting the state labels of each state to the
        // concatenation of the state labels of its equivalence class.

        if (aut.has_state_info())   /* If there are no state labels
                                                    this step can be ignored */
        {
            /* Create a vector for the new labels */
            typename std::remove_reference<decltype(aut.state_labels())>::type
                                                  new_labels(num_eq_classes());

            state_type i(0);                                                    assert(i < aut.num_states());
            do
            {
                const state_type new_index(get_eq_class(i));
                new_labels[new_index]=new_labels[new_index]+aut.state_label(i);
            }
            while (++i < aut.num_states());

            aut.set_num_states(num_eq_classes(), false);                        assert(0 == aut.num_state_labels());
            //m_aut.clear_state_labels();
            new_labels.swap(aut.state_labels());
        }
        else
        {
            aut.set_num_states(num_eq_classes(), false);
        }

        aut.set_initial_state(get_eq_class(aut.initial_state()));
    }


    /// \brief Check whether two states are in the same equivalence class.
    /// \param s first state that needs to be compared.
    /// \param t second state that needs to be compared.
    /// \returns true iff the two states are in the same equivalence class.
    bool in_same_class(state_type const s, state_type const t) const
    {
        return part_st.block(s) == part_st.block(t);
    }
  private:

    /*--------------------------- main algorithm ----------------------------*/

    /// \brief Create a partition satisfying the main invariant
    /// \details Before the actual bisimulation minimisation can start, this
    /// function needs to be called to create a partition that satisfies the
    /// main invariant of the efficient O(m log n) branching bisimulation
    /// minimisation.
    ///
    /// It puts all non-inert transitions into a single bunch, containing one
    /// action_block-slice for each action label.  It creates a single block
    /// (or possibly two, if there are states that never will do any visible
    /// action).  As a side effect, it deletes all transitions from the LTS
    /// that is stored with the partitioner;  information about the transitions
    /// is kept in data structures that are suitable for the efficient
    /// algorithm.
    ///
    /// For divergence-preserving branching bisimulation, we only need to treat
    /// tau-self-loops as non-inert transitions.  In other texts, this is
    /// sometimes described as temporarily renaming the tau-self-loops to
    /// self-loops with a special label.  However, as there are no other
    /// non-inert tau transitions, we can simply put them in their own
    /// action_block-slice, separate from the inert tau transitions.  (It would
    /// be an error to mix the inert transitions with the self-loops in the
    /// same slice.)
    void create_initial_partition()
    {
        mCRL2log(log::verbose) << "An O(m log n) "
             << (branching ? (preserve_divergence
                                           ? "divergence-preserving branching "
                                           : "branching ")
                           : "")
             << "bisimulation partitioner created for " << part_st.state_size()
             << " states and " << aut.num_transitions() << " transitions.\n";

        if (part_st.state_size() > 2 * aut.num_transitions() + 1)
        {
            mCRL2log(log::warning) << "There are several isolated states "
                "without incoming or outgoing transition. It is not "
                "guaranteed that branching bisimulation minimisation runs in "
                "time O(m log n).\n";
        }

        sort_transitions(aut.get_transitions(), tgt_lbl_src);
        mCRL2log(log::verbose) << "Carried out sorting\n";
        // create one block for all states
        bisim_dnj::block_t* B(
                #ifdef USE_POOL_ALLOCATOR
                    simple_list<bisim_dnj::block_bunch_slice_t>::get_pool().
                    template construct<bisim_dnj::block_t>
                #else
                    new bisim_dnj::block_t
                #endif
                     (part_st.permutation.data(),
                      part_st.permutation.data_end(), part_st.nr_of_blocks++));

        // Iterate over the transitions to count how to order them in
        // part_trans_t

        // counters for the non-inert outgoing and incoming transitions per
        // state are provided in part_st.state_info.  These counters have been
        // initialised to zero in the constructor of part_state_t.
        // counters for the non-inert transition per label are stored in
        // action_label.
                                                                                assert(action_label.size() == aut.num_action_labels());
        // counter for the total number of inert transitions:
        trans_type inert_transitions(0);
        for (const transition& t: aut.get_transitions())
        {
            if (branching&&aut.is_tau(aut.apply_hidden_label_map(t.label()))&&  ((
                                                            t.from() != t.to())   || (assert(preserve_divergence), false)))
            {
                // The transition is inert.
                ++part_st.state_info[t.from()].succ_inert.count;
                ++inert_transitions;

                // The source state should become non-bottom:
                if (part_st.state_info[t.from()].pos < B->nonbottom_begin)
                {
                    std::swap(*part_st.state_info[t.from()].pos,
                                                        *--B->nonbottom_begin);
                    // we do not yet update the marked_bottom_begin pointer
                }
            }
            else
            {
                // The transition is non-inert.  (It may be a self-loop).
                ++part_st.state_info[t.from()].untested_to_U_eqv.count;
                ++action_label[aut.apply_hidden_label_map(t.label())].count;
            }
            ++part_st.state_info[t.to()].pred_inert.count;
        }
        // Now we update the marked_bottom_begin pointer:
        B->marked_bottom_begin = B->nonbottom_begin;

        // set the pointers to transition slices in the state info entries

        // We set them all to the end of the respective slice here.  Then, with
        // every transition, the pointer will be reduced by one, so that after
        // placing all transitions it will point to the beginning of the slice.

        bisim_dnj::pred_entry* next_pred_begin(&part_tr.pred.begin()[1]);
        bisim_dnj::succ_entry* next_succ_begin(&part_tr.succ.begin()[1]);
        bisim_dnj::state_info_entry* state_iter(part_st.state_info.data());     assert(state_iter < part_st.state_info.data_end());
        do
        {
            state_iter->bl.ed_noninert_end = next_pred_begin;
            next_pred_begin += state_iter->pred_inert.count;
            state_iter->pred_inert.convert_to_iterator(next_pred_begin);

            // create slice descriptors in part_tr.succ for each state with
            /* outgoing transitions.                                         */ assert(nullptr != next_succ_begin);
            state_iter->untested_to_U_eqv.convert_to_iterator(
                        next_succ_begin + state_iter->untested_to_U_eqv.count);
            if (next_succ_begin < state_iter->untested_to_U_eqv.begin)
            {                                                                   assert(nullptr != state_iter->untested_to_U_eqv.begin);
                next_succ_begin->begin_or_before_end =
                                       state_iter->untested_to_U_eqv.begin - 1;
                for (bisim_dnj::succ_entry* const
                                              out_slice_begin(next_succ_begin);
                     ++next_succ_begin < state_iter->untested_to_U_eqv.begin; )
                {
                    next_succ_begin->begin_or_before_end = out_slice_begin;     // mCRL2complexity(next_succ_begin->block_bunch->pred, ...) -- subsumed in the
                }                                                               // call below

                B->mark(state_iter->pos);
            }
            state_iter->succ_inert.convert_to_iterator(next_succ_begin +
                                                 state_iter->succ_inert.count);
                                                                                #ifndef NDEBUG
                                                                                    while (next_succ_begin < state_iter->succ_inert.begin)
                                                                                    {   assert(nullptr == next_succ_begin->begin_or_before_end);
                                                                                        ++next_succ_begin;
                                                                                    }
                                                                                #endif
            next_succ_begin = state_iter->succ_inert.begin;                     // mCRL2complexity(*state_iter, ...) -- subsumed in the call at the end
        }
        while (++state_iter < part_st.state_info.data_end());

        // Line 2.4: Pi_t := { { all non-inert transitions } }
        part_tr.action_block_inert_begin =
                           part_tr.action_block.data_end() - inert_transitions; assert(part_tr.action_block.data() <= part_tr.action_block_inert_begin);
        part_tr.block_bunch_inert_begin =
                            part_tr.block_bunch.data_end() - inert_transitions; assert(part_tr.block_bunch.data() < part_tr.block_bunch_inert_begin);
        bisim_dnj::bunch_t* bunch(nullptr);

        if (1 + part_tr.block_bunch.data() < part_tr.block_bunch_inert_begin)
        {                                                                       assert(part_tr.action_block.data() < part_tr.action_block_inert_begin);
            // create a single bunch containing all non-inert transitions
            bunch =
                #ifdef USE_POOL_ALLOCATOR
                    simple_list<bisim_dnj::block_bunch_slice_t>::get_pool().
                    template construct<bisim_dnj::bunch_t>
                #else
                    new bisim_dnj::bunch_t
                #endif
                                            (part_tr.action_block.data(),
                                             part_tr.action_block_inert_begin); assert(nullptr != bunch);  assert(part_tr.splitter_list.empty());
            ++part_tr.nr_of_bunches;                                            assert(1 == part_tr.nr_of_bunches);

            // create a single block_bunch entry for all non-inert transitions
            part_tr.splitter_list.emplace_front(
                                part_tr.block_bunch_inert_begin, bunch, false); assert(!part_tr.splitter_list.empty());
            ++part_tr.nr_of_block_bunch_slices;                                 assert(1 == part_tr.nr_of_block_bunch_slices);
        }

        // create slice descriptors in part_tr.action_block for each label

        // The action_block array shall have the tau transitions at the end:
        // first the non-inert tau transitions (during initialisation, that are
        // only the tau self-loops), then the tau transitions that have become
        // non-inert and finally the inert transitions.
        // Transitions with other labels are placed from beginning to end.
        // Every such transition block except the last one ends with a dummy
        /* entry.  If there are transition labels without transitions,       */ assert(part_tr.action_block.size() ==
        /* multiple dummy entries will be placed side-by-side.               */                               aut.num_transitions() + action_label.size() - 1);
        bisim_dnj::action_block_entry*
                          next_action_label_begin(part_tr.action_block.data());
        trans_type const n_square(part_st.state_size() * part_st.state_size());
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    trans_type max_transitions = n_square;
                                                                                #endif
        label_type label(action_label.size());                                  assert(0 < label);
        do
        {
            --label;
            if (0 < action_label[label].count)
            {                                                                   assert(nullptr != bunch);
                if (++part_tr.nr_of_action_block_slices == 2)
                {
                    // This is the second action_block-slice, so the bunch is
                    // not yet marked as nontrivial but it should be.
                    part_tr.make_nontrivial(bunch);
                }
                if (n_square < action_label[label].count)
                {
                    mCRL2log(log::warning) << "There are "
                        << action_label[label].count << ' '
                        << pp(aut.action_label(label)) << "-transitions.  "
                        "This is more than n^2 (= " << n_square << "). It is "
                        "not guaranteed that branching bisimulation "
                        "minimisation runs in time O(m log n).\n";
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    if (max_transitions < action_label[label].count)
                                                                                    {   max_transitions = action_label[label].count;   }
                                                                                #endif
                }
                // initialise begin_or_before_end pointers for this
                // action_block-slice
                action_label[label].convert_to_iterator(
                          next_action_label_begin + action_label[label].count);
                next_action_label_begin->begin_or_before_end =
                                                 action_label[label].begin - 1; assert(nullptr != next_action_label_begin->begin_or_before_end);
                bisim_dnj::action_block_entry* const
                             action_block_slice_begin(next_action_label_begin); assert(nullptr != action_block_slice_begin);
                while (++next_action_label_begin < action_label[label].begin)
                {
                    next_action_label_begin->begin_or_before_end =
                                                      action_block_slice_begin; // mCRL2complexity(next_action_label_begin->succ->block_bunch->pred, ...) --
                }                                                               // subsumed in the call at the end
            }
            else
            {
                action_label[label].convert_to_iterator(
                                                      next_action_label_begin);
                if (0 != label && aut.num_transitions() < action_label.size())
                {
                    mCRL2log(log::warning) << "Action label "
                        << pp(aut.action_label(label)) << " has no "
                        "transitions, and the number of action labels exceeds "
                        "the number of transitions. It is not guaranteed that "
                        "branching bisimulation minimisation runs in time "
                        "O(m log n).\n";
                }
            }
        }
        while (0 < label && (/* insert a dummy entry                         */ assert(next_action_label_begin < part_tr.action_block_inert_begin),
                     next_action_label_begin->succ = nullptr,
                     next_action_label_begin->begin_or_before_end = nullptr,
                                             ++next_action_label_begin, true)); assert(next_action_label_begin == part_tr.action_block_inert_begin);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
        /* distribute the transitions over the data structures               */     check_complexity::init(2 * max_transitions);
                                                                                #endif
        bisim_dnj::block_bunch_entry*
                              next_block_bunch(1 + part_tr.block_bunch.data());
        for (const transition& t: aut.get_transitions())
        {
            bisim_dnj::state_info_entry* const
                                  source(part_st.state_info.data() + t.from());
            bisim_dnj::state_info_entry* const
                                    target(part_st.state_info.data() + t.to());
            bisim_dnj::succ_entry* succ_pos;
            bisim_dnj::block_bunch_entry* block_bunch_pos;
            bisim_dnj::pred_entry* pred_pos;
            bisim_dnj::action_block_entry* action_block_pos;

            if (branching&&aut.is_tau(aut.apply_hidden_label_map(t.label()))&&  ((
                                                            t.from() != t.to())   || (assert(preserve_divergence), false)))
            {
                // It is a (normal) inert transition: place near the end of the
                // respective pred/succ slices, just before the other inert
                // transitions.
                succ_pos = --source->succ_inert.begin;                          assert(nullptr == succ_pos->begin_or_before_end);
                block_bunch_pos = part_tr.block_bunch.data_end() -
                                                             inert_transitions; assert(block_bunch_pos >= part_tr.block_bunch_inert_begin);
                pred_pos = --target->pred_inert.begin;                          assert(block_bunch_pos->slice.is_null());
                action_block_pos = part_tr.action_block.data_end() -
                                                             inert_transitions; assert(action_block_pos >= part_tr.action_block_inert_begin);
                action_block_pos->begin_or_before_end = nullptr;
                --inert_transitions;
            }
            else
            {
                // It is a non-inert transition (possibly a self-loop): place
                // at the end of the respective succ slice and at the beginning
                // of the respective pred slice.
                succ_pos =
                        --part_st.state_info[t.from()].untested_to_U_eqv.begin; assert(nullptr != succ_pos->begin_or_before_end);
                                                                                assert(nullptr != succ_pos->begin_or_before_end->begin_or_before_end);
                                                                                assert(succ_pos->begin_or_before_end <= succ_pos ||
                                                                                               succ_pos->begin_or_before_end->begin_or_before_end == succ_pos);
                block_bunch_pos = next_block_bunch++;                           assert(block_bunch_pos < part_tr.block_bunch_inert_begin);
                pred_pos = target->bl.ed_noninert_end++;
                action_block_pos =
                   --action_label[aut.apply_hidden_label_map(t.label())].begin; assert(nullptr != action_block_pos->begin_or_before_end);
                                                                                assert(nullptr != action_block_pos->begin_or_before_end->begin_or_before_end);
                                                                                assert(action_block_pos->begin_or_before_end <= action_block_pos ||
                                                                                                action_block_pos->begin_or_before_end->
                                                                                                                      begin_or_before_end == action_block_pos);
                                                                                assert(!part_tr.splitter_list.empty());
                block_bunch_pos->slice = part_tr.splitter_list.begin();         assert(action_block_pos < part_tr.action_block_inert_begin);
            }                                                                   assert(target->bl.ed_noninert_end <= target->pred_inert.begin);
            succ_pos->block_bunch = block_bunch_pos;
            block_bunch_pos->pred = pred_pos;
            pred_pos->action_block = action_block_pos;
            pred_pos->source = source;
            pred_pos->target = target;                                          assert(nullptr != succ_pos);
            action_block_pos->succ = succ_pos;                                  // mCRL2complexity(pred_pos, ...) -- subsumed in the call at the end
        }                                                                       assert(0 == inert_transitions);
        /* delete transitions already -- they are no longer needed.  We will */ assert(next_block_bunch == part_tr.block_bunch_inert_begin);
        // add new transitions at the end of minimisation.
        aut.clear_transitions();

        state_iter = part_st.state_info.data();                                 assert(state_iter < part_st.state_info.data_end());
        do
        {
            state_iter->bl.ock = B;
        }
        while (++state_iter < part_st.state_info.data_end());

        if (nullptr != bunch)
        {
            while (nullptr == bunch->begin->succ)
            {                                                                   assert(nullptr == bunch->begin->begin_or_before_end);
                ++bunch->begin;                                                 assert(bunch->begin < bunch->end);
            }                                                                   assert(nullptr != bunch->begin->begin_or_before_end);
            while (nullptr == bunch->end[-1].succ)
            {                                                                   assert(nullptr == bunch->end[-1].begin_or_before_end);
                --bunch->end;                                                   assert(bunch->begin < bunch->end);
            }                                                                   assert(nullptr != bunch->end[-1].begin_or_before_end);

            /* Line 2.2: B_vis := { s in S | there exists a visible          */ mCRL2complexity(B, add_work(check_complexity::
            /*                               transition that is reachable    */                                          create_initial_partition, 1U), *this);
            //                               from s }
            //           B_invis := S \ B_vis
            // Line 2.3: Pi_s := { B_vis, B_invis } \ { emptyset }
                // At this point, all states with a visible transition are
                // marked.
            if (0 < B->marked_size())
            {                                                                   ONLY_IF_DEBUG( part_st.print_part(*this);
                                                                                               part_tr.print_trans(*this); )
                bisim_dnj::block_bunch_slice_iter_t slice =
                                                 part_tr.splitter_list.begin();
                if (1 < B->size())
                {
                    B = split(B, /* splitter block_bunch */ slice,
                      extend_from_marked_states__add_new_noninert_to_splitter);
                    // We can ignore possible new non-inert transitions, as
                    // every R-bottom state already has a transition in bunch.
                    B->marked_nonbottom_begin = B->end;
                }
                else
                {                                                               assert(B->nonbottom_begin == B->end);
                    /* A block with 1 state will not be split.  However, we  */ assert(B->marked_nonbottom_begin == B->end);
                    // still have to make the splitter stable.
                    B->stable_block_bunch.splice(B->stable_block_bunch.end(),
                                                 part_tr.splitter_list, slice);
                    slice->make_stable();
                }                                                               assert(!B->stable_block_bunch.empty());  assert(part_tr.splitter_list.empty());
                                                                                assert(B->stable_block_bunch.front().end <= part_tr.block_bunch_inert_begin);
                                                                                assert(1 + part_tr.block_bunch.data() < B->stable_block_bunch.front().end);
                B->marked_bottom_begin = B->nonbottom_begin;                    assert(!B->stable_block_bunch.front().empty());
            }
        }                                                                       else  assert(0 == B->marked_size());
    }
                                                                                #ifndef NDEBUG
                                                                                    /// \brief assert that the data structure is consistent and stable
                                                                                    /// \details The data structure is tested against a large number of
                                                                                    /// assertions to ensure that everything is consistent, e. g. pointers that
                                                                                    /// should point to successors of state s actually point to a transition
                                                                                    /// that starts in s.
                                                                                    ///
                                                                                    /// Additionally, it is asserted that the partition is stable. i. e. every
                                                                                    /// bottom state in every block can reach exactly every bunch in the list
                                                                                    /// of bunches that should be reachable from it, and every nonbottom state
                                                                                    /// can reach a subset of them.
                                                                                    void assert_stability() const
                                                                                    {
                                                                                        part_st.assert_consistency(*this);

                                                                                        assert(part_tr.succ.size() == part_tr.block_bunch.size() + 1);
                                                                                        assert(part_tr.pred.size() == part_tr.block_bunch.size() + 1);
                                                                                        assert(part_tr.action_block.size() ==
                                                                                                         part_tr.block_bunch.size() + action_label.size() - 2);
                                                                                        if (part_tr.block_bunch.empty())  return;

                                                                                        assert(part_tr.splitter_list.empty());
                                                                                        /* for (const block_bunch_slice_t& block_bunch : part_tr.splitter_list)
                                                                                        {
                                                                                            assert(!block_bunch.is_stable());
                                                                                        } */

                                                                                        trans_type true_nr_of_block_bunch_slices(0);
                                                                                        // for all blocks
                                                                                        const bisim_dnj::permutation_entry*
                                                                                                                         perm_iter(part_st.permutation.data());
                                                                                        assert(perm_iter < part_st.permutation.data_end());
                                                                                        do
                                                                                        {
                                                                                            const bisim_dnj::block_t* const block(perm_iter->st->bl.ock);
                                                                                            unsigned const max_block(check_complexity::log_n -
                                                                                                                       check_complexity::ilog2(block->size()));
                                                                                            // iterators have no predefined hash, so we store pointers:
                                                                                            std::unordered_set<const bisim_dnj::block_bunch_slice_t*>
                                                                                                                                         block_bunch_check_set;
                                                                                            #ifndef USE_SIMPLE_LIST
                                                                                                block_bunch_check_set.reserve(
                                                                                                                             block->stable_block_bunch.size());
                                                                                            #endif

                                                                                            // for all stable block_bunch-slices of the block
                                                                                            for (const bisim_dnj::block_bunch_slice_t& block_bunch :
                                                                                                                                     block->stable_block_bunch)
                                                                                            {
                                                                                                assert(block_bunch.source_block() == block);
                                                                                                assert(block_bunch.is_stable());
                                                                                                block_bunch_check_set.insert(&block_bunch);
                                                                                                mCRL2complexity(&block_bunch, no_temporary_work(
                                                                                                           block_bunch.bunch->max_work_counter(*this)), *this);
                                                                                                ++true_nr_of_block_bunch_slices;
                                                                                            }

                                                                                            // for all states in the block
                                                                                            do
                                                                                            {
                                                                                                trans_type block_bunch_count(0);
                                                                                                const bisim_dnj::state_info_entry* const state(perm_iter->st);
                                                                                                assert(state!=part_tr.succ.front().block_bunch->pred->source);
                                                                                                // for all out-slices of the state
                                                                                                for (const bisim_dnj::succ_entry*
                                                                                                                        out_slice_end(state->succ_inert.begin);
                                                                                                        state == out_slice_end[-1].block_bunch->pred->source; )
                                                                                                {   assert(!out_slice_end[-1].block_bunch->slice.is_null());
                                                                                                    bisim_dnj::block_bunch_slice_const_iter_t const
                                                                                                       block_bunch_slice(out_slice_end[-1].block_bunch->slice);
                                                                                                    const bisim_dnj::bunch_t* const bunch(
                                                                                                                                     block_bunch_slice->bunch);
                                                                                                    assert(block == block_bunch_slice->source_block());
                                                                                                    if (block_bunch_slice->is_stable())
                                                                                                    {
                                                                                                        assert(1 == block_bunch_check_set.count(
                                                                                                                                         &*block_bunch_slice));
                                                                                                        ++block_bunch_count;
                                                                                                    }
                                                                                                    else  assert(0); // i. e. all block_bunch-slices should
                                                                                                                         // be stable
                                                                                                    const bisim_dnj::succ_entry* const out_slice_begin(
                                                                                                                        out_slice_end[-1].begin_or_before_end);
                                                                                                    assert(nullptr != out_slice_begin);
                                                                                                    assert(out_slice_begin < out_slice_end);
                                                                                                    assert(nullptr != out_slice_begin->begin_or_before_end);
                                                                                                    assert(out_slice_begin->begin_or_before_end + 1 ==
                                                                                                                                                out_slice_end);

                                                                                                    // for all transitions in the out-slice
                                                                                                    do
                                                                                                    {
                                                                                                        --out_slice_end;
                                                                                                        assert(bunch->begin <=
                                                                                                               out_slice_end->block_bunch->pred->action_block);
                                                                                                        assert(out_slice_end->block_bunch->pred->
                                                                                                                                    action_block < bunch->end);
                                                                                                        assert(out_slice_end->block_bunch->slice ==
                                                                                                                                            block_bunch_slice);
                                                                                                        assert(nullptr != out_slice_end->begin_or_before_end);
                                                                                                        if (out_slice_end->block_bunch + 1 !=
                                                                                                                                        block_bunch_slice->end)
                                                                                                        {
                                                                                                            assert(out_slice_end->block_bunch + 1 <
                                                                                                                                       block_bunch_slice->end);
                                                                                                            assert(out_slice_end->block_bunch[1].slice ==
                                                                                                                                            block_bunch_slice);
                                                                                                        }
                                                                                                        mCRL2complexity(out_slice_end->block_bunch->pred,
                                                                                                            no_temporary_work(max_block,
                                                                                                                    check_complexity::log_n -
                                                                                                                    check_complexity::ilog2(out_slice_end->
                                                                                                                    block_bunch->pred->target->bl.ock->size()),
                                                                                                                    perm_iter < block->nonbottom_begin),*this);
                                                                                                    }
                                                                                                    while (out_slice_begin < out_slice_end &&
                                                                                                           (assert(out_slice_begin ==
                                                                                                                   out_slice_end->begin_or_before_end), true));
                                                                                                }
                                                                                                if (perm_iter < block->nonbottom_begin)
                                                                                                {
                                                                                                    assert(block_bunch_check_set.size() == block_bunch_count);
                                                                                                }
                                                                                            }
                                                                                            while (++perm_iter < block->end);
                                                                                        }
                                                                                        while (perm_iter < part_st.permutation.data_end());
                                                                                        assert(part_tr.nr_of_block_bunch_slices ==
                                                                                                                                true_nr_of_block_bunch_slices);
                                                                                        assert(part_tr.action_block.data()<=part_tr.action_block_inert_begin);
                                                                                        assert(part_tr.block_bunch.data() < part_tr.block_bunch_inert_begin);
                                                                                        if (branching)
                                                                                        {   assert(part_tr.action_block_inert_begin <=
                                                                                                                              part_tr.action_block.data_end());
                                                                                            assert(part_tr.block_bunch_inert_begin <=
                                                                                                                               part_tr.block_bunch.data_end());
                                                                                            assert(part_tr.block_bunch.data_end() -
                                                                                                                             part_tr.block_bunch_inert_begin ==
                                                                                             part_tr.action_block.data_end()-part_tr.action_block_inert_begin);

                                                                                            // for all inert transitions
                                                                                            for (const bisim_dnj::action_block_entry* action_block(
                                                                                                                             part_tr.action_block_inert_begin);
                                                                                                                action_block < part_tr.action_block.data_end();
                                                                                                                                                ++action_block)
                                                                                            {   assert(nullptr == action_block->begin_or_before_end);
                                                                                                const bisim_dnj::succ_entry* const
                                                                                                                                 succ_iter(action_block->succ);
                                                                                                assert(nullptr != succ_iter);
                                                                                                assert(succ_iter->block_bunch->slice.is_null());
                                                                                                const bisim_dnj::pred_entry* const
                                                                                                                       pred_iter(succ_iter->block_bunch->pred);
                                                                                                assert(action_block == pred_iter->action_block);
                                                                                                assert(part_tr.block_bunch_inert_begin <=
                                                                                                                                       succ_iter->block_bunch);
                                                                                                assert(pred_iter->source != pred_iter->target);
                                                                                                assert(pred_iter->source->bl.ock == pred_iter->target->bl.ock);
                                                                                                assert(pred_iter->source->succ_inert.begin <= succ_iter);
                                                                                                assert(pred_iter->source->succ_inert.begin == succ_iter ||
                                                                                                   succ_iter[-1].block_bunch->pred->source==pred_iter->source);
                                                                                                assert(pred_iter->target->pred_inert.begin <= pred_iter);
                                                                                                assert(pred_iter->target->pred_inert.begin == pred_iter ||
                                                                                                                    pred_iter[-1].target == pred_iter->target);
                                                                                                unsigned const max_block(check_complexity::log_n -
                                                                                                   check_complexity::ilog2(pred_iter->target->bl.ock->size()));
                                                                                                mCRL2complexity(pred_iter, no_temporary_work(max_block,
                                                                                                                                     max_block, false), *this);
                                                                                            }
                                                                                        }
                                                                                        else
                                                                                        {
                                                                                            assert(!preserve_divergence);
                                                                                            assert(part_tr.action_block_inert_begin ==
                                                                                                                              part_tr.action_block.data_end());
                                                                                            assert(part_tr.block_bunch_inert_begin ==
                                                                                                                               part_tr.block_bunch.data_end());
                                                                                        }
                                                                                        const bisim_dnj::action_block_entry*
                                                                                                            action_slice_end(part_tr.action_block_inert_begin);
                                                                                        trans_type true_nr_of_bunches(0);
                                                                                        trans_type true_nr_of_nontrivial_bunches(0);
                                                                                        trans_type true_nr_of_action_block_slices(0);
                                                                                        // for all action labels and bunches
                                                                                        label_type label(0);
                                                                                        assert(label < action_label.size());
                                                                                        const bisim_dnj::bunch_t* previous_bunch(nullptr);
                                                                                        do
                                                                                        {
                                                                                            assert(part_tr.action_block.data() <= action_label[label].begin);
                                                                                            assert(action_label[label].begin <= action_slice_end);
                                                                                            assert(action_slice_end <= part_tr.action_block_inert_begin);
                                                                                            // for all action_block slices
                                                                                            for (const bisim_dnj::action_block_entry*
                                                                                                                     action_block_slice_end(action_slice_end);
                                                                                                          action_label[label].begin < action_block_slice_end; )
                                                                                            {
                                                                                                const bisim_dnj::action_block_entry* const
                                                                                                        action_block_slice_begin(
                                                                                                               action_block_slice_end[-1].begin_or_before_end);
                                                                                                assert(nullptr != action_block_slice_begin);
                                                                                                assert(action_block_slice_begin < action_block_slice_end);
                                                                                                assert(action_block_slice_end <= action_slice_end);
                                                                                                assert(nullptr != action_block_slice_begin->succ);
                                                                                                const bisim_dnj::block_t* const
                                                                                                        target_block(action_block_slice_begin->
                                                                                                                      succ->block_bunch->pred->target->bl.ock);
                                                                                                const bisim_dnj::bunch_t* const
                                                                                                                bunch(action_block_slice_begin->succ->bunch());
                                                                                                if (previous_bunch != bunch)
                                                                                                {
                                                                                                    assert(nullptr == previous_bunch);
                                                                                                    previous_bunch = bunch;
                                                                                                    assert(bunch->end == action_block_slice_end);
                                                                                                    if (bunch->begin == action_block_slice_begin)
                                                                                                    {
                                                                                                        // Perhaps this does not always hold; sometimes, an
                                                                                                        // action_block slice disappears but the bunch cannot
                                                                                                        // be made trivial.
                                                                                                        assert(bunch->is_trivial());
                                                                                                    }
                                                                                                    else
                                                                                                    {
                                                                                                        assert(!bunch->is_trivial());
                                                                                                        ++true_nr_of_nontrivial_bunches;
                                                                                                    }
                                                                                                    mCRL2complexity(bunch, no_temporary_work(
                                                                                                                       bunch->max_work_counter(*this)), *this);
                                                                                                    ++true_nr_of_bunches;
                                                                                                }
                                                                                                if(bunch->begin == action_block_slice_begin)
                                                                                                {
                                                                                                    previous_bunch = nullptr;
                                                                                                }
                                                                                                else  assert(bunch->begin < action_block_slice_begin);

                                                                                                assert(action_block_slice_begin->begin_or_before_end + 1 ==
                                                                                                                                       action_block_slice_end);
                                                                                                // for all transitions in the action_block slice
                                                                                                const bisim_dnj::action_block_entry*
                                                                                                                          action_block(action_block_slice_end);
                                                                                                do
                                                                                                {
                                                                                                    --action_block;
                                                                                                    const bisim_dnj::succ_entry* const
                                                                                                                                 succ_iter(action_block->succ);
                                                                                                    assert(nullptr != succ_iter);
                                                                                                    const bisim_dnj::pred_entry* const
                                                                                                                       pred_iter(succ_iter->block_bunch->pred);
                                                                                                    assert(action_block == pred_iter->action_block);
                                                                                                    assert(succ_iter->block_bunch <
                                                                                                                              part_tr.block_bunch_inert_begin);
                                                                                                    assert(!branching || !aut.is_tau(label) ||
                                                                                                        pred_iter->source->bl.ock!=pred_iter->target->bl.ock ||
                                                                                                                     (preserve_divergence &&
                                                                                                                      pred_iter->source == pred_iter->target));
                                                                                                    assert(succ_iter < pred_iter->source->succ_inert.begin);
                                                                                                    assert(succ_iter+1==pred_iter->source->succ_inert.begin ||
                                                                                                                succ_iter[1].block_bunch->pred->source ==
                                                                                                                                            pred_iter->source);
                                                                                                    assert(pred_iter < pred_iter->target->pred_inert.begin);
                                                                                                    assert(pred_iter+1==pred_iter->target->pred_inert.begin ||
                                                                                                                     pred_iter[1].target == pred_iter->target);
                                                                                                    assert(target_block == pred_iter->target->bl.ock);
                                                                                                    assert(bunch == succ_iter->bunch());
                                                                                                }
                                                                                                while (action_block_slice_begin < action_block &&
                                                                                                   (// some properties only need to be checked for states that
                                                                                                    // are not the first one:
                                                                                                    assert(action_block->begin_or_before_end ==
                                                                                                                             action_block_slice_begin), true));
                                                                                                action_block_slice_end = action_block_slice_begin;
                                                                                                ++true_nr_of_action_block_slices;
                                                                                            }
                                                                                            if (action_slice_end < part_tr.action_block_inert_begin)
                                                                                            {
                                                                                                // there is a dummy transition between action labels
                                                                                                assert(nullptr == action_slice_end->succ);
                                                                                                assert(nullptr == action_slice_end->begin_or_before_end);
                                                                                            }
                                                                                        }
                                                                                        while (++label < action_label.size() &&
                                                                                                 (action_slice_end = action_label[label - 1].begin - 1, true));
                                                                                        assert(nullptr == previous_bunch);
                                                                                        assert(part_tr.nr_of_bunches == true_nr_of_bunches);
                                                                                        assert(part_tr.nr_of_nontrivial_bunches ==
                                                                                                                                true_nr_of_nontrivial_bunches);
                                                                                        assert(part_tr.nr_of_action_block_slices ==
                                                                                                                               true_nr_of_action_block_slices);
                                                                                    }
                                                                                #endif
    /// \brief Run (branching) bisimulation minimisation in time O(m log n)
    /// \details This function assumes that the partitioner object stores a LTS
    /// with a partition satisfying the invariant:
    ///
    /// > If a state contains a transition in a bunch, then every bottom state
    /// > in the same block contains a transition in that bunch.
    ///
    /// The function runs the efficient O(m log n) algorithm for branching
    /// bisimulation minimisation on the LTS that has been stored in the
    /// partitioner:  As long as there are nontrivial bunches, it selects one,
    /// subdivides it into two bunches and then stabilises the partition for
    /// these bunches.  As a result, the partition stored in the partitioner
    /// will become stable.
    ///
    /// Parameters and return value are implicit with this function:  the LTS,
    /// the partition and the flags of the bisimulation algorithm are all
    /// stored in the partitioner object.
    void refine_partition_until_it_becomes_stable()
    {
        // Line 2.5: for all non-trivial bunches bunch_T in Pi_t do
        std::clock_t next_print_time = std::clock();
        const std::clock_t rounded_start_time=next_print_time-CLOCKS_PER_SEC/2;
        // const double log_initial_nr_of_action_block_slices =
        //                   100 / std::log(part_tr.nr_of_action_block_slices);
        for (;;)
        {                                                                       // mCRL2complexity(...) -- this loop will be ascribed to (the transitions in)
                                                                                // the new bunch below.
            /*------------------ find a non-trivial bunch -------------------*/ ONLY_IF_DEBUG( part_st.print_part(*this);  part_tr.print_trans(*this);
                                                                                                                                          assert_stability(); )
            /* Line 2.6: Select some a in Act and B' in Pi_s such that       */ assert(part_tr.nr_of_bunches + part_tr.nr_of_nontrivial_bunches <=
            /*                           |bunch_T_a_Bprime| <= 1/2 |bunch_T| */                                             part_tr.nr_of_action_block_slices);
            bisim_dnj::bunch_t* const bunch_T(part_tr.get_some_nontrivial());
            if (mCRL2logEnabled(log::verbose))
            {
                if (std::clock_t now = std::clock(); next_print_time <= now ||
                                                            nullptr == bunch_T)
                {

                    /* -  -  -  -  -print progress information-  -  -  -  - */

                    // The formula below should ensure that `next_print_time`
                    // increases by a whole number of minutes, so that the
                    // progress information is printed every minute (or, if
                    // one iteration takes more than one minute, after a whole
                    // number of minutes).
                    next_print_time+=((now-next_print_time)/(60*CLOCKS_PER_SEC)
                                                    + 1) * (60*CLOCKS_PER_SEC);
                    now = (now - rounded_start_time) / CLOCKS_PER_SEC;
                    if (0 != now)
                    {
                        if (60 <= now)
                        {
                            if (3600 <= now)
                            {
                                mCRL2log(log::verbose)
                                    << now / 3600 << " h ";
                                now %= 3600;
                            }
                            mCRL2log(log::verbose)
                                << now / 60 << " min ";
                            now %= 60;
                        }
                        mCRL2log(log::verbose) << now
                              << " sec passed since starting the main loop.\n";
                    }
                    #define PRINT_SG_PL(counter, sg_string, pl_string)        \
                            (counter)                                         \
                            << (1 == (counter) ? (sg_string) : (pl_string))
                    mCRL2log(log::verbose)
                        << (nullptr == bunch_T ? "The reduced LTS contains "
                                        : "The reduced LTS contains at least ")
                        << PRINT_SG_PL(part_st.nr_of_blocks,
                                                 " state and ", " states and ")
                        << PRINT_SG_PL(part_tr.nr_of_block_bunch_slices,
                                              " transition.", " transitions.");
                    if (1 < part_tr.nr_of_action_block_slices)
                    {
                        #define PRINT_INT_PERCENTAGE(num,denom)               \
                                (((num) * 200 + (denom)) / (denom) / 2)
                        mCRL2log(log::verbose) << " Estimated "
                            << PRINT_INT_PERCENTAGE(part_tr.nr_of_bunches - 1,
                                         part_tr.nr_of_action_block_slices - 1)
                            << "% done.";
                        #undef PRINT_INT_PERCENTAGE
                    }
                    mCRL2log(log::verbose)
                    //  << " Logarithmic estimate: "
                    //  << (int)(100.5+std::log((double) part_tr.nr_of_bunches/
                    //                      part_tr.nr_of_action_block_slices)
                    //                  *log_initial_nr_of_action_block_slices)
                    //  << "% done."
                        << "\nThe current partition contains ";
                    if (branching)
                    {
                        mCRL2log(log::verbose)
                            << PRINT_SG_PL(part_tr.nr_of_new_bottom_states,
                                " new bottom state, ", " new bottom states, ");
                    }                                                           else  assert(0 == part_tr.nr_of_new_bottom_states);
                    mCRL2log(log::verbose)
                        << PRINT_SG_PL(part_tr.nr_of_bunches,
                                    " bunch (of which ", " bunches (of which ")
                        << PRINT_SG_PL(part_tr.nr_of_nontrivial_bunches,
                             " is nontrivial), and ", " are nontrivial), and ")
                        << PRINT_SG_PL(part_tr.nr_of_action_block_slices,
                          " action-block-slice.\n", " action-block-slices.\n");
                    #undef PRINT_SG_PL
                }
            }
            if (nullptr == bunch_T)  { break; }                                 ONLY_IF_DEBUG( mCRL2log(log::debug) << "Refining "
            /* Line 2.7: Pi_t := Pi_t \ { bunch_T } union                    */                                          << bunch_T->debug_id(*this) << '\n'; )
            /*              { bunch_T_a_Bprime, bunch_T \ bunch_T_a_Bprime } */ assert(part_tr.nr_of_bunches < part_tr.nr_of_action_block_slices);
            bisim_dnj::bunch_t* const bunch_T_a_Bprime(
                         bunch_T->split_off_small_action_block_slice(part_tr));
                                                                                ONLY_IF_DEBUG( mCRL2log(log::debug) << "Splitting off "
            /*------------ find predecessors of bunch_T_a_Bprime ------------*/                                 << bunch_T_a_Bprime->debug_id(*this) << '\n'; )
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
            /* Line 2.8: for all B in splittableBlocks(bunch_T_a_Bprime) do  */     unsigned const max_splitter_counter(
                /* we actually run through the transitions in T--a-->B'      */                                     bunch_T_a_Bprime->max_work_counter(*this));
                                                                                #endif
            bisim_dnj::action_block_entry* splitter_iter(
                                                      bunch_T_a_Bprime->begin); assert(splitter_iter < bunch_T_a_Bprime->end);
            do
            {                                                                   assert(nullptr != splitter_iter->succ);
                bisim_dnj::state_info_entry* const
                        source(splitter_iter->succ->block_bunch->pred->source); assert(splitter_iter->succ->block_bunch->pred->action_block == splitter_iter);
                // Line 2.11: Mark all transitions in bunch_T_a_Bprime
                    // actually we mark the source state (i.e. register it's in
                    // R)
                bool const first_transition_of_state(
                                            source->bl.ock->mark(source->pos));
                // Line 2.9:  Add bunch_T_a_Bprime as primary to the splitter
                //            list
                // Line 2.10: Add T_B--> \ bunch_T_a_Bprime as secondary to the
                //            splitter list
                part_tr.first_move_transition_to_new_bunch(splitter_iter,
                                  bunch_T_a_Bprime, first_transition_of_state); // mCRL2complexity(splitter_iter->succ->block_bunch->pred, ...) -- subsumed
            // Line 2.13: end for                                               // in the call below
            }
            while (++splitter_iter < bunch_T_a_Bprime->end);

            // We cannot join the loop above with the loop below!

            // Line 2.8: for all B in splittableBlocks(T--a-->B') do
            splitter_iter = bunch_T_a_Bprime->begin;                            assert(splitter_iter < bunch_T_a_Bprime->end);
            do
            {
                // Line 2.12: For every state with both marked outgoing
                //            transitions and outgoing transitions in
                //            T_B--> \ bunch_T_a_Bprime, mark one such
                //            transition
                part_tr.second_move_transition_to_new_bunch(splitter_iter,      ONLY_IF_DEBUG( *this, bunch_T_a_Bprime, )
                                                                      bunch_T); // mCRL2complexity(splitter_iter->succ->block_bunch->pred, ...) -- subsumed
            // Line 2.13: end for                                               // in the call below
            }
            while (++splitter_iter < bunch_T_a_Bprime->end);                    mCRL2complexity(bunch_T_a_Bprime,
                                                                                         add_work(check_complexity::refine_partition_until_stable__find_pred,
            /*----------------- stabilise the partition again ---------------*/                                                  max_splitter_counter), *this);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
            /* Line 2.14: for all T'_B--> in the splitter list (in order) do */     bisim_dnj::block_bunch_slice_iter_or_null_t bbslice_T_a_Bprime_B(nullptr);
                                                                                #endif
            while (!part_tr.splitter_list.empty())
            {
                bisim_dnj::block_bunch_slice_iter_t splitter_Tprime_B(          // We have to call mCRL2complexity here because `splitter_Tprime_B` may be
                                                part_tr.splitter_list.begin()); // split up later.
                bisim_dnj::block_t* block_B(splitter_Tprime_B->source_block()); assert(!splitter_Tprime_B->is_stable());
                bool const is_primary_splitter = 0 < block_B->marked_size();    assert(!splitter_Tprime_B->empty());
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    bool add_stabilize_to_bottom_transns_succeeded = true;
                                                                                    if (is_primary_splitter)
                                                                                    {
                                                                                        assert(bbslice_T_a_Bprime_B.is_null());
                                                                                        // assign work to this splitter bunch
                                                                                        mCRL2complexity(splitter_Tprime_B, add_work(
                                                                                                check_complexity::refine_partition_until_stable__stabilize,
                                                                                                                                 max_splitter_counter), *this);
                                                                                    }
                                                                                    else if (!bbslice_T_a_Bprime_B.is_null())
                                                                                    {
                                                                                        // assign work to this the corresponding block_bunch-slice of
                                                                                        // bunch_T_a_Bprime
                                                                                        mCRL2complexity(bbslice_T_a_Bprime_B,
                                                                                        add_work(check_complexity::
                                                                                                refine_partition_until_stable__stabilize_for_large_splitter,
                                                                                                                                 max_splitter_counter), *this);
                                                                                    }
                                                                                    else
                                                                                    {
                                                                                        // This must be a refinement to stabilize for new bottom states.
                                                                                        // assign work to the new bottom states in this block_bunch-slice
                                                                                        add_stabilize_to_bottom_transns_succeeded = splitter_Tprime_B->
                                                                                            add_work_to_bottom_transns(check_complexity::
                                                                                                refine_partition_until_stable__stabilize_new_noninert_a_priori,
                                                                                                                                                    1U, *this);
                                                                                    }
                                                                                #endif
                if (1 < block_B->size())
                {
                    bisim_dnj::permutation_entry* const
                                                 block_B_begin(block_B->begin); assert(block_B_begin->st->pos == block_B_begin);
                    // Line 2.15: (R, U) := split(B, T'_B-->)
                    // Line 2.16: Remove T'_B--> from the splitter list
                    // Line 2.17: Pi_s := Pi_s \ { B } union { R, U } \ { {} }
                    bisim_dnj::block_t*block_R=split(block_B,splitter_Tprime_B,
                                is_primary_splitter ? extend_from_marked_states
                                                    : extend_from_splitter);
                    if (block_B_begin < block_R->begin)
                    {
                        // The refinement was non-trivial.

                        // Line 2.18: if T'_B--> was a primary splitter then
                        if (is_primary_splitter)
                        {                                                       assert(splitter_Tprime_B->bunch == bunch_T_a_Bprime);
                            // Line 2.19: Remove T_U--> \ T_U--a-->B' from the
                            //            splitter list
                            bisim_dnj::block_t* const
                                            block_U(block_B_begin->st->bl.ock); assert(block_U->end == block_R->begin);
                            bisim_dnj::block_bunch_slice_iter_t U_splitter(
                                                part_tr.splitter_list.begin()); assert(0 == block_U->marked_size());
                            if (part_tr.splitter_list.end() != U_splitter &&
                               (U_splitter->source_block() == block_U ||
                                (++U_splitter != part_tr.splitter_list.end() &&
                                       U_splitter->source_block() == block_U)))
                            {                                                   assert(!U_splitter->is_stable());
                                                                                assert(U_splitter->bunch == bunch_T);
                                block_U->stable_block_bunch.splice(
                                            block_U->stable_block_bunch.end(),
                                            part_tr.splitter_list, U_splitter);
                                U_splitter->make_stable();
                            }
                                                                                #ifndef NDEBUG
                                                                                    // There should be no block-bunch-slice for the co-splitter that is
                                                                                    // still unstable.
                                                                                    for (bisim_dnj::block_bunch_slice_const_iter_t
                                                                                                                          iter(part_tr.splitter_list.cbegin());
                                                                                                                  part_tr.splitter_list.cend() != iter; ++iter)
                                                                                    {   assert(!iter->is_stable());
                                                                                        assert(iter->source_block() != block_U);
                        /* Line 2.20: end if                                 */     }
                                                                                #endif
                        }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    else
                                                                                    {
                                                                                        // account for work that couldn't be accounted for earlier (because we
                                                                                        // didn't know yet which state would become a new bottom state)
                                                                                        if (!add_stabilize_to_bottom_transns_succeeded)
                                                                                        {   assert(splitter_Tprime_B->add_work_to_bottom_transns(
                                                                                            check_complexity::
                                                                                            refine_partition_until_stable__stabilize_new_noninert_a_posteriori,
                                                                                                                                                   1U, *this));
                                                                                        }
                                                                                        if (splitter_Tprime_B->work_counter.has_temporary_work())
                                                                                        {   assert(splitter_Tprime_B->add_work_to_bottom_transns(
                                                                                                    check_complexity::
                                                                                                    handle_new_noninert_transns__make_unstable_a_posteriori,
                                                                                                                                                   1U, *this));
                                                                                            splitter_Tprime_B->work_counter.reset_temporary_work();
                        /* Line 2.21: if R--tau-->U is not empty (i. e. R    */         }
                        /*               has new non-inert transitions) then */     }
                                                                                #endif
                        if (0 < block_R->marked_size())
                        {                                                       ONLY_IF_DEBUG( const bisim_dnj::block_bunch_entry* const splitter_end =
                            /* Line 2.22: Create a new bunch containing      */                                                       splitter_Tprime_B->end; )
                            //            exactly R--tau-->U, add R--tau-->U to
                            //            the splitter list, and mark all its
                            //            transitions
                            // to
                            // Line 2.28: For each bottom state, mark one of
                            //            its outgoing transitions in every
                            //            T_N--> where it has one
                            block_R = handle_new_noninert_transns(
                                                   block_R, splitter_Tprime_B);
                                                                                #ifndef NDEBUG
                                                                                    if (splitter_end[-1].pred->source->bl.ock == block_R)
                                                                                    {   assert(!splitter_end[-1].slice.is_null());
                                                                                        splitter_Tprime_B =
                                                                                                  (bisim_dnj::block_bunch_slice_iter_t) splitter_end[-1].slice;
                                                                                    }
                        /* Line 2.29: end if                                 */     assert(nullptr == block_R || splitter_Tprime_B->source_block() == block_R);
                                                                                #endif
                        }
                    }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    else
                                                                                    {   assert(0 == block_R->marked_size());
                                                                                        assert(add_stabilize_to_bottom_transns_succeeded);
                                                                                        // now splitter must have some transitions that start in bottom states:
                                                                                        if (splitter_Tprime_B->work_counter.has_temporary_work())
                                                                                        {   assert(!is_primary_splitter);
                                                                                            if (!splitter_Tprime_B->add_work_to_bottom_transns(check_complexity
                                                                                                  ::handle_new_noninert_transns__make_unstable_a_posteriori,
                                                                                                                                                    1U, *this))
                                                                                            {  assert(0);  }
                                                                                            splitter_Tprime_B->work_counter.reset_temporary_work();
                                                                                        }
                                                                                    }
                                                                                    block_B = block_R;
                                                                                #endif
                }
                else
                {                                                               assert(block_B->nonbottom_begin == block_B->end);
                    /* A block with 1 state will not be split.  However, we  */ assert(block_B->marked_nonbottom_begin == block_B->end);
                    // may have to unmark all states.
                    block_B->marked_bottom_begin = block_B->end;
                    block_B->stable_block_bunch.splice(
                                     block_B->stable_block_bunch.end(),
                                     part_tr.splitter_list, splitter_Tprime_B);
                    splitter_Tprime_B->make_stable();                           assert(add_stabilize_to_bottom_transns_succeeded);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    // now splitter must have some transitions that start in bottom states:
                                                                                    if (splitter_Tprime_B->work_counter.has_temporary_work())
                                                                                    {   assert(!is_primary_splitter);
                                                                                        if (!splitter_Tprime_B->add_work_to_bottom_transns(check_complexity::
                                                                                                    handle_new_noninert_transns__make_unstable_a_posteriori,
                                                                                                                                                   1U, *this))
                                                                                        {  assert(0);  }
                                                                                        splitter_Tprime_B->work_counter.reset_temporary_work();
                                                                                    }
                                                                                #endif
                }
                                                                                #ifndef NDEBUG
                                                                                    if (is_primary_splitter && !part_tr.splitter_list.empty() &&
                                                                                                       part_tr.splitter_list.front().bunch == bunch_T &&
                                                                                                       part_tr.splitter_list.front().source_block() == block_B)
                                                                                    {   // The next block_bunch-slice to be handled is the one in the large
                                                                                        // splitter corresponding to the current splitter.  In that iteration,
                                                                                        // we will need the current splitter block_bunch-slice.
                                                                                        assert(nullptr != block_B);
                                                                                        assert(splitter_Tprime_B->source_block() == block_B);
                                                                                        assert(splitter_Tprime_B->bunch == bunch_T_a_Bprime);
                                                                                        bbslice_T_a_Bprime_B = splitter_Tprime_B;
                                                                                    }
            /* Line 2.30: end for                                            */     else  bbslice_T_a_Bprime_B = nullptr;
                                                                                #endif
            }
        // Line 2.31: end for
        }                                                                       assert(part_tr.nr_of_bunches == part_tr.nr_of_action_block_slices);
                                                                                assert(0 == part_tr.nr_of_nontrivial_bunches);

        // store the labels with the action_block-slices
        // As every action_block-slice is a (trivial) bunch at the same time,
        // we can reuse the field next_nontrivial_and_label.label (instead of
        // next_nontrivial_and_label.next_nontrivial) to store the label.
        const bisim_dnj::action_block_entry*
                       action_block_iter_end(part_tr.action_block_inert_begin);
        label_type label(0);                                                    assert(label < action_label.size());
        do
        {
            for (bisim_dnj::action_block_entry*
                                 action_block_iter(action_label[label].begin);
              action_block_iter < action_block_iter_end;
              action_block_iter = action_block_iter->begin_or_before_end + 1)
            {                                                                   assert(nullptr != action_block_iter->succ);
                                                                                assert(action_block_iter->succ->block_bunch->pred->action_block ==
                                                                                                                                            action_block_iter);
                                                                                assert(action_block_iter->succ->bunch()->is_trivial());
                action_block_iter->succ->bunch()->
                                       next_nontrivial_and_label.label = label; assert(nullptr != action_block_iter->begin_or_before_end);
                                                                                assert(action_block_iter <= action_block_iter->begin_or_before_end);
            }
        }
        while (++label < action_label.size() &&
            (action_block_iter_end = action_label[label - 1].begin - 1, true));
    }

    /*----------------- Split -- Algorithm 3 of [JGKW 2020] -----------------*/

    /// \brief Split a block according to a splitter
    /// \details The function splits `block_B` into the R-subblock (states
    /// with a transition in `splitter_T`) and the U-subblock (states without a
    /// transition in `splitter_T`).  Depending on `mode`, the states are
    /// primed as follows:
    ///
    /// - If `mode == extend_from_marked_states`, then all states with strong
    ///   transition(s) must have been marked already.
    /// - If `mode == extend_from_marked_states__add_new_noninert_to_splitter`,
    ///   states are marked as above.  The only difference is the handling of
    ///   new non-inert transitions.
    /// - If `mode == extend_from_splitter`, then no states must be marked;
    ///   the initial states with a transition in `splitter_T` are searched by
    ///   `split()` itself.  Every bottom state with strong transition(s)
    ///   needs to have at least one marked strong transition.
    ///
    /// The  function  will  also  adapt  all  data  structures  and  determine
    /// which  transitions  have  changed  from  inert  to  non-inert.   States
    /// with  a  new  non-inert  transition  will  be  marked  upon  returning.
    /// Normally,  the  new  non-inert  transitions  are  moved  to  a  new
    /// bunch,  which  will  be  specially  created.   However,  if  `mode ==
    /// extend_from_marked_states__add_new_noninert_to_splitter`, then the new
    /// non-inert transitions will be added to `splitter_T` (which must hold
    /// transitions that have just become non-inert before this call to
    /// `split()`).  If the resulting block contains marked states, the caller
    /// has to call `handle_new_noninert_transns()` to stabilise the block
    /// because the new bunch may make the block unstable.
    /// \param block_B     block that needs to be refined
    /// \param splitter_T  transition set that makes the block unstable
    /// \param mode        indicates how to find states with a transition in
    ///                    `splitter_T`, as described above
    /// \returns (a pointer to) the R-subblock.  It is an error to call the
    /// function with settings that lead to an empty R-subblock.  (An empty
    /// U-subblock is ok.)
    bisim_dnj::block_t* split(bisim_dnj::block_t* const block_B,
              const bisim_dnj::block_bunch_slice_iter_t splitter_T,
                                                       enum refine_mode_t mode)
    {                                                                           assert(block_B == splitter_T->source_block());
                                                                                #ifndef NDEBUG
                                                                                    mCRL2log(log::debug) << "split("
                                                                                        << block_B->debug_id(*this)
                                                                                        << ',' << splitter_T->debug_id(*this)
                                                                                        << (extend_from_marked_states__add_new_noninert_to_splitter == mode
                                                                                           ? ",extend_from_marked_states__add_new_noninert_to_splitter)\n"
                                                                                           : (extend_from_marked_states == mode
                                                                                             ? ",extend_from_marked_states)\n"
                                                                                             : (extend_from_splitter == mode
                                                                                               ? ",extend_from_splitter)\n"
                                                                                               : ",UNKNOWN MODE)\n")));
                                                                                #endif
        bisim_dnj::block_t* block_R=nullptr;                                   assert(!splitter_T->is_stable());  assert(1 < block_B->size());
        union R_s_iter_t
        {
            bisim_dnj::block_bunch_entry* splitter_iter;
            bisim_dnj::permutation_entry* block;
        } R_s_iter;

        if (extend_from_splitter == mode)
        {                                                                       assert(0 == block_B->marked_size());
            // Line 3.2: R := B--Marked(T)--> ; U := Bottom(B) \ R
            R_s_iter.splitter_iter = splitter_T->end;                           assert(splitter_T->marked_begin <= R_s_iter.splitter_iter);
            while (splitter_T->marked_begin < R_s_iter.splitter_iter)
            {                                                                   assert(&part_tr.block_bunch.cbegin()[1] < R_s_iter.splitter_iter);
                --R_s_iter.splitter_iter;
                bisim_dnj::state_info_entry* const
                                       s(R_s_iter.splitter_iter->pred->source); assert(s->bl.ock == block_B);  assert(s->pos->st == s);
                block_B->mark(s->pos);
                // We cannot stop, even if the R-subblock becomes too large,
                // because we need to mark all bottom states that are not in U.
            }
        }                                                                       else  {  assert(0 < block_B->marked_size());
                                                                                         assert(splitter_T->marked_begin == splitter_T->end);  }
        block_B->stable_block_bunch.splice(block_B->stable_block_bunch.end(),
                                            part_tr.splitter_list, splitter_T);
        splitter_T->make_stable();

        COROUTINES_SECTION
            // shared variables of both coroutines
            bisim_dnj::permutation_entry*
                           untested_to_U_defined_end(block_B->nonbottom_begin);
            bisim_dnj::permutation_entry*
                                    U_nonbottom_end(untested_to_U_defined_end);

            // variable declarations of the U-coroutine
            bisim_dnj::permutation_entry* U_s_iter;
            bisim_dnj::pred_entry* U_t_iter;
            bisim_dnj::state_info_entry* U_t;
            const bisim_dnj::succ_entry* U_u_iter;

            // variable declarations of the R-coroutine
            bisim_dnj::pred_entry* R_t_iter;

            COROUTINE_LABELS(   (SPLIT_R_PREDECESSOR_HANDLED)
                                (SPLIT_U_PREDECESSOR_HANDLED)
                                (SPLIT_R_STATE_HANDLED)
                                (SPLIT_U_STATE_HANDLED)
                                (SPLIT_U_TESTING)
                                (SPLIT_R_COLLECT_SPLITTER))

            /*------------------------ find U-states ------------------------*/

            COROUTINE
                // Line 3.21l: if |U| > |B|/2 then
                if(block_B->size() / 2 < block_B->unmarked_bottom_size())
                {
                    // Line 3.22l: Abort this coroutine
                    ABORT_THIS_COROUTINE();
                // Line 3.23l: end if
                }
                if (0 == block_B->unmarked_bottom_size())
                {
                    // all bottom states are in R, so there cannot be any
                    // U-states.  Unmark all states, as there are no
                    // transitions that have become non-inert.
                    block_B->marked_nonbottom_begin = block_B->end;
                    block_B->marked_bottom_begin = block_B->nonbottom_begin;
                    block_R = block_B;
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    finalise_U_is_smaller(nullptr, block_R, *this);
                                                                                #endif
                    TERMINATE_COROUTINE_SUCCESSFULLY();
                }

                 /* -  -  -  -  -  -  -  visit U-states  -  -  -  -  -  -  - */

                if (U_nonbottom_end < block_B->marked_nonbottom_begin)
                {
                    // Line 3.5l: for all s in U while |U| < |B|/2 do
                    U_s_iter = block_B->begin;
                    COROUTINE_DO_WHILE (SPLIT_U_STATE_HANDLED,
                                                    U_s_iter < U_nonbottom_end)
                    {
                        /* Line 3.6l: for all inert transitions t--tau-->s do*/ assert(part_tr.pred.front().target != U_s_iter->st);
                        COROUTINE_FOR(SPLIT_U_PREDECESSOR_HANDLED,
                                  U_t_iter = U_s_iter->st->pred_inert.begin,
                                  U_t_iter->target == U_s_iter->st, ++U_t_iter)
                        {
                            U_t = U_t_iter->source;                             assert(block_B->nonbottom_begin <= U_t->pos);
                            /* Line 3.7l: if t in R then  Skip state t       */ assert(U_t->pos < block_B->end);
                            if (block_B->marked_nonbottom_begin <= U_t->pos)
                            {
                                goto continuation;
                            }
                            // Line 3.8l: if untested[t] is undefined then
                            if (untested_to_U_defined_end <= U_t->pos)
                            {
                                // Line 3.9l: untested[t] :=
                                //                    |{ t--tau-->u | u in B }|
                                U_t->untested_to_U_eqv.begin =
                                                         U_t->succ_inert.begin;
                                std::swap(*U_t->pos,
                                                 *untested_to_U_defined_end++);
                            // Line 3.10l: end if
                            }                                                   assert(U_t != part_tr.succ.back().block_bunch->pred->source);
                            // Line 3.11l: untested[t] := untested[t] − 1
                            ++U_t->untested_to_U_eqv.begin;
                            // Line 3.12l: if untested[t]>0 then  Skip state t
                            if (U_t == U_t->untested_to_U_eqv.
                                              begin->block_bunch->pred->source)
                            {
                                goto continuation;
                            }
                            // Line 3.13l: if not (B--T--> subset R) then
                            if (extend_from_splitter == mode)
                            {                                                   assert(U_t != part_tr.succ.front().block_bunch->pred->source);
                                // Line 3.14l: for all non-inert
                                //             t --alpha--> u do
                                U_u_iter = U_t->succ_inert.begin;               assert(part_tr.succ.data() < U_u_iter);
                                COROUTINE_WHILE(SPLIT_U_TESTING, U_t ==
                                        U_u_iter[-1].block_bunch->pred->source)
                                {
                                    U_u_iter=U_u_iter[-1].begin_or_before_end;  assert(nullptr != U_u_iter);
                                    /* Line 3.15l: if t --alpha--> u in T    */ assert(U_u_iter->block_bunch->pred->source == U_t);
                                    /*             then  Skip t              */ assert(!U_u_iter->block_bunch->slice.is_null());
                                    bisim_dnj::block_bunch_slice_const_iter_t
                                            const block_bunch(
                                                 U_u_iter->block_bunch->slice);
                                    if (&*block_bunch == &*splitter_T)
                                    {
                                        goto continuation;
                                               // i. e. break and then continue
                                    }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    bisim_dnj::succ_entry::add_work_to_out_slice(*this, U_u_iter,
                                /* Line 3.16l: end for                       */                      check_complexity::split_U__test_noninert_transitions, 1U);
                                                                                #endif
                                }
                                END_COROUTINE_WHILE;
                            // Line 3.17l: end if
                            }                                                   assert(U_nonbottom_end <= U_t->pos);
                            /* Line 3.18l: Add t to U                        */ assert(U_t->pos < untested_to_U_defined_end);
                            std::swap(*U_t->pos, *U_nonbottom_end++);
                            // Line 3.21l: if |U| > |B|/2 then
                            if (block_B->size() / 2 <
                                    U_nonbottom_end-block_B->nonbottom_begin +
                                               block_B->unmarked_bottom_size())
                            {
                                // Line 3.22l: Abort this coroutine
                                // As the U-coroutine is now aborted, the
                                // untested_to_U values are no longer relevant.
                                // The assignment tells the R-coroutine that it
                                // doesn't need to make complicated swaps any
                                // more to keep untested properly initialized.
                                untested_to_U_defined_end = U_nonbottom_end;
                                ABORT_THIS_COROUTINE();
                            // Line 3.23l: end if
                            }
                        // Line 3.19l: end for
                    continuation:                                               mCRL2complexity(U_t_iter, add_work(
                                                                                          check_complexity::split_U__handle_transition_to_U_state, 1U), *this);
                        }
                        END_COROUTINE_FOR;                                      mCRL2complexity(U_s_iter->st, add_work(
                    /* Line 3.20l: end for                                   */           check_complexity::split_U__find_predecessors_of_U_state, 1U), *this);
                        ++U_s_iter;
                        if(block_B->marked_bottom_begin == U_s_iter)
                        {
                            U_s_iter = block_B->nonbottom_begin;
                        }
                    }
                    END_COROUTINE_DO_WHILE;
                }

                /*-  -  -  -  -  -  -  split off U-block  -  -  -  -  -  -  -*/

                // Line 3.24l: Abort the other coroutine
                ABORT_OTHER_COROUTINE();
                // Line 2.17: Pi_s := Pi_s \ { B } union ({ R, U } \ { {} })
                    // All non-U states are in R.
                block_B->marked_nonbottom_begin = U_nonbottom_end;
                block_R = block_B;
                bisim_dnj::block_t* const block_U(
                    block_R->split_off_block(bisim_dnj::new_block_is_U,         ONLY_IF_DEBUG( *this, )
                                                      part_st.nr_of_blocks++));
                // Line 2.16: Remove Tprime_B--> = Tprime_R--> from the
                //            splitter list
                /* and the remainder of Line 2.17                            */ assert(0 == block_U->marked_size());  assert(0 == block_R->marked_size());
                part_tr.adapt_transitions_for_new_block(block_U, block_R,       ONLY_IF_DEBUG( *this, )
                    extend_from_marked_states__add_new_noninert_to_splitter ==
                                  mode, splitter_T, bisim_dnj::new_block_is_U);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    finalise_U_is_smaller(block_U, block_R, *this);
                                                                                #endif
            END_COROUTINE

            /*------------------------ find R-states ------------------------*/

            COROUTINE
                // Line 3.21r: if |R| > |B|/2 then
                if (block_B->size() / 2 < block_B->marked_size())
                {
                    // Line 3.22r: Abort this coroutine
                    ABORT_THIS_COROUTINE();
                // Line 3.23r: end if
                }

                /* -  -  -  -  -  collect states from B--T-->  -  -  -  -  - */

                if (extend_from_splitter == mode)
                {
                    // Line 3.4r: R := R union B--(T \ Marked(T))-->
                    if (U_nonbottom_end < block_B->marked_nonbottom_begin)
                    {                                                           assert(part_tr.block_bunch.front().slice != splitter_T);
                        COROUTINE_WHILE (SPLIT_R_COLLECT_SPLITTER,
                                R_s_iter.splitter_iter[-1].slice == splitter_T)
                        {                                                       assert(&part_tr.block_bunch.cbegin()[1] < R_s_iter.splitter_iter);
                            --R_s_iter.splitter_iter;
                            bisim_dnj::state_info_entry* const s(
                                         R_s_iter.splitter_iter->pred->source); assert(s->bl.ock == block_B);  assert(s->pos->st == s);
                            if (block_B->nonbottom_begin <= s->pos)
                            {                                                   assert(U_nonbottom_end <= s->pos);
                                if (s->pos < untested_to_U_defined_end)
                                {
                                    // The non-bottom state has a transition
                                    // to a visited U-state, so untested is
                                    // initialised; however, now it is
                                    // discovered to be in R anyway.
                                    std::swap(*s->pos,
                                                 *--untested_to_U_defined_end);
                                }
                                if (block_B->mark_nonbottom(s->pos) &&
                                // Line 3.21r: if |R| > |B|/2 then
                                    block_B->size()/2 < block_B->marked_size())
                                {
                                    // Line 3.22r: Abort this coroutine
                                    ABORT_THIS_COROUTINE();
                                // Line 3.23r: end if
                                }
                            }                                                   else  assert(block_B->marked_bottom_begin <= s->pos);
                                                                                mCRL2complexity(R_s_iter.splitter_iter->pred, add_work(
                                                                                        check_complexity::split_R__handle_transition_from_R_state, 1U), *this);
                        }
                        END_COROUTINE_WHILE;

                        // Indicate to the U-coroutine that all states in
                        // B--T--> are now in R.
                            // The shared variable `mode` is used
                            // instead of a separate shared variable.
                        mode = extend_from_marked_states;
                    }
                                                                                #ifndef NDEBUG
                                                                                    else
                                                                                    {
                                                                                        // assert that all non-marked transitions in `splitter_T` start in
                                                                                        // marked states
                                                                                        assert(part_tr.block_bunch.front().slice != splitter_T);
                                                                                        while (R_s_iter.splitter_iter[-1].slice == splitter_T)
                                                                                        {
                                                                                            assert(&part_tr.block_bunch.cbegin()[1] < R_s_iter.splitter_iter);
                                                                                            --R_s_iter.splitter_iter;
                                                                                            bisim_dnj::state_info_entry* const
                                                                                                                       s(R_s_iter.splitter_iter->pred->source);
                                                                                            assert(s->bl.ock == block_B);  assert(s->pos->st == s);
                                                                                            assert(s->pos < block_B->nonbottom_begin ||
                                                                                                                    block_B->marked_nonbottom_begin <= s->pos);
                                                                                            assert(block_B->marked_bottom_begin <= s->pos);
                                                                                        }
                                                                                    }
                                                                                #endif
                }

                /* -  -  -  -  -  -  -  visit R-states  -  -  -  -  -  -  - */  assert(0 != block_B->marked_size());

                if (U_nonbottom_end < block_B->marked_nonbottom_begin)
                {
                    // Line 3.5r: for all s in R while |R| < |B|/2 do
                    R_s_iter.block = block_B->nonbottom_begin;
                    if (block_B->marked_bottom_begin == R_s_iter.block)
                    {
                        // It may happen that all found states are non-bottom
                        // states.  (In that case, some of these states will
                        // become new bottom states.)
                        R_s_iter.block = block_B->end;
                    }                                                           assert(block_B->marked_nonbottom_begin != R_s_iter.block);
                    COROUTINE_DO_WHILE(SPLIT_R_STATE_HANDLED,
                             block_B->marked_nonbottom_begin != R_s_iter.block)
                    {
                        --R_s_iter.block;                                       assert(part_tr.pred.back().target != R_s_iter.block->st);
                        // Line 3.7r: for all inert transitions t--tau-->s do
                        COROUTINE_FOR (SPLIT_R_PREDECESSOR_HANDLED,
                            R_t_iter = R_s_iter.block->st->pred_inert.begin,
                            R_t_iter->target == R_s_iter.block->st, ++R_t_iter)
                        {
                            bisim_dnj::state_info_entry* const
                                                           t(R_t_iter->source); assert(U_nonbottom_end <= t->pos);
                            /* Line 3.18r: Add t to R                        */ assert(t->pos->st == t);  assert(t->pos < block_B->end);
                            if (t->pos < untested_to_U_defined_end)
                            {
                                // The state has a transition to a U-state, so
                                // untested is initialised; however, now it is
                                // discovered to be in R anyway.
                                std::swap(*t->pos,
                                                 *--untested_to_U_defined_end);
                            }
                            if (block_B->mark_nonbottom(t->pos) &&
                            // Line 3.21r: if |R| > |B|/2 then
                                  block_B->size() / 2 < block_B->marked_size())
                            {
                                // Line 3.22r: Abort this coroutine
                                ABORT_THIS_COROUTINE();
                            // Line 3.23r: end if
                            }                                                   mCRL2complexity(R_t_iter, add_work(
                        /* Line 3.19r: end for                               */           check_complexity::split_R__handle_transition_to_R_state, 1U), *this);
                        }
                        END_COROUTINE_FOR;                                      mCRL2complexity(R_s_iter.block->st, add_work(
                                                                                                check_complexity::split_R__find_predecessors_of_R_state,
                    /* Line 3.20r: end for                                   */                                                                    1U), *this);
                        if (block_B->marked_bottom_begin == R_s_iter.block &&
                                     R_s_iter.block < block_B->nonbottom_begin)
                        {
                            R_s_iter.block = block_B->end;
                        }
                    }
                    END_COROUTINE_DO_WHILE;
                }

                /*-  -  -  -  -  -  -  split off R-block  -  -  -  -  -  -  -*/

                // Line 3.24r: Abort the other coroutine
                ABORT_OTHER_COROUTINE();
                // Line 2.17: Pi_s := Pi_s \ { B } union ({ R, U } \ { {} })
                    // All non-R states are in U.
                block_R = block_B->split_off_block(bisim_dnj::new_block_is_R,   ONLY_IF_DEBUG( *this, )
                                                       part_st.nr_of_blocks++);
                // Line 2.16: Remove Tprime_B--> = Tprime_R--> from the
                //            splitter list
                /* and the remainder of Line 2.17                            */ assert(0 == block_B->marked_size());  assert(0 == block_R->marked_size());
                part_tr.adapt_transitions_for_new_block(block_R, block_B,       ONLY_IF_DEBUG( *this, )
                    extend_from_marked_states__add_new_noninert_to_splitter ==
                                  mode, splitter_T, bisim_dnj::new_block_is_R);
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    finalise_R_is_smaller(block_B, block_R, *this);
                                                                                #endif
            END_COROUTINE
        END_COROUTINES_SECTION
        return block_R;
    }

    /*-- Handle new non-inert transitions -- Lines 2.22-2.28 in [JGKW2020] --*/

    /// \brief Handle a block with new non-inert transitions
    /// \details When this function starts, it assumes that the states with a
    /// new non-inert transition in `block_R` are marked.  It is an error if it
    /// does not contain any marked states.
    ///
    /// The function separates the states with new non-inert transitions from
    /// those without;  as a result, the N-subblock (which contains states
    /// with new non-inert transitions) will contain at least one new bottom
    /// state (and no old bottom states).  It then registers almost all
    /// block_bunch-slices of this N-subblock as unstable and marks one
    /// transition per block_bunch-slice and new bottom state.  Only the
    /// block_bunch-slice containing the new non-inert transitions and, if
    /// possible, `bbslice_Tprime_R` are not registered as unstable.
    /// \param block_R           block containing states with new non-inert
    ///                          transitions that need to be stabilised
    /// \param bbslice_Tprime_R  splitter of the last separation before, i.e.
    ///                          the splitter that made these transitions
    ///                          non-inert (`block_R` should already be stable
    ///                          w.r.t. `bbslice_Tprime_R`).
    /// \returns the block containing the old bottom states (and every state in
    ///          `block_R` that cannot reach any new non-inert transition),
    ///          i.e. the U-subblock of the separation
    bisim_dnj::block_t* handle_new_noninert_transns(
                    bisim_dnj::block_t* const block_R,
                    bisim_dnj::block_bunch_slice_const_iter_t bbslice_Tprime_R)
    {                                                                           assert(block_R == bbslice_Tprime_R->source_block());
        bisim_dnj::block_t* block_Rprime;                                       assert(&part_tr.block_bunch.cbegin()[1] < part_tr.block_bunch_inert_begin);
        bisim_dnj::block_t* block_N;                                            assert(!part_tr.block_bunch_inert_begin[-1].slice.is_null());
        bisim_dnj::block_bunch_slice_iter_t const bbslice_R_tau_U(
                                    part_tr.block_bunch_inert_begin[-1].slice); assert(bbslice_Tprime_R->is_stable());
        /* Line 2.23: (N, R') := split(R, R--tau-->U)                        */ assert(!bbslice_R_tau_U->is_stable());
        /* Line 2.24: Remove R--tau-->U from the splitter list               */ assert(block_R == bbslice_R_tau_U->source_block());
        /* Line 2.25: Pi_s := Pi_s \ { R } union { N, R' } \ { emptyset }    */ assert(0 < block_R->marked_size());
        // Line 2.26: Add N--tau-->R' to the bunch containing R--tau-->U
        if (0 < block_R->unmarked_bottom_size())
        {                                                                       assert(part_tr.splitter_list.begin() != part_tr.splitter_list.end());
                                                                                #ifndef NDEBUG
                                                                                    bool const next_splitter_is_of_same_block =
                                                                                            part_tr.splitter_list.begin() != bbslice_R_tau_U &&
                                                                                                       part_tr.splitter_list.front().source_block() == block_R;
                                                                                #endif
            block_N = split(block_R, bbslice_R_tau_U,
                      extend_from_marked_states__add_new_noninert_to_splitter); assert(part_st.permutation.data() < block_N->begin);
            block_Rprime = block_N->begin[-1].st->bl.ock;
                                                                                #ifndef NDEBUG
                                                                                    // If the first element of the splitter list was a block_bunch-slice of
                                                                                    // block_N, it was split up.  The condition below checks whether the
                                                                                    // N-subblock's (= the block with new bottom states) slice is placed before
                                                                                    // the R'-subblock's (= the block with old bottom states).
                                                                                    if (next_splitter_is_of_same_block &&
                                                                                        (assert(part_tr.splitter_list.begin() != part_tr.splitter_list.end()),
                                                                                                        part_tr.splitter_list.front().source_block()==block_N))
                                                                                    {
                                                                                        bisim_dnj::block_bunch_slice_iter_t const bbslice_T_Rprime(
                                                                                                                     std::next(part_tr.splitter_list.begin()));
                                                                                        if (part_tr.splitter_list.end() != bbslice_T_Rprime &&
                                                                                                              bbslice_T_Rprime->source_block() == block_Rprime)
                                                                                        {
                                                                                            // The R'-subblock's slice must be the first in the splitter list.
                                                                                            // This is necessary in Debug-mode to ensure that the cost of
                                                                                            // refining R' is accounted for correctly.
                                                                                            part_tr.splitter_list.splice(part_tr.splitter_list.begin(),
            /* If more new noninert transitions are found, we do not need to */                                       part_tr.splitter_list, bbslice_T_Rprime);
            /* separate them further, as every bottom state already has a    */         }
            /* transition in bbslice_R_tau_U->bunch.                         */     }
                                                                                #endif
            if (0 < block_N->marked_bottom_size())
            {
                // Not only new non-inert transitions, but also new bottom
                // states have been found.  In that case, we also have to
                // refine w.r.t. the last splitter, as the new bottom states in
                // block_N may be unstable under it.  We set the variable
                // `bbslice_Tprime_R` to `bbslice_R_tau_U` so it won't disturb
                // in the test below.
                bbslice_Tprime_R = bbslice_R_tau_U;
                block_N->marked_bottom_begin = block_N->nonbottom_begin;
            }
            else if (bbslice_Tprime_R->source_block() != block_N)
            {                                                                   assert(bbslice_Tprime_R->source_block() == block_Rprime);
                // bbslice_Tprime_R contained transitions from every (old and
                // new) bottom state in the block.  It has been split, and now
                // it contains transitions from the block with old bottom
                /* states; however, we need the block_bunch-slice with       */ assert(!bbslice_Tprime_R->end->slice.is_null());
                /* transitions from the block with new bottom states.        */ assert(bbslice_Tprime_R->end < part_tr.block_bunch_inert_begin);
                bbslice_Tprime_R = (bisim_dnj::block_bunch_slice_const_iter_t)
                                                  bbslice_Tprime_R->end->slice; assert(bbslice_Tprime_R->source_block() == block_N);
            }
        }
        else
        {
            block_N = block_R;
            // make bbslice_R_tau_U stable
            block_N->stable_block_bunch.splice(
                                       block_N->stable_block_bunch.end(),
                                       part_tr.splitter_list, bbslice_R_tau_U);
            bbslice_R_tau_U->make_stable();
            block_N->marked_bottom_begin = block_N->nonbottom_begin;
            block_Rprime = nullptr;
        }
        block_N->marked_nonbottom_begin = block_N->end;

        if (1 >= block_N->size())  return block_Rprime;

        // Line 2.27: Insert all T_N--> as secondary into the splitter list
            // However, the bunch of new noninert transitions and the bunch
            // that was the last splitter do not need to be handled (as long
            // as there are no further new bottom states).
            // We cannot do this in time O(1) because we need to call
            // `make_unstable()` for each block_bunch-slice individually.
        for (bisim_dnj::block_bunch_slice_iter_t bbslice_T_N(
                                          block_N->stable_block_bunch.begin());
                            block_N->stable_block_bunch.end() != bbslice_T_N; )
        {                                                                       assert(bbslice_T_N->is_stable());
            bisim_dnj::block_bunch_slice_iter_t const
                                      next_bbslice_T_N(std::next(bbslice_T_N));
            if (&*bbslice_T_N != &*bbslice_Tprime_R &&
                                            &*bbslice_T_N != &*bbslice_R_tau_U)
            {
                // In Debug mode, we have to place the new splitters at the end
                // of the splitter list -- after a refinement with a primary
                // splitter, the corresponding refinement with the large
                // splitter should follow immediately, to ensure that the cost
                // for refining `block_Rprime` is accounted for correctly.
                part_tr.splitter_list.splice(part_tr.splitter_list.end(),
                                     block_N->stable_block_bunch, bbslice_T_N);
                bbslice_T_N->make_unstable();
            }
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    // Try to assign this work to a transition from a bottom state in
                                                                                    // bbslice_T_N.
                                                                                    // If that does not succeed, temporarily assign it to the block_bunch
                                                                                    // itself.  Later, we shall find a bottom state to which this work can be
                                                                                    // assigned.
                                                                                    assert(!bbslice_T_N->work_counter.has_temporary_work());
                                                                                    if (!bbslice_T_N->add_work_to_bottom_transns(check_complexity::
                                                                                               handle_new_noninert_transns__make_unstable_a_priori, 1U, *this))
                                                                                    {   mCRL2complexity(bbslice_T_N, add_work(check_complexity::
                                                                                                  handle_new_noninert_transns__make_unstable_temp, 1U), *this);
                                                                                        assert(bbslice_T_N->work_counter.has_temporary_work());
                                                                                        assert(!bbslice_T_N->is_stable());
                                                                                    }
                                                                                #endif
            bbslice_T_N = next_bbslice_T_N;
        }

        // Line 2.28: For each bottom state, mark one of its outgoing
        //            transitions in each in every T_N--> where it has one
        bisim_dnj::permutation_entry* s_iter(block_N->begin);                   assert(s_iter < block_N->nonbottom_begin);
        do
        {
            bisim_dnj::state_info_entry* const s(s_iter->st);                   assert(s->pos == s_iter);
            // for all out-slices of s do
            for (bisim_dnj::succ_entry* succ_iter(s->succ_inert.begin);
                                s == succ_iter[-1].block_bunch->pred->source; )
            {                                                                   assert(succ_iter[-1].begin_or_before_end < succ_iter);
                succ_iter = succ_iter[-1].begin_or_before_end;                  assert(nullptr != succ_iter);
                /* Mark the first transition in the out-slice in its         */ assert(succ_iter->block_bunch->pred->action_block->succ == succ_iter);
                /* block_bunch-slice                                         */ assert(s == succ_iter->block_bunch->pred->source);
                bisim_dnj::block_bunch_entry* const
                                   old_block_bunch_pos(succ_iter->block_bunch); assert(!old_block_bunch_pos->slice.is_null());
                bisim_dnj::block_bunch_slice_iter_t const
                           bbslice_T_N((bisim_dnj::block_bunch_slice_iter_t)
                                                   old_block_bunch_pos->slice);
                if (!bbslice_T_N->is_stable())
                {                                                               assert(&*bbslice_T_N != &*bbslice_Tprime_R && bbslice_T_N != bbslice_R_tau_U);
                    bisim_dnj::block_bunch_entry* const
                            new_block_bunch_pos(bbslice_T_N->marked_begin - 1);
                    // It may happen that the transition was already paid
                    // for earlier, namely if it once was in bunch_T
                    if (old_block_bunch_pos <= new_block_bunch_pos)
                    {
                        bbslice_T_N->marked_begin = new_block_bunch_pos;        assert(new_block_bunch_pos->slice == bbslice_T_N);
                        std::swap(old_block_bunch_pos->pred,
                                                    new_block_bunch_pos->pred); assert(nullptr != old_block_bunch_pos->pred->action_block->succ);
                        old_block_bunch_pos->pred->action_block->succ->
                                             block_bunch = old_block_bunch_pos; assert(new_block_bunch_pos->pred->action_block->succ == succ_iter);
                        succ_iter->block_bunch = new_block_bunch_pos;           // add_work(succ_iter->block_bunch->pred, ...) -- subsumed in the call below
                    }
                }                                                               else assert(&*bbslice_T_N==&*bbslice_Tprime_R || bbslice_T_N==bbslice_R_tau_U);
            }                                                                   mCRL2complexity(s, add_work(
                                                                                                    check_complexity::handle_new_noninert_transns, 1U), *this);
        }
        while (++s_iter < block_N->nonbottom_begin);

        return block_Rprime;
    }
};

///@} (end of group part_refine)
                                                                                #if !defined(NDEBUG) || defined(COUNT_WORK_BALANCE)
                                                                                    namespace bisim_dnj {

                                                                                    /// \brief moves temporary work counters to normal ones if the U-block is
                                                                                    /// smaller
                                                                                    /// \param block_U      the smaller subblock containing the U-states (can
                                                                                    ///                     be nullptr)
                                                                                    /// \param block_R      the larger subblock containing the R-states
                                                                                    /// \param partitioner  the partitioner data structure, used to write
                                                                                    ///                     diagnostic messages
                                                                                    template <class LTS_TYPE>
                                                                                    static void finalise_U_is_smaller(const block_t* const block_U,
                                                                                                            const block_t* const block_R,
                                                                                                            const bisim_partitioner_dnj<LTS_TYPE>& partitioner)
                                                                                    {
                                                                                        if (nullptr != block_U)
                                                                                        {
                                                                                            unsigned const max_U_block(check_complexity::log_n -
                                                                                                                     check_complexity::ilog2(block_U->size()));
                                                                                            // finalise work counters for the U-states and their transitions
                                                                                            const permutation_entry* s_iter(block_U->begin);
                                                                                            assert(s_iter < block_U->end);
                                                                                            do
                                                                                            {
                                                                                                const state_info_entry* const s(s_iter->st);
                                                                                                mCRL2complexity(s, finalise_work(
                                                                                                    check_complexity::split_U__find_predecessors_of_U_state,
                                                                                                    check_complexity::split__find_predecessors_of_R_or_U_state,
                                                                                                                                    max_U_block), partitioner);
                                                                                                assert(s != partitioner.part_tr.pred.back().target);
                                                                                                for (const pred_entry* pred_iter(s->pred_inert.begin);
                                                                                                                           s == pred_iter->target; ++pred_iter)
                                                                                                {
                                                                                                    mCRL2complexity(pred_iter, finalise_work(check_complexity::
                                                                                                            split_U__handle_transition_to_U_state,
                                                                                                            check_complexity::
                                                                                                            split__handle_transition_to_R_or_U_state,
                                                                                                                                    max_U_block), partitioner);
                                                                                                }
                                                                                                // Sometimes, inert transitions become transitions from R- to
                                                                                                // U-states; therefore, we also have to walk through the
                                                                                                // noninert predecessors of U-states:
                                                                                                assert(s != partitioner.part_tr.pred.front().target);
                                                                                                for (const pred_entry* pred_iter(s->pred_inert.begin);
                                                                                                                                  s == (--pred_iter)->target; )
                                                                                                {
                                                                                                    mCRL2complexity(pred_iter, finalise_work(check_complexity::
                                                                                                            split_U__handle_transition_to_U_state,
                                                                                                            check_complexity::
                                                                                                            split__handle_transition_to_R_or_U_state,
                                                                                                                                    max_U_block), partitioner);
                                                                                                }
                                                                                                assert(s != partitioner.part_tr.succ.front().
                                                                                                                                    block_bunch->pred->source);
                                                                                                for (const succ_entry* succ_iter(s->succ_inert.begin);
                                                                                                               s == (--succ_iter)->block_bunch->pred->source; )
                                                                                                {
                                                                                                    mCRL2complexity(succ_iter->block_bunch->pred,finalise_work(
                                                                                                          check_complexity::
                                                                                                          split_U__test_noninert_transitions,
                                                                                                          check_complexity::
                                                                                                          split__handle_transition_from_R_or_U_state,
                                                                                                                                    max_U_block), partitioner);
                                                                                                }
                                                                                            }
                                                                                            while (++s_iter < block_U->end);
                                                                                        }
                                                                                        // cancel work counters for the R-states and their transitions, and
                                                                                        // also account for work done in the U-coroutine on R-states
                                                                                        const permutation_entry* s_iter(block_R->begin);
                                                                                        assert(s_iter < block_R->end);
                                                                                        do
                                                                                        {
                                                                                            const state_info_entry* const s(s_iter->st);
                                                                                            mCRL2complexity(s, cancel_work(check_complexity::
                                                                                                          split_R__find_predecessors_of_R_state), partitioner);
                                                                                            assert(s != partitioner.part_tr.pred.back().target);
                                                                                            for (const pred_entry* pred_iter(s->pred_inert.begin);
                                                                                                                           s == pred_iter->target; ++pred_iter)
                                                                                            {
                                                                                                mCRL2complexity(pred_iter, cancel_work(check_complexity::
                                                                                                          split_R__handle_transition_to_R_state), partitioner);
                                                                                            }
                                                                                            assert(s !=
                                                                                                   partitioner.part_tr.succ.front().block_bunch->pred->source);
                                                                                            for (const succ_entry* succ_iter(s->succ_inert.begin);
                                                                                                               s == (--succ_iter)->block_bunch->pred->source; )
                                                                                            {
                                                                                                mCRL2complexity(succ_iter->block_bunch->pred, cancel_work(
                                                                                                        check_complexity::
                                                                                                        split_R__handle_transition_from_R_state), partitioner);
                                                                                                // the following counter measures work done in the
                                                                                                // U-coroutine that found R-states.
                                                                                                mCRL2complexity(succ_iter->block_bunch->pred,finalise_work(
                                                                                                       check_complexity::
                                                                                                       split_U__test_noninert_transitions,
                                                                                                       check_complexity::
                                                                                                       split__test_noninert_transitions_found_new_bottom_state,
                                                                                                                                             1U), partitioner);
                                                                                            }
                                                                                        }
                                                                                        while (++s_iter < block_R->end);

                                                                                        check_complexity::check_temporary_work();
                                                                                        (void) partitioner; // avoid unused variable warning
                                                                                    }

                                                                                    /// \brief moves temporary work counters to normal ones if the R-block is
                                                                                    /// smaller
                                                                                    /// \param block_U      the larger subblock containing the U-states
                                                                                    /// \param block_R      the smaller but non-empty subblock containing the
                                                                                    ///                     R-states
                                                                                    /// \param partitioner  the partitioner data structure, used to write
                                                                                    ///                     diagnostic messages
                                                                                    template <class LTS_TYPE>
                                                                                    static void finalise_R_is_smaller(const block_t* const block_U,
                                                                                                            const block_t* const block_R,
                                                                                                            const bisim_partitioner_dnj<LTS_TYPE>& partitioner)
                                                                                    {
                                                                                        unsigned const max_R_block(check_complexity::log_n -
                                                                                                                     check_complexity::ilog2(block_R->size()));
                                                                                        // cancel work counters for the U-states and their transitions
                                                                                        const permutation_entry* s_iter(block_U->begin);
                                                                                        assert(s_iter < block_U->end);
                                                                                        do
                                                                                        {
                                                                                            const state_info_entry* const s(s_iter->st);
                                                                                            mCRL2complexity(s, cancel_work(check_complexity::
                                                                                                          split_U__find_predecessors_of_U_state), partitioner);
                                                                                            assert(s != partitioner.part_tr.pred.back().target);
                                                                                            for (const pred_entry* pred_iter(s->pred_inert.begin);
                                                                                                                           s == pred_iter->target; ++pred_iter)
                                                                                            {
                                                                                                mCRL2complexity(pred_iter, cancel_work(check_complexity::
                                                                                                          split_U__handle_transition_to_U_state), partitioner);
                                                                                            }
                                                                                            // Sometimes, inert transitions become transitions from R- to
                                                                                            // U-states; therefore, we also have to walk through the
                                                                                            // noninert predecessors of U-states:
                                                                                            assert(s != partitioner.part_tr.pred.front().target);
                                                                                            for (const pred_entry* pred_iter(s->pred_inert.begin);
                                                                                                                                  s == (--pred_iter)->target; )
                                                                                            {
                                                                                                mCRL2complexity(pred_iter, cancel_work(check_complexity::
                                                                                                          split_U__handle_transition_to_U_state), partitioner);
                                                                                            }
                                                                                            assert(s !=
                                                                                                 partitioner.part_tr.succ.front().block_bunch->pred->source);
                                                                                            for (const succ_entry* succ_iter(s->succ_inert.begin);
                                                                                                               s == (--succ_iter)->block_bunch->pred->source; )
                                                                                            {
                                                                                                mCRL2complexity(succ_iter->block_bunch->pred, cancel_work(
                                                                                                             check_complexity::
                                                                                                             split_U__test_noninert_transitions), partitioner);
                                                                                            }
                                                                                        }
                                                                                        while (++s_iter < block_U->end);
                                                                                        // finalise work counters for the R-states and their transitions
                                                                                        s_iter = block_R->begin;
                                                                                        assert(s_iter < block_R->end);
                                                                                        do
                                                                                        {
                                                                                            const state_info_entry* const s(s_iter->st);
                                                                                            mCRL2complexity(s, finalise_work(
                                                                                                check_complexity::split_R__find_predecessors_of_R_state,
                                                                                                check_complexity::split__find_predecessors_of_R_or_U_state,
                                                                                                                                    max_R_block), partitioner);
                                                                                            assert(s != partitioner.part_tr.pred.back().target);
                                                                                            for (const pred_entry* pred_iter(s->pred_inert.begin);
                                                                                                                           s == pred_iter->target; ++pred_iter)
                                                                                            {
                                                                                                mCRL2complexity(pred_iter, finalise_work(
                                                                                                    check_complexity::split_R__handle_transition_to_R_state,
                                                                                                    check_complexity::split__handle_transition_to_R_or_U_state,
                                                                                                                                    max_R_block), partitioner);
                                                                                            }
                                                                                            assert(s !=
                                                                                                   partitioner.part_tr.succ.front().block_bunch->pred->source);
                                                                                            for (const succ_entry* succ_iter(s->succ_inert.begin);
                                                                                                               s == (--succ_iter)->block_bunch->pred->source; )
                                                                                            {
                                                                                                mCRL2complexity(succ_iter->block_bunch->pred, finalise_work(
                                                                                                        check_complexity::
                                                                                                        split_R__handle_transition_from_R_state,
                                                                                                        check_complexity::
                                                                                                        split__handle_transition_from_R_or_U_state,
                                                                                                                                    max_R_block), partitioner);
                                                                                                // the following counter actually is work done in the
                                                                                                // U-coroutine that found R-states.
                                                                                                mCRL2complexity(succ_iter->block_bunch->pred, cancel_work(
                                                                                                             check_complexity::
                                                                                                             split_U__test_noninert_transitions), partitioner);
                                                                                            }
                                                                                        }
                                                                                        while (++s_iter < block_R->end);
                                                                                        check_complexity::check_temporary_work();
                                                                                        (void) partitioner; // avoid unused variable warning
                                                                                    }

                                                                                    } // end namespace bisim_dnj
                                                                                #endif





/* ************************************************************************* */
/*                                                                           */
/*                             I N T E R F A C E                             */
/*                                                                           */
/* ************************************************************************* */





/// \defgroup part_interface
/// \brief nonmember functions serving as interface with the rest of mCRL2
/// \details These functions are copied, almost without changes, from
/// liblts_bisim_gw.h, which was written by Anton Wijs.
///@{

/// \brief Reduce transition system l with respect to strong or
/// (divergence-preserving) branching bisimulation.
/// \param[in,out] l                   The transition system that is reduced.
/// \param         branching           If true branching bisimulation is
///                                    applied, otherwise strong bisimulation.
/// \param         preserve_divergence Indicates whether loops of internal
///                                    actions on states must be preserved.  If
///                                    false these are removed.  If true these
///                                    are preserved.
template <class LTS_TYPE>
void bisimulation_reduce_dnj(LTS_TYPE& l, bool const branching = false,
                                        bool const preserve_divergence = false)
{
    if (1 >= l.num_states())
    {
        // LTSs with 1 state also need to be reduced because some users call
        // bisimulation minimisation just to remove duplicated transitions.
        mCRL2log(log::warning) << "There is only 1 state in the LTS. It is "
                "not guaranteed that branching bisimulation minimisation runs "
                "in time O(m log n).\n";
    }
    // Line 2.1: Find tau-SCCs and contract each of them to a single state
    const std::clock_t start_SCC=std::clock();
    mCRL2log(log::verbose) << "Start SCC\n";
    if (branching)
    {
        timer.start("scc_reduce");
        scc_reduce(l, preserve_divergence);
        // If only 1 state remains after this contraction, we are already
        // finished because scc_reduce() also removes duplicated transitions.
        timer.finish("scc_reduce");
        if (1 >= l.num_states())  return;
    }
    // Now apply the branching bisimulation reduction algorithm.  If there
    // are no taus, this will automatically yield strong bisimulation.
    const std::clock_t start_part=std::clock();
    mCRL2log(log::verbose) << "Start Partitioning\n";
    bisim_partitioner_dnj<LTS_TYPE> bisim_part(l, branching,
                                                          preserve_divergence);
    timer.finish("reduction");

    // Assign the reduced LTS
    const std::clock_t end_part=std::clock();
    mCRL2log(log::verbose) << "Start finalizing\n";
    bisim_part.finalize_minimized_LTS();

    if (mCRL2logEnabled(log::verbose))
    {
        const std::clock_t end_finalizing=std::clock();
        const int prec=std::lrint(std::log10(CLOCKS_PER_SEC)+0.19897000433602);
            // For example, if CLOCKS_PER_SEC>=     20: >=2 digits
            //              If CLOCKS_PER_SEC>=    200: >=3 digits
            //              If CLOCKS_PER_SEC>=2000000: >=7 digits

        double runtime[5];
        runtime[0]=(double) (end_finalizing                        -                        start_SCC)/CLOCKS_PER_SEC; // total time
        runtime[1]=(double) (                                                    start_part-start_SCC)/CLOCKS_PER_SEC;
        runtime[2]=(double) (                        bisim_part.end_initial_part-start_part          )/CLOCKS_PER_SEC;
        runtime[3]=(double) (               end_part-bisim_part.end_initial_part                     )/CLOCKS_PER_SEC;
        runtime[4]=(double) (end_finalizing-end_part                                                 )/CLOCKS_PER_SEC;
        if (runtime[0]>=60.0)
        {
            int min[sizeof(runtime)/sizeof(runtime[0])];
            for (unsigned i = 0; i < sizeof(runtime)/sizeof(runtime[0]); ++i)
            {
                min[i] = trunc(runtime[i] / 60.0);
                runtime[i] -= 60 * min[i];
            }
            if (min[0]>=60)
            {
                int h[sizeof(runtime)/sizeof(runtime[0])];
                for (unsigned i=0; i < sizeof(runtime)/sizeof(runtime[0]); ++i)
                {
                    h[i] = min[i] / 60;
                    min[i] %= 60;
                }
                int width = trunc(log10(h[0])) + 1;

                mCRL2log(log::verbose) << std::fixed << std::setprecision(prec)
                    << "Time spent on contracting SCCs: " << std::setw(width) << h[1] << "h " << std::setw(2) << min[1] << "min " << std::setw(prec+3) << runtime[1] << "s\n"
                       "Time spent on initial partition:" << std::setw(width) << h[2] << "h " << std::setw(2) << min[2] << "min " << std::setw(prec+3) << runtime[2] << "s\n"
                       "Time spent on refining:         " << std::setw(width) << h[3] << "h " << std::setw(2) << min[3] << "min " << std::setw(prec+3) << runtime[3] << "s\n"
                       "Time spent on finalizing:       " << std::setw(width) << h[4] << "h " << std::setw(2) << min[4] << "min " << std::setw(prec+3) << runtime[4] << "s\n"
                       "Total CPU time:                 " << std::setw(width) << h[0] << "h " << std::setw(2) << min[0] << "min " << std::setw(prec+3) << runtime[0] << "s\n"
                    << std::defaultfloat;
            }
            else
            {
                mCRL2log(log::verbose) << std::fixed << std::setprecision(prec)
                    << "Time spent on contracting SCCs: " << std::setw(2) << min[1] << "min " << std::setw(prec+3) << runtime[1] << "s\n"
                       "Time spent on initial partition:" << std::setw(2) << min[2] << "min " << std::setw(prec+3) << runtime[2] << "s\n"
                       "Time spent on refining:         " << std::setw(2) << min[3] << "min " << std::setw(prec+3) << runtime[3] << "s\n"
                       "Time spent on finalizing:       " << std::setw(2) << min[4] << "min " << std::setw(prec+3) << runtime[4] << "s\n"
                       "Total CPU time:                 " << std::setw(2) << min[0] << "min " << std::setw(prec+3) << runtime[0] << "s\n"
                    << std::defaultfloat;
            }
        }
        else
        {
            mCRL2log(log::verbose) << std::fixed << std::setprecision(prec)
                << "Time spent on contracting SCCs: " << std::setw(prec+3) << runtime[1] << "s\n"
                   "Time spent on initial partition:" << std::setw(prec+3) << runtime[2] << "s\n"
                   "Time spent on refining:         " << std::setw(prec+3) << runtime[3] << "s\n"
                   "Time spent on finalizing:       " << std::setw(prec+3) << runtime[4] << "s\n"
                   "Total CPU time:                 " << std::setw(prec+3) << runtime[0] << "s\n"
                << std::defaultfloat;
        }
    }
}


/// \brief Checks whether the two initial states of two LTSs are strong or
/// (divergence-preserving) branching bisimilar.
/// \details This routine uses the O(m log n) branching bisimulation algorithm
/// developed in 2018 by David N. Jansen.  It runs in O(m log n) time and uses
/// O(m) memory, where n is the number of states and m is the number of
/// transitions.
///
/// The LTSs l1 and l2 are not usable anymore after this call.
/// \param[in,out] l1                  A first transition system.
/// \param[in,out] l2                  A second transistion system.
/// \param         branching           If true branching bisimulation is used,
///                                    otherwise strong bisimulation is
///                                    applied.
/// \param         preserve_divergence If true and branching is true, preserve
///                                    tau loops on states.
/// \param         generate_counter_examples  (non-functional, only in the
///                                    interface for historical reasons)
/// \returns True iff the initial states of the transition systems l1 and l2
/// are ((divergence-preserving) branching) bisimilar.
template <class LTS_TYPE>
bool destructive_bisimulation_compare_dnj(LTS_TYPE& l1, LTS_TYPE& l2,
        bool const branching = false, bool const preserve_divergence = false,
        bool const generate_counter_examples = false,
        const std::string& /*counter_example_file*/ = "",
        bool /*structured_output*/ = false)
{
    if (generate_counter_examples)
    {
        mCRL2log(log::warning) << "The JGKW20 branching bisimulation "
                              "algorithm does not generate counterexamples.\n";
    }
    std::size_t init_l2(l2.initial_state() + l1.num_states());
    detail::merge(l1, std::move(l2));
    l2.clear(); // No use for l2 anymore.

    // Line 2.1: Find tau-SCCs and contract each of them to a single state
    if (branching)
    {
        scc_partitioner<LTS_TYPE> scc_part(l1);
        scc_part.replace_transition_system(preserve_divergence);
        init_l2 = scc_part.get_eq_class(init_l2);
    }                                                                           else  assert(!preserve_divergence);
                                                                                assert(1 < l1.num_states());
    bisim_partitioner_dnj<LTS_TYPE> bisim_part(l1, branching,
                                                          preserve_divergence);

    return bisim_part.in_same_class(l1.initial_state(), init_l2);
}


/// \brief Checks whether the two initial states of two LTSs are strong or
/// (divergence-preserving) branching bisimilar.
/// \details The LTSs l1 and l2 are first duplicated and subsequently reduced
/// modulo bisimulation.  If memory is a concern, one could consider to use
/// destructive_bisimulation_compare().  This routine uses the O(m log n)
/// branching bisimulation algorithm developed in 2018 by David N. Jansen.  It
/// runs in O(m log n) time and uses O(m) memory, where n is the number of
/// states and m is the number of transitions.
/// \param l1                  A first transition system.
/// \param l2                  A second transistion system.
/// \param branching           If true branching bisimulation is used,
///                            otherwise strong bisimulation is applied.
/// \param preserve_divergence If true and branching is true, preserve tau
///                            loops on states.
/// \retval True iff the initial states of the transition systems l1 and l2
/// are ((divergence-preserving) branching) bisimilar.
template <class LTS_TYPE>
inline bool bisimulation_compare_dnj(const LTS_TYPE& l1, const LTS_TYPE& l2,
          bool const branching = false, bool const preserve_divergence = false)
{
    LTS_TYPE l1_copy(l1);
    LTS_TYPE l2_copy(l2);
    return destructive_bisimulation_compare_dnj(l1_copy, l2_copy, branching,
                                                          preserve_divergence);
}

///@} (end of group part_interface)

} // end namespace detail
} // end namespace lts
} // end namespace mcrl2

#endif // ifndef LIBLTS_BISIM_DNJ_H
