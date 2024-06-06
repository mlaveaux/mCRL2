// Author(s): Jan Friso Groote
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/// \file lts/detail/liblts_bisim_gj.h
///
/// \brief O(m log n)-time branching bisimulation algorithm similar to liblts_bisim_dnj.h
///        which does not use bunches, i.e., partitions of transitions. This algorithm
///        should be slightly faster, but in particular use less memory than liblts_bisim_dnj.h.
///        Otherwise the functionality is exactly the same. 
///

#ifndef LIBLTS_BISIM_GJ_H
#define LIBLTS_BISIM_GJ_H

#include <coroutine>

#include "mcrl2/lts/detail/liblts_scc.h"
#include "mcrl2/lts/detail/liblts_merge.h"
// #include "mcrl2/lts/detail/coroutine.h"
// #include "mcrl2/lts/detail/check_complexity.h"
// #include "mcrl2/lts/detail/fixed_vector.h"


namespace mcrl2
{
namespace lts
{
namespace detail
{
namespace bisimulation_gj
{

// Forward declaration.
struct block_type;
struct transition_type;

typedef std::size_t block_index;
typedef std::size_t transition_index;
typedef std::size_t state_index;
constexpr transition_index null_transition=-1;

// Below the four main data structures are listed.
struct state_type
{
  block_index block=0;
  std::vector<transition_index>::iterator start_incoming_transitions;
  std::vector<transition_index>::iterator start_outgoing_transitions;
  std::vector<transition_index>::iterator start_outgoing_non_inert_transitions;
};

struct transition_type
{
  // The position of the transition type corresponds to m_aut.transitions(). 
  // std::size_t from, label, to are found in m_aut.transitions.
  std::forward_list<transition_index > :: const_iterator transitions_per_block_to_constellation;
  transition_index next_L_B_C_element;
  transition_index previous_L_B_C_element;
  std::vector<std::size_t>::iterator trans_count;
};

struct block_type
{
  state_index start_bottom_states;
  state_index start_non_bottom_states;
  std::forward_list< transition_index > block_to_constellation;

  block_type(state_index i)
    : start_bottom_states(i),
      start_non_bottom_states(i)
  {}
};


struct constellation_type
{
  block_index start_block;
  block_index end_block;

  constellation_type(const block_index start, const block_index end)
   : start_block(start), 
     end_block(end)
  {}
};


} // end namespace bisimulation_gj
} // end namespace detail


/*=============================================================================
=                                 main class                                  =
=============================================================================*/


using namespace mcrl2::lts::detail::bisimulation_gj;



/// \class bisim_partitioner_gj
/// \brief implements the main algorithm for the branching bisimulation
/// quotient
template <class LTS_TYPE>
class bisim_partitioner_gj
{
  protected:
    typedef typename LTS_TYPE::labels_sizes_type label_index;
    typedef typename LTS_TYPE::states_sizes_type state_index;
    typedef std::unordered_set<state_index> set_of_states_type;
    typedef std::unordered_set<size_t> set_of_constellations_type;
    typedef std::unordered_map<label_index, set_of_states_type > states_per_action_label_type;
    typedef std::unordered_map<block_index, set_of_states_type > states_per_block_type;

    /// \brief automaton that is being reduced
    LTS_TYPE& m_aut;
    
    // Generic data structures.
    std::vector<state_type> m_states;
    std::vector<std::reference_wrapper<transition_type>> m_incoming_transitions;
    std::vector<std::reference_wrapper<transition_type>> m_outgoing_transitions;
    std::vector<transition_type> m_transitions;
    std::deque<std::size_t> m_state_to_constellation_count;
    std::vector<std::reference_wrapper<state_type>> m_states_in_blocks;
    std::vector<block_type> m_blocks;
    std::vector<constellation_type> m_constellations;

    /// \brief true iff branching (not strong) bisimulation has been requested
    const bool m_branching;
  
    /// \brief true iff divergence-preserving branching bisimulation has been
    /// requested
    /// \details Note that this field must be false if strong bisimulation has
    /// been requested.  There is no such thing as divergence-preserving strong
    /// bisimulation.
    const bool m_preserve_divergence;

    /// The following variable contains all non trivial constellations.
    set_of_constellations_type non_trivial_constellations;

    std::size_t number_of_states_in_block(std::size_t block_index) const
    {
      if (m_blocks.size()==block_index+1) // This is the last block.
      {
        return m_states_in_blocks.size()-m_blocks[block_index].start_bottom_states; 
      }
      return m_blocks[block_index+1].start_bottom_states-m_blocks[block_index].start_bottom_states; 
    }

  public:
    /// \brief constructor
    /// \details The constructor constructs the data structures and immediately
    /// calculates the partition corresponding with the bisimulation quotient.
    /// It destroys the transitions on the LTS (to save memory) but does not
    /// adapt the LTS to represent the quotient's transitions.
    /// It is assumed that there are no tau-loops in aut.
    /// \param aut                 LTS that needs to be reduced
    /// \param branching           If true branching bisimulation is used,
    ///                                otherwise strong bisimulation is
    ///                                applied.
    /// \param preserve_divergence If true and branching is true, preserve
    ///                                tau loops on states.
    bisim_partitioner_gj(LTS_TYPE& aut, 
                         const bool branching = false,
                         const bool preserve_divergence = false)
      : m_aut(aut),
        m_states(aut.num_states()),
        m_transitions(aut.num_transitions()),
        m_blocks(1),
        m_constellations(1,constellation_type(0,1)),   // Algorithm 1, line 1.2.
        m_branching(branching),
        m_preserve_divergence(preserve_divergence)
    {                                                                           
      assert(m_branching || !m_preserve_divergence);
      create_initial_partition();        
      refine_partition_until_it_becomes_stable();
    }


    /// \brief Calculate the number of equivalence classes
    /// \details The number of equivalence classes (which is valid after the
    /// partition has been constructed) is equal to the number of states in the
    /// bisimulation quotient.
    state_type num_eq_classes() const
    {
      return m_blocks.size();
    }


    /// \brief Get the equivalence class of a state
    /// \details After running the minimisation algorithm, this function
    /// produces the number of the equivalence class of a state.  This number
    /// is the same as the number of the state in the minimised LTS to which
    /// the original state is mapped.
    /// \param s state whose equivalence class needs to be found
    /// \returns sequence number of the equivalence class of state s
    state_type get_eq_class(const state_type s) const
    {
      assert(s<m_blocks.size());
      return m_states[s].block;
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
      std::unordered_set<transition> T;
      for(const transition& t: m_aut.transitions())
      {
        T.insert(transition(get_eq_class(t.from()), t.label(), t.to()));
      }
      for (const transition t: T)
      {
        m_aut.add_transition(t);
      }

      // Merge the states, by setting the state labels of each state to the
      // concatenation of the state labels of its equivalence class.

      if (m_aut.has_state_info())   /* If there are no state labels this step is not needed */
      {
        /* Create a vector for the new labels */
        std::vector<typename LTS_TYPE::state_label_t> new_labels(num_eq_classes());

    
        for(std::size_t i=0; i<m_aut.num_states(); ++i)
        {
          const state_type new_index(get_eq_class(i));
          new_labels[new_index]=new_labels[new_index]+m_aut.state_label(i);
        }

        m_aut.set_num_states(num_eq_classes());
        for (std::size_t i=0; i<num_eq_classes(); ++i)
        {
          m_aut.set_state_label(i, new_labels[i]);
        }
      }
      else
      {
        m_aut.set_num_states(num_eq_classes());
      }

      m_aut.set_initial_state(get_eq_class(m_aut.initial_state()));
    }


    /// \brief Check whether two states are in the same equivalence class.
    /// \param s first state that needs to be compared.
    /// \param t second state that needs to be compared.
    /// \returns true iff the two states are in the same equivalence class.
    bool in_same_class(state_type const s, state_type const t) const
    {
      return get_eq_class(s) == get_eq_class(t);
    }
  protected:

    /*--------------------------- main algorithm ----------------------------*/


    /*----------------- SplitB -- Algorithm 3 of [GJ 2024] -----------------*/

    void swap_states_in_states_in_block(block_index pos1, block_index pos2, block_index pos3)
    {
      block_index temp=m_states_in_block[pos3]
      m_states_in_block[pos3]=m_states_in_block[pos2]
      m_blocks[m_states_in_block[pos3]].ref_states_in_block=pos3;
      m_states_in_block[pos2]=m_states_in_block[pos1]
      m_blocks[m_states_in_block[pos2]].ref_states_in_block=pos2;
      m_states_in_block[pos1]temp;
      m_blocks[m_states_in_block[pos1]].ref_states_in_block=pos1;
    }
    
    void simpleSplitB(const block_index B, const std::unordered_set<state_index>& M)
    {
      std::unordered_map<size_t> count;
      std::vector<state_index> U, U_todo;
      std::vector<state_index> R, R_todo;
      typedef enum { initializing, state_checking, aborted } status_type;
      status_type U_status=initializing, R_status=state_checking;
      block_index R_i=m_blocks[B].start_bottom_states;

      if (2*M.size()<=number_states_in_block())
      {
        R_todo:=M;
      }
      else R_status=aborted;

      // start coroutines.
      while (true)
      {
        if (R_status=state_checking) 
        {
          if (R_todo.empty())
          {
            // split_block B into R and B\R.
            
            m_blocks.emplace_back(B.bottom_states);
            const block_index new_block_index=m_blocks.size()-1;
            for(state_index s: R)
            {
              m_states[s].block=new_block_index;
              std::size_t pos=m_states[s].ref_states_in_block;
              if (pos>=m_blocks[B].non_bottom_states) // the state is a non bottom state.
              {
                swap_states_in_states_in_block(pos,m_blocks[B].bottom_states,m_blocks[B].non_bottom_states);
              }
              else // the state is a non bottom state
              {
                swap_states_in_states_in_block(pos,m_blocks[new_block_index].bottom_states,m_blocks[B].non_bottom_states);
              }
              m_blocks[B].bottom_states++;
              m_blocks[B].non_bottom_states++;
            }
          }
        }
      }
    }

    template <class SET_OF_STATES>
    void splitB(const block_type& B, const SET_OF_STATES& M, const detail::label_type a, const constellation_type& C)
    {
      
    }


    void create_initial_partition()
    {
      mCRL2log(log::verbose) << "An O(m log n) "
           << (m_branching ? (m_preserve_divergence
                                         ? "divergence-preserving branching "
                                         : "branching ")
                         : "")
           << "bisimulation partitioner created for " << m_aut.num_states()
           << " states and " << m_aut.num_transitions() << " transitions.\n";
      // Initialisation.

      // Initialise m_incoming_transitions and m_transitions.trans_count and m_transitions.transitions_per_block_to_constellation.
      typedef std::unordered_multimap<typename std::pair<typename LTS_TYPE::states_sizes_type, typename LTS_TYPE::labels_size_type>, 
                                      transition&> temporary_store_type;
      temporary_store_type temporary_store;
      for(const transition& t: m_aut.transitions())
      {
        temporary_store[std::pair(t.from(),t.label())]=t;
      }
      state_index current_from_state=-1;
      label_index current_label=-1;
      m_incoming_transitions.reserve(m_aut.num_transitions());
      for(auto [_,t]: temporary_store)
      {
        m_incoming_transitions.emplace(t);
        if (t.from()!=current_from_state || t.label()!=current_label)
        {
          m_state_to_constellation_count.emplace_back(1);
          m_blocks[0].block_to_constellation.push_front(std::list(1,t));
          current_label=t.label();
        }
        else
        {
          assert(m_state_to_constellation_count.size()>0);
          m_state_to_constellation_count.back()++;
          m_blocks[0].block_to_constellation.front().push_front(t);
        }
        if (t.from()!=current_from_state)
        {
          m_states[t.from()].start_incoming_transitions=m_incoming_transitions.end()-1;
          current_from_state=t.from();
        }
        std::size_t transition_index=std::distance(t,m_aut.transitions().begin());
        m_transitions[transition_index].trans_count=m_state_to_constellation_count.end()-1;
        m_transitions[transition_index].transitions_per_block_to_constellation=m_blocks[0].block_to_constellation.front().begin();
      }
      temporary_store.clear();
      
      // Initialise m_outgoing_transitions and
      // initialise m_states_in_blocks, together with start_bottom_states start_non_bottom_states in m_blocks.
      std::vector<bool> state_has_outgoing_tau(m_aut.num_states(),false);
      for(const transition& t: m_aut.transitions())
      {
        temporary_store[std::pair(t.to(),t.label())]=t;
        if (m_aut.is_tau(t))
        {
          state_has_outgoing_tau[t.from()]=true;
        }
      }
      m_outgoing_transitions.reserve(m_aut.num_transitions());
      typename LTS_TYPE::states_sizes_type current_to_state=-1;
      for(auto [_,t]: temporary_store)
      {
        m_outgoing_transitions.emplace(t);
        if (t.to()!=current_to_state)
        {
          m_states[t.to()].start_outgoing_transitions=m_outgoing_transitions.end()-1;
          current_to_state=t.to();
        }
      }
      temporary_store=temporary_store_type(); // release memory. 

      m_states_in_blocks.reserve(m_aut.num_states());
      std::size_t i=0;
      for(bool b: state_has_outgoing_tau)
      {
        if (b)
        {
          m_states_in_blocks.emplace_back(i);
        }
        i++;
      }
      m_blocks[0].start_bottom_states=0;
      m_blocks[0].start_non_bottom_states=i;
      i=0;
      for(bool b: state_has_outgoing_tau)
      {
        if (!b)
        {
          m_states_in_blocks.emplace_back(i);
        }
        i++;
      }

      // The data structures are now initialized.
      // The following implements line 1.3 of Algorithm 1. 
      states_per_action_label_type states_per_action_label;
      for(const transition& t: m_aut.transitions())
      {
        states_per_action_label[t.label()]=states_per_action_label[t.label()].insert(t.from());
      }

      for(const set_of_states_type& stateset: states_per_action_label)
      {
        for(const set_of_states_type& M: stateset)
        {
          states_per_block_type Bprime;
          for(const state_index s: M)
          {
            Bprime[m_states[s].block].insert(s);
          }
          
          for(auto [block_index, split_states]: Bprime)
          {
            // Check whether the block B, indexed by block_index, can be split.
            // This means that the bottom states of B are not all in the split_states.
            const block_type& B=m_blocks[block_index];
            for(state_index i=B.start_bottom_states; i<B.start_non_bottom_states; ++i)
            {
              if (!split_states.contains(i))
              {
                simpleSplitB(block_index, split_states);
                i=B.start_non_bottom_states; // This means break the loop.
              }
            }
          }
        }
      }
    }
 
    void refine_partition_until_it_becomes_stable()
    {
      // This represents the while loop in Algorithm 1 from line 1.6 to 1.25.

      // Algorithm 1, line 1.6.
      while (!non_trivial_constellations.empty())
      {
        const set_of_constellations_type::const_iterator i=non_trivial_constellations.begin();
        std::size_t constellation_index= *i;
        non_trivial_constellations.extract(i);

        // Algorithm 1, line 1.7.
        std::size_t index_block_B=m_constellations[constellation_index].start_block;
        if (number_of_states_in_block(index_block_B)<number_of_states_in_block(index_block_B+1))
        {
          m_states_in_blocks[index_block_B].swap(m_states_in_blocks[index_block_B+1]);
        }
        
        // Algorithm 1, line 1.8.
        m_constellations[constellation_index].start_block=index_block_B+1;
        if (m_constellations[constellation_index].start_block!=m_constellations[constellation_index].end_block) // Constellation is not trivial.
        {
          non_trivial_constellations.insert(constellation_index);
        }
        m_constellations.emplace_back(index_block_B, index_block_B);
        // Here the variables block_to_constellation and the doubly linked list L_B->C in blocks must be still be updated.
        // Moreover, the variable m_state_to_constellation_count in transitions requires updating.
        // This happens below.

        // Algorithm 1, line 1.9.
        states_per_action_label_type calM;
        typedef std::unordered_map<std::pair<state_index, label_index>, std::size_t> state_label_to_size_t_map;
        state_label_to_size_t_map newly_created_state_to_constellation_count_entry;
        
 
        // Walk through all states in block B
        for(std::vector<std::reference_wrapper<state_type>>::iterator i=m_states_in_blocks[index_block_B].start_bottom_states;
               i!=m_states_in_blocks[index_block_B+1].start_bottom_states; ++i)
        {
          // and visit the incoming transitions. 
          for(std::vector<transition_index>::iterator j=i->get().start_incoming_transitions;
              j!=(i+1)->get().start_incoming_transitions; ++j)
          {
            const transition& t=m_aut.transitions[*j];
            calM[t.label()].insert(t.from());
            // Update m_state_to_constellation_count.
            if (m_states[t.to()].block==index_block_B)
            {
              const std::size_t new_position=m_state_to_constellation_count.size();
              const std::size_t found_position=
                                   newly_created_state_to_constellation_count_entry.try_emplace(
                                                                      std::pair(t.from(),t.label()),
                                                                      new_position)->second;
              if (new_position==found_position)
              {
                m_state_to_constellation_count.push_back(1);
                (m_transitions[t].trans_count)--;
                m_transitions[t].trans_count=m_state_to_constellation_count.end()-1;
              }
            }
            // Update the doubly linked list L_B->C in blocks.
            if (m_states[t.from()].block==index_block_B)
            {
              std::forward_list<transition_index > :: iterator this_block_to_constellation=
                                            m_transitions[*j].transitions_per_block_to_constellation;
              std::forward_list<transition_index > :: iterator next_block_to_constellation=
                                            ++std::forward_list<transition_index > :: iterator(this_block_to_constellation);
              if (next_block_to_constellation==m_blocks[m_states[t.from()].block].end() ||
                  *next_block_to_constellation==null_transition ||
                  m_blocks[m_aut.transitions()[*next_block_to_constellation].to()]!=index_block_B ||
                  m_aut.transitions()[*next_block_to_constellation].label()!=t.label())
              { 
                // Make a new entry in the list next_block_to_constellation;
                m_blocks[m_states[m_transitions[*j].from()].block].block_to_constellation.insert_after(this_block_to_constellation, *j);
                // Move the current transition to the next list.
                // First check whether this_block_to_constellation contains exactly transition *j.
                // It must be replaced by a later or earlier element from the L_B_C_list.
                if (*this_block_to_constellation==*j)
                {
                  if (m_transitions[*j].next_L_B_C_element!=null_transition)
                  { // TODO: CHECK FOR INERTNESS
                    *this_block_to_constellation= *m_transitions[*j].next_L_B_C_element;
                  }
                  else if (m_transitions[*j].previous_L_B_C_element!=null_transition)
                  { 
                    *this_block_to_constellation= *m_transitions[*j].previous_L_B_C_element;
                  }
                  else
                  {
                    // This is the last element of this L_B_C_list. 
                    //
                    *this_block_to_constellation=null_transition; // TODO: move next list to this list. 
                  }
                }
                // Rewire the connections to insert in the new list. 
                next_block_to_constellation= ++std::forward_list<transition_index > :: iterator(this_block_to_constellation);
                if (*next_block_to_constellationXXX)
                {}

              }
            }
          }
        }
        
        // Algorithm 1, line 1.10.
        for(const set_of_states_type& M: calM)
        {
          states_per_block_type Bprime;
          for(const state_index s: M)
          {
            Bprime[]
          }
          // XXXXX Split existing blocks by stateset.
        }
      }
      
      
      
    }

};






/* ************************************************************************* */
/*                                                                           */
/*                             I N T E R F A C E                             */
/*                                                                           */
/* ************************************************************************* */





/// \defgroup part_interface
/// \brief nonmember functions serving as interface with the rest of mCRL2
/// \details These functions are copied, almost without changes, from
/// liblts_bisim_gw.h, which was written by Anton Wijs.


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
void bisimulation_reduce_gj(LTS_TYPE& l, const bool branching = false,
                                        const bool preserve_divergence = false)
{
    if (1 >= l.num_states())
    {
        mCRL2log(log::warning) << "There is only 1 state in the LTS. It is not "
                "guaranteed that branching bisimulation minimisation runs in "
                "time O(m log n).\n";
    }
    // Line 2.1: Find tau-SCCs and contract each of them to a single state
    if (branching)
    {
        scc_reduce(l, preserve_divergence);
    }

    // Now apply the branching bisimulation reduction algorithm.  If there
    // are no taus, this will automatically yield strong bisimulation.
    bisim_partitioner_gj<LTS_TYPE> bisim_part(l, branching,
                                                          preserve_divergence);

    // Assign the reduced LTS
    bisim_part.finalize_minimized_LTS();
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
bool destructive_bisimulation_compare_gj(LTS_TYPE& l1, LTS_TYPE& l2,
        const bool branching = false, const bool preserve_divergence = false,
        const bool generate_counter_examples = false,
        const std::string& /*counter_example_file*/ = "", bool /*structured_output*/ = false)
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
        detail::scc_partitioner<LTS_TYPE> scc_part(l1);
        scc_part.replace_transition_system(preserve_divergence);
        init_l2 = scc_part.get_eq_class(init_l2);
    }                                                                           else  assert(!preserve_divergence);
                                                                                assert(1 < l1.num_states());
    bisim_partitioner_gj<LTS_TYPE> bisim_part(l1, branching,
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
inline bool bisimulation_compare_gj(const LTS_TYPE& l1, const LTS_TYPE& l2,
          const bool branching = false, const bool preserve_divergence = false)
{
    LTS_TYPE l1_copy(l1);
    LTS_TYPE l2_copy(l2);
    return destructive_bisimulation_compare_gj(l1_copy, l2_copy, branching,
                                                          preserve_divergence);
}


} // end namespace lts
} // end namespace mcrl2

#endif // ifndef LIBLTS_BISIM_GJ_H