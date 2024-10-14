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

// TODO:
// Merge identifying whether there is a splitter and actual splitting (Done. No performance effect).
// Use BLC lists for the main split.
// Maintain a co-splitter and a splitter to enable a co-split.
// Eliminate two pointers in transition_type (done).
// JFG: Optimise swap. 

#ifndef LIBLTS_BISIM_GJ_H
#define LIBLTS_BISIM_GJ_H

#include <forward_list>
#include <deque>
#include "mcrl2/utilities/hash_utility.h"
#include "mcrl2/lts/detail/liblts_scc.h"
#include "mcrl2/lts/detail/liblts_merge.h"
#include "mcrl2/utilities/execution_timer.h"

#ifndef NDEBUG
#define CHECK_COMPLEXITY_GJ // Check whether coroutines etc. satisfy the O(m log n) time complexity constraint for the concrete input.
                            // Outcomment to disable. Works only in debug mode.
#endif

#ifdef CHECK_COMPLEXITY_GJ
  #include "mcrl2/lts/detail/check_complexity.h"
  // Using __VA_ARGS__ is not handled appropriately by MSVC.
  #define mCRL2complexity_gj(unit, call, info_for_debug) mCRL2complexity(unit, call, info_for_debug)
#else
  #define mCRL2complexity_gj(unit, call, info_for_debug)  do{}while(0)
#endif

namespace mcrl2
{
namespace lts
{
namespace detail
{

template <class LTS_TYPE> class bisim_partitioner_gj;

namespace bisimulation_gj
{

// Forward declaration.
struct block_type;
struct transition_type;
struct transition_pointer_pair;

typedef std::size_t state_index;
typedef std::size_t transition_index;
typedef std::size_t block_index;


typedef std::size_t label_index;
typedef std::size_t constellation_index;
typedef std::vector<transition_pointer_pair>::iterator outgoing_transitions_it;
typedef std::vector<transition_pointer_pair>::reverse_iterator outgoing_transitions_reverse_it;
typedef std::vector<transition_pointer_pair>::const_iterator outgoing_transitions_const_it;

constexpr transition_index null_transition=-1;
constexpr label_index null_action=-1;
constexpr state_index null_state=-1;
constexpr block_index null_block=-1;
constexpr transition_index undefined=-1;
constexpr transition_index Rmarked=-2;

// The function clear takes care that a container frees memory when it is cleared and it is large.
template <class CONTAINER>
void clear(CONTAINER& c)
{
  if (c.size()>1000) { c=CONTAINER(); } else { c.clear(); }
}

// Private linked list that uses less memory.

template <class T>
struct linked_list_node;

template <class T>
struct linked_list_iterator
{
  linked_list_node<T>* m_iterator=nullptr;

  linked_list_iterator()=default;

  linked_list_iterator(linked_list_node<T>* t)
    : m_iterator(t)
  {}

  linked_list_iterator operator ++()
  {
    assert(nullptr != m_iterator);
    *this=m_iterator->next();
    return m_iterator;
  }

  linked_list_node<T>& operator *() const
  {
    assert(nullptr != m_iterator);
    return *m_iterator;
  }

  linked_list_node<T>* operator ->() const
  {
    assert(nullptr != m_iterator);
    return m_iterator;
  }

  bool operator !=(linked_list_iterator other) const
  {
    return m_iterator!=other.m_iterator;
  }

  bool operator ==(linked_list_iterator other) const
  {
    return m_iterator==other.m_iterator;
  }

};

template <class T>
struct linked_list_node: public T
{
  typedef  linked_list_iterator<T> iterator;
  iterator m_next;
  iterator m_prev;

  linked_list_node(const T& t, iterator next, iterator prev)
   : T(t),
     m_next(next),
     m_prev(prev)
  {}

  iterator& next()
  {
    return m_next;
  }

  iterator& prev()
  {
    return m_prev;
  }

  T& content()
  {
    return static_cast<T&>(*this);
  }
};

template <class T>
struct global_linked_list_administration
{
  std::deque<linked_list_node<T> > m_content;
  linked_list_iterator<T> m_free_list=nullptr;
};

template <class T>
struct linked_list
{
  typedef linked_list_iterator<T> iterator;

  global_linked_list_administration<T>& glla() { static global_linked_list_administration<T> glla; return glla; }
  iterator m_initial_node=nullptr;

  iterator begin() const
  {
    return m_initial_node;
  }

  static iterator end()
  {
    return nullptr;
  }

  bool empty() const
  {
    return m_initial_node==nullptr;
  }

#ifndef NDEBUG
  bool check_linked_list() const
  {
    if (empty())
    {
      return true;
    }
    iterator i=m_initial_node;
    if (i->prev()!=nullptr)
    {
      return false;
    }
    while (i->next()!=nullptr)
    {
      if (i->next()->prev()!=i)
      {
        return false;
      }
      i=i->next();
      if (i->prev()->next()!=i)
      {
        return false;
      }
    }
    return true;
  }
#endif

  // Puts a new element after the current element indicated by pos, unless
  // pos==end(), in which it is put in front of the list.
  template <class... Args>
  iterator emplace_after(iterator pos, Args&&... args)
  {
    if (pos==nullptr)
    {
      return emplace_front(args...);
    }

    iterator new_position;
    if (glla().m_free_list==nullptr)
    {
      new_position=&glla().m_content.emplace_back(T(args...), pos->next(), pos);
    }
    else
    {
      // Take an element from the free list.
      new_position=glla().m_free_list;
      glla().m_free_list=glla().m_free_list->next();
      new_position->prev()=pos;
      new_position->next()=pos->next();;
      new_position->content()=T(args...);
    }
    if (pos->next()!=nullptr)
    {
      pos->next()->prev()=new_position;
    }
    pos->next()=new_position;

    return new_position;
  }

  template <class... Args>
  iterator emplace_front(Args&&... args)
  {
    iterator new_position;
    if (glla().m_free_list==nullptr)
    {
      new_position=&glla().m_content.emplace_back(T(args...), m_initial_node, nullptr);  // Deliver the address to new position.
    }
    else
    {
      // Take an element from the free list.
      new_position=glla().m_free_list;
      glla().m_free_list=glla().m_free_list->next();
      new_position->content()=T(args...);
      new_position->next()=m_initial_node;
      new_position->prev()=nullptr;
    }
    if (m_initial_node!=nullptr)
    {
      m_initial_node->prev()=new_position;
    }
    m_initial_node=new_position;

    return new_position;
  }

  //iterator push_front(const T& t)
  //{
  //  return emplace_front(t);
  //}

  // The function moves the element pointed at by from_pos (that is in the list
  // indicated by the 2nd parameter) just after position to_pos (that is in
  // this list). If to_pos == nullptr, move the element to the beginning of this
  // list.
  void splice_to_after(iterator const to_pos, linked_list<T>& other_list, iterator const from_pos)
  {
    // remove element from_pos from its original list
    if (from_pos->next() != nullptr)
    {
      assert(from_pos == from_pos->next()->prev());
      from_pos->next()->prev() = from_pos->prev();
    }
    if (from_pos->prev() != nullptr)
    {
      assert(from_pos == from_pos->prev()->next());
      from_pos->prev()->next() = from_pos->next();
    }
    else
    {
      assert(from_pos == other_list.m_initial_node);
      other_list.m_initial_node = from_pos->next();
    }
    // update the pointers of from_pos and insert from_pos into this list
    from_pos->prev() = to_pos;
    if (to_pos != nullptr)
    {
      from_pos->next() = to_pos->next();
      to_pos->next() = from_pos;
    }
    else
    {
      from_pos->next() = m_initial_node;
      m_initial_node = from_pos;
    }
    if (from_pos->next() != nullptr)
    {
      assert(to_pos == from_pos->next()->prev());
      from_pos->next()->prev() = from_pos;
    }
  }

  void erase(iterator const pos)
  {
    if (pos->next()!=nullptr)
    {
      assert(pos == pos->next()->prev());
      pos->next()->prev()=pos->prev();
    }
    if (pos->prev()!=nullptr)
    {
      assert(pos == pos->prev()->next());
      pos->prev()->next()=pos->next();
    }
    else
    {
      assert(pos == m_initial_node);
      m_initial_node=pos->next();
    }
    pos->next()=glla().m_free_list;
    glla().m_free_list=pos;
#ifndef NDEBUG
    pos->prev()=nullptr;
#endif
  }
};

struct transition_pointer_pair
{
  // David suggests to rename this struct to outgoing_transition_type, so it describes more clearly its function.

  // David suggests to change the following field to BLC_list_iterator ref_BLC_transitions to make it more clear what it is.
  transition_index transition;
  outgoing_transitions_it start_same_saC; // Refers to the last state with the same state, action and constellation,
                                          // unless it is the last, which refers to the first state.

  // The default initialiser does not initialize the fields of this struct.
  transition_pointer_pair()
  {}

  transition_pointer_pair(const transition_index t, const outgoing_transitions_it sssaC)
   : transition(t),
     start_same_saC(sssaC)
  {}
};

struct label_count_sum_triple
// David suggests to call this a pair, not a triple. The triple is not really used.
{
  transition_index label_counter=0;
  transition_index not_investigated=0;

  // The default initialiser does not initialize the fields of this struct.
  label_count_sum_triple()
  {}
};

/* struct label_count_sum_triple: label_count_sum_triple
{
  transition_index cumulative_label_counter=0;

  // The default initialiser does not initialize the fields of this struct.
  label_count_sum_triple()
  {}
}; */

class todo_state_vector
{
  std::size_t m_todo_indicator=0;
  std::vector<state_index> m_vec;

  public:
    void add_todo(const state_index s)
    {
      assert(!find(s));
      m_vec.push_back(s);
    }

    // Move a state from the todo part to the definitive vector.
    state_index move_from_todo()
    {
      assert(!todo_is_empty());
      assert(m_todo_indicator<m_vec.size());
      state_index result=m_vec[m_todo_indicator];
      m_todo_indicator++;
      return result;
    }

    std::size_t size() const
    {
      return m_vec.size();
    }

    std::size_t todo_is_empty() const
    {
      return m_vec.size()==m_todo_indicator;
    }

    std::size_t empty() const
    {
      return m_vec.empty();
    }

    bool find(const state_index s) const
    {
      return std::find(m_vec.begin(), m_vec.end(), s)!=m_vec.end();
    }

    std::vector<state_index>::const_iterator begin() const
    {
      return m_vec.begin();
    }

    std::vector<state_index>::const_iterator end() const
    {
      return m_vec.end();
    }

    void clear()
    {
      m_todo_indicator=0;
      bisimulation_gj::clear(m_vec);
    }

    void clear_todo()
    {
      m_todo_indicator=m_vec.size();
    }
};



// Below the four main data structures are listed.
struct state_type_gj
{
  block_index block=0;
  // std::vector<transition_index>::iterator start_incoming_inert_transitions;
  std::vector<transition>::iterator start_incoming_transitions;
  // std::vector<transition_index>::iterator start_incoming_non_inert_transitions;
  outgoing_transitions_it start_outgoing_transitions;
  std::vector<state_index>::iterator ref_states_in_blocks;
  transition_index no_of_outgoing_inert_transitions=0;
  transition_index counter=undefined; // This field is used to store local information while splitting. While set to -1 (undefined)
                                 // it is considered to be undefined.
                                 // When set to -2 (Rmarked) it is considered to be marked for being in R or R_todo.
  #ifdef CHECK_COMPLEXITY_GJ
    /// \brief print a short state identification for debugging
    template<class LTS_TYPE>
    std::string debug_id_short(const bisim_partitioner_gj<LTS_TYPE>& partitioner) const
    {
        assert(&partitioner.m_states.front() <= this);
        assert(this <= &partitioner.m_states.back());
        return std::to_string(this - &partitioner.m_states.front());
    }

    /// \brief print a state identification for debugging
    template<class LTS_TYPE>
    std::string debug_id(const bisim_partitioner_gj<LTS_TYPE>& partitioner) const
    {
        return "state " + debug_id_short(partitioner);
    }

    mutable check_complexity::state_gj_counter_t work_counter;
  #endif
};

// The struct below facilitates to walk through a LBC_list starting from an arbitrary transition.
typedef std::vector<transition_index>::iterator BLC_list_iterator;
typedef std::vector<transition_index>::const_iterator BLC_list_const_iterator;

// The following type gives the start and end indications of the transitions for the same block, label and constellation
// in the array m_BLC_transitions.
struct BLC_indicators
{
  BLC_list_iterator start_same_BLC;
  BLC_list_iterator start_marked_BLC;
// David suggests: To save a pointer, one could first store the UNmarked transitions and then the marked transitions in a BLC slice.
// i.e. one replaces the two pointers end_marked_BLC and end_same_BLC by one start_marked_BLC.
// If a BLC set does not contain any marked transitions, then start_marked_BLC points to its end;
// if a BLC set does contain some marked transitions, one is allowed to go through them all
// and find the actual end of the BLC set by checking whether transitions actually have the
// correct start block, action label, and end constellation. (To find these, one has to assume
// that a BLC set is not empty, so *start_same_BLC always points at a transition in the BLC
// set, even if all transitions are marked.)
// Practically, not_all_bottom_states_are_marked() should also unmark all transitions
// and splitB() receives the old value of start_marked_BLC as one of its parameters,
// similar to its current parameter splitter_end_unmarked_BLC.
  #ifndef SAVE_BLC_MEMORY
    BLC_list_iterator end_same_BLC;
  #endif

  BLC_indicators(BLC_list_iterator start, BLC_list_iterator end)
   : start_same_BLC(start),
     start_marked_BLC(end)
     #ifndef SAVE_BLC_MEMORY
       , end_same_BLC(end)
     #endif
  {}

  #ifdef CHECK_COMPLEXITY_GJ
    /// \brief print a B_to_C slice identification for debugging
    /// \details This function is only available if compiled in Debug mode.
    template<class LTS_TYPE>
    std::string debug_id(const bisim_partitioner_gj<LTS_TYPE>& partitioner) const
    {
        assert(partitioner.m_BLC_transitions.begin() <= start_same_BLC);
        BLC_list_const_iterator my_end_same_BLC = partitioner.calculate_end_same_BLC(*this);
        assert(start_same_BLC < my_end_same_BLC);
        assert(my_end_same_BLC <= partitioner.m_BLC_transitions.end());
        std::string result("BLC slice ");
        result += partitioner.m_blocks[partitioner.m_states[partitioner.m_aut.get_transitions()[*start_same_BLC].from()].block].debug_id(partitioner);
        result += " -> ";
        result += partitioner.m_constellations[partitioner.m_blocks[partitioner.m_states[partitioner.m_aut.get_transitions()[*start_same_BLC].to()].block].constellation].debug_id(partitioner);
        result += " containing the ";
        if (std::distance<BLC_list_const_iterator>(start_same_BLC, my_end_same_BLC) > 1)
        {
            result += std::to_string(std::distance<BLC_list_const_iterator>(start_same_BLC, my_end_same_BLC));
            result += " transitions ";
        }
        else
            result += "transition ";
        BLC_list_const_iterator iter = start_same_BLC;
        if (start_marked_BLC == iter)
        {
            result += "| ";
        }
        result += partitioner.m_transitions[*iter].debug_id_short(partitioner);
        if (std::distance<BLC_list_const_iterator>(start_same_BLC, my_end_same_BLC) > 4)
        {
            ++iter;
            result += start_marked_BLC == iter ? " | " : ", ";
            result += partitioner.m_transitions[*iter].debug_id_short(partitioner);
            result += std::next(iter) == start_marked_BLC ? " | ..."
                      : (std::next(iter) < start_marked_BLC && start_marked_BLC <= my_end_same_BLC - 3 ? ", ..|.." : ", ...");
            iter = my_end_same_BLC - 3;
        }
        while (++iter != my_end_same_BLC)
        {
            result += start_marked_BLC == iter ? " | " : ", ";
            result += partitioner.m_transitions[*iter].debug_id_short(partitioner);
        }
        if (start_marked_BLC == iter)
        {
            result += " |";
        }
        return result;
    }

    mutable check_complexity::BLC_gj_counter_t work_counter;
  #endif
};

struct transition_type
{
  // The position of the transition type corresponds to m_aut.get_transitions().
  // std::size_t from, label, to are found in m_aut.get_transitions().
  linked_list<BLC_indicators>::iterator transitions_per_block_to_constellation;
  // std::vector<transition_index>::iterator ref_incoming_transitions;
  outgoing_transitions_it ref_outgoing_transitions;  // This refers to the position of this transition in m_outgoing_transitions.
                                                     // During initialisation m_outgoing_transitions contains the indices of this
                                                     // transition. After initialisation m_outgoing_transitions refers to the corresponding
                                                     // entry in m_BLC_transitions, of which the field transition contains the index
                                                     // of this transition.
  // BLC_list_iterator ref_BLC_list;  Access through ref_outgoing_transitions.

  #ifdef CHECK_COMPLEXITY_GJ
    /// \brief print a short transition identification for debugging
    /// \details This function is only available if compiled in Debug mode.
    template<class LTS_TYPE>
    std::string debug_id_short(const bisim_partitioner_gj<LTS_TYPE>& partitioner) const
    {
        assert(&partitioner.m_transitions.front() <= this);
        assert(this <= &partitioner.m_transitions.back());
        const transition& t = partitioner.m_aut.get_transitions()[this - &partitioner.m_transitions.front()];
        return partitioner.m_states[t.from()].debug_id_short(partitioner) + " -" +
               pp(partitioner.m_aut.action_label(t.label())) + "-> " +
               partitioner.m_states[t.to()].debug_id_short(partitioner);
    }

    /// \brief print a transition identification for debugging
    /// \details This function is only available if compiled in Debug mode.
    template<class LTS_TYPE>
    std::string debug_id(const bisim_partitioner_gj<LTS_TYPE>& partitioner) const
    {
        return "transition " + debug_id_short(partitioner);
    }

    mutable check_complexity::trans_gj_counter_t work_counter;
#endif
};

struct block_type
{
  constellation_index constellation : (sizeof(constellation_index) * CHAR_BIT - 1);
  unsigned contains_new_bottom_states : 1;
  std::vector<state_index>::iterator start_bottom_states;
  std::vector<state_index>::iterator start_non_bottom_states;
  std::vector<state_index>::iterator end_states;
// David thinks that end_states may perhaps be suppressed.
// We need the size of a block in two cases: to choose a small block in a constellation,
// and to decide whether to abort a coroutine early in simple_splitB().
// In both cases it is enough to get an upper bound on the size of the block,
// and if a constellation is a contiguous slice in m_states_in_blocks (as I
// suggested elsewhere), the constellation can provide this upper bound.
  // The list below seems too expensive. Maybe a cheaper construction is possible. Certainly the size of the list is not important.
  linked_list< BLC_indicators > block_to_constellation;

  block_type(const std::vector<state_index>::iterator beginning_of_states, constellation_index c)
    : constellation(c),
      contains_new_bottom_states(false),
      start_bottom_states(beginning_of_states),
      start_non_bottom_states(beginning_of_states),
      end_states(beginning_of_states)
  {}

  #ifdef CHECK_COMPLEXITY_GJ
    /// \brief print a block identification for debugging
    template<class LTS_TYPE>
    inline std::string debug_id(const bisim_partitioner_gj<LTS_TYPE>& partitioner) const
    {
        assert(!partitioner.m_blocks.empty());
        assert(&partitioner.m_blocks.front() <= this);
        assert(this <= &partitioner.m_blocks.back());
        assert(partitioner.m_states_in_blocks.begin() <= start_bottom_states);
        assert(start_bottom_states <= start_non_bottom_states);
        assert(start_non_bottom_states <= end_states);
        assert(end_states <= partitioner.m_states_in_blocks.end());
        return "block [" + std::to_string(std::distance<std::vector<state_index>::const_iterator>(partitioner.m_states_in_blocks.begin(), start_bottom_states)) + "," + std::to_string(std::distance<std::vector<state_index>::const_iterator>(partitioner.m_states_in_blocks.begin(),end_states)) + ")"
                    " (#" + std::to_string(std::distance(&partitioner.m_blocks.front(), this)) + ")";
    }

    mutable check_complexity::block_gj_counter_t work_counter;
  #endif
};

struct constellation_type
{
// David suggests: group the states not only per block but also per constellation in m_states_in_blocks.
// Then, a constellation can be a contiguous area in m_states_in_blocks,
// and we do not need to store a list of block indices to describe a constellation
// (but only the index of the first and last state in m_states_in_blocks that belongs to the constellation).
// With that grouping, it is also not necessary to store constellations in a vector
// or a list; one can allocate each constellation as its own data structure
// and use a pointer to constellation_type as identifier of a constellation.
// JFG answers: I do not see this. Moving a block out of the constellation is tricky.
// David answers: One can still move out the first or the last block, whichever is smaller.
  std::vector<state_index>::iterator start_const_states;
  std::vector<state_index>::iterator end_const_states;
  constellation_type(const std::vector<state_index>::iterator new_start, const std::vector<state_index>::iterator new_end)
    : start_const_states(new_start),
      end_const_states(new_end)
  {}

  #ifndef NDEBUG
    /// \brief print a constellation identification for debugging
    template<class LTS_TYPE>
    inline std::string debug_id(const bisim_partitioner_gj<LTS_TYPE>& partitioner) const
    {
        assert(&partitioner.m_constellations.front() <= this);
        assert(this <= &partitioner.m_constellations.back());
        return "constellation " + std::to_string(this - &partitioner.m_constellations.front());
    }
  #endif
};

} // end namespace bisimulation_gj


/*=============================================================================
=                                 main class                                  =
=============================================================================*/


using namespace mcrl2::lts::detail::bisimulation_gj;

/// \class bisim_partitioner_gj
/// \brief implements the main algorithm for the branching bisimulation quotient
template <class LTS_TYPE>
class bisim_partitioner_gj
{
  protected:

    typedef std::unordered_set<state_index> set_of_states_type;
    typedef std::unordered_set<transition_index> set_of_transitions_type;
    typedef std::vector<constellation_index> set_of_constellations_type;

    typedef std::unordered_map<std::pair<label_index, constellation_index>, set_of_states_type>
                      label_constellation_to_set_of_states_map;
    typedef std::unordered_map<std::pair<block_index, label_index>, transition_index> block_label_to_size_t_map;

    #ifndef NDEBUG
      public: // needed for the debugging functions, e.g. debug_id().
    #endif
    /// \brief automaton that is being reduced
    LTS_TYPE& m_aut;

    // Generic data structures.
    std::vector<state_type_gj> m_states;
    // std::vector<transition_index> m_incoming_transitions;
    std::vector<transition_pointer_pair> m_outgoing_transitions;
                                                                  // During refining this contains the index in m_BLC_transition, of which
                                                                  // the transition field contains the index of the transition.
    std::vector<transition_type> m_transitions;
    std::vector<state_index> m_states_in_blocks;
    // David suggests to change the type of m_states_in_blocks to std::vector<std::vector<state_type>::iterator>.
    // This will allow to store pointers-to-states in m_states_in_blocks, and will reduce the complexity of address calculations.
    std::vector<block_type> m_blocks;
    std::vector<constellation_type> m_constellations;
    // David suggests to allocate blocks and constellations in global_linked_list_administration.
    // (There is a pool allocator that can be used like this in liblts_bisim_dnj.h.)
    // Then, one can store pointers to block_type and constellation_type instead of numbers;
    // the type checking of pointers is more strict than the type checking of integer types,
    // so it becomes impossible to assign a constellation number to a block or v.v.
    // Also, it will reduce the complexity of address calculations.
    std::vector<transition_index> m_BLC_transitions;
  protected:
    std::vector<block_index> m_blocks_with_new_bottom_states;
    // Below are the two vectors that contain the marked and unmarked states, which
    // are internally split in a part for states to be investigated, and a part for
    // states that belong definitively to this set.
    todo_state_vector m_R, m_U;
    std::vector<state_index> m_U_counter_reset_vector;
    // The following variable contains all non trivial constellations.
    set_of_constellations_type m_non_trivial_constellations;


    /// \brief true iff branching (not strong) bisimulation has been requested
    const bool m_branching;

    /// \brief true iff divergence-preserving branching bisimulation has been
    /// requested
    /// \details Note that this field must be false if strong bisimulation has
    /// been requested.  There is no such thing as divergence-preserving strong
    /// bisimulation.
    const bool m_preserve_divergence;

    // The function tests whether transition t is inert during initialisation, i.e. when there is only one source/target block.
    bool is_inert_during_init(const transition& t) const
    {
      return m_branching && m_aut.is_tau(m_aut.apply_hidden_label_map(t.label())) && (!m_preserve_divergence || t.from() != t.to());
    }

    // The function assumes that m_branching is true and tests whether transition t is inert during initialisation under that condition
    bool is_inert_during_init_if_branching(const transition& t) const
    {
      assert(m_branching);
      return m_aut.is_tau(m_aut.apply_hidden_label_map(t.label())) && (!m_preserve_divergence || t.from() != t.to());
    }

    // The function calculates the label index of transition t, where tau-self-loops get the special index m_aut.num_action_labels() if divergence needs to be preserved
    label_index label_or_divergence(const transition& t) const
    {
      label_index result = m_aut.apply_hidden_label_map(t.label());
      if (m_preserve_divergence && (assert(m_branching), t.from() == t.to()) && m_aut.is_tau(result))
      {
        return m_aut.num_action_labels();
      }
      return result;
    }

    // The function decides whether iterator it points to a transition that is part
    // of the BLC set ind. (This is non-trivial if SAVE_BLC_MEMORY.)
    bool points_into_BLC_set(BLC_list_const_iterator it, const BLC_indicators& ind) const
    {
      assert(ind.start_same_BLC <= ind.start_marked_BLC);
      assert(ind.start_same_BLC <= it);
      #ifdef SAVE_BLC_MEMORY
        if (m_BLC_transitions.end() <= it)
        {
          return false;
        }
        const transition& first_t = m_aut.get_transitions()[*ind.start_same_BLC];
        const transition& it_t = m_aut.get_transitions()[*it];
        return m_states[first_t.from()].block == m_states[it_t.from()].block &&
               label_or_divergence(first_t) == label_or_divergence(it_t) &&
               m_blocks[m_states[first_t.to()].block].constellation == m_blocks[m_states[it_t.to()].block].constellation;
      #else
        assert(ind.start_marked_BLC <= ind.end_same_BLC);
        return it < ind.end_same_BLC;
      #endif
    }

    // This function returns true iff the BLC set ind contains at least one
    // marked transition.
    bool has_marked_transitions(const BLC_indicators& ind) const
    {
      return points_into_BLC_set(ind.start_marked_BLC, ind);
    }

    // This function unmarks all transitions in BLC set ind.
    // If SAVE_BLC_MEMORY, it is linear in the number of marked transitions, but
    // that is ok, as transition markings come with the right to visit the
    // transitions.
    void unmark_all_transitions(BLC_indicators& ind)
    {
      assert(ind.start_same_BLC <= ind.start_marked_BLC);
      #ifdef SAVE_BLC_MEMORY
        BLC_list_iterator end_same_BLC = ind.start_marked_BLC;
        for (; points_into_BLC_set(end_same_BLC, ind); ++it)
        {
          // mCRL2complexity(m_aut.get_transitions()[*ind], add_work(...), *this);
              // not needed, as we unmark transitions and we are allowed to visit
              // every marked transition a fixed number of times
        }
        ind.start_marked_BLC = end_same_BLC;
      #else
        ind.start_marked_BLC = ind.end_same_BLC;
      #endif
    }

  public:
    BLC_list_const_iterator calculate_end_same_BLC(const BLC_indicators& ind) const
    {
      assert(ind.start_same_BLC <= ind.start_marked_BLC);
      #ifdef SAVE_BLC_MEMORY
        BLC_list_const_iterator it = ind.start_marked_BLC;
        for (; points_into_BLC_set(it, ind); ++it)
        {}
        return it;
      #else
        return ind.end_same_BLC;
      #endif
    }

#ifndef NDEBUG 
    // This suppresses many unused variable warnings.
    // This function calculates the end of the transitions in BLC set ind.
    // Because it is slow if SAVE_BLC_MEMORY, it is only defined in debug mode
    // and should only be called in debugging functions.

    // This function counts the number of marked transitions. Because it is slow
    // if SAVE_BLC_MEMORY, it is only defined in debug mode and should only be
    // called in debugging functions.
    transition_index count_marked_transitions(const BLC_indicators& ind) const
    {
      return std::distance<BLC_list_const_iterator>(ind.start_marked_BLC, calculate_end_same_BLC(ind));
    }

  protected:
    void check_transitions(const bool check_temporary_complexity_counters, const bool check_block_to_constellation = true) const
    {
      // This routine can only be used after initialisation.
      for(std::size_t ti=0; ti<m_transitions.size(); ++ti)
      {
        const BLC_list_const_iterator btc_ti=m_BLC_transitions.begin() + m_transitions[ti].ref_outgoing_transitions->transition;
        assert(*btc_ti==ti);

        const transition& t=m_aut.get_transitions()[ti];
        assert(&*m_states[t.to()].start_incoming_transitions <= &t);
        assert(&t <= &*std::prev(m_aut.get_transitions().end()));
        if (t.to() + 1 != m_states.size())
        {
          assert(&t <= &*std::prev(m_states[t.to() + 1].start_incoming_transitions));
        }

        assert(m_states[t.from()].start_outgoing_transitions <= m_transitions[ti].ref_outgoing_transitions);
        if (t.from() + 1 == m_states.size())
          assert(m_transitions[ti].ref_outgoing_transitions < m_outgoing_transitions.end());
        else
          assert(m_transitions[ti].ref_outgoing_transitions < m_states[t.from() + 1].start_outgoing_transitions);

        assert(m_transitions[ti].transitions_per_block_to_constellation->start_same_BLC <= btc_ti);
        assert(points_into_BLC_set(btc_ti, *m_transitions[ti].transitions_per_block_to_constellation));

        if (!check_block_to_constellation)
          return;

        const block_index b=m_states[t.from()].block;

        const label_index t_label = label_or_divergence(t);
        bool found=false;
        for(const BLC_indicators& blc: m_blocks[b].block_to_constellation)
        {
          if (blc.start_same_BLC != blc.start_marked_BLC)
          {
            assert(blc.start_same_BLC < blc.start_marked_BLC);
            assert(points_into_BLC_set(std::prev(blc.start_marked_BLC), blc));
          }
          #ifndef SAVE_BLC_MEMORY
            assert(blc.start_same_BLC < blc.end_same_BLC);
          #endif
          transition& first_t = m_aut.get_transitions()[*blc.start_same_BLC];
          assert(b == m_states[first_t.from()].block);
          if (t_label == label_or_divergence(first_t) &&
              m_blocks[m_states[first_t.to()].block].constellation == m_blocks[m_states[t.to()].block].constellation)
          {
// if (found) { std::cerr << "Found multiple BLC sets with transitions (block " << b << " -" << m_aut.action_label(t.label()) << "-> constellation " << m_blocks[m_states[t.to()].block].constellation << ")\n"; }
            assert(!found);
            assert(blc.start_same_BLC <= btc_ti);
            assert(points_into_BLC_set(btc_ti, blc));
            assert(&blc == &*m_transitions[ti].transitions_per_block_to_constellation);
            found = true;
          }
        }
        assert(found);
        if (check_temporary_complexity_counters)
        {
          #ifdef CHECK_COMPLEXITY_GJ
            const block_index targetb = m_states[t.to()].block;
            const unsigned max_sourceB = check_complexity::log_n - check_complexity::ilog2(number_of_states_in_block(b));
            const unsigned max_targetC = check_complexity::log_n - check_complexity::ilog2(number_of_states_in_constellation(m_blocks[targetb].constellation));
            const unsigned max_targetB = check_complexity::log_n - check_complexity::ilog2(number_of_states_in_block(targetb));
            mCRL2complexity_gj(&m_transitions[ti], no_temporary_work(max_sourceB, max_targetC, max_targetB,
                  0 == m_states[t.from()].no_of_outgoing_inert_transitions), *this);
          #endif
        }
      }
    }

    bool check_data_structures(const std::string& tag, const bool initialisation=false, const bool check_temporary_complexity_counters=true) const
    {
assert(!initialisation);
      mCRL2log(log::debug) << "Check data structures: " << tag << ".\n";
      assert(m_states.size()==m_aut.num_states());
      assert(m_outgoing_transitions.size()==m_aut.num_transitions());

      // Check that the elements in m_states are well formed.
      for(state_index si=0; si< m_states.size(); si++)
      {
        const state_type_gj& s=m_states[si];

        assert(s.counter==undefined);
        assert(m_blocks[s.block].start_bottom_states<m_blocks[s.block].start_non_bottom_states);
        assert(m_blocks[s.block].start_non_bottom_states<=m_blocks[s.block].end_states);

        assert(std::find(m_blocks[s.block].start_bottom_states, m_blocks[s.block].end_states,si)!=m_blocks[s.block].end_states);

        // The construction below is added to enable compilation on Windows. 
        outgoing_transitions_const_it end_it1;
        if (si+1>=m_states.size()) { end_it1=m_outgoing_transitions.cend(); } else { end_it1=m_states[si+1].start_outgoing_transitions; }
        for(outgoing_transitions_const_it it=s.start_outgoing_transitions; it!=end_it1; ++it)
        {
          const transition& t=m_aut.get_transitions()
                                [initialisation ?it->transition :m_BLC_transitions[it->transition]];
// if (t.from() != si) { std::cerr << m_transitions[m_BLC_transitions[it->transition]].debug_id(*this) << " is an outgoing transition of state " << si << "!\n"; }
          assert(t.from()==si);
          assert(!initialisation || m_transitions[it->transition].ref_outgoing_transitions==it);
          assert(initialisation || m_transitions[m_BLC_transitions[it->transition]].ref_outgoing_transitions==it);
          assert((it->start_same_saC>it && it->start_same_saC<m_outgoing_transitions.end() &&
                        ((it+1)->start_same_saC==it->start_same_saC || (it+1)->start_same_saC<=it)) ||
                 (it->start_same_saC<=it && (it+1==m_outgoing_transitions.end() || (it+1)->start_same_saC>it)));
// if (it->start_same_saC < it->start_same_saC->start_same_saC) { std::cerr << "Now checking transitions " << m_transitions[m_BLC_transitions[it->start_same_saC->transition]].debug_id_short(*this) << " ... " << m_transitions[m_BLC_transitions[it->start_same_saC->start_same_saC->transition]].debug_id_short(*this) << '\n'; }
          const label_index t_label = label_or_divergence(t);
          // The following for loop is only executed if it is the last transition in the saC-slice.
          for(outgoing_transitions_const_it itt=it->start_same_saC; itt<it->start_same_saC->start_same_saC; ++itt)
          {
            const transition& t1=m_aut.get_transitions()
                                 [initialisation?itt->transition:m_BLC_transitions[itt->transition]];
// if (t1.from()!=si) { assert(!initialisation); std::cerr << m_transitions[m_BLC_transitions[itt->transition]].debug_id(*this) << " does not start in state " << si << '\n'; }
            assert(t1.from()==si);
            assert(label_or_divergence(t1) == t_label);
            assert(m_blocks[m_states[t.to()].block].constellation==m_blocks[m_states[t1.to()].block].constellation);
          }
        }
        assert(*(s.ref_states_in_blocks)==si);

        // Check that for each state the outgoing transitions satisfy the following invariant.
        // First there are inert transitions. Then there are other transitions sorted per label
        // and constellation.
        std::unordered_set<std::pair<label_index, constellation_index> > constellations_seen;

        // The construction below is to enable translation on Windows. 
        outgoing_transitions_const_it end_it2;
        if (si+1>=m_states.size()) { end_it2=m_outgoing_transitions.cend(); } else { end_it2=m_states[si+1].start_outgoing_transitions; }
        for(outgoing_transitions_const_it it=s.start_outgoing_transitions; it!=end_it2; ++it)
        {
          const transition& t=m_aut.get_transitions()[/*initialisation?it->transition:*/m_BLC_transitions[it->transition]];
          const label_index label = label_or_divergence(t);
          // Check that if the target constellation, if not new, is equal to the target constellation of the previous outgoing transition.
          const constellation_index t_to_constellation = m_blocks[m_states[t.to()].block].constellation;
          if (constellations_seen.count(std::pair(label,t_to_constellation))>0)
          {
            assert(it!=s.start_outgoing_transitions);
            const transition& old_t=m_aut.get_transitions()[/*initialisation?std::prev(it)->transition:*/m_BLC_transitions[std::prev(it)->transition]];
            assert(label_or_divergence(old_t) == label);
            assert(t_to_constellation==m_blocks[m_states[old_t.to()].block].constellation);
          }
          constellations_seen.emplace(label,t_to_constellation);
        }
      }
      // Check that the elements in m_transitions are well formed.
      if (!initialisation)
      {
        check_transitions(check_temporary_complexity_counters);
      }
      // Check that the elements in m_blocks are well formed.
      {
        set_of_transitions_type all_transitions;
        for(block_index bi=0; bi<m_blocks.size(); ++bi)
        {
          const block_type& b=m_blocks[bi];
          const constellation_type& c=m_constellations[b.constellation];
          assert(b.start_bottom_states<m_states_in_blocks.end());
          assert(b.start_non_bottom_states<=m_states_in_blocks.end());
          assert(b.start_non_bottom_states>=m_states_in_blocks.begin());

          assert(m_states_in_blocks.begin() <= c.start_const_states);
          assert(c.start_const_states <= b.start_bottom_states);
          assert(b.start_bottom_states < b.start_non_bottom_states);
          assert(b.start_non_bottom_states <= b.end_states);
          assert(b.end_states <= c.end_const_states);
          assert(c.end_const_states <= m_states_in_blocks.end());

          #ifdef CHECK_COMPLEXITY_GJ
            unsigned max_B = check_complexity::log_n - check_complexity::ilog2(number_of_states_in_block(bi));
            unsigned max_C = check_complexity::log_n - check_complexity::ilog2(number_of_states_in_constellation(b.constellation));
          #endif
          for(typename std::vector<state_index>::const_iterator is=b.start_bottom_states;
                   is!=b.start_non_bottom_states; ++is)
          {
            const state_type_gj& s=m_states[*is];
            assert(s.block==bi);
            assert(s.no_of_outgoing_inert_transitions==0);
            if (check_temporary_complexity_counters)
            {
              mCRL2complexity_gj(&s, no_temporary_work(max_B, true), *this);
            }
          }
          for(typename std::vector<state_index>::const_iterator is=b.start_non_bottom_states;
                   is!=b.end_states; ++is)
          {
            const state_type_gj& s=m_states[*is];
            assert(s.block==bi);
            assert(s.no_of_outgoing_inert_transitions>0);
            // Because there cannot be new bottom states among the non-bottom states,
            // we can always check the temporary work of non-bottom states:
            mCRL2complexity_gj(&s, no_temporary_work(max_B, false), *this);
          }
          // Because a block has no temporary or new-bottom-state-related counters,
          // we can always check its temporary work:
          mCRL2complexity_gj(&b, no_temporary_work(max_C, max_B), *this);

          assert(b.block_to_constellation.check_linked_list());
          for(linked_list< BLC_indicators >::iterator ind=b.block_to_constellation.begin();
                     ind!=b.block_to_constellation.end(); ++ind)
          {
            const transition& first_transition=m_aut.get_transitions()[*(ind->start_same_BLC)];
            const label_index first_transition_label = label_or_divergence(first_transition);
            for(BLC_list_const_iterator i=ind->start_same_BLC; points_into_BLC_set(i, *ind); ++i)
            {
              const transition& t=m_aut.get_transitions()[*i];
              assert(m_transitions[*i].transitions_per_block_to_constellation == ind);
              all_transitions.emplace(*i);
              assert(m_states[t.from()].block==bi);
              assert(m_blocks[m_states[t.to()].block].constellation==
                               m_blocks[m_states[first_transition.to()].block].constellation);
              assert(label_or_divergence(t) == first_transition_label);
              if (is_inert_during_init(t) && m_blocks[m_states[t.to()].block].constellation==m_blocks[bi].constellation)
              {
                // The inert transitions should be in the first element of block_to_constellation:
                assert(b.block_to_constellation.begin()==ind);
              }
            }
            if (check_temporary_complexity_counters)
            {
              mCRL2complexity_gj(ind, no_temporary_work(max_C, check_complexity::log_n - check_complexity::ilog2(number_of_states_in_constellation(m_blocks[m_states[first_transition.to()].block].constellation))), *this);
            }
          }
        }
        assert(initialisation || all_transitions.size()==m_transitions.size());
        // destruct all_transitions here
      }

      // TODO Check that the elements in m_constellations are well formed.
      {
        std::unordered_set<block_index> all_blocks;
        for(constellation_index ci=0; ci<m_constellations.size(); ci++)
        {
          for (std::vector<state_index>::const_iterator constln_it = m_constellations[ci].start_const_states; constln_it < m_constellations[ci].end_const_states; )
          {
            const block_index bi = m_states[*constln_it].block;
            assert(bi<m_blocks.size());
            assert(all_blocks.emplace(bi).second);  // Block is not already present. Otherwise a block occurs in two constellations.
            constln_it = m_blocks[bi].end_states;
          }
        }
        assert(all_blocks.size()==m_blocks.size());
        // destruct all_blocks here
      }

      // Check that the states in m_states_in_blocks refer to with ref_states_in_block to the right position.
      // and that a state is correctly designated as a (non-)bottom state.
      for(typename std::vector<state_index>::const_iterator si=m_states_in_blocks.begin(); si!=m_states_in_blocks.end(); ++si)
      {
        assert(si==m_states[*si].ref_states_in_blocks);
      }

      // Check that the blocks in m_blocks_with_new_bottom_states are bottom states.
      for(const block_index bi: m_blocks_with_new_bottom_states)
      {
        assert(m_blocks[bi].contains_new_bottom_states);
      }

      // Check that the non-trivial constellations are non trivial.
      for(const constellation_index ci: m_non_trivial_constellations)
      {
        // There are at least two blocks in a non-trivial constellation.
        const block_index first_bi = m_states[*m_constellations[ci].start_const_states].block;
        const block_index last_bi = m_states[*std::prev(m_constellations[ci].end_const_states)].block;
        assert(first_bi != last_bi);
      }
      return true;
    }
#endif //#ifndef NDEBUG

    /// Checks the following invariant:
    ///     If a block has a constellation-non-inert transition, then every
    ///     bottom state has a constellation-non-inert transition with the same
    ///     label to the same target constellation.
    /// It is assumed that the BLC data structure is correct, so we conveniently
    /// use that to verify the invariant.
    ///
    /// The function can also check a partial invariant while stabilisation has
    /// not yet finished. If calM != nullptr, then we have:
    ///     The above invariant may be violated for BLC sets that are still to
    ///     be stabilized, as given by the main splitters in calM.
    ///     (calM_elt indicates how far stabilization has handled calM already.)
    ///     (block_label_to_cotransition indicates the co-splitters that belong
    ///     to the main splitters in calM.)
    ///     It may also be violated for blocks that contain new bottom states,
    ///     as indicated by m_blocks_with_new_bottom_states.
    ///
    /// Additionally, the function ensures that only transitions in BLC sets
    /// satisfying the above conditions are marked:
    ///     Transitions may only be marked in BLC sets that are still to be
    ///     stabilized, as given by calM (including co-splitters); they may also
    ///     be marked if they start in new bottom states, as indicated by
    ///     m_blocks_with_new_bottom_states.
    bool check_stability(const std::string& tag,
                         const std::vector<std::pair<BLC_list_iterator, BLC_list_iterator> >* calM = nullptr,
                         const std::pair<BLC_list_iterator, BLC_list_iterator>* calM_elt = nullptr,
                         const block_label_to_size_t_map* const block_label_to_cotransition = nullptr) const
    {
#ifndef NDEBUG
      mCRL2log(log::debug) << "Check stability: " << tag << ".\n";
      for(block_index bi=0; bi<m_blocks.size(); ++bi)
      {
        const block_type& b=m_blocks[bi];
        for(linked_list< BLC_indicators >::iterator ind=b.block_to_constellation.begin();
                     ind!=b.block_to_constellation.end(); ++ind)
        {
          set_of_states_type all_source_bottom_states;

          #ifndef SAVE_BLC_MEMORY
            assert(ind->start_same_BLC < ind->end_same_BLC);
          #endif
          const transition& first_t = m_aut.get_transitions()[*ind->start_same_BLC];
          const label_index first_t_label = label_or_divergence(first_t);
          const bool all_transitions_in_BLC_are_inert = is_inert_during_init(first_t) &&
                                                        m_blocks[m_states[first_t.to()].block].constellation == b.constellation;
          assert(!all_transitions_in_BLC_are_inert || b.block_to_constellation.begin() == ind);
          for(BLC_list_const_iterator i=ind->start_same_BLC; points_into_BLC_set(i, *ind); ++i)
          {
            assert(m_BLC_transitions.begin() <= i);
            assert(i < m_BLC_transitions.end());
            const transition& t=m_aut.get_transitions()[*i];
// if (m_states[t.from()].block != bi) { std::cerr << m_transitions[*ind->start_same_BLC].debug_id(*this) << " should start in block " << bi << '\n'; }
            assert(m_states[t.from()].block == bi);
            assert(label_or_divergence(t) == first_t_label);
            assert(m_blocks[m_states[t.to()].block].constellation == m_blocks[m_states[first_t.to()].block].constellation);
            if (is_inert_during_init(t) && m_blocks[m_states[t.to()].block].constellation==b.constellation)
            {
              assert(all_transitions_in_BLC_are_inert);
            }
            else
            {
              // This is a constellation-non-inert transition.
              assert(!all_transitions_in_BLC_are_inert);
              if (0 == m_states[t.from()].no_of_outgoing_inert_transitions)
              {
                assert(b.start_bottom_states <= m_states[t.from()].ref_states_in_blocks);
                assert(m_states[t.from()].ref_states_in_blocks < b.start_non_bottom_states);
                all_source_bottom_states.emplace(t.from());
              }
              else
              {
                assert(b.start_non_bottom_states <= m_states[t.from()].ref_states_in_blocks);
                assert(m_states[t.from()].ref_states_in_blocks < b.end_states);
              }
            }
          }
          assert(all_source_bottom_states.size() <= static_cast<std::size_t>(std::distance(b.start_bottom_states, b.start_non_bottom_states)));
          // check that every bottom state has a transition in this BLC entry:
          bool eventual_instability_is_ok = true;
          bool eventual_marking_is_ok = true;
          if (!all_transitions_in_BLC_are_inert &&
              all_source_bottom_states.size()!=static_cast<std::size_t>(std::distance(b.start_bottom_states, b.start_non_bottom_states)))
          {
            // only splitters should be instable.
            mCRL2log(log::debug) << "Not all " << std::distance(b.start_bottom_states, b.start_non_bottom_states)<< " bottom states have a transition in the " << ind->debug_id(*this) << ": transitions found from states";
            for (set_of_states_type::iterator asbc_it = all_source_bottom_states.begin() ; asbc_it != all_source_bottom_states.end() ; ++asbc_it) { mCRL2log(log::debug) << ' ' << *asbc_it; }
            mCRL2log(log::debug) << '\n';
            eventual_instability_is_ok = false;
          }
          if (has_marked_transitions(*ind))
          {
            // only splitters should contain marked transitions.
            mCRL2log(log::debug) << ind->debug_id(*this) << " contains " << count_marked_transitions(*ind) << " marked transitions.\n";
            eventual_marking_is_ok = false;
          }
          if (m_blocks[bi].contains_new_bottom_states)
          {
            /* I would like the following to check more closely because in a
               block with new bottom states, one should have...
            if (!eventual_marking_is_ok)
            {
              eventual_marking_is_ok = true;
              for (BLC_list_const_iterator i=ind->start_marked_BLC; points_into_BLC_set(i, *ind); ++i)
              {
                const state_index from = m_aut.get_transitions()[*i].from();
                // assert(m_states[from].block == bi); -- already checked earlier
                if (0 != m_states[from].no_of_outgoing_inert_transitions)
                {
                  // the state is a non-bottom state
                  eventual_marking_is_ok = false;
                  break;
                }
              }
              if (eventual_marking_is_ok)
              {
                mCRL2log(log::debug) << "  This is ok because all marked transitions begin in new bottom states of block " << bi << ".\n";
                eventual_instability_is_ok = true;
              }
            } */
            if (!(eventual_instability_is_ok && eventual_marking_is_ok))
            {
              mCRL2log(log::debug) << "  This is ok because block " << bi << " contains new bottom states.\n";
              eventual_instability_is_ok = true;
              eventual_marking_is_ok = true;
            }
          }
          if (!(eventual_instability_is_ok && eventual_marking_is_ok) && nullptr != calM && calM->begin() != calM->end())
          {
            std::vector<std::pair<BLC_list_iterator, BLC_list_iterator> >::const_iterator calM_iter = calM->begin();
            if (nullptr != calM_elt)
            {
              for(;;)
              {
                assert(calM->end() != calM_iter);
                if (calM_iter->first <= calM_elt->first && calM_elt->second <= calM_iter->second)
                {
                  break;
                }
                ++calM_iter;
              }
              if (calM_elt->first <= ind->start_same_BLC && ind->start_same_BLC <= calM_iter->second && !points_into_BLC_set(calM_elt->second, *ind))
              {
                mCRL2log(log::debug) << "  This is ok because the BLC set (block " << bi << " -" << pp(m_aut.action_label(first_t.label())) << "-> constellation " << m_blocks[m_states[first_t.to()].block].constellation << ") is soon going to be a main splitter.\n";
                eventual_instability_is_ok = true;
                eventual_marking_is_ok = true;
              }
              else if (nullptr != block_label_to_cotransition)
              {
                for (BLC_list_const_iterator ind_iter = calM_elt->second; ind_iter > calM_elt->first; )
                {
                  ind_iter = m_transitions[*std::prev(ind_iter)].transitions_per_block_to_constellation->start_same_BLC;
                  const transition& t = m_aut.get_transitions()[*ind_iter];
                  block_label_to_size_t_map::const_iterator co_iter = block_label_to_cotransition->find(std::pair(m_states[t.from()].block, label_or_divergence(t)));
                  if (block_label_to_cotransition->end() != co_iter && null_transition != co_iter->second && m_transitions[co_iter->second].transitions_per_block_to_constellation == ind)
                  {
                    const transition& co_t = m_aut.get_transitions()[co_iter->second];
                    mCRL2log(log::debug) << "  This is ok because the BLC set (block " << m_states[co_t.from()].block << " -" << pp(m_aut.action_label(co_t.label())) << "-> constellation " << m_blocks[m_states[co_t.to()].block].constellation << ") is soon going to be a co-splitter.\n";
                    eventual_instability_is_ok = true;
                    eventual_marking_is_ok = true;
                    break;
                  }
                }
              }
              ++calM_iter;
            }
            for(; !(eventual_instability_is_ok && eventual_marking_is_ok) && calM->end() != calM_iter; ++calM_iter)
            {
              if (calM_iter->first <= ind->start_same_BLC && ind->start_same_BLC <= calM_iter->second && !points_into_BLC_set(calM_iter->second, *ind))
              {
                mCRL2log(log::debug) << "  This is ok because the BLC set (block " << bi << " -" << pp(m_aut.action_label(first_t.label())) << "-> constellation " << m_blocks[m_states[first_t.to()].block].constellation << ") is going to be a main splitter later.\n";
                eventual_instability_is_ok = true;
                eventual_marking_is_ok = true;
              }
              else if (nullptr != block_label_to_cotransition)
              {
                for (BLC_list_const_iterator ind_iter = calM_iter->second; ind_iter > calM_iter->first; )
                {
                  ind_iter = m_transitions[*std::prev(ind_iter)].transitions_per_block_to_constellation->start_same_BLC;
                  const transition& t = m_aut.get_transitions()[*ind_iter];
                  block_label_to_size_t_map::const_iterator co_iter = block_label_to_cotransition->find(std::pair(m_states[t.from()].block, label_or_divergence(t)));
                  if (block_label_to_cotransition->end() != co_iter && null_transition != co_iter->second && m_transitions[co_iter->second].transitions_per_block_to_constellation == ind)
                  {
                    const transition& co_t = m_aut.get_transitions()[co_iter->second];
                    mCRL2log(log::debug) << "  This is ok because the BLC set (block " << m_states[co_t.from()].block << " -" << pp(m_aut.action_label(co_t.label())) << "-> constellation " << m_blocks[m_states[co_t.to()].block].constellation << ") is going to be a co-splitter later.\n";
                    eventual_instability_is_ok = true;
                    eventual_marking_is_ok = true;
                    break;
                  }
                }
              }
            }
          }
          assert(eventual_marking_is_ok);
          assert(eventual_instability_is_ok);
        }
      }
      mCRL2log(log::debug) << "Check stability finished: " << tag << ".\n";
#endif //#ifndef NDEBUG
      return true;
    }

  #ifdef CHECK_COMPLEXITY_GJ
    // assign work to all transitions from state si with label a to constellation C.
    // It is ensured that there will be some work assigned.
    void add_work_to_same_saC(const bool initialisation, const state_index si, const label_index a, const constellation_index C, const enum check_complexity::counter_type ctr, const unsigned max_value) const
    {
      assert(DONT_COUNT_TEMPORARY != max_value);
      if (m_aut.num_action_labels() == a)
      {
        assert(m_preserve_divergence);
        assert(C == m_blocks[m_states[si].block].constellation);
      }
      else
      {
        assert(a < m_aut.num_action_labels());
        assert(m_aut.apply_hidden_label_map(a) == a);
      }
      bool work_assigned = false;
      for (outgoing_transitions_const_it outtrans = (si + 1 >= m_states.size() ? m_outgoing_transitions.end() : m_states[si +  1].start_outgoing_transitions);
                  outtrans != m_states[si].start_outgoing_transitions; )
      {
        --outtrans;
        const transition& t = m_aut.get_transitions()[initialisation?outtrans->transition:m_BLC_transitions[outtrans->transition]];
        if (label_or_divergence(t) == a && m_blocks[m_states[t.to()].block].constellation == C)
        {
          mCRL2complexity_gj(&m_transitions[initialisation?outtrans->transition:m_BLC_transitions[outtrans->transition]], add_work(ctr, max_value), *this);
          work_assigned = true;
        }
      }
      if (!work_assigned)
      {
        mCRL2log(log::error) << "No suitable transition " << m_states[si].debug_id_short(*this)
            << " -" << pp(m_aut.action_label(m_aut.num_action_labels() == a ? m_aut.tau_label_index() : a)) << "-> "
            << (m_aut.num_action_labels() == a ? m_states[si].debug_id_short(*this) : m_constellations[C].debug_id(*this))
            << " found to assign work for counter \""
                << check_complexity::work_names[ctr - check_complexity::BLOCK_MIN] << "\"\n";
        exit(EXIT_FAILURE);
      }
    }
  #else
    #define add_work_to_same_saC(initialisation, si, a, C, ctr, max_value) do{}while (0)
  #endif //#ifdef CHECK_COMPLEXITY_GJ

    void display_BLC_list(const block_index bi) const
    {
      mCRL2log(log::debug) << "\n  BLC_List\n";
      for(const BLC_indicators& blc_it: m_blocks[bi].block_to_constellation)
      {
        mCRL2log(log::debug) << "\n    BLC_sublist:  " << std::distance<BLC_list_const_iterator>(m_BLC_transitions.begin(),blc_it.start_same_BLC) << " -- "
                             << std::distance<BLC_list_const_iterator>(m_BLC_transitions.begin(),calculate_end_same_BLC(blc_it)) << "\n";
        for(BLC_list_const_iterator i=blc_it.start_same_BLC; points_into_BLC_set(i, blc_it); ++i)
        {
          if (i == blc_it.start_marked_BLC)
          {
            mCRL2log(log::debug) << "        (The following transitions are marked.)\n";
          }
          const transition& t=m_aut.get_transitions()[*i];
          mCRL2log(log::debug) << "        " << t.from() << " -" << pp(m_aut.action_label(t.label())) << "-> " << t.to();
          if (is_inert_during_init(t) && m_states[t.from()].block == m_states[t.to()].block)
          {
            mCRL2log(log::debug) << " (block-inert)\n";
          }
          else if (is_inert_during_init(t) && m_blocks[m_states[t.from()].block].constellation == m_blocks[m_states[t.to()].block].constellation)
          {
            mCRL2log(log::debug) << " (constellation-inert)\n";
          }
          else if (m_preserve_divergence && t.from() == t.to() && m_aut.is_tau(m_aut.apply_hidden_label_map(t.label())))
          {
            mCRL2log(log::debug) << " (divergent self-loop)\n";
          }
          else
          {
            mCRL2log(log::debug) << " (in constellation " << m_blocks[m_states[t.to()].block].constellation << ")\n";
          }
        }
      }
      mCRL2log(log::debug) << "  BLC_List end\n";
    }


    void print_data_structures(const std::string& header, const bool initialisation=false) const
    {
      if (!mCRL2logEnabled(log::debug))  return;
      mCRL2log(log::debug) << "========= PRINT DATASTRUCTURE: " << header << " =======================================\n";
      mCRL2log(log::debug) << "++++++++++++++++++++  States     ++++++++++++++++++++++++++++\n";
      for(state_index si=0; si<m_states.size(); ++si)
      {
        mCRL2log(log::debug) << "State " << si <<" (Block: " << m_states[si].block <<"):\n";
        mCRL2log(log::debug) << "  #Inert outgoing transitions: " << m_states[si].no_of_outgoing_inert_transitions << '\n';

        mCRL2log(log::debug) << "  Incoming transitions:\n";
        std::vector<transition>::const_iterator end=(si+1==m_states.size()?m_aut.get_transitions().end():m_states[si+1].start_incoming_transitions);
        for(std::vector<transition>::const_iterator it=m_states[si].start_incoming_transitions; it!=end; ++it)
        {
           mCRL2log(log::debug) << "  " << ptr(*it) << "\n";
        }

        mCRL2log(log::debug) << "  Outgoing transitions:\n";
        for(outgoing_transitions_const_it it=m_states[si].start_outgoing_transitions;
                        it!=m_outgoing_transitions.end() &&
                        (si+1>=m_states.size() || it!=m_states[si+1].start_outgoing_transitions);
                     ++it)
        {
           const transition& t=m_aut.get_transitions()[initialisation?it->transition:m_BLC_transitions[it->transition]];
           mCRL2log(log::debug) << "  " << t.from() << " -" << m_aut.action_label(t.label()) << "-> " << t.to() << "\n";;
        }
        mCRL2log(log::debug) << "  Ref states in blocks: " << *(m_states[si].ref_states_in_blocks) << ". Must be " << si <<".\n";
        mCRL2log(log::debug) << "---------------------------------------------------\n";
      }
      mCRL2log(log::debug) << "++++++++++++++++++++ Transitions ++++++++++++++++++++++++++++\n";
      for(state_index ti=0; ti<m_transitions.size(); ++ti)
      {
        const transition& t=m_aut.get_transitions()[ti];
        mCRL2log(log::debug) << "Transition " << ti <<": " << t.from()
                                              << " -" << m_aut.action_label(t.label()) << "-> "
                                              << t.to() << "\n";
      }

      mCRL2log(log::debug) << "++++++++++++++++++++ Blocks ++++++++++++++++++++++++++++\n";
      for(block_index bi=0; bi<m_blocks.size(); ++bi)
      {
        mCRL2log(log::debug) << "  Block " << bi << " (const: " << m_blocks[bi].constellation <<"):\n";
        mCRL2log(log::debug) << "  Bottom states: ";
        for(typename std::vector<state_index>::const_iterator sit=m_blocks[bi].start_bottom_states;
                        sit!=m_blocks[bi].start_non_bottom_states; ++sit)
        {
          mCRL2log(log::debug) << *sit << "  ";
        }
        if (m_branching)
        {
          mCRL2log(log::debug) << "\n  Non-bottom states: ";
          for(typename std::vector<state_index>::const_iterator sit=m_blocks[bi].start_non_bottom_states;
                                 sit!=m_blocks[bi].end_states; ++sit)
          {
            mCRL2log(log::debug) << *sit << "  ";
          }
        }
        else
        {
          assert(m_blocks[bi].start_non_bottom_states == m_blocks[bi].end_states);
        }
        if (!initialisation)
        {
          display_BLC_list(bi);
        }
        mCRL2log(log::debug) << "\n";
      }

      mCRL2log(log::debug) << "++++++++++++++++++++ Constellations ++++++++++++++++++++++++++++\n";
      for(state_index ci=0; ci<m_constellations.size(); ++ci)
      {
        mCRL2log(log::debug) << "  Constellation " << ci << ":\n";
        mCRL2log(log::debug) << "    Blocks in constellation:";
        for (std::vector<state_index>::const_iterator constln_it = m_constellations[ci].start_const_states; constln_it < m_constellations[ci].end_const_states; )
        {
          const block_index bi = m_states[*constln_it].block;
          mCRL2log(log::debug) << " " << bi;
          constln_it = m_blocks[bi].end_states;
        }
        mCRL2log(log::debug) << "\n";
      }
      mCRL2log(log::debug) << "Non-trivial constellations:";
      for(const constellation_index ci: m_non_trivial_constellations)
      {
        mCRL2log(log::debug) << " " << ci;
      }

      mCRL2log(log::debug) << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      mCRL2log(log::debug) << "Outgoing transitions:\n";

      for(outgoing_transitions_const_it pi = m_outgoing_transitions.begin(); pi < m_outgoing_transitions.end(); ++pi)
      {
        const transition& t=m_aut.get_transitions()[initialisation?pi->transition:m_BLC_transitions[pi->transition]];
        mCRL2log(log::debug) << "  " << t.from() << " -" << m_aut.action_label(t.label()) << "-> " << t.to();
        if (m_outgoing_transitions.begin() <= pi->start_same_saC && pi->start_same_saC < m_outgoing_transitions.end())
        {
          const transition& t1=m_aut.get_transitions()[initialisation?pi->start_same_saC->transition:m_BLC_transitions[pi->start_same_saC->transition]];
          mCRL2log(log::debug) << "  \t(same saC: " << t1.from() << " -" << m_aut.action_label(t1.label()) << "-> " << t1.to();
          const label_index t_label = label_or_divergence(t);
          if (pi->start_same_saC->start_same_saC == pi)
          {
            // Transition t must be the beginning and/or the end of a saC-slice
            if (pi->start_same_saC >= pi && pi > m_outgoing_transitions.begin())
            {
              // Transition t must be the beginning of a saC-slice
              const transition& prev_t = m_aut.get_transitions()[initialisation ? std::prev(pi)->transition : m_BLC_transitions[std::prev(pi)->transition]];
              if (prev_t.from() == t.from() &&
                  label_or_divergence(prev_t) == t_label &&
                  m_blocks[m_states[prev_t.to()].block].constellation == m_blocks[m_states[t.to()].block].constellation)
              {
                mCRL2log(log::debug) << " -- error: not the beginning of a saC-slice";
              }
            }
            if (pi->start_same_saC <= pi && std::next(pi) < m_outgoing_transitions.end())
            {
              // Transition t must be the end of a saC-slice
              const transition& next_t = m_aut.get_transitions()[initialisation ? std::next(pi)->transition : m_BLC_transitions[std::next(pi)->transition]];
              if (next_t.from() == t.from() &&
                  label_or_divergence(next_t) == t_label &&
                  m_blocks[m_states[next_t.to()].block].constellation == m_blocks[m_states[t.to()].block].constellation)
              {
                mCRL2log(log::debug) << " -- error: not the end of a saC-slice";
              }
            }
          }
          else if (pi->start_same_saC > pi ? pi->start_same_saC->start_same_saC > pi : pi->start_same_saC->start_same_saC < pi)
          {
            mCRL2log(log::debug) << " -- error: not in its own saC-slice";
          }
          mCRL2log(log::debug) << ')';
        }
        mCRL2log(log::debug) << '\n';
      }

      mCRL2log(log::debug) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      mCRL2log(log::debug) << "New bottom blocks to be investigated:";

      for(block_index bi: m_blocks_with_new_bottom_states)
      {
        mCRL2log(log::debug) << " " << bi;
      }

      mCRL2log(log::debug) << "\n========= END PRINT DATASTRUCTURE: " << header << " =======================================\n";
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
        m_outgoing_transitions(aut.num_transitions()),
        m_transitions(aut.num_transitions()),
        m_states_in_blocks(aut.num_states()),
        m_blocks(1,{m_states_in_blocks.begin(),0}),
        m_constellations(1,constellation_type(m_states_in_blocks.begin(), m_states_in_blocks.end())),   // Algorithm 1, line 1.2.
        m_BLC_transitions(aut.num_transitions()),
        m_branching(branching),
        m_preserve_divergence(preserve_divergence)
    {
      assert(m_branching || !m_preserve_divergence);
      mCRL2log(log::verbose) << "Start initialisation.\n";
      create_initial_partition();
      mCRL2log(log::verbose) << "After initialisation there are " << m_blocks.size() << " equivalence classes. Start refining. \n";
      refine_partition_until_it_becomes_stable();
      assert(check_data_structures("READY"));
    }


    /// \brief Calculate the number of equivalence classes
    /// \details The number of equivalence classes (which is valid after the
    /// partition has been constructed) is equal to the number of states in the
    /// bisimulation quotient.
    std::size_t num_eq_classes() const
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
    state_index get_eq_class(const state_index si) const
    {
      assert(si<m_states.size());
      return m_states[si].block;
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
      // The transitions are most efficiently directly extracted from the block_to_constellation lists in blocks.
      std::vector<transition> T;
      for(block_index bi=0; bi<m_blocks.size(); ++bi)
      {
        const block_type& B=m_blocks[bi];
        //mCRL2complexity_gj(&B, add_work(..., 1), *this);
            // Because every block is touched exactly once, we do not store a physical counter for this.
        for(const BLC_indicators blc_ind: B.block_to_constellation)
        {
          // mCRL2complexity_gj(&blc_ind, add_work(..., 1), *this);
              // Because every BLC set is touched exactly once, we do not store a physical counter for this.
          const transition& t= m_aut.get_transitions()[*blc_ind.start_same_BLC];
          const transition_index new_to=get_eq_class(t.to());
          if (!is_inert_during_init(t) || bi!=new_to)
          {
            T.emplace_back(bi, t.label(), new_to);
          }
        }
      }
      m_aut.clear_transitions();
      for (const transition& t: T)
      {
        //mCRL2complexity_gj(..., add_work(..., 1), *this);
            // we do not add a counter because every transition has been generated by one of the above iterations.
        m_aut.add_transition(t);
      }
      //
      // Merge the states, by setting the state labels of each state to the
      // concatenation of the state labels of its equivalence class.

      if (m_aut.has_state_info())   /* If there are no state labels this step is not needed */
      {
        /* Create a vector for the new labels */
        std::vector<typename LTS_TYPE::state_label_t> new_labels(num_eq_classes());


        for(std::size_t i=0; i<m_aut.num_states(); ++i)
        {
          //mCRL2complexity_gj(&m_states[i], add_work(..., 1), *this);
              // Because every state is touched exactly once, we do not store a physical counter for this.
          const state_index new_index(get_eq_class(i));
          new_labels[new_index]=new_labels[new_index]+m_aut.state_label(i);
        }

        m_aut.set_num_states(num_eq_classes());
        for (std::size_t i=0; i<num_eq_classes(); ++i)
        {
          // mCRL2complexity_gj(&m_blocks[i], add_work(check_complexity::finalize_minimized_LTS__set_labels_of_block, 1), *this);
              // Because every block is touched exactly once, we do not store a physical counter for this.
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
    bool in_same_class(state_index const s, state_index const t) const
    {
      return get_eq_class(s) == get_eq_class(t);
    }
  protected:

    std::string ptr(const transition& t) const
    {
      return std::to_string(t.from()) + " -" + pp(m_aut.action_label(t.label())) + "-> " + std::to_string(t.to());
    }

    std::string ptr(const transition_index ti) const
    {
      const transition& t=m_aut.get_transitions()[ti];
      return ptr(t);
    }

    /*--------------------------- main algorithm ----------------------------*/

    /*----------------- splitB -- Algorithm 3 of [GJ 2024] -----------------*/

    state_index number_of_states_in_block(const block_index B) const
    {
      assert(m_blocks[B].start_bottom_states < m_blocks[B].end_states);
      return std::distance(m_blocks[B].start_bottom_states, m_blocks[B].end_states);
    }

    state_index number_of_states_in_constellation(const constellation_index C) const
    {
      assert(m_constellations[C].start_const_states < m_constellations[C].end_const_states);
      return std::distance(m_constellations[C].start_const_states, m_constellations[C].end_const_states);
    }

    void swap_states_in_states_in_block(
              typename std::vector<state_index>::iterator pos1,
              typename std::vector<state_index>::iterator pos2)
    {
      assert(m_states_in_blocks.begin() <= pos1);
      assert(pos1 < m_states_in_blocks.end());
      assert(m_states_in_blocks.begin() <= pos2);
      assert(pos2 < m_states_in_blocks.end());
      if (pos1 != pos2)
      {
        std::swap(*pos1,*pos2);
        m_states[*pos1].ref_states_in_blocks=pos1;
        m_states[*pos2].ref_states_in_blocks=pos2;
      }
    }

    // Move pos1 to pos2, pos2 to pos3 and pos3 to pos1;
    // The function requires that pos3 lies in between pos1 and pos2.
    void swap_states_in_states_in_block(
              typename std::vector<state_index>::iterator pos1,
              typename std::vector<state_index>::iterator pos2,
              typename std::vector<state_index>::iterator pos3)
    {
      assert(m_states_in_blocks.begin() <= pos2);
      assert(pos2 <= pos3); assert(pos3 <= pos1);
      assert(pos2 <= pos3); assert(pos3 <= pos1);
      assert(pos1 < m_states_in_blocks.end());
      if (pos1 == pos3 || pos2 == pos3)
      {
        swap_states_in_states_in_block(pos1, pos2);
        return;
      }
      const state_index temp=*pos1;
      *pos1=*pos3;
      *pos3=*pos2;
      *pos2=temp;

      m_states[*pos1].ref_states_in_blocks=pos1;
      m_states[*pos2].ref_states_in_blocks=pos2;
      m_states[*pos3].ref_states_in_blocks=pos3;
    }

    // Swap the range [pos1, pos1 + count) with the range [pos2, pos2 + count).
    // pos1 must come before pos2.
    // (If the ranges overlap, only swap the non-overlapping part.)
    // The function requires count > 0 and pos1 < pos2
    // (this is sufficient for how it's used below: to swap new bottom states into their proper places;
    // also, the work counters assume that [pos2, pos2 + count) is assigned the work.)
    void multiple_swap_states_in_states_in_block(
              typename std::vector<state_index>::iterator pos1,
              typename std::vector<state_index>::iterator pos2,
              state_index count
              #ifdef CHECK_COMPLEXITY_GJ
                , const state_index max_B
              #endif
              )
    {
      assert(0 < count);
      assert(count <= m_states_in_blocks.size());
      assert(m_states_in_blocks.begin() <= pos1);
      if (pos1 > pos2)  std::swap(pos1, pos2);
      assert(pos2 <= m_states_in_blocks.end() - count);
      if (std::distance(pos2, pos1 + count) > 0)
      {
        // state_index is unsigned, so we cannot test "if (overlap > 0)"
        state_index overlap = std::distance(pos2, pos1 + count);
        count -= overlap;
        pos2 += overlap;
      }
      state_index temp = *pos1;
      while (--count > 0)
      {
        mCRL2complexity_gj(&m_states[*pos2], add_work(check_complexity::multiple_swap_states_in_block__swap_state_in_small_block, max_B), *this);
        *pos1 = *pos2;
        m_states[*pos1].ref_states_in_blocks = pos1;
        ++pos1;
        *pos2 = *pos1;
        m_states[*pos2].ref_states_in_blocks = pos2;
        ++pos2;
      }
      *pos1 = *pos2;
      m_states[*pos1].ref_states_in_blocks = pos1;
      *pos2 = temp;
      m_states[*pos2].ref_states_in_blocks = pos2;

      #ifndef NDEBUG
        for (std::vector<state_type_gj>::const_iterator si = m_states.cbegin(); si < m_states.cend(); ++si)
        {
          assert(static_cast<std::size_t>(std::distance(m_states.cbegin(), si)) == *si->ref_states_in_blocks);
        }
      #endif
    }

    // marks the transition indicated by out_pos.
    // (We use an outgoing_transitions_it because it points to the m_BLC_transitions entry that needs to be updated.)
    void mark_BLC_transition(const outgoing_transitions_it out_pos)
    {
      BLC_list_iterator old_pos = m_BLC_transitions.begin() + out_pos->transition;
      linked_list<BLC_indicators>::iterator ind = m_transitions[*old_pos].transitions_per_block_to_constellation;
      assert(ind->start_same_BLC <= old_pos);
      assert(old_pos < m_BLC_transitions.end());
      assert(points_into_BLC_set(old_pos, *ind));
      if (old_pos < ind->start_marked_BLC)
      {
        // The transition is not marked
        BLC_list_iterator new_pos = std::prev(ind->start_marked_BLC);
        assert(ind->start_same_BLC <= new_pos);
        assert(new_pos < m_BLC_transitions.end());
        assert(points_into_BLC_set(new_pos, *ind));
        if (old_pos < new_pos)
        {
          std::swap(*old_pos, *new_pos);
          m_transitions[*old_pos].ref_outgoing_transitions->transition = std::distance(m_BLC_transitions.begin(), old_pos);
          assert(out_pos == m_transitions[*new_pos].ref_outgoing_transitions);
          out_pos->transition = std::distance(m_BLC_transitions.begin(), new_pos);
        }
        ind->start_marked_BLC--;
      }

      #ifndef NDEBUG
        for (BLC_list_const_iterator it = m_BLC_transitions.cbegin(); it < m_BLC_transitions.cend(); ++it)
        {
          assert(m_transitions[*it].ref_outgoing_transitions->transition == static_cast<std::size_t>(std::distance(m_BLC_transitions.cbegin(), it)));
          assert(m_transitions[*it].transitions_per_block_to_constellation->start_same_BLC <= it);
          assert(points_into_BLC_set(it, *m_transitions[*it].transitions_per_block_to_constellation));
        }
      #endif
    }

    // Split the block B by moving the elements in R to the front in m_states, and add a
    // new element B_new at the end of m_blocks referring to R. Adapt B.start_bottom_states,
    // B.start_non_bottom_states and B.end_states, and do the same for B_new.
    block_index split_block_B_into_R_and_BminR(
                     const block_index B,
                     std::vector<state_index>::iterator first_bottom_state_in_R,
                     std::vector<state_index>::iterator last_bottom_state_in_R,
                     const todo_state_vector& R)
    {
//std::cerr << "SPLIT BLOCK " << B << " by removing"; for(auto s = first_bottom_state; s < last_bottom_state; ++s){ std::cerr << ' ' << *s;} for(auto s:R){ std::cerr << ' ' << s; } std::cerr << "\n";
      assert(m_blocks[B].start_bottom_states <= first_bottom_state_in_R);
      assert(first_bottom_state_in_R <= last_bottom_state_in_R);
      assert(last_bottom_state_in_R <= m_blocks[B].start_non_bottom_states);
      // Basic administration. Make a new block and add it to the current constellation.
      const block_index B_new=m_blocks.size();
      const constellation_index ci = m_blocks[B].constellation;
      m_blocks.emplace_back(m_blocks[B].start_bottom_states,ci);
      #ifdef CHECK_COMPLEXITY_GJ
        m_blocks[B_new].work_counter = m_blocks[B].work_counter;
      #endif
      if (m_states[*m_constellations[ci].start_const_states].block == m_states[*std::prev(m_constellations[ci].end_const_states)].block) // This constellation is trivial.
      {
        // This constellation is trivial, as it will be split add it to the non trivial constellations.
        assert(std::find(m_non_trivial_constellations.begin(),
                         m_non_trivial_constellations.end(),
                         ci)==m_non_trivial_constellations.end());
        m_non_trivial_constellations.emplace_back(ci);
      }

      // Carry out the split.
      #ifdef CHECK_COMPLEXITY_GJ
        // The size of the new block is not yet fixed.
        const state_index max_B = check_complexity::log_n - check_complexity::ilog2(std::distance(first_bottom_state_in_R, last_bottom_state_in_R) + R.size());
      #endif
      if (m_blocks[B].start_bottom_states < first_bottom_state_in_R)
      {
        multiple_swap_states_in_states_in_block(m_blocks[B].start_bottom_states, first_bottom_state_in_R,
                    std::distance(first_bottom_state_in_R, last_bottom_state_in_R)
                    #ifdef CHECK_COMPLEXITY_GJ
                      , max_B
                    #endif
                    );
        last_bottom_state_in_R = m_blocks[B].start_bottom_states + std::distance(first_bottom_state_in_R, last_bottom_state_in_R);
        first_bottom_state_in_R = m_blocks[B].start_bottom_states;
      }
      m_blocks[B_new].start_non_bottom_states = last_bottom_state_in_R;
      m_blocks[B_new].end_states = last_bottom_state_in_R;
      m_blocks[B].start_bottom_states = last_bottom_state_in_R;
      for(state_index s: R)
      {
        mCRL2complexity_gj(&m_states[s], add_work(check_complexity::split_block_B_into_R_and_BminR__carry_out_split,
                max_B), *this);
//std::cerr << "MOVE STATE TO NEW BLOCK: " << s << "\n";
        m_states[s].block=B_new;
        typename std::vector<state_index>::iterator pos=m_states[s].ref_states_in_blocks;
        assert(pos>=m_blocks[B].start_non_bottom_states); // the state is a non bottom state.
            // pos --> B.start_bottom_states --> B.start_non_bottom_states --> pos.
          swap_states_in_states_in_block(pos, m_blocks[B].start_bottom_states, m_blocks[B].start_non_bottom_states);
          m_blocks[B].start_non_bottom_states++;
          m_blocks[B].start_bottom_states++;
          m_blocks[B_new].end_states++;
          assert(m_blocks[B].start_bottom_states<=m_blocks[B].start_non_bottom_states);
          assert(m_blocks[B_new].start_bottom_states<m_blocks[B_new].end_states);
      }
      return B_new;
    }

    // It is assumed that the new block is located precisely before the old_block in m_BLC_transitions.
    // This routine can not be used in the initialisation phase. It can only be used during refinement.
    // The routine returns true if the last element of old_BLC_block has been removed.
    bool swap_in_the_doubly_linked_list_LBC_in_blocks(
               const transition_index ti,
               linked_list<BLC_indicators>::iterator new_BLC_block,
               linked_list<BLC_indicators>::iterator old_BLC_block
               #ifdef SAVE_BLC_MEMORY
                 , const block_index old_from_block,
                   const constellation_index old_to_constellation
               #endif
               )
    {
      assert(new_BLC_block->start_same_BLC <= new_BLC_block->start_marked_BLC);
      assert(old_BLC_block->start_same_BLC <= old_BLC_block->start_marked_BLC);
      #ifndef SAVE_BLC_MEMORY
        assert(new_BLC_block->start_marked_BLC <= new_BLC_block->end_same_BLC);
        assert(new_BLC_block->end_same_BLC==old_BLC_block->start_same_BLC);
        assert(old_BLC_block->start_marked_BLC <= old_BLC_block->end_same_BLC);
      #endif
      BLC_list_iterator old_position=m_BLC_transitions.begin()+m_transitions[ti].ref_outgoing_transitions->transition;
      assert(old_BLC_block->start_same_BLC <= old_position);
      assert(points_into_BLC_set(old_position, *old_BLC_block));
      assert(m_transitions[ti].transitions_per_block_to_constellation == old_BLC_block);
      assert(ti == *old_position);
      if (old_position < old_BLC_block->start_marked_BLC)
      {
// std::cerr << "Moving unmarked " << m_transitions[*old_position].debug_id(*this);
        assert(old_BLC_block->start_same_BLC <= old_position);
        assert(new_BLC_block->start_marked_BLC <= old_BLC_block->start_same_BLC);
        if (old_position != new_BLC_block->start_marked_BLC)
        {
          transition_index temp = *old_position;
          *old_position = *old_BLC_block->start_same_BLC;
          *old_BLC_block->start_same_BLC = *new_BLC_block->start_marked_BLC;
          *new_BLC_block->start_marked_BLC = temp;
          m_transitions[*old_position].ref_outgoing_transitions->transition = std::distance(m_BLC_transitions.begin(),old_position);
          m_transitions[*old_BLC_block->start_same_BLC].ref_outgoing_transitions->transition = std::distance(m_BLC_transitions.begin(), old_BLC_block->start_same_BLC);
          m_transitions[*new_BLC_block->start_marked_BLC].ref_outgoing_transitions->transition = std::distance(m_BLC_transitions.begin(), new_BLC_block->start_marked_BLC);
        }
        else
        {
          assert(old_position == old_BLC_block->start_same_BLC);
        }
        new_BLC_block->start_marked_BLC++;
      }
      else
      {
// std::cerr << "Moving marked " << m_transitions[*old_position].debug_id(*this);
        assert(old_BLC_block->start_marked_BLC <= old_position);
        assert(old_BLC_block->start_same_BLC <= old_BLC_block->start_marked_BLC);
        if (old_position != old_BLC_block->start_same_BLC)
        {
          transition_index temp = *old_position;
          *old_position = *old_BLC_block->start_marked_BLC;
          *old_BLC_block->start_marked_BLC = *old_BLC_block->start_same_BLC;
          *old_BLC_block->start_same_BLC = temp;
          m_transitions[*old_position].ref_outgoing_transitions->transition = std::distance(m_BLC_transitions.begin(),old_position);
          m_transitions[*old_BLC_block->start_marked_BLC].ref_outgoing_transitions->transition = std::distance(m_BLC_transitions.begin(), old_BLC_block->start_marked_BLC);
          m_transitions[*old_BLC_block->start_same_BLC].ref_outgoing_transitions->transition = std::distance(m_BLC_transitions.begin(), old_BLC_block->start_same_BLC);
        }
        else
        {
          assert(old_position == old_BLC_block->start_marked_BLC);
        }
        old_BLC_block->start_marked_BLC++;
      }
      m_transitions[ti].transitions_per_block_to_constellation=new_BLC_block;
      #ifndef SAVE_BLC_MEMORY
        new_BLC_block->end_same_BLC++;
      #endif
      old_BLC_block->start_same_BLC++;
// std::cerr << " to new " << new_BLC_block->debug_id(*this) << '\n';
      #ifdef SAVE_BLC_MEMORY
        if (old_BLC_block->start_same_BLC >= m_BLC_transitions.end())
        {
          return true;
        }
        const transition& old_t = m_aut.get_transitions()[*old_BLC_block->start_same_BLC];
        return m_states[old_t.from()].block != old_from_block ||
               label_or_divergence(old_t) != label_or_divergence(m_aut.get_transitions()[ti]) ||
               m_blocks[m_states[old_t.to()].block].constellation != old_to_constellation;
      #else
        return old_BLC_block->start_same_BLC==old_BLC_block->end_same_BLC;
      #endif
    }

    // Move the transition t with transition index ti to a new
    // LBC list as the target state switches to a new constellation.
    // Returns true if a new BLC-entry has been created.
    bool update_the_doubly_linked_list_LBC_new_constellation(
               const block_index index_block_B,
               const transition& t,
               const transition_index ti)
    {
      assert(m_states[t.to()].block==index_block_B);
      assert(&m_aut.get_transitions()[ti] == &t);
      bool last_element_removed;
      bool new_block_created = false;
      linked_list<BLC_indicators>::iterator this_block_to_constellation=
                           m_transitions[ti].transitions_per_block_to_constellation;
      assert(this_block_to_constellation!= m_blocks[m_states[t.from()].block].block_to_constellation.end());
      // if this transition is inert, it is inserted in a block in front. Otherwise, it is inserted after
      // the current element in the list.
      if (is_inert_during_init(t) && m_states[t.from()].block==index_block_B)
      {
        linked_list<BLC_indicators>::iterator first_block_to_constellation=m_blocks[m_states[t.from()].block].block_to_constellation.begin();
        #ifndef SAVE_BLC_MEMORY
          assert(first_block_to_constellation->start_same_BLC < first_block_to_constellation->end_same_BLC);
        #endif
        assert(m_states[m_aut.get_transitions()[*(first_block_to_constellation->start_same_BLC)].from()].block==index_block_B);
        assert(m_aut.is_tau(m_aut.apply_hidden_label_map(m_aut.get_transitions()[*(first_block_to_constellation->start_same_BLC)].label())));
        if (first_block_to_constellation==this_block_to_constellation)
        {
          // Make a new entry in the list block_to_constellation, at the beginning;

          first_block_to_constellation=
                  m_blocks[m_states[t.from()].block].block_to_constellation.
                           emplace_front(//first_block_to_constellation,
                                   this_block_to_constellation->start_same_BLC,
                                   this_block_to_constellation->start_same_BLC);
          #ifdef CHECK_COMPLEXITY_GJ
            first_block_to_constellation->work_counter = this_block_to_constellation->work_counter;
          #endif
          new_block_created = true;
        }
        else  assert(m_states[m_aut.get_transitions()[*(first_block_to_constellation->start_same_BLC)].to()].block==index_block_B);
        last_element_removed=swap_in_the_doubly_linked_list_LBC_in_blocks(ti,
                    first_block_to_constellation, this_block_to_constellation
                    #ifdef SAVE_BLC_MEMORY
                      , /*m_states[t.from()].block*/ index_block_B, m_blocks[index_block_B].constellation
                    #endif
                    );
      }
      else
      {
        linked_list<BLC_indicators>::iterator next_block_to_constellation=this_block_to_constellation;
        ++next_block_to_constellation;
        const label_index t_label = label_or_divergence(t);
        const transition* first_t;
        if (next_block_to_constellation==m_blocks[m_states[t.from()].block].block_to_constellation.end() ||
            (first_t = &m_aut.get_transitions()[*(next_block_to_constellation->start_same_BLC)],
             assert(m_states[first_t->from()].block==m_states[t.from()].block),
             m_states[first_t->to()].block!=index_block_B) ||
            (label_or_divergence(*first_t) != t_label))
        {
// std::cerr << "Creating new BLC set for " << m_transitions[ti].debug_id(*this) << ": ";
          // Make a new entry in the list next_block_to_constellation, after the current list element.
          next_block_to_constellation=
                  m_blocks[m_states[t.from()].block].block_to_constellation.
                           emplace_after(this_block_to_constellation,
                                         this_block_to_constellation->start_same_BLC,
                                         this_block_to_constellation->start_same_BLC);
          #ifdef CHECK_COMPLEXITY_GJ
            next_block_to_constellation->work_counter = this_block_to_constellation->work_counter;
          #endif
          new_block_created = true;
        }
// else { std::cerr << "Extending existing BLC set for " << m_transitions[ti].debug_id(*this) << ": "; }
        last_element_removed=swap_in_the_doubly_linked_list_LBC_in_blocks(ti,
                    next_block_to_constellation, this_block_to_constellation
                    #ifdef SAVE_BLC_MEMORY
                      , m_states[t.from()].block, m_blocks[index_block_B].constellation
                    #endif
                    );
      }

      if (last_element_removed)
      {
        m_blocks[m_states[t.from()].block].block_to_constellation.erase(this_block_to_constellation);
      }
      #ifndef NDEBUG
        check_transitions(false, false);
      #endif
      return new_block_created;
    }

    // Update the LBC list of a transition, when the from state of the transition moves
    // from block old_bi to new_bi.
    transition_index update_the_doubly_linked_list_LBC_new_block(
               const block_index old_bi,
               const block_index new_bi,
               const transition_index ti)
    {
// std::cerr << "update_the_doubly_linked_list_LBC_new_block(old_bi = " << old_bi << ", new_bi = " << new_bi << ", transition_index = " << ti << ")\n";
      const transition& t=m_aut.get_transitions()[ti];
      transition_index remaining_transition=null_transition;

      assert(m_states[t.from()].block==new_bi);

      linked_list<BLC_indicators>::iterator this_block_to_constellation=
                           m_transitions[ti].transitions_per_block_to_constellation;
      transition_index co_transition=null_transition;
      bool co_block_found=false;
      if (this_block_to_constellation->start_same_BLC!=m_BLC_transitions.begin())
      {
        co_transition=*(this_block_to_constellation->start_same_BLC-1);
        const transition& co_t=m_aut.get_transitions()[co_transition];
        co_block_found=m_states[co_t.from()].block==new_bi &&
                       label_or_divergence(co_t) == label_or_divergence(t) &&
                       m_blocks[m_states[co_t.to()].block].constellation==m_blocks[m_states[t.to()].block].constellation;
      }

      bool last_element_removed;

      if (!co_block_found)
      {
        // Make a new entry in the list next_block_to_constellation;

        // Put inert tau's to the front. Otherwise, the new block is put after the current block.
        linked_list<BLC_indicators>::iterator new_position;
        BLC_list_iterator old_BLC_start=this_block_to_constellation->start_same_BLC;
        if (m_blocks[new_bi].block_to_constellation.empty() ||
            (is_inert_during_init(t) &&
             m_blocks[new_bi].constellation==m_blocks[m_states[t.to()].block].constellation))
        {
          m_blocks[new_bi].block_to_constellation.emplace_front(old_BLC_start, old_BLC_start);
          new_position=m_blocks[new_bi].block_to_constellation.begin();
        }
        else
        {
          new_position=m_blocks[new_bi].block_to_constellation.begin();
          new_position= m_blocks[new_bi].block_to_constellation.emplace_after(new_position,old_BLC_start, old_BLC_start);
        }
        #ifdef CHECK_COMPLEXITY_GJ
          new_position->work_counter = this_block_to_constellation->work_counter;
        #endif
        last_element_removed=swap_in_the_doubly_linked_list_LBC_in_blocks(ti,
                    new_position, this_block_to_constellation
                    #ifdef SAVE_BLC_MEMORY
                      , old_bi, m_blocks[m_states[t.to()].block].constellation
                    #endif
                    );
      }
      else
      {
        // Move the current transition to the next list indicated by the iterator it.
        linked_list<BLC_indicators>::iterator new_BLC_block= m_transitions[co_transition].transitions_per_block_to_constellation;
        last_element_removed=swap_in_the_doubly_linked_list_LBC_in_blocks(ti,
                    new_BLC_block, this_block_to_constellation
                    #ifdef SAVE_BLC_MEMORY
                      , old_bi, m_blocks[m_states[t.to()].block].constellation
                    #endif
                    );
      }

      if (last_element_removed)
      {
        // Remove this element.
        m_blocks[old_bi].block_to_constellation.erase(this_block_to_constellation);
      }
      else
      {
        remaining_transition= *(this_block_to_constellation->start_same_BLC);
      }
      #ifndef NDEBUG
        check_transitions(false, false);
      #endif

      return remaining_transition;
    }

    // Set m_states[s].counter:=undefined for all s in m_R and m_U.
    void clear_state_counters(bool restrict_to_R=false)
    {
      for(const state_index si: m_R)
      {
        assert(Rmarked == m_states[si].counter); // this allows us to charge the work in this loop to setting the counter to Rmarked
        m_states[si].counter=undefined;
      }
      if (restrict_to_R)
      {
        return;
      }
      for(const state_index si: m_U_counter_reset_vector)
      {
        // this work is charged to adding a value to m_U_counter_reset_vector
        assert(undefined != m_states[si].counter || m_R.find(si));
        assert(Rmarked != m_states[si].counter);
        m_states[si].counter=undefined;
      }
      clear(m_U_counter_reset_vector);
    }

    // Calculate the states R in block B that can inertly reach a state with a transition in splitter
    // and split B in R and B\R. The complexity is conform the smallest block of R and B\R.
    // The LBC_list and bottom states are not updated.
    // Provide the index of the newly created block as a result. This block is the smallest of R and B\R.
    // Return M_in_bi=true iff the new block bi is the one with the transitions in the splitter.

    // The function assumes that R has already been pre-filled with the sources of
    // the marked transitions in the splitter, and therefore a split is always
    // needed. Bottom states in R are those in the range
    // [m_blocks[B].start_bottom_states, first_unmarked_bottom_state),
    // non-bottom states in R are already in m_R.

    // Occasionally it is necessary to split using only a subset of a splitter.
    // Then, splitter_end_unmarked_BLC can be used to indicate that this is the case.
    // The function assumes that the *marked* transitions in the splitter have
    // already been used to separate the bottom states in block B. This is indicated
    // using first_unmarked_bottom_state, an iterator into m_states_in_blocks.

    // Return value: the index of subblock R
    block_index simple_splitB(const block_index B,
                              linked_list<BLC_indicators>::iterator splitter,
                              const std::vector<state_index>::iterator first_unmarked_bottom_state,
                              // const bool initialisation,
                              const BLC_list_iterator splitter_end_unmarked_BLC)
// David suggests to change this function as follows:
// - Generally, the function takes a BLC-set as splitter.
//   Some transitions in the BLC-set may be *marked.*
//   The function first moves all sources of marked transitions to R;
//   it assumes that all bottom states that should go to R have a marked transition.
//   (Some non-bottom states may also have such a transition.)
//   Then, it moves the remaining bottom states to U.
// - If at this point it appears that R contains all bottom states, no split is needed.
//   (i.e. it is not necessary to check before calling whether there will be an actual split.)
// - After that, the coroutines are started to separate the non-bottom states
//   as in the current code.
// - The bottom states are not actually inserted into m_R or m_U; the information
//   which states are in R or U is stored through an iterator into m_states_in_blocks,
//   known elsewhere as first_unmarked_bottom_state. This will allow to create
//   the bottom part of U quickly.
//
// - There is one special case: when a block has new bottom states, its old and
//   new bottom states need to be separated. This is done by using the tau-
//   transitions from R to U as splitter. These tau-transitions are
//   constellation-inert; the same BLC-set may contain other constellation-inert
//   tau-transitions out of R but they should be disregarded. In this special
//   case, we do not take the full BLC-set but only the marked transitions in it
//   as splitter.
//   This can be handled by giving an upper bound to the end of the BLC-set.
//
// Always using a BLC-set as splitter would require to initialise the BLC-sets
// before the initial stabilisation. That simplifies some calculations slightly,
// as one does not need to distinguish between initialisation and later phases.
    {
// std::cerr << "simple_splitB(block " << B << ",...)\n";
      assert(1 < number_of_states_in_block(B));
      assert(!m_blocks[B].contains_new_bottom_states);
      assert(m_U.empty());
      assert(m_U_counter_reset_vector.empty());
      // assert(splitter_end_unmarked_BLC == old value of splitter->start_marked_BLC || splitter_end_unmarked_BLC == splitter->start_same_BLC);
      typedef enum { initializing, state_checking, aborted, aborted_after_initialisation,
                     incoming_inert_transition_checking, outgoing_action_constellation_check } status_type;
      status_type U_status=state_checking;
      status_type R_status=initializing;
      BLC_list_iterator M_it = splitter->start_same_BLC;
      assert(m_BLC_transitions.begin() <= splitter->start_same_BLC);
      assert(splitter->start_same_BLC <= M_it);
      assert(M_it <= splitter_end_unmarked_BLC);
      assert(splitter_end_unmarked_BLC <= splitter->start_marked_BLC);
      #ifndef SAVE_BLC_MEMORY
        assert(splitter->start_marked_BLC <= splitter->end_same_BLC);
        assert(splitter->start_same_BLC < splitter->end_same_BLC);
        assert(splitter->end_same_BLC <= m_BLC_transitions.end());
      #endif
      assert(!has_marked_transitions(*splitter));
      assert(m_states[m_aut.get_transitions()[*splitter->start_same_BLC].from()].block == B);
      std::vector<transition>::iterator current_U_incoming_transition_iterator;
      std::vector<transition>::iterator current_U_incoming_transition_iterator_end;
      state_index current_U_outgoing_state = null_state;
      outgoing_transitions_it current_U_outgoing_transition_iterator;
      outgoing_transitions_it current_U_outgoing_transition_iterator_end;
      std::vector<transition>::iterator current_R_incoming_transition_iterator;
      std::vector<transition>::iterator current_R_incoming_transition_iterator_end;
      std::vector<state_index>::iterator current_R_incoming_bottom_state_iterator = m_blocks[B].start_bottom_states;
      assert(current_R_incoming_bottom_state_iterator <= first_unmarked_bottom_state);
      assert(current_R_incoming_bottom_state_iterator < first_unmarked_bottom_state || !m_R.empty() || M_it < splitter_end_unmarked_BLC);
      std::vector<state_index>::iterator current_U_incoming_bottom_state_iterator = first_unmarked_bottom_state;
      assert(current_U_incoming_bottom_state_iterator < m_blocks[B].start_non_bottom_states);

      const state_index max_R_nonbottom_size = number_of_states_in_block(B)/2 - std::distance(m_blocks[B].start_bottom_states, first_unmarked_bottom_state); // can underflow
      if (number_of_states_in_block(B)/2 < std::distance(m_blocks[B].start_bottom_states, first_unmarked_bottom_state) + m_R.size())
      {
        R_status = (M_it == splitter_end_unmarked_BLC) ? aborted_after_initialisation : aborted;
      }
      else if (m_blocks[B].start_non_bottom_states == m_blocks[B].end_states)
      {
        // There are no non-bottom states, hence we do not need to carry out the tau closure.
        // Also, there cannot be any new bottom states.
        assert(m_blocks[B].start_bottom_states < first_unmarked_bottom_state);
        assert(std::distance(m_blocks[B].start_bottom_states, first_unmarked_bottom_state) <= std::distance(first_unmarked_bottom_state, m_blocks[B].start_non_bottom_states));
        assert(m_R.empty()); // m_R.clear(); is superfluous
        assert(m_U.empty()); // m_U.clear(); is superfluous;
        assert(m_U_counter_reset_vector.empty()); // clear_state_counters(); is superfluous
        // split_block B into R and B\R.
        return split_block_B_into_R_and_BminR(B, m_blocks[B].start_bottom_states, first_unmarked_bottom_state, m_R);
      }
      else if (M_it == splitter_end_unmarked_BLC)
      {
        // There are no more transitions in the splitter whose source states need to be added to R.
        assert(m_blocks[B].start_bottom_states < first_unmarked_bottom_state || !m_R.empty());
        R_status=state_checking;
      }
      const state_index max_U_nonbottom_size = number_of_states_in_block(B)/2 - std::distance(first_unmarked_bottom_state, m_blocks[B].start_non_bottom_states); // can underflow
      if (static_cast<std::ptrdiff_t>(number_of_states_in_block(B)/2) < std::distance(first_unmarked_bottom_state, m_blocks[B].start_non_bottom_states))
      {
        assert(aborted != R_status); assert(aborted_after_initialisation != R_status);
        U_status=aborted_after_initialisation;
      }
      else if (m_blocks[B].start_non_bottom_states == m_blocks[B].end_states)
      {
        // There are no non-bottom states, hence we do not need to carry out the tau closure.
        // Also, there cannot be any new bottom states.
        assert(first_unmarked_bottom_state < m_blocks[B].start_non_bottom_states);
        assert(std::distance(first_unmarked_bottom_state, m_blocks[B].start_non_bottom_states) <= std::distance(m_blocks[B].start_bottom_states, first_unmarked_bottom_state));
        assert(m_R.empty()); // m_R.clear(); is superfluous
        assert(m_U.empty()); // m_U.clear(); is superfluous;
        assert(m_U_counter_reset_vector.empty()); // clear_state_counters(); is superfluous
        // split_block B into U and B\U.
        split_block_B_into_R_and_BminR(B, first_unmarked_bottom_state, m_blocks[B].end_states, m_U);
        return B;
      }
      assert(m_blocks[B].start_non_bottom_states < m_blocks[B].end_states);
      const transition& first_t = m_aut.get_transitions()[*splitter->start_same_BLC];
      const label_index a = label_or_divergence(first_t);
      const constellation_index C = m_blocks[m_states[first_t.to()].block].constellation;

      // Algorithm 3, line 3.2 left.

      // start coroutines. Each co-routine handles one state, and then gives control
      // to the other co-routine. The coroutines can be found sequentially below surrounded
      // by a while loop.

//std::cerr << "simple_splitB() before while\n";
      while (true)
      {
        assert(U_status!=aborted_after_initialisation || (R_status!=aborted && R_status!=aborted_after_initialisation));
#ifndef NDEBUG
        for(state_index si=0; si<m_states.size(); ++si)
        {
          if (m_states[si].block != B || 0 == m_states[si].no_of_outgoing_inert_transitions)
          {
            assert(undefined == m_states[si].counter);
            assert(!m_R.find(si));
            assert(!m_U.find(si));
            assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), si) == m_U_counter_reset_vector.end());
          }
          else
          {
            switch(m_states[si].counter)
            {
            case undefined:  assert(!m_U.find(si)); assert(!m_R.find(si));
                             assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), si) == m_U_counter_reset_vector.end());
                             break;
            case Rmarked:    assert( m_R.find(si)); assert(!m_U.find(si));
                             // It can happen that the state has a tau-transition to a U-state, then it will be in the m_U_counter_reset_vector.
                             break;
            case 0:          assert(!m_R.find(si)); // It can happen that the state is in U or is not in U
                             assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), si) != m_U_counter_reset_vector.end());
                             break;
            default:         assert(!m_R.find(si)); assert(!m_U.find(si));
                             assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), si) != m_U_counter_reset_vector.end());
                             break;
            }
          }
        }
#endif
        // The code for the right co-routine.
        switch (R_status)
        {
          case initializing:
          {
            // Algorithm 3, line 3.3, right.
              assert(M_it < splitter_end_unmarked_BLC);
              const state_index si= m_aut.get_transitions()[*M_it].from();
              mCRL2complexity_gj(&m_transitions[*M_it], add_work(check_complexity::simple_splitB_R__handle_transition_from_R_state, 1), *this);
              assert(m_states[si].block == B);
              assert(!is_inert_during_init(m_aut.get_transitions()[*M_it]) || m_blocks[B].constellation != m_blocks[m_states[m_aut.get_transitions()[*M_it].to()].block].constellation);
              ++M_it;
              if (0 == m_states[si].no_of_outgoing_inert_transitions)
              {
                // The state is a bottom state, it should be in R already
                assert(m_blocks[B].start_bottom_states <= m_states[si].ref_states_in_blocks);
                assert(m_states[si].ref_states_in_blocks < first_unmarked_bottom_state);
                assert(!m_R.find(si));
              }
              else if (m_states[si].counter != Rmarked)
              {
                // The state is a nonbottom state that is not yet in R
                assert(m_blocks[B].start_non_bottom_states <= m_states[si].ref_states_in_blocks);
                assert(m_states[si].ref_states_in_blocks < m_blocks[B].end_states);
                assert(!m_R.find(si));
                m_R.add_todo(si);
                m_states[si].counter=Rmarked;
//std::cerr << "R_todo1 insert: " << si << "\n";
                if (m_R.size() > max_R_nonbottom_size)
                {
                  assert(aborted_after_initialisation != U_status);
                  R_status=aborted;
                  break;
                }
              }
              else assert(m_R.find(si));
              assert(!m_U.find(si));
              if (M_it == splitter_end_unmarked_BLC)
              {
                assert(m_blocks[B].start_non_bottom_states < m_blocks[B].end_states);
                R_status=state_checking;
              }
            break;
          }
          case state_checking:
          {
              const state_index s = current_R_incoming_bottom_state_iterator < first_unmarked_bottom_state
                                    ? *current_R_incoming_bottom_state_iterator++
                                    : m_R.move_from_todo();
              mCRL2complexity_gj(&m_states[s], add_work(check_complexity::simple_splitB_R__find_predecessors, 1), *this);
//std::cerr << "R insert: " << s << "\n";
              assert(m_states[s].block == B);
              if (s+1==m_states.size())
              {
                current_R_incoming_transition_iterator_end=m_aut.get_transitions().end();
              }
              else
              {
                current_R_incoming_transition_iterator_end=m_states[s+1].start_incoming_transitions;
              }
              current_R_incoming_transition_iterator=m_states[s].start_incoming_transitions;
              if (current_R_incoming_transition_iterator!=current_R_incoming_transition_iterator_end &&
                  m_aut.is_tau(m_aut.apply_hidden_label_map(current_R_incoming_transition_iterator->label())))
              {
                R_status=incoming_inert_transition_checking;
              }
            break;
          }
          case incoming_inert_transition_checking:
          {
              assert(current_R_incoming_transition_iterator < current_R_incoming_transition_iterator_end);
              const transition& tr= *current_R_incoming_transition_iterator;
              assert(m_aut.is_tau(m_aut.apply_hidden_label_map(tr.label())));
              mCRL2complexity_gj(&m_transitions[std::distance(m_aut.get_transitions().begin(), current_R_incoming_transition_iterator)],
                        add_work(check_complexity::simple_splitB_R__handle_transition_to_R_state, 1), *this);
              assert(m_states[tr.to()].block == B);
              if (m_states[tr.from()].block==B && !(m_preserve_divergence && tr.from() == tr.to()))
              {
                if (m_states[tr.from()].counter!=Rmarked)
                {
                  assert(!m_R.find(tr.from()));
//std::cerr << "R_todo2 insert: " << tr.from() << "\n";
                  m_R.add_todo(tr.from());
                  m_states[tr.from()].counter=Rmarked;

                  // Algorithm 3, line 3.10 and line 3.11, right.
                  if (m_R.size() > max_R_nonbottom_size)
                  {
                    assert(aborted_after_initialisation != U_status);
                    R_status=aborted_after_initialisation;
                    break;
                  }
                }
                else assert(m_R.find(tr.from()));
                assert(!m_U.find(tr.from()));
              }
              ++current_R_incoming_transition_iterator;
              if (current_R_incoming_transition_iterator==current_R_incoming_transition_iterator_end ||
                  !m_aut.is_tau(m_aut.apply_hidden_label_map(current_R_incoming_transition_iterator->label())))
              {
                R_status=state_checking;
              }
              break;
          }
          default:
          {
            assert(aborted == R_status || aborted_after_initialisation == R_status);
            break;
          }
        }
        if (state_checking == R_status && current_R_incoming_bottom_state_iterator == first_unmarked_bottom_state && m_R.todo_is_empty())
        {
//std::cerr << "R empty: " << "\n";
          // split_block B into R and B\R.
          assert(0 < std::distance(m_blocks[B].start_bottom_states, first_unmarked_bottom_state) + m_R.size());
          assert(std::distance(m_blocks[B].start_bottom_states, first_unmarked_bottom_state) + m_R.size() <= number_of_states_in_block(B)/2);
          clear_state_counters();
          m_U.clear();
          const block_index block_index_of_R = split_block_B_into_R_and_BminR(B, m_blocks[B].start_bottom_states, first_unmarked_bottom_state, m_R);
          m_R.clear();
          return block_index_of_R;
        }

#ifndef NDEBUG
        for(state_index si=0; si<m_states.size(); ++si)
        {
          if (m_states[si].block != B || 0 == m_states[si].no_of_outgoing_inert_transitions)
          {
            assert(undefined == m_states[si].counter);
            assert(!m_R.find(si));
            assert(!m_U.find(si));
            assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), si) == m_U_counter_reset_vector.end());
          }
          else
          {
            switch(m_states[si].counter)
            {
            case undefined:  assert(!m_U.find(si)); assert(!m_R.find(si));
                             assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), si) == m_U_counter_reset_vector.end());
                             break;
            case Rmarked:    assert( m_R.find(si)); assert(!m_U.find(si));
                             // It can happen that the state has a tau-transition to a U-state, then it will be in the m_U_counter_reset_vector.
                             break;
            case 0:          assert(!m_R.find(si)); // It can happen that the state is in U or is not in U
                             assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), si) != m_U_counter_reset_vector.end());
                             break;
            default:         assert(!m_R.find(si)); assert(!m_U.find(si));
                             assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), si) != m_U_counter_reset_vector.end());
                             break;
            }
          }
        }
#endif
        // The code for the left co-routine.
        switch (U_status)
        {
          case state_checking:
          {
//std::cerr << "U_state_checking\n";

            // Algorithm 3, line 3.23 and line 3.24, left.
            const state_index s = current_U_incoming_bottom_state_iterator < m_blocks[B].start_non_bottom_states
                                  ? *current_U_incoming_bottom_state_iterator++
                                  : m_U.move_from_todo();
              assert(!m_R.find(s));
              mCRL2complexity_gj(&m_states[s], add_work(check_complexity::simple_splitB_U__find_predecessors, 1), *this);
//std::cerr << "U insert/ U_todo_remove: " << s << "\n";
              current_U_incoming_transition_iterator=m_states[s].start_incoming_transitions;
              current_U_incoming_transition_iterator_end = s+1 >= m_states.size() ? m_aut.get_transitions().end() : m_states[s+1].start_incoming_transitions;
              if (current_U_incoming_transition_iterator != current_U_incoming_transition_iterator_end &&
                  m_aut.is_tau(m_aut.apply_hidden_label_map(current_U_incoming_transition_iterator->label())))
              {
                U_status=incoming_inert_transition_checking;
              }
              break;
          }
          case incoming_inert_transition_checking:
          {
//std::cerr << "U_incoming_inert_transition_checking\n";
            // Algorithm 3, line 3.8, left.
            assert(current_U_incoming_transition_iterator < current_U_incoming_transition_iterator_end);
            assert(m_aut.is_tau(m_aut.apply_hidden_label_map(current_U_incoming_transition_iterator->label())));
            // Check one incoming transition.
            // Algorithm 3, line 3.12, left.
            mCRL2complexity_gj(&m_transitions[std::distance(m_aut.get_transitions().begin(), current_U_incoming_transition_iterator)], add_work(check_complexity::simple_splitB_U__handle_transition_to_U_state, 1), *this);
            current_U_outgoing_state = current_U_incoming_transition_iterator->from();
            assert(m_states[current_U_incoming_transition_iterator->to()].block == B);
            current_U_incoming_transition_iterator++;
//std::cerr << "FROM " << current_U_outgoing_state << "\n";
            if (m_states[current_U_outgoing_state].block==B && !(m_preserve_divergence && current_U_outgoing_state == std::prev(current_U_incoming_transition_iterator)->to()))
            {
              assert(!m_U.find(current_U_outgoing_state));
              if (m_states[current_U_outgoing_state].counter != Rmarked)
              {
                if (m_states[current_U_outgoing_state].counter==undefined) // count(current_U_outgoing_state) is undefined;
                {
                  // Algorithm 3, line 3.13, left.
                  // Algorithm 3, line 3.15 and 3.18, left.
                  m_states[current_U_outgoing_state].counter=m_states[current_U_outgoing_state].no_of_outgoing_inert_transitions-1;
                  m_U_counter_reset_vector.push_back(current_U_outgoing_state);
                }
                else
                {
                  // Algorithm 3, line 3.18, left.
                  assert(std::find(m_U_counter_reset_vector.begin(), m_U_counter_reset_vector.end(), current_U_outgoing_state) != m_U_counter_reset_vector.end());
                  assert(m_states[current_U_outgoing_state].counter>0);
                  m_states[current_U_outgoing_state].counter--;
                }
//std::cerr << "COUNTER " << m_states[current_U_outgoing_state].counter << "\n";
                // Algorithm 3, line 3.19, left.
                if (m_states[current_U_outgoing_state].counter==0)
                {
                  if (!(R_status==state_checking || R_status==incoming_inert_transition_checking || R_status==aborted_after_initialisation))
                  {
                    // Start searching for an outgoing transition with action a to constellation C.
                    current_U_outgoing_transition_iterator = m_states[current_U_outgoing_state].start_outgoing_transitions;
                    assert(m_outgoing_transitions.begin() <= current_U_outgoing_transition_iterator);
                    assert(current_U_outgoing_transition_iterator < m_outgoing_transitions.end());
                    current_U_outgoing_transition_iterator_end = (current_U_outgoing_state+1 >= m_states.size() ? m_outgoing_transitions.end() : m_states[current_U_outgoing_state+1].start_outgoing_transitions);
                    assert(current_U_outgoing_transition_iterator < current_U_outgoing_transition_iterator_end);
                    assert(m_aut.get_transitions()[m_BLC_transitions[current_U_outgoing_transition_iterator->transition]].from() == current_U_outgoing_state);
                    U_status=outgoing_action_constellation_check;
                    break;
                  }
                  else
                  {
                    // The state can be added to U_todo immediately.
                    #ifndef NDEBUG
                      // check that the state has no transition in the splitter
                      for (outgoing_transitions_it out_it = m_states[current_U_outgoing_state].start_outgoing_transitions;
                           out_it < m_outgoing_transitions.end() && (current_U_outgoing_state+1 >= m_states.size() || out_it < m_states[current_U_outgoing_state+1].start_outgoing_transitions); ++out_it)
                      {
                        assert(m_outgoing_transitions.begin() <= out_it);
                        assert(out_it < m_outgoing_transitions.end());
                        assert(0 <= out_it->transition);
                        assert(out_it->transition < m_aut.num_transitions());
                        assert(0 <= m_BLC_transitions[out_it->transition]);
                        assert(m_BLC_transitions[out_it->transition] < m_aut.num_transitions());
                        assert(m_transitions[m_BLC_transitions[out_it->transition]].ref_outgoing_transitions == out_it);
                        const transition& t = m_aut.get_transitions()[/*initialisation?current_U_outgoing_transition_iterator->transition:*/m_BLC_transitions[out_it->transition]];
                        assert(t.from() == current_U_outgoing_state);
                        if (a == label_or_divergence(t) && C == m_blocks[m_states[t.to()].block].constellation)
                        {
                          // The transition is in the splitter, so it must be in the part of the splitter that is disregarded.
// std::cerr << "State " << current_U_outgoing_state << " has a transition in the splitter, namely " << m_transitions[m_BLC_transitions[out_it->transition]].debug_id_short(*this) << '\n';
                          assert(m_BLC_transitions.begin() + out_it->transition >= splitter_end_unmarked_BLC);
                          assert(points_into_BLC_set(m_BLC_transitions.begin() + out_it->transition, *splitter));
                        }
                      }
                    #endif
                    m_U.add_todo(current_U_outgoing_state);
                    // Algorithm 3, line 3.10 and line 3.11 left.
                    if (m_U.size() > max_U_nonbottom_size)
                    {
                      assert(aborted != R_status); assert(aborted_after_initialisation != R_status);
                      U_status=aborted_after_initialisation;
                      break;
                    }
                  }
                }
              }
              else assert(m_R.find(current_U_outgoing_state));
            }
            if (current_U_incoming_transition_iterator == current_U_incoming_transition_iterator_end ||
                !m_aut.is_tau(m_aut.apply_hidden_label_map(current_U_incoming_transition_iterator->label())))
            {
              U_status = state_checking;
            }
            break;
          }
          case outgoing_action_constellation_check:
          {
//std::cerr << "U_outgoing_action_constellation_check\n";
            assert(current_U_outgoing_transition_iterator != current_U_outgoing_transition_iterator_end);
            // assert(splitter_end_unmarked_BLC == old value of splitter->start_marked_BLC); -- can no longer be checked
            #ifdef CHECK_COMPLEXITY_GJ
              mCRL2complexity_gj((&m_transitions[/*initialisation?current_U_outgoing_transition_iterator->transition:*/m_BLC_transitions[current_U_outgoing_transition_iterator->transition]]), add_work(check_complexity::simple_splitB_U__handle_transition_from_potential_U_state, 1), *this);
              // This is one step in the coroutine, so we should assign the work to exactly one transition.
              // But to make sure, we also mark the other transitions that we skipped in the optimisation.
              for (outgoing_transitions_it out_it = current_U_outgoing_transition_iterator; out_it < current_U_outgoing_transition_iterator->start_same_saC; )
              {
                ++out_it;
                mCRL2complexity_gj(&m_transitions[/*initialisation?out_it->transition:*/m_BLC_transitions[out_it->transition]], add_work_notemporary(check_complexity::simple_splitB_U__handle_transition_from_potential_U_state, 1), *this);
              }
            #endif
            const transition& t_local=m_aut.get_transitions()
                          [/*initialisation
                               ?current_U_outgoing_transition_iterator->transition
                               :*/m_BLC_transitions[current_U_outgoing_transition_iterator->transition]];
            current_U_outgoing_transition_iterator=current_U_outgoing_transition_iterator->start_same_saC; // This is an optimisation.
            ++current_U_outgoing_transition_iterator;

            assert(t_local.from() == current_U_outgoing_state);
            assert(m_branching);
            const label_index t_local_label = label_or_divergence(t_local);
            if (m_blocks[m_states[t_local.to()].block].constellation==C &&
                t_local_label == a &&
                (// transition t_local is not constellation-inert:
                 m_blocks[B].constellation != C || m_aut.num_action_labels() == t_local_label || !m_aut.is_tau(t_local_label)))
            {
                // This state must be blocked.
            }
            else if (current_U_outgoing_transition_iterator == current_U_outgoing_transition_iterator_end)
            {
              // assert(U.find(current_U_outgoing_state)==U.end());
              assert(!m_U.find(current_U_outgoing_state));
//std::cerr << "U_todo4 insert: " << current_U_outgoing_state << "   " << m_U.size() << "    " << number_of_states_in_block(B) << "\n";
              m_U.add_todo(current_U_outgoing_state);
              // Algorithm 3, line 3.10 and line 3.11 left.
              if (m_U.size() > max_U_nonbottom_size)
              {
                assert(aborted != R_status); assert(aborted_after_initialisation != R_status);
                U_status=aborted_after_initialisation;
                break;
              }
            }
            else  break;

            U_status = incoming_inert_transition_checking;
            if (current_U_incoming_transition_iterator == current_U_incoming_transition_iterator_end ||
                !m_aut.is_tau(m_aut.apply_hidden_label_map(current_U_incoming_transition_iterator->label())))
            {
              U_status = state_checking;
            }
            break;
          }
          default:
          {
            assert(aborted_after_initialisation == U_status);
            break;
          }
        }
        if (state_checking == U_status && current_U_incoming_bottom_state_iterator == m_blocks[B].start_non_bottom_states && m_U.todo_is_empty())
        {
          // split_block B into U and B\U.
//std::cerr << "U_todo empty:\n";
          assert(0 < std::distance(first_unmarked_bottom_state, m_blocks[B].start_non_bottom_states) + m_U.size());
          assert(std::distance(first_unmarked_bottom_state, m_blocks[B].start_non_bottom_states) + m_U.size() <= number_of_states_in_block(B)/2);
          clear_state_counters();
          m_R.clear();
          split_block_B_into_R_and_BminR(B, first_unmarked_bottom_state, m_blocks[B].start_non_bottom_states, m_U);
          m_U.clear();
          return B;
        }
      }
      assert(0);
    }

    // Make this transition non-inert.
    // The transition must go from one block to another but it cannot be constellation-inert yet.
    void make_transition_non_inert(const transition& t)
    {
// std::cerr << "Transition " << t.from() << " -" << m_aut.action_label(t.label()) << "-> " << t.to() << " becomes non-inert.\n";
      assert(is_inert_during_init(t));
      assert(m_states[t.to()].block!=m_states[t.from()].block);
      assert(m_blocks[m_states[t.to()].block].constellation == m_blocks[m_states[t.from()].block].constellation);
      m_states[t.from()].no_of_outgoing_inert_transitions--;
    }

    void change_non_bottom_state_to_bottom_state(const state_index si)
    {
      // Move this former non-bottom state to the bottom states.
      // The block of si is not yet inserted into the set of blocks with new bottom states.
      block_index bi = m_states[si].block;
// std::cerr << m_states[si].debug_id(*this) << " becomes a new bottom state of " << m_blocks[bi].debug_id(*this) << ".\n";
      assert(0 == m_states[si].no_of_outgoing_inert_transitions);
      assert(!m_blocks[bi].contains_new_bottom_states);
      swap_states_in_states_in_block(m_states[si].ref_states_in_blocks, m_blocks[bi].start_non_bottom_states);
      m_blocks[bi].start_non_bottom_states++;
    }

    // Split block B in R, being the inert-tau transitive closure of M contains
    // states that must be in block, and M\R. M_nonmarked, minus those in unmarked_blocker, are those in the other block.
    // The splitting is done in time O(min(|R|,|B\R|). Returns the block index of the R-block.
    block_index splitB(linked_list<BLC_indicators>::iterator splitter,
                       std::vector<state_index>::iterator first_unmarked_bottom_state,
                       const BLC_list_iterator splitter_end_unmarked_BLC /* = splitter.start_marked_BLC -- but this default argument is not allowed */,
                       std::function<void(const block_index, const block_index, const transition_index, const transition_index)>
                                                 update_block_label_to_cotransition = [](const block_index, const block_index, const transition_index, const transition_index){},
                       // const bool initialisation=false,
                       //std::function<void(const transition_index, const transition_index, const block_index)> process_transition=
                       //                                 [](const transition_index, const transition_index, const block_index){},
                       const bool split_off_new_bottom_states = true)
    {
// std::cerr << "splitB(splitter = " << splitter->debug_id(*this) << ", first_unmarked_bottom_state = " << m_states[*first_unmarked_bottom_state].debug_id(*this) << ", splitter_end_unmarked_BLC = "
//<< (split_off_new_bottom_states && splitter_end_unmarked_BLC == splitter->start_marked_BLC ? "start_marked_BLC" : (splitter_end_unmarked_BLC == splitter->start_same_BLC ? "start_same_BLC" : "?")) << ", ..., split_off_new_bottom_states = " << split_off_new_bottom_states << ")\n";
      const block_index B = m_states[m_aut.get_transitions()[*splitter->start_same_BLC].from()].block;
//std::cerr << "Marked bottom states:"; for (std::vector<state_index>::iterator it = m_blocks[B].start_bottom_states; it < first_unmarked_bottom_state; ++it) { std::cerr << ' ' << *it; }
//std::cerr << "\nUnmarked bottom states:"; for (std::vector<state_index>::iterator it = first_unmarked_bottom_state; it < m_blocks[B].start_non_bottom_states; ++it) { std::cerr << ' ' << *it; } std::cerr << "\nAdditionally, " << m_R.size() << " non-bottom states have been marked.\n";
      assert(!has_marked_transitions(*splitter));
      if (1 >= number_of_states_in_block(B))
      {
        mCRL2log(log::debug) << "Trying to split up singleton block " << B << '\n';
        clear_state_counters();
        return null_block;
      }
      block_index R_block = simple_splitB(B, splitter, first_unmarked_bottom_state, splitter_end_unmarked_BLC);
      const block_index bi = m_blocks.size() - 1;
// std::cerr << "Split block of size " << number_of_states_in_block(B) + number_of_states_in_block(bi) << " taking away " << number_of_states_in_block(bi) << " states\n";
      assert(bi == R_block || B == R_block);
      assert(number_of_states_in_block(B) >= number_of_states_in_block(bi));

      // Because we visit all states of block bi and almost all their incoming and outgoing transitions,
      // we subsume all this bookkeeping in a single block counter:
      mCRL2complexity_gj(&m_blocks[bi], add_work(check_complexity::splitB__update_BLC_of_smaller_subblock, check_complexity::log_n - check_complexity::ilog2(number_of_states_in_block(bi))), *this);
      // Update the LBC_list, and bottom states, and invariant on inert transitions.
      const std::vector<state_index>::iterator start_new_bottom_states = m_blocks[R_block].start_non_bottom_states;
      linked_list<BLC_indicators>::iterator R_to_U_tau_splitter = m_blocks[R_block].block_to_constellation.end();
      bool skip_transitions_in_splitter = false;
      assert(m_blocks[bi].block_to_constellation.begin() == m_blocks[bi].block_to_constellation.end());
      if (split_off_new_bottom_states && bi == R_block)
      {
        // move splitter as a whole from its current list to the new list where it will be part of:
        // (This only works if the whole splitter BLC set has moved, i.e. only if split_off_new_bottom_states == true.)
        m_blocks[R_block].block_to_constellation.splice_to_after(m_blocks[R_block].block_to_constellation.end(),
                                m_blocks[B].block_to_constellation, splitter);
        skip_transitions_in_splitter = true;
      }
      // Recall new LBC positions.
      for(typename std::vector<state_index>::iterator ssi=m_blocks[bi].start_bottom_states;
                                                      ssi!=m_blocks[bi].end_states;
                                                      ++ssi)
      {
        const state_index si=*ssi;
        state_type_gj& s= m_states[si];
        assert(s.ref_states_in_blocks == ssi);
        // mCRL2complexity_gj(s, add_work(..., max_bi_counter), *this);
            // is subsumed in the above call
        s.block=bi;

        // if (!initialisation)  // update the BLC_lists.
        {
          const outgoing_transitions_it end_it=((si+1)==m_states.size())?m_outgoing_transitions.end():m_states[si+1].start_outgoing_transitions;
          for(outgoing_transitions_it ti=s.start_outgoing_transitions; ti!=end_it; ti++)
          {
            assert(m_aut.get_transitions()[m_BLC_transitions[ti->transition]].from() == si);
            // mCRL2complexity_gj(&m_transitions[ti->transition], add_work(..., max_bi_counter), *this);
                // is subsumed in the above call
            // transition_index old_remaining_transition=update_the_doubly_linked_list_LBC_new_block(B, bi, ti->transition, new_BLC_list_entries);
            transition_index old_remaining_transition;
            if (!skip_transitions_in_splitter || m_transitions[m_BLC_transitions[ti->transition]].transitions_per_block_to_constellation != splitter)
            {
              old_remaining_transition=update_the_doubly_linked_list_LBC_new_block(B, bi, m_BLC_transitions[ti->transition]);
            }
            else
            {
              // No need to move the splitter transitions individually.
              old_remaining_transition = null_transition;
            }
            //process_transition(m_BLC_transitions[ti->transition], old_remaining_transition, B);
            update_block_label_to_cotransition(B, bi, m_BLC_transitions[ti->transition], old_remaining_transition);
// TODO: check what happens if the U-subblock is the new block;
// then perhaps update_block_label_to_cotransition() is not needed? Or is it still needed?
          }
        }

        // Situation below is only relevant if the new block contains the R-states:
        if (bi == R_block && 0 < s.no_of_outgoing_inert_transitions)
        {
          assert(ssi >= m_blocks[R_block].start_non_bottom_states);
          // si is a non_bottom_state in the smallest block containing M..
          bool non_bottom_state_becomes_bottom_state = true;

          const outgoing_transitions_it end_it=((*ssi)+1>=m_states.size())?m_outgoing_transitions.end():m_states[(*ssi)+1].start_outgoing_transitions;
          for(outgoing_transitions_it ti=s.start_outgoing_transitions; ti!=end_it; ti++)
          {
            // mCRL2complexity_gj(&m_transitions[ti->transition], add_work(..., max_bi_counter), *this);
                // is subsumed in the above call
            const transition& t=m_aut.get_transitions()[/*initialisation?ti->transition:*/m_BLC_transitions[ti->transition]];
            assert(t.from() == *ssi);
            if (is_inert_during_init_if_branching(t))
            {
              if (m_states[t.to()].block==B)
              {
                // This is a transition that has become non-inert.
                make_transition_non_inert(t);
                const linked_list<BLC_indicators>::iterator new_splitter = m_transitions[m_BLC_transitions[ti->transition]].transitions_per_block_to_constellation;
                if (R_to_U_tau_splitter == m_blocks[R_block].block_to_constellation.end())
                {
                  assert(!has_marked_transitions(*new_splitter));
                }
                else
                {
                  assert(R_to_U_tau_splitter == new_splitter);
                  assert(has_marked_transitions(*new_splitter));
                }
                R_to_U_tau_splitter = new_splitter;
                // immediately mark this transition, in case we get new bottom states:
                mark_BLC_transition(ti);

                // The LBC-list of this transition will be updated below, as it is now non-inert.
// David suggests: But the transition is not yet constellation-inert so nothing should change in BLC.
              }
              else if (m_states[t.to()].block==R_block)
              {
                non_bottom_state_becomes_bottom_state=false; // There is an outgoing inert tau. State remains non-bottom.
              }
            }
          }
          if (non_bottom_state_becomes_bottom_state)
          {
            // The state at si has become a bottom_state.
              // Note that the call below damages the value of *ssi. Here it is not anymore equal to si.
            assert(m_blocks[R_block].block_to_constellation.end() != R_to_U_tau_splitter);
            change_non_bottom_state_to_bottom_state(si);
          }
        }

        // Investigate the incoming formerly inert tau transitions.
        if (bi != R_block && m_blocks[R_block].start_non_bottom_states < m_blocks[R_block].end_states)
        {
          const std::vector<transition>::iterator it_end = si+1>=m_states.size() ? m_aut.get_transitions().end() : m_states[si+1].start_incoming_transitions;
          for(std::vector<transition>::iterator it=s.start_incoming_transitions;
                        it!=it_end; it++)
          {
            const transition& t=*it;
            // mCRL2complexity_gj(&m_transitions[std::distance(m_aut.get_transitions().begin(), it)], add_work(..., max_bi_counter), *this);
                // is subsumed in the above call
            assert(t.to() == si);
            if (!m_aut.is_tau(m_aut.apply_hidden_label_map(t.label())))
            {
              break; // All tau transitions have been investigated.
            }

            const state_index from=t.from();
            if (m_states[from].block==R_block && !(m_preserve_divergence && from == si))
            {
              // This transition did become non-inert.
              make_transition_non_inert(t);
              const linked_list<BLC_indicators>::iterator new_splitter = m_transitions[std::distance(m_aut.get_transitions().begin(), it)].transitions_per_block_to_constellation;
              if (R_to_U_tau_splitter == m_blocks[R_block].block_to_constellation.end())
              {
                assert(!has_marked_transitions(*new_splitter));
              }
              else
              {
                assert(R_to_U_tau_splitter == new_splitter);
                assert(has_marked_transitions(*new_splitter));
              }
              R_to_U_tau_splitter = new_splitter;
              // immediately mark this transition, in case we get new bottom states:
              mark_BLC_transition(m_transitions[std::distance(m_aut.get_transitions().begin(), it)].ref_outgoing_transitions);

              // Check whether from is a new bottom state.
              if (m_states[from].no_of_outgoing_inert_transitions==0)
              {
                // This state has no more outgoing inert transitions. It becomes a bottom state.
                change_non_bottom_state_to_bottom_state(from);
              }
            }
          }
        }
      }
      assert(m_blocks[R_block].start_bottom_states < m_blocks[R_block].start_non_bottom_states);

      #ifdef CHECK_COMPLEXITY_GJ
        unsigned const max_block(check_complexity::log_n - check_complexity::ilog2(number_of_states_in_block(bi)));
        if (bi == R_block)
        {
          // account for the work in R
          for (typename std::vector<state_index>::iterator s = m_blocks[bi].start_bottom_states ;
                              s != m_blocks[bi].end_states ; ++s)
          {
            mCRL2complexity_gj(&m_states[*s], finalise_work(check_complexity::simple_splitB_R__find_predecessors, check_complexity::simple_splitB__find_predecessors_of_R_or_U_state, max_block), *this);
            // incoming tau-transitions of s
            std::vector<transition>::iterator ti_end = *s + 1 >= m_states.size() ? m_aut.get_transitions().end() : m_states[*s+1].start_incoming_transitions;
            for (std::vector<transition>::iterator ti = m_states[*s].start_incoming_transitions; ti != ti_end; ++ti)
            {
              if (!m_aut.is_tau(m_aut.apply_hidden_label_map(ti->label())))  break;
              mCRL2complexity_gj(&m_transitions[std::distance(m_aut.get_transitions().begin(), ti)], finalise_work(check_complexity::simple_splitB_R__handle_transition_to_R_state, check_complexity::simple_splitB__handle_transition_to_R_or_U_state, max_block), *this);
            }
            // outgoing transitions of s
            for (outgoing_transitions_it ti = (*s+1 >= m_states.size() ? m_outgoing_transitions.end() : m_states[*s+1].start_outgoing_transitions);
              ti != m_states[*s].start_outgoing_transitions; )
            {
              ti--;
              mCRL2complexity_gj(&m_transitions[/*initialisation?ti->transition:*/m_BLC_transitions[ti->transition]], finalise_work(check_complexity::simple_splitB_R__handle_transition_from_R_state, check_complexity::simple_splitB__handle_transition_from_R_or_U_state, max_block), *this);
              // We also need to cancel the work on outgoing transitions of U-state candidates that turned out to be new bottom states:
              mCRL2complexity_gj(&m_transitions[/*initialisation?ti->transition:*/m_BLC_transitions[ti->transition]], cancel_work(check_complexity::simple_splitB_U__handle_transition_from_potential_U_state), *this);
            }
          }
          // ensure not too much work has been done on U
          for (typename std::vector<state_index>::iterator s = m_blocks[B].start_bottom_states ;
                              s != m_blocks[B].end_states ; ++s )
          {
            mCRL2complexity_gj(&m_states[*s], cancel_work(check_complexity::simple_splitB_U__find_bottom_state), *this);
            mCRL2complexity_gj(&m_states[*s], cancel_work(check_complexity::simple_splitB_U__find_predecessors), *this);
            // incoming tau-transitions of s
            std::vector<transition>::iterator ti_end = *s + 1 >= m_states.size() ? m_aut.get_transitions().end() : m_states[*s+1].start_incoming_transitions;
            for (std::vector<transition>::iterator ti = m_states[*s].start_incoming_transitions; ti != ti_end; ++ti)
            {
              if (!m_aut.is_tau(m_aut.apply_hidden_label_map(ti->label())))  break;
              mCRL2complexity_gj(&m_transitions[std::distance(m_aut.get_transitions().begin(), ti)], cancel_work(check_complexity::simple_splitB_U__handle_transition_to_U_state), *this);
            }
            // outgoing transitions of s
            for (outgoing_transitions_it ti = (*s + 1 >= m_states.size() ? m_outgoing_transitions.end() : m_states[*s+1].start_outgoing_transitions);
              ti != m_states[*s].start_outgoing_transitions; )
            {
              ti--;
              mCRL2complexity_gj(&m_transitions[/*initialisation?ti->transition:*/m_BLC_transitions[ti->transition]], cancel_work(check_complexity::simple_splitB_U__handle_transition_from_potential_U_state), *this);
            }
          }
        }
        else
        {
          // account for the work in U
          for (typename std::vector<state_index>::iterator s = m_blocks[bi].start_bottom_states ;
                              s != m_blocks[bi].end_states ; ++s)
          {
            mCRL2complexity_gj(&m_states[*s], finalise_work(check_complexity::simple_splitB_U__find_bottom_state, check_complexity::simple_splitB__find_bottom_state, max_block), *this);
            mCRL2complexity_gj(&m_states[*s], finalise_work(check_complexity::simple_splitB_U__find_predecessors, check_complexity::simple_splitB__find_predecessors_of_R_or_U_state, max_block), *this);
            // incoming tau-transitions of s
            std::vector<transition>::iterator ti_end = *s + 1 >= m_states.size() ? m_aut.get_transitions().end() : m_states[*s+1].start_incoming_transitions;
            for (std::vector<transition>::iterator ti = m_states[*s].start_incoming_transitions; ti != ti_end; ++ti)
            {
              if (!m_aut.is_tau(m_aut.apply_hidden_label_map(ti->label())))  break;
              mCRL2complexity_gj(&m_transitions[std::distance(m_aut.get_transitions().begin(), ti)], finalise_work(check_complexity::simple_splitB_U__handle_transition_to_U_state, check_complexity::simple_splitB__handle_transition_to_R_or_U_state, max_block), *this);
            }
            // outgoing transitions of s
            for (outgoing_transitions_it ti = *s + 1 >= m_states.size() ? m_outgoing_transitions.end() : m_states[*s+1].start_outgoing_transitions;
              ti != m_states[*s].start_outgoing_transitions; )
            {
              ti--;
              mCRL2complexity_gj(&m_transitions[/*initialisation?ti->transition:*/m_BLC_transitions[ti->transition]], finalise_work(check_complexity::simple_splitB_U__handle_transition_from_potential_U_state, check_complexity::simple_splitB__handle_transition_from_R_or_U_state, max_block), *this);
            }
          }
          // ensure not too much work has been done on R
          for (typename std::vector<state_index>::iterator s = m_blocks[B].start_bottom_states ;
                              s != m_blocks[B].end_states ; ++s )
          {
            mCRL2complexity_gj(&m_states[*s], cancel_work(check_complexity::simple_splitB_R__find_predecessors), *this);
            // incoming tau-transitions of s
            std::vector<transition>::iterator ti_end = *s + 1 >= m_states.size() ? m_aut.get_transitions().end() : m_states[*s+1].start_incoming_transitions;
            for (std::vector<transition>::iterator ti = m_states[*s].start_incoming_transitions; ti != ti_end; ++ti)
            {
              if (!m_aut.is_tau(m_aut.apply_hidden_label_map(ti->label())))  break;
              mCRL2complexity_gj(&m_transitions[std::distance(m_aut.get_transitions().begin(), ti)], cancel_work(check_complexity::simple_splitB_R__handle_transition_to_R_state), *this);
            }
            // outgoing transitions of s
            for (outgoing_transitions_it ti = (*s + 1 >= m_states.size() ? m_outgoing_transitions.end() : m_states[*s+1].start_outgoing_transitions);
              ti != m_states[*s].start_outgoing_transitions; )
            {
              ti--;
              mCRL2complexity_gj(&m_transitions[/*initialisation?ti->transition:*/m_BLC_transitions[ti->transition]], cancel_work(check_complexity::simple_splitB_R__handle_transition_from_R_state), *this);
              // We also need to move the work on outgoing transitions of U-state candidates that turned out to be new bottom states:
              mCRL2complexity_gj(&m_transitions[/*initialisation?ti->transition:*/m_BLC_transitions[ti->transition]], finalise_work(check_complexity::simple_splitB_U__handle_transition_from_potential_U_state, check_complexity::simple_splitB__test_outgoing_transitions_found_new_bottom_state,
                      0 == m_states[*s].no_of_outgoing_inert_transitions ? 1 : 0), *this);
            }
          }
        }
        check_complexity::check_temporary_work();
      #endif // ifdef CHECK_COMPLEXITY_GJ

      if (split_off_new_bottom_states && start_new_bottom_states < m_blocks[R_block].start_non_bottom_states)
      {
        // There are new bottom states, and we have to separate them immediately from the old bottom states.
        // This will make sure that the next call of stabilizeB() will not encounter blocks with both old and new bottom states.
        assert(m_blocks[R_block].block_to_constellation.end() != R_to_U_tau_splitter);
        assert(has_marked_transitions(*R_to_U_tau_splitter));

        block_index new_bottom_block;
        if (start_new_bottom_states == m_blocks[R_block].start_bottom_states)
        {
          // all bottom states in this block are new
          // unmark_all_transitions(R_to_U_tau_splitter); -- will be done below
          new_bottom_block = R_block;
          R_block = null_block;
        }
        else
        {
          // Some tau-transitions from R to U may come out of states that are not (yet) new bottom states.
          // Therefore we still have to go through the movement of constructing m_R:
          const BLC_list_iterator splitter_end_unmarked_BLC = R_to_U_tau_splitter->start_same_BLC;
          first_unmarked_bottom_state = not_all_bottom_states_are_touched(R_to_U_tau_splitter
                      #ifndef NDEBUG
                        , splitter_end_unmarked_BLC
                      #endif
                      );
          assert(std::distance(start_new_bottom_states, m_blocks[R_block].start_non_bottom_states) ==
                 std::distance(m_blocks[R_block].start_bottom_states, first_unmarked_bottom_state));
          assert(m_blocks[R_block].start_bottom_states < first_unmarked_bottom_state);
          assert(first_unmarked_bottom_state < m_blocks[R_block].start_non_bottom_states);
          new_bottom_block = splitB(R_to_U_tau_splitter, first_unmarked_bottom_state,
                                    splitter_end_unmarked_BLC,
                                    update_block_label_to_cotransition,
                                    //process_transition,
                                    false);
          if (R_block == new_bottom_block)
          {
            R_block = m_blocks.size() - 1;
          }
        }
        assert(0 <= new_bottom_block);  assert(new_bottom_block < m_blocks.size());
        assert(!m_blocks[new_bottom_block].contains_new_bottom_states);
// std::cerr << "new_bottom_block = " << new_bottom_block << ", R_block = " << static_cast<std::make_signed<block_index>::type>(R_block) << '\n';
        m_blocks[new_bottom_block].contains_new_bottom_states = true;
        m_blocks_with_new_bottom_states.push_back(new_bottom_block);

        // unmark all transitions out of the new bottom state block:
        // (This is in particular needed if the same block still needs to be
        // stabilized under other action labels.)
        for (BLC_indicators& ind: m_blocks[new_bottom_block].block_to_constellation)
        {
          // TODO: mCRL2complexity(&ind, add_work(..., 1), *this);
            // if there are transitions from new bottom states in this BLC set,
            // assign the work to these transitions;
            // otherwise, assign the work later.
          // if (has_marked_transitions(ind)) { assert(to_constellation(ind) == new_constellation || to_constellation(ind) == old_constellation); }
          unmark_all_transitions(ind);
        }

// David suggests: Instead of block_label_to_cotransition,
// one could add one pointer to the block_type.
// This pointer points at a BLC set (or perhaps at some transition).
// This BLC set is the co-splitter that has been found most recently.
// The pointer needs to be updated when a block is split.

// Possibly one can completely get away without an additional pointer,
// by stating: the second element of the list of BLC indicators should be the
// co-splitter. (The first element is the constellation-inert tau transitions.)
// This would require more thought about how to construct a new BLC list
// when a block is split; one probably would want to distinguish emplace_front,
// emplace_second and emplace_third or so. (emplace_second should test whether
// the first entry in the list of BLC indicators contains the constellation-inert
// transitions. If not, it actually changes its behaviour and continues like
// emplace_front.)

      }
      else if (m_blocks[R_block].block_to_constellation.end() != R_to_U_tau_splitter)
      {
        // unmark the transitions from R to U.
        // We do not need to stabilize the (formerly inert) tau-transitions from R to U,
        // because they are still in the same constellation.
        unmark_all_transitions(*R_to_U_tau_splitter);
      }

      return R_block;
    }

    void accumulate_entries_into_not_investigated(std::vector<label_count_sum_triple>& action_counter,
                            const std::vector<block_index>& todo_stack)
    {
      transition_index sum=0;
      for(block_index index: todo_stack)
      {
        // The caller has to account for this work.
        // Typically, it can be accounted for because todo_stack has been filled
        // by actions that are counted.
        action_counter[index].not_investigated=sum;
        sum=sum+action_counter[index].label_counter;
      }
    }

    void accumulate_entries(std::vector<transition_index>& counter)
    {
      transition_index sum=0;
      for(transition_index& index: counter)
      {
        transition_index n=index;
        index=sum;
        sum=sum+n;
      }
    }

    void reset_entries(std::vector<label_count_sum_triple>& action_counter,
                       std::vector<block_index>& todo_stack)
    {
      for(block_index index: todo_stack)
      {
        // To account for this work, we need to ensure that todo_stack has only
        // been filled by actions that have been accounted for already.
        // it is not necessary to reset the cumulative_label_counter;
        action_counter[index].label_counter=0;
      }
      todo_stack.clear();
    }

    transition_index accumulate_entries(std::vector<transition_index>& action_counter,
                                   const std::vector<label_index>& todo_stack) const
    {
      transition_index sum=0;
      for(label_index index: todo_stack)
      {
        transition_index n=sum;
        sum=sum+action_counter[index];
        action_counter[index]=n;
      }
      return sum;
    }

/*    // Group the elements from begin up to end, using a range from [0,number of blocks),
    // where each transition pinpointed by the iterator has as value its source block.

    void group_in_situ(const std::vector<transition_index>::iterator begin,
                       const std::vector<transition_index>::iterator end,
                       std::vector<block_index>& todo_stack,
                       std::vector<label_count_sum_triple>& value_counter)
    {
      // Initialise the action counter.
      // The work in reset_entries() can be subsumed under mCRL2complexity calls
      // that have been issued in earlier executions of group_in_situ().
      reset_entries(value_counter, todo_stack);
      #ifdef CHECK_COMPLEXITY_GJ
        const constellation_index new_constellation = m_constellations.size()-1;
        const unsigned max_C = check_complexity::log_n - check_complexity::ilog2(number_of_states_in_constellation(new_constellation));
      #endif
      for(std::vector<transition_index>::iterator i=begin; i!=end; ++i)
      {
        const transition& t = m_aut.get_transitions()[*i];
        #ifdef CHECK_COMPLEXITY_GJ
          assert(m_blocks[m_states[t.to()].block].constellation == new_constellation);
          mCRL2complexity_gj(&m_transitions[*i], add_work(check_complexity::group_in_situ__count_transitions_per_block, max_C), *this);
        #endif
        const block_index n=m_states[t.from()].block;

        // mCRL2complexity_gj(..., add_work(..., ...), *this);
        if (value_counter[n].label_counter==0)
        {
          todo_stack.push_back(n);
        }
        value_counter[n].label_counter++;
      }

      // The work in accumulate_entries_into_not_investigated() can be subsumed
      // under the above call to mCRL2complexity, because an entry in todo_stack
      // is only made if there is at least one transition from that block.
      accumulate_entries_into_not_investigated(value_counter, todo_stack);

      std::vector<block_index>::iterator current_value=todo_stack.begin();
      for(std::vector<transition_index>::iterator i=begin; i!=end; )
      {
        mCRL2complexity_gj(&m_transitions[*i], add_work(check_complexity::group_in_situ__swap_transition, max_C), *this);
        block_index n=m_states[m_aut.get_transitions()[*i].from()].block;
        if (n==*current_value)
        {
          value_counter[n].label_counter--;
          value_counter[n].not_investigated++;
          ++i;
          while (assert(current_value!=todo_stack.end()), value_counter[n].label_counter==0)
          {
            #ifdef CHECK_COMPLEXITY_GJ
              // This work needs to be assigned to some transition from block n.
              // Just to make sure we assign it to all such transitions:
              std::vector<transition_index>::iterator work_i = i;
              assert(begin != work_i);
              --work_i;
              assert(m_states[m_aut.get_transitions()[*work_i].from()].block == n);
              do
              {
                mCRL2complexity_gj(&m_transitions[*work_i], add_work(check_complexity::group_in_situ__skip_to_next_block, max_C), *this);
              }
              while (begin != work_i && m_states[m_aut.get_transitions()[*--work_i].from()].block == n);
            #endif
            current_value++;
            if (current_value!=todo_stack.end())
            {
              n = *current_value;
              i=begin+value_counter[n].not_investigated; // Jump to the first non investigated action.
            }
            else
            {
              assert(i == end);
              break; // exit the while and the for loop.
            }
          }
        }
        else
        {
          // Find the first transition with a different label than t.label to swap with.
          std::vector<transition_index>::iterator new_position=begin+value_counter[n].not_investigated;
          while (m_states[m_aut.get_transitions()[*new_position].from()].block==n)
          {
            mCRL2complexity_gj(&m_transitions[*new_position], add_work(check_complexity::group_in_situ__swap_transition, max_C), *this);
            value_counter[n].not_investigated++;
            value_counter[n].label_counter--;
            new_position++;
            assert(new_position!=end);
          }
          assert(value_counter[n].label_counter>0);
          std::swap(*i, *new_position);
          value_counter[n].not_investigated++;
          value_counter[n].label_counter--;
        }
      }
    } */

//================================================= Create initial partition ========================================================
    void create_initial_partition()
    {
      mCRL2log(log::verbose) << "An O(m log n) "
           << (m_branching ? (m_preserve_divergence
                                         ? "divergence-preserving branching "
                                         : "branching ")
                         : "")
           << "bisimulation partitioner created for " << m_aut.num_states()
           << " states and " << m_aut.num_transitions() << " transitions (using the experimental algorithm GJ2024).\n";
      // Initialisation.
      #ifdef CHECK_COMPLEXITY_GJ
        check_complexity::init(2 * m_aut.num_states()); // we need ``2*'' because there is one additional call to splitB during initialisation
      #endif

      // Initialise m_incoming_(non-)inert-transitions, m_outgoing_transitions, and m_states[si].no_of_outgoing_transitions
      //group_transitions_on_label(m_aut.get_transitions(),
      //                          [](const transition& t){ return m_aut.apply_hidden_label_map(t.label()); },
      //                          m_aut.num_action_labels(), m_aut.tau_label_index()); // sort on label. Tau transitions come first.
      // group_transitions_on_label(m_aut.get_transitions(),
      //                           [](const transition& t){ return t.from(); },
      //                           m_aut.num_states(), 0); // sort on label. Tau transitions come first.
      // group_transitions_on_label_tgt(m_aut.get_transitions(), m_aut.num_action_labels(), m_aut.tau_label_index(), m_aut.num_states()); // sort on label. Tau transitions come first.
      // group_transitions_on_tgt(m_aut.get_transitions(), m_aut.num_action_labels(), m_aut.tau_label_index(), m_aut.num_states()); // sort on label. Tau transitions come first.
      // sort_transitions(m_aut.get_transitions(), lbl_tgt_src);
// David suggests: I think it is enough to sort according to tgt_lbl.
// JFG answers: Agreed. But I believe this will cost performance. For 1394-fin-vvlarge this saves 1 second to sort, but
//                      requires five more seconds to carry out the splitting. Apparently, there is benefit to have src together.
//                      This may have been measured on an older version of the code.
      sort_transitions(m_aut.get_transitions(), m_aut.hidden_label_set(), tgt_lbl_src); // THIS IS NOW ESSENTIAL.
      // sort_transitions(m_aut.get_transitions(), src_lbl_tgt);
      // sort_transitions(m_aut.get_transitions(), tgt_lbl);
      // sort_transitions(m_aut.get_transitions(), target);

      // Count the number of occurring action labels.
      assert((unsigned) m_preserve_divergence <= 1);

      mCRL2log(log::verbose) << "Start initialisation of the BLC list in the initialisation\n";
      {
        std::vector<label_index> todo_stack_actions;
        std::vector<transition_index> count_transitions_per_action(m_aut.num_action_labels() + (unsigned) m_preserve_divergence, 0);
// David suggests: The above allocation may take time up to O(|Act|).
// This is a place where the number of actions plays a role.
// JFG answers: I think we should accept that the algorithm has a complexity of .... + O(|Act|). Act can be assumed to be smaller
// than m, and n can be assumed to be bigger than 1. In that case O(|Act|) will be subsumed. 
        if (m_branching)
        {
          todo_stack_actions.push_back(m_aut.tau_label_index()); // ensure that inert transitions come first
          count_transitions_per_action[m_aut.tau_label_index()] = 1; // set the number of transitions to a nonzero value so it doesn't trigger todo_stack_actions.push_back(...) in the loop
        }
        for(transition_index ti=0; ti<m_aut.num_transitions(); ++ti)
        {
          const transition& t=m_aut.get_transitions()[ti];
          // mCRL2complexity_gj(&m_transitions[ti], add_work(..., 1), *this);
            // Because every transition is touched exactly once, we do not store a physical counter for this.
          const label_index label = label_or_divergence(t);
          transition_index& c=count_transitions_per_action[label];
          if (c==0)
          {
            todo_stack_actions.push_back(label);
          }
          c++;
        }
        if (m_branching)
        {
          assert(m_aut.is_tau(todo_stack_actions.front()));
          --count_transitions_per_action[m_aut.tau_label_index()];
        }
//std::cerr << "COUNT_TRANSITIONS PER ACT1    "; for(auto s: count_transitions_per_action){ std::cerr << s << "  "; } std::cerr << "\n";
        accumulate_entries(count_transitions_per_action, todo_stack_actions);
//std::cerr << "COUNT_TRANSITIONS PER ACT2    "; for(auto s: count_transitions_per_action){ std::cerr << s << "  "; } std::cerr << "\n";
        for(transition_index ti=0; ti<m_aut.num_transitions(); ++ti)
        {
          // mCRL2complexity_gj(&m_transitions[ti], add_work(..., 1), *this);
            // Because every transition is touched exactly once, we do not store a physical counter for this.
          const transition& t=m_aut.get_transitions()[ti];
          const label_index label = label_or_divergence(t);
          transition_index& c=count_transitions_per_action[label];
          assert(0 <= c); assert(c < m_aut.num_transitions());
          m_BLC_transitions[c]=ti;
          c++;
        }

        // create BLC_indicators for every action label:
        BLC_list_iterator start_index = m_BLC_transitions.begin();
        for(label_index a: todo_stack_actions)
        {
          // mCRL2complexity_gj(..., add_work(..., 1), *this);
              // not needed because the inner loop is always executed (except possibly for 1 iteration)
//std::cerr << "  Initialising m_BLC_transitions for action " << (m_aut.num_action_labels() == a ? "(tau-self-loops)" : pp(m_aut.action_label(a))) << '\n';
          const BLC_list_iterator end_index = m_BLC_transitions.begin() + count_transitions_per_action[a];
          assert(end_index <= m_BLC_transitions.end());
          if (start_index == end_index)
          {
            assert(m_branching); assert(m_aut.is_tau(a));
          }
          else
          {
            assert(start_index < end_index);
            // create a BLC_indicator and insert it into the list...
            typename linked_list<BLC_indicators>::iterator current_BLC = m_blocks[0].block_to_constellation.emplace_after(m_blocks[0].block_to_constellation.begin(), start_index, end_index);
            if (!is_inert_during_init(m_aut.get_transitions()[*start_index]))
            {
              current_BLC->start_marked_BLC = start_index; // mark all states in this BLC_indicator for the initial stabilization
            }
            do
            {
              // mCRL2complexity_gj(&m_transitions[*start_index], add_work(..., 1), *this);
                  // Because every transition is touched exactly once, we do not store a physical counter for this.
              m_transitions[*start_index].transitions_per_block_to_constellation = current_BLC;
            }
            while (++start_index < end_index);
          }
        }
        assert(start_index == m_BLC_transitions.end());

        // destroy and deallocate todo_stack_actions and count_transitions_per_action here.
      }

      // Group transitions per outgoing state.
      mCRL2log(log::verbose) << "Start setting outgoing transitions\n";
      {

        std::vector<transition_index> count_outgoing_transitions_per_state(m_aut.num_states(), 0);
        for(const transition& t: m_aut.get_transitions())
        {
          // mCRL2complexity_gj(&m_transitions[std::distance(&*m_aut.get_transitions.begin(), &t)], add_work(..., 1), *this);
            // Because every transition is touched exactly once, we do not store a physical counter for this.
          count_outgoing_transitions_per_state[t.from()]++;
          if (is_inert_during_init(t))
          {
            m_states[t.from()].no_of_outgoing_inert_transitions++;
          }
        }

        // We now set the outgoing transition per state pointer to the first non-inert transition.
        // The counters for outgoing transitions calculated above are reset to 0
        // and will later contain the number of transitions already stored.
        // Every time an inert transition is stored, the outgoing transition per state pointer is reduced by one.
        outgoing_transitions_it current_outgoing_transitions = m_outgoing_transitions.begin();

        // place transitions and set the pointers to incoming/outgoing transitions
        for (state_index s = 0; s < m_aut.num_states(); ++s)
        {
          // mCRL2complexity_gj(&m_states[s], add_work(..., 1), *this);
            // Because every state is touched exactly once, we do not store a physical counter for this.
          m_states[s].start_outgoing_transitions = current_outgoing_transitions + m_states[s].no_of_outgoing_inert_transitions;
          current_outgoing_transitions += count_outgoing_transitions_per_state[s];
          count_outgoing_transitions_per_state[s] = 0; // meaning of this counter changes to: number of outgoing transitions already stored
        }
        assert(current_outgoing_transitions == m_outgoing_transitions.end());

        mCRL2log(log::verbose) << "Moving incoming and outgoing transitions\n";

        for(BLC_list_iterator ti = m_BLC_transitions.begin(); ti < m_BLC_transitions.end(); ++ti)
        {
          // mCRL2complexity_gj(&m_transitions[*ti], add_work(..., 1), *this);
            // Because every transition is touched exactly once, we do not store a physical counter for this.
          const transition& t=m_aut.get_transitions()[*ti];
          if (is_inert_during_init(t))
          {
            m_transitions[*ti].ref_outgoing_transitions = --m_states[t.from()].start_outgoing_transitions;
          }
          else
          {
            m_transitions[*ti].ref_outgoing_transitions = m_states[t.from()].start_outgoing_transitions + count_outgoing_transitions_per_state[t.from()];
          }
          m_transitions[*ti].ref_outgoing_transitions->transition = std::distance(m_BLC_transitions.begin(), ti);
          ++count_outgoing_transitions_per_state[t.from()];
        }
        // destroy and deallocate count_outgoing_transitions_per_state here.
      }

      state_index current_state=null_state;
      assert(current_state + 1 == 0);
      // bool tau_transitions_passed=true;
      // TODO: This should be combined with another pass through all transitions.
      for(std::vector<transition>::iterator it=m_aut.get_transitions().begin(); it!=m_aut.get_transitions().end(); it++)
      {
        // mCRL2complexity_gj(&m_transitions[std::distance(m_aut.get_transitions().begin(), it)], add_work(..., 1), *this);
            // Because every transition is touched exactly once, we do not store a physical counter for this.
        const transition& t=*it;
        if (t.to()!=current_state)
        {
          for (state_index i=current_state+1; i<=t.to(); ++i)
          {
            // ensure that every state is visited at most once:
            mCRL2complexity_gj(&m_states[i], add_work(check_complexity::create_initial_partition__set_start_incoming_transitions, 1), *this);
//std::cerr << "SET start_incoming_transitions for state " << i << "\n";
            m_states[i].start_incoming_transitions=it;
          }
          current_state=t.to();
        }
      }
      for (state_index i=current_state+1; i<m_states.size(); ++i)
      {
        mCRL2complexity_gj(&m_states[i], add_work(check_complexity::create_initial_partition__set_start_incoming_transitions, 1), *this);
//std::cerr << "SET residual start_incoming_transitions for state " << i << "\n";
        m_states[i].start_incoming_transitions=m_aut.get_transitions().end();
      }

      // Set the start_same_saC fields in m_outgoing_transitions.
      outgoing_transitions_it it = m_outgoing_transitions.end();
      if (m_outgoing_transitions.begin() < it)
      {
        --it;
        const transition& t = m_aut.get_transitions()[m_BLC_transitions[it->transition]];
        state_index current_state = t.from();
        label_index current_label = label_or_divergence(t);
        outgoing_transitions_it current_end_same_saC = it;
        while (m_outgoing_transitions.begin() < it)
        {
          --it;
          // mCRL2complexity_gj(&m_transitions[it->transition], add_work(..., 1), *this);
              // Because every transition is touched exactly once, we do not store a physical counter for this.
          const transition& t = m_aut.get_transitions()[m_BLC_transitions[it->transition]];
          const label_index new_label = label_or_divergence(t);
          if (current_state == t.from() && current_label == new_label)
          {
            // We encounter a transition with the same saC. Let it refer to the end.
            it->start_same_saC = current_end_same_saC;
          }
          else
          {
            // We encounter a transition with a different saC.
            current_state = t.from();
            current_label = new_label;
            current_end_same_saC->start_same_saC = std::next(it);
            current_end_same_saC = it;
          }
        }
        current_end_same_saC->start_same_saC = m_outgoing_transitions.begin();
      }

//std::cerr << "Start filling states_in_blocks\n";
      m_states_in_blocks.resize(m_aut.num_states());
      typename std::vector<state_index>::iterator lower_i=m_states_in_blocks.begin(), upper_i=m_states_in_blocks.end();
      for (state_index i=0; i < m_aut.num_states(); ++i)
      {
        // mCRL2complexity_gj(&m_states[i], add_work(..., 1), *this);
            // Because every state is touched exactly once, we do not store a physical counter for this.
        if (0 < m_states[i].no_of_outgoing_inert_transitions)
        {
          --upper_i;
          *upper_i = i;
          m_states[i].ref_states_in_blocks = upper_i;
        }
        else
        {
          *lower_i = i;
          m_states[i].ref_states_in_blocks = lower_i;
          ++lower_i;
        }
      }
      assert(lower_i == upper_i);
      m_blocks[0].start_bottom_states=m_states_in_blocks.begin();
      m_blocks[0].start_non_bottom_states = lower_i;
      m_blocks[0].end_states=m_states_in_blocks.end();

      // The data structures have now been completely initialized.

      // print_data_structures("After initial reading before splitting in the initialisation");
      assert(check_data_structures("After initial reading before splitting in the initialisation"));

      // The initial partition has been constructed. Continue with the initialisation.
      mCRL2log(log::verbose) << "Start refining in the initialisation\n";

//std::cerr << "COUNT_STATES PER ACT     "; for(auto s: count_transitions_per_action){ std::cerr << s << "  "; } std::cerr << "\n";
//std::cerr << "STATES PER ACTION LABEL  "; for(transition_index ti: m_BLC_transitions){ std::cerr << ti << "  "; } std::cerr << "\n";
//std::cerr << "STATES PER ACTION LABELB "; for(transition_index ti: m_BLC_transitions){ std::cerr << m_states[m_aut.get_transitions()[ti].from()].block << "  "; } std::cerr << "\n";

      // In principle, we now stabilize for every BLC_indicator.
      // All transitions have already been marked in every BLC_indicator, so we can just run normal stabilisations.
      std::pair<BLC_list_iterator, BLC_list_iterator> stabilize_pair(m_BLC_transitions.begin(), m_BLC_transitions.end());

      #ifndef NDEBUG
        const std::vector<std::pair<BLC_list_iterator, BLC_list_iterator> > transition_array_for_check_stability(1, stabilize_pair);
      #endif
      while (stabilize_pair.first < stabilize_pair.second)
      {
        #ifndef NDEBUG
          check_stability("Initialisation loop", &transition_array_for_check_stability, &stabilize_pair);
        #endif
        linked_list<BLC_indicators>::iterator splitter = m_transitions[*std::prev(stabilize_pair.second)].transitions_per_block_to_constellation;
        // mCRL2complexity_gj(...);
            // not needed, as the splitter has marked transitions and we are allowed to visit each marked transition a fixed number of times.
        stabilize_pair.second = splitter->start_same_BLC;
        assert(stabilize_pair.first <= stabilize_pair.second);
        assert(points_into_BLC_set(stabilize_pair.second, *splitter));

        const transition& first_t = m_aut.get_transitions()[*splitter->start_same_BLC];
        if (m_blocks[m_states[first_t.from()].block].contains_new_bottom_states)
        {
          // This slice contains new bottom states and therefore should not be stabilized
          assert(!has_marked_transitions(*splitter));
// std::cerr << splitter->debug_id(*this) << " is skipped because it starts in a block with new bottom states\n";
        }
        else if (is_inert_during_init(first_t))
        {
          // This slice contains constellation-inert transitions;
          // we do not need to stabilize under it.
          assert(!has_marked_transitions(*splitter));
// std::cerr << splitter->debug_id(*this) << " is skipped because it contains constellation-inert transitions\n";
        }
        else
        {
// std::cerr << "Now stabilizing under " << splitter->debug_id(*this) << '\n';
          assert(splitter->start_marked_BLC == splitter->start_same_BLC);
          const block_index source_block = m_states[first_t.from()].block;
          const BLC_list_iterator splitter_end_unmarked_BLC = splitter->start_marked_BLC;
          const std::vector<state_index>::iterator first_unmarked_bottom_state = not_all_bottom_states_are_touched(splitter
                    #ifndef NDEBUG
                      , splitter_end_unmarked_BLC
                    #endif
                    );

          if (first_unmarked_bottom_state < m_blocks[source_block].start_non_bottom_states)
          {
            splitB(splitter, first_unmarked_bottom_state, splitter_end_unmarked_BLC);
          }
        }
      }

      // Algorithm 1, line 1.4 is implicitly done in the call to splitB above.

      // Algorithm 1, line 1.5.
      //print_data_structures("End initialisation");
      check_stability("End initialisation");
      mCRL2log(log::verbose) << "Start stabilizing in the initialisation\n";
      stabilizeB();
      assert(check_data_structures("End initialisation"));
    }

    // Algorithm 4. Stabilize the current partition with respect to the current constellation
    // given that the blocks in m_blocks_with_new_bottom_states do contain new bottom states.
    // Stabilisation is always called after initialisation, i.e., m_incoming_transitions[ti].transition refers
    // to a position in m_BLC_transitions, where the transition index of this transition can be found.

    void stabilizeB()
    {
      // 0. We make a special split: the new bottom states are separated from the old bottom states,
      //    by using the block-non-inert tau-transitions just mentioned as a splitter.
      // -- This special split already has been performed in splitB().

      while (!m_blocks_with_new_bottom_states.empty())
      {
        // mCRL2complexity_gj(all bottom states, add_work(..., 1), *this);
            // not necessary, as the inner loop is always executed

        // 1. Administration: all transitions out of new bottom states are marked.
        #ifdef CHECK_COMPLEXITY_GJ
          std::vector<std::pair<BLC_list_const_iterator, BLC_list_const_iterator> > work_to_assign_later;
        #endif
        // Qhat contains the slices of BLC transitions that still need stabilization
        std::vector<std::pair<BLC_list_iterator, BLC_list_iterator> > Qhat;
        for(const block_index bi: m_blocks_with_new_bottom_states)
        {
          // mCRL2complexity_gj(bottom states in m_blocks[bi], add_work(..., 1), *this);
              // not necessary, as the inner loop (mark all transitions out of bottom states) is always executed
          assert(m_blocks[bi].contains_new_bottom_states);
          m_blocks[bi].contains_new_bottom_states = false;

// 2. Register all other sets with transitions out of blocks containing new bottom states as splitters.
//    (If such a set contains transitions out of new bottom states, the work can be ascribed to them;
//    otherwise more new bottom states will be found later, to which the work can be ascribed in hindsight.)
//    This registration is done by registering slices in m_BLC_transitions.
          typename linked_list<BLC_indicators>::iterator ind = m_blocks[bi].block_to_constellation.begin();
          if (m_blocks[bi].block_to_constellation.end() != ind)
          {
            const transition& first_t = m_aut.get_transitions()[*ind->start_same_BLC];
            assert(m_states[first_t.from()].block == bi);
            if (is_inert_during_init_if_branching(first_t) &&
                m_blocks[bi].constellation == m_blocks[m_states[first_t.to()].block].constellation)
            {
              // The first BLC-set is constellation-inert, so skip it
              assert(!has_marked_transitions(*ind));
              ++ind;
            }
            for (; m_blocks[bi].block_to_constellation.end() != ind; ++ind)
            {
              assert(!has_marked_transitions(*ind));
#ifndef NDEBUG
              const transition& first_t = m_aut.get_transitions()[*ind->start_same_BLC];
              assert(m_states[first_t.from()].block == bi);
              assert(!is_inert_during_init(first_t) ||
                     m_blocks[bi].constellation != m_blocks[m_states[first_t.to()].block].constellation);
#endif
              // BLC set transitions are not constellation-inert, so we need to stabilize under them
              Qhat.emplace_back(ind->start_same_BLC, ind->start_marked_BLC);
              // In the above line we need that !has_marked_transitions(*ind).

              #ifdef CHECK_COMPLEXITY_GJ
                // The work is assigned to the transitions out of new bottom states in ind.
                // Try to find a new bottom state to which to assign it.
                bool work_assigned = false;
                // assign the work to the transitions out of bottom states in this BLC-set
                for (BLC_list_const_iterator work_it = ind->start_same_BLC; work_it < ind->start_marked_BLC; ++work_it)
                {
                  // assign the work to this transition
                  if (0 == m_states[m_aut.get_transitions()[*work_it].from()].no_of_outgoing_inert_transitions)
                  {
                    mCRL2complexity_gj(&m_transitions[*work_it], add_work(check_complexity::stabilizeB__initialize_Qhat, 1), *this);
                    work_assigned = true;
                  }
                }
                if (!work_assigned)
                {
                  // We register that we still have to find a transition from a new bottom state in this slice.
                  // we should do the following:
                  // - store the range [blc_it->start_same_BLC, blc_it->end_same_BLC) in a temporary variable
                  //   (We cannot use the blc_it pointer here because the source block
                  //   might be split up later, in a part that contains only old
                  //   bottom states and a part that does contain new bottom states.
                  //   Then, we should decide which part[s] contain new bottom states...)
                  // - expect that one of the source states of this range will become
                  //   a new bottom state in m_blocks_with_new_bottom_states
                  //   later in the main loop of stabilizeB().
                  // - towards the end of the main loop, go through the ranges stored in the temporary variable
                  //   and mark the transitions that start in a bottm state in m_blocks_with_new_bottom_states.
                  //   For every range, we must mark at least one transition.
mCRL2log(log::debug) << "Haven't yet found a transition from a new bottom state in " << ind->debug_id(*this) << " to assign the work to\n";
                  work_to_assign_later.emplace_back(ind->start_same_BLC, ind->start_marked_BLC);
                  // In the above line we need that !has_marked_transitions(*ind).
                }
              #endif
            }
          }

          // mark all transitions out of (new) bottom states
          assert(m_blocks[bi].start_bottom_states < m_blocks[bi].start_non_bottom_states);
          for (std::vector<state_index>::iterator si = m_blocks[bi].start_bottom_states; si < m_blocks[bi].start_non_bottom_states; ++si)
          {
            mCRL2complexity_gj(&m_states[*si], add_work(check_complexity::stabilizeB__distribute_states_over_Phat, 1), *this);
            assert(m_states[*si].block == bi);
            outgoing_transitions_it end_it = 1+*si >= m_states.size() ? m_outgoing_transitions.end() : m_states[1+*si].start_outgoing_transitions;
            for (outgoing_transitions_it ti = m_states[*si].start_outgoing_transitions; ti < end_it; ++ti)
            {
              // mCRL2complexity_gj(&m_transitions[m_BLC_transitions[ti->transition]], add_work(..., 1), *this);
                  // subsumed under the above counter
              const transition& t = m_aut.get_transitions()[m_BLC_transitions[ti->transition]];
              assert(t.from() == *si);
              if (!is_inert_during_init_if_branching(t) ||
                  m_blocks[bi].constellation != m_blocks[m_states[t.to()].block].constellation)
              {
                // the transition is not constellation-inert, so mark it
                mark_BLC_transition(ti);
              }
              else
              {
                // skip all other constellation-inert transitions (this is an optimization)
                assert(ti <= ti->start_same_saC);
                ti = ti->start_same_saC;
              }
            }
          }
        }
        clear(m_blocks_with_new_bottom_states);

// 3. As long as there are registered slices in m_BLC_transitions, select any one of them.
//    Take the first BLC_indicator that has transitions in this slice; remove it from the slice;
//    if the slice is now empty remove it from the register.
//    Do a normal splitB() under this splitter.
//    If more new bottom states are created, store them in the new m_blocks_with_new_bottom_states.

        // Algorithm 4, line 4.8.
        for (std::pair<BLC_list_iterator, BLC_list_iterator> Qhat_elt: Qhat)
        {
          // Algorithm 4, line 4.9.
          // mCRL2complexity_gj(..., add_work(..., max_C), *this);
              // not needed as the inner loop is always executed at least once.
          //print_data_structures("Main stabilize loop");
          assert(check_data_structures("New bottom state loop", false, false));
          check_stability("New bottom state loop", &Qhat, &Qhat_elt);

          assert(Qhat_elt.first < Qhat_elt.second);
          do
          {
            linked_list<BLC_indicators>::iterator splitter = m_transitions[*std::prev(Qhat_elt.second)].transitions_per_block_to_constellation;
            // TODO: mCRL2complexity_gj(splitter, add_work(..., max_C), *this);
                // work should be added to the current (and possibly future) new bottom states of splitter.
            Qhat_elt.second = splitter->start_same_BLC;

            const transition& first_t = m_aut.get_transitions()[*splitter->start_same_BLC];
            const block_type& from_block = m_blocks[m_states[first_t.from()].block];
            if (from_block.contains_new_bottom_states)
            {
              // The block contains even more new bottom states and does not need to be stabilized in this round.
              assert(!has_marked_transitions(*splitter));
            }
            else if (is_inert_during_init_if_branching(first_t) &&
                     from_block.constellation == m_blocks[m_states[first_t.to()].block].constellation)
            {
              // the BLC slice only contains constellation-inert transitions.
              // We do not need to do stabilize this.
              assert(!has_marked_transitions(*splitter));
              // In fact, the BLC set should not have been added at all!
              assert(0);
            }
            else
            {
              // Check whether the bottom states of Bpp are not all included in Mleft.

          // Algorithm 4, line 4.10.
              const BLC_list_iterator splitter_end_unmarked_BLC = splitter->start_marked_BLC;
              std::vector<state_index>::iterator first_unmarked_bottom_state = not_all_bottom_states_are_touched(splitter
                        #ifndef NDEBUG
                          , splitter_end_unmarked_BLC
                        #endif
                        );
              if (first_unmarked_bottom_state < from_block.start_non_bottom_states)
              {
// std::cerr << "PERFORM A NEW BOTTOM STATE SPLIT\n";
            // Algorithm 4, line 4.11, and implicitly 4.12, 4.13 and 4.18.
                splitB(splitter, first_unmarked_bottom_state, splitter_end_unmarked_BLC);
              }
            }
          }
          while (Qhat_elt.first < Qhat_elt.second);
        }

        #ifdef CHECK_COMPLEXITY_GJ
          // We now have to try and find further new bottom states to which to assign
          // the initialization of Qhat that had not yet been assigned earlier.
          for (std::vector<std::pair<BLC_list_const_iterator, BLC_list_const_iterator> >::iterator qhat_it = work_to_assign_later.begin(); qhat_it != work_to_assign_later.end(); ++qhat_it)
          {
            bool new_bottom_state_with_transition_found = false;
            for (BLC_list_const_iterator work_it = qhat_it->first; work_it < qhat_it->second; ++work_it)
            {
              const state_index t_from = m_aut.get_transitions()[*work_it].from();
              if (0 == m_states[t_from].no_of_outgoing_inert_transitions)
              {
                // t_from is a new bottom state, so we can assign the work to this transition
// std::cerr << m_transitions[*work_it].debug_id(*this) << " is assigned work on bottom states afterwards\n";
                assert(m_blocks[m_states[t_from].block].contains_new_bottom_states);
                mCRL2complexity_gj(&m_transitions[*work_it], add_work(check_complexity::stabilizeB__initialize_Qhat_afterwards, 1), *this);
                new_bottom_state_with_transition_found = true;
              }
            }
            assert(new_bottom_state_with_transition_found);
          }
        #endif
        // destroy and deallocate work_to_assign_later here
// 4. After all registered slices have been handled: If more new bottom states have been created,
//    restart from 0. above. Otherwise finish the procedure.
        // Algorithm 4, line 4.17.
      }

      return;
    }

    void maintain_block_label_to_cotransition(
                   const block_index old_block,
                   const block_index new_block,
                   const transition_index moved_transition,
                   const transition_index alternative_transition,
                   block_label_to_size_t_map& block_label_to_cotransition,
                   const constellation_index ci) const
    {
      const transition& t_move=m_aut.get_transitions()[moved_transition];
      const label_index lbl = label_or_divergence(t_move);
      assert(m_states[t_move.from()].block == new_block);
      if (null_transition != alternative_transition)
      {
        assert(old_block == m_states[m_aut.get_transitions()[alternative_transition].from()].block ||
               new_block == m_states[m_aut.get_transitions()[alternative_transition].from()].block);
        assert(lbl == label_or_divergence(m_aut.get_transitions()[alternative_transition]));
        assert(m_blocks[m_states[t_move.to()].block].constellation == m_blocks[m_states[m_aut.get_transitions()[alternative_transition].to()].block].constellation);
      }
      if (m_blocks[m_states[t_move.to()].block].constellation==ci &&
          (!is_inert_during_init(t_move) || m_blocks[m_states[t_move.from()].block].constellation!=ci))
      {
        // This is a non-constellation-inert transition to the current co-constellation.

        typename block_label_to_size_t_map::iterator bltc_it=
                       block_label_to_cotransition.find(std::pair(old_block,lbl));
        if (bltc_it!=block_label_to_cotransition.end())
        {
          if (bltc_it->second==moved_transition)
          {
            // This transition is being moved. Find a replacement in block_label_to_cotransition.
            bltc_it->second=alternative_transition;
          }
        }

        // Check whether there is a representation for the new_block in block_label_to_cotransition.
        bltc_it=block_label_to_cotransition.find(std::pair(new_block,lbl));
        if (bltc_it==block_label_to_cotransition.end())
        {
          // No such transition exists as yet. Give moved transition this purpose.
          block_label_to_cotransition[std::pair(new_block,lbl)]=moved_transition;
        }
      }
    }

    // finds a splitter for the tau-transitions from block B (= the new constellation) to old_constellation,
    // if such a splitter exists.
    // It uses the fact that the first element block_to_constellation contains the inert transitions (if there are any),
    // and just after splitting the new constellation off from the old one,
    // the element immediately after that the tau-transitions from the new to the old constellation.
    linked_list<BLC_indicators>::iterator find_inert_co_transition_for_block(const block_index index_block_B, const constellation_index old_constellation)
    {
      linked_list< BLC_indicators >::iterator btc_it= m_blocks[index_block_B].block_to_constellation.begin();
      if (btc_it == m_blocks[index_block_B].block_to_constellation.end())
      {
        // The new constellation has no outgoing transitions at all.
        return m_blocks[index_block_B].block_to_constellation.end();
      }
      const transition& btc_t=m_aut.get_transitions()[*(btc_it->start_same_BLC)];
      if (!is_inert_during_init_if_branching(btc_t))
      {
        // The new constellation has no outgoing tau-transitions at all (except possibly tau-self-loops, for divergence-preserving branching bisimulation).
        return m_blocks[index_block_B].block_to_constellation.end();
      }
      if (m_blocks[m_states[btc_t.to()].block].constellation==old_constellation)
      {
        // The new constellation has no inert transitions but it does have tau-transitions to the old constellation
        // (which were inert before).
        return &*btc_it;
      }
      if (m_blocks[m_states[btc_t.to()].block].constellation != m_constellations.size() - 1)
      {
        // The new constellation, before it was separated from the old one,
        // had no constellation-inert outgoing transitions.
        return m_blocks[index_block_B].block_to_constellation.end();
      }
      // *btc_it is the BLC_indicator for the inert transitions of the new constellation.
      // Try the second element in the list:
      ++btc_it;
      if (btc_it == m_blocks[index_block_B].block_to_constellation.end())
      {
        // The new constellation has no other outgoing transitions.
        return m_blocks[index_block_B].block_to_constellation.end();
      }
      const transition& btc2_t=m_aut.get_transitions()[*(btc_it->start_same_BLC)];
      if (!is_inert_during_init_if_branching(btc2_t) ||
          m_blocks[m_states[btc2_t.to()].block].constellation != old_constellation) // The new constellation has no tau-transitions to the old constellation.
      {
        return m_blocks[index_block_B].block_to_constellation.end();
      }
      return &*btc_it;
    }

/*
    // This routine can only be called after initialisation, i.e., when m_outgoing_transitions[t].transition refers to
    // a position in m_BLC_transitions.
    bool state_has_outgoing_co_transition(const transition_index transition_to_bi, const constellation_index old_constellation)
    {
      // i1 refers to the position of the transition_to_bi;
      outgoing_transitions_it i1=m_transitions[transition_to_bi].ref_outgoing_transitions;
      // i2 refers to the last position, unless i1 is the last, then i2 is the first of all transitions with the same s, a and C.
      outgoing_transitions_it i2=i1->start_same_saC;
      if (i2<=i1)
      {
        // So, i1 is the last element with the same s, a, C.
        i2=i1+1;
      }
      else
      {
        // In this case i2 refers to the last element with the same s, a, C.
        i2++;
      }
      if (i2==m_outgoing_transitions.end())
      {
        return false;
      }
      const transition& t1=m_aut.get_transitions()[transition_to_bi];
      const transition& t2=m_aut.get_transitions()[m_BLC_transitions[i2->transition]];
      return t1.from()==t2.from() && m_aut.apply_hidden_label_map(t1.label())==m_aut.apply_hidden_label_map(t2.label()) && m_blocks[m_states[t2.to()].block].constellation==old_constellation;
    }

    // This function determines whether all bottom states in B have outgoing co-transitions. If yes false is reported.
    // If no, true is reported and the source states with outgoing co-transitions are added to m_R and those without outgoing
    // co-transitions are added to m_u. The counters of these states are set to Rmarked or to 0. This already initialises
    // these states for the splitting process.
    // This routine is called after initialisation.
    bool some_bottom_state_has_no_outgoing_co_transition(block_index B,
                                                         std::vector<transition_index>::iterator transitions_begin,
                                                         std::vector<transition_index>::iterator transitions_end,
                                                         const constellation_index old_constellation,
                                                         std::vector<state_index>::iterator& first_unmarked_bottom_state)
    {
      first_unmarked_bottom_state = m_blocks[B].start_bottom_states;
      #ifdef CHECK_COMPLEXITY_GJ
        const constellation_index new_constellation = m_constellations.size()-1;
        const unsigned max_C = check_complexity::log_n - check_complexity::ilog2(number_of_states_in_constellation(new_constellation));
      #endif
      for(std::vector<transition_index>::iterator ti=transitions_begin; ti!=transitions_end; ++ti)
      {
        const transition& t=m_aut.get_transitions()[*ti];
//std::cerr << "INSPECT TRANSITION  " << t.from() << " -" << m_aut.action_label(t.label()) << "-> " << t.to() << "\n";
        #ifdef CHECK_COMPLEXITY_GJ
          assert(new_constellation==m_blocks[m_states[t.to()].block].constellation);
          mCRL2complexity_gj(&m_transitions[*ti], add_work(check_complexity::some_bottom_state_has_no_outgoing_co_transition__handle_transition, max_C), *this);
        #endif
        const state_index s = t.from();
        assert(m_states[s].ref_states_in_blocks>=m_blocks[B].start_bottom_states);
        assert(m_states[s].ref_states_in_blocks<m_blocks[B].end_states);
        if (m_states[s].ref_states_in_blocks<m_blocks[B].start_non_bottom_states &&
            m_states[s].counter==undefined)
        {
            if (state_has_outgoing_co_transition(*ti,old_constellation))
            {
              m_states[s].counter=Rmarked;
              m_R.add_todo(s);
              const std::vector<state_index>::iterator pos_s=m_states[s].ref_states_in_blocks;
              assert(pos_s < m_blocks[B].start_non_bottom_states);
              assert(first_unmarked_bottom_state <= pos_s);
              swap_states_in_states_in_block(first_unmarked_bottom_state, pos_s); // Move marked states to the front.
              first_unmarked_bottom_state++;
            }
            else
            {
              // We need to register that we added this state to U
              // because it might have other transitions to new_constellation,
              // so it will be revisited
              m_U_counter_reset_vector.push_back(t.from());
              m_states[t.from()].counter=0;
              m_U.add_todo(t.from());

            }
        }
      }

      assert(m_R.size() + m_U.size() == static_cast<std::size_t>(std::distance(m_blocks[B].start_bottom_states, m_blocks[B].start_non_bottom_states)));
      if (first_unmarked_bottom_state != m_blocks[B].start_non_bottom_states)
      {
        return true; // A split can commence.
      }
      // Otherwise, reset the marks.
      clear_state_counters();
      m_R.clear();
      m_U.clear();
      return false;
    }

    // Check if there is a state in W that has no outgoing a transition to some constellation.
    // If so, return false, but set in m_R and m_U whether those states in W have or have no
    // outgoing transitions. Set m_states[s].counter accordingly.
    // If all states in W have outgoing transitions with the label and constellation, leave
    // m_R, m_U, m_states[s].counters and m_U_counter_reset vector untouched.

    bool W_empty(const set_of_states_type& W, const set_of_states_type& aux,
                 #ifdef CHECK_COMPLEXITY_GJ
                   const label_index a,
                   const constellation_index C,
                 #endif
                 std::vector<state_index>::iterator& first_unmarked_bottom_state
                 )
    {
      bool W_empty=true;
//std::cerr << "W: "; for(auto s: W) { std::cerr << s << " "; } std::cerr << "\n";
      #ifndef NDEBUG
        #ifdef CHECK_COMPLEXITY_GJ
          assert(m_aut.apply_hidden_label_map(a) == a);
        #endif
        const block_type& B = m_blocks[m_states[*first_unmarked_bottom_state].block];
        assert(B.start_bottom_states == first_unmarked_bottom_state);
        // assert(static_cast<std::ptrdiff_t>(std::distance(B.start_bottom_states, B.start_non_bottom_states)) == W.size());
      #endif
      for(const state_index si: W)
      {
        assert(0 == m_states[si].no_of_outgoing_inert_transitions);
        assert(m_states[*first_unmarked_bottom_state].block == m_states[si].block);
        if (aux.count(si)==0)
        {
          // I do not know how to account for the work in this case.
          W_empty=false;
          m_U.add_todo(si);
        }
        else
        {
          add_work_to_same_saC(false, si, a, C, check_complexity::W_empty__find_new_bottom_state_in_R, 1);
          assert(undefined == m_states[si].counter);
          m_states[si].counter=Rmarked;
          m_R.add_todo(si);
          const std::vector<state_index>::iterator pos_s=m_states[si].ref_states_in_blocks;
          assert(first_unmarked_bottom_state <= pos_s);
          assert(pos_s < B.start_non_bottom_states);
          swap_states_in_states_in_block(first_unmarked_bottom_state, pos_s); // Move marked states to the front.
          first_unmarked_bottom_state++;
        }
      }
      if (!W_empty)
      {
        return false; // A split can commence.
      }
      // Otherwise, reset the marks.
      clear_state_counters();
      m_R.clear();
      m_U.clear();
      return true;
    } */

    // Mark all states that are sources of the marked transitions in splitter
    // and move them to the beginning.
    // (Marking of bottom states is done by moving them to the beginning of the
    // bottom states. Marking of non-bottom states is done by adding them to m_R
    // and setting m_states[*].counter = Rmarked.)
    // If all bottom states are marked, reset markers and m_R. Otherwise, if not all
    // bottom states are touched, leave the marked non-bottom states in m_R and leave
    // the markers in m_states[*].counter in place.
    // The marked bottom states are moved to the front in m_states_in_blocks and
    // the return value indicates the position (in m_states_in_blocks)
    // of the first non-marked bottom state.
    std::vector<state_index>::iterator not_all_bottom_states_are_touched(linked_list<BLC_indicators>::iterator splitter
            #ifndef NDEBUG
              , const BLC_list_const_iterator splitter_end_unmarked_BLC /* = splitter->start_marked_BLC -- but this default argument is not allowed */
            #endif
            )
    {
// std::cerr << "not_all_bottom_states_are_touched(" << splitter->debug_id(*this) << ")\n";
      const block_index bi = m_states[m_aut.get_transitions()[*splitter->start_same_BLC].from()].block;
      const block_type& B=m_blocks[bi];
      std::vector<state_index>::iterator first_unmarked_bottom_state = B.start_bottom_states;
      BLC_list_iterator marked_t_it = splitter->start_marked_BLC;
      for(; points_into_BLC_set(marked_t_it, *splitter); ++marked_t_it)
      {
        const transition& t = m_aut.get_transitions()[*marked_t_it];
        const state_index s=t.from();
        assert(m_states[s].block == bi);
        // mCRL2complexity_gj(&m_transitions[*i], add_work(...), *this);
            // not needed because this work can be attributed to the marking of the transition
        const std::vector<state_index>::iterator pos_s=m_states[s].ref_states_in_blocks;
        assert(B.start_bottom_states <= pos_s);
        assert(pos_s < B.end_states);
        if (first_unmarked_bottom_state <= pos_s)
        {
          if (pos_s<B.start_non_bottom_states)
          {
            assert(0 == m_states[s].no_of_outgoing_inert_transitions);
            assert(undefined == m_states[s].counter);
            swap_states_in_states_in_block(first_unmarked_bottom_state, pos_s); // Move marked states to the front.
            first_unmarked_bottom_state++;
          }
          else if (Rmarked != m_states[s].counter)
          {
            assert(0 < m_states[s].no_of_outgoing_inert_transitions);
            assert(undefined == m_states[s].counter);
            m_R.add_todo(s);
            m_states[s].counter=Rmarked;
          }
        }
      }
      #ifndef NDEBUG
        // ensure that the unmarked transitions do not add any bottom states.
        assert(splitter->start_same_BLC <= splitter_end_unmarked_BLC);
        for(BLC_list_const_iterator i = splitter->start_same_BLC; i < splitter_end_unmarked_BLC; ++i)
        {
          assert(splitter_end_unmarked_BLC == splitter->start_marked_BLC);
          const transition& t = m_aut.get_transitions()[*i];
          const state_index s=t.from();
          assert(m_states[s].block == bi);
          const std::vector<state_index>::const_iterator pos_s=m_states[s].ref_states_in_blocks;
          assert(*pos_s == s);
          assert(B.start_bottom_states <= pos_s);
          assert(pos_s < B.end_states);
          if (0 == m_states[s].no_of_outgoing_inert_transitions)
          {
            // State s is a bottom state. It should already have been marked.
            assert(pos_s < first_unmarked_bottom_state);
          }
          else
          {
            assert(B.start_non_bottom_states <= pos_s);
          }
        }
      #endif
      // also unmark all transitions:
      splitter->start_marked_BLC = marked_t_it;
      if (first_unmarked_bottom_state==B.start_non_bottom_states)
      {
        // All bottom states are marked. No splitting is possible. Reset m_R, m_states[s].counter for s in m_R.
        clear_state_counters(false);
        m_R.clear();
      }
      return first_unmarked_bottom_state;
    }

/*
    // This routine can only be called after initialisation.
    bool hatU_does_not_cover_B_bottom(const block_index index_block_B,
                                      const constellation_index old_constellation,
                                      std::vector<state_index>::iterator first_unmarked_bottom_state)
    {
      mCRL2complexity_gj(&m_blocks[index_block_B], add_work(check_complexity::hatU_does_not_cover_B_bottom__handle_bottom_states_and_their_outgoing_transitions_in_splitter, check_complexity::log_n - check_complexity::ilog2(number_of_states_in_block(index_block_B))), *this);
      assert(m_branching);
      first_unmarked_bottom_state = m_blocks[index_block_B].start_bottom_states;
      for(typename std::vector<state_index>::iterator si=m_blocks[index_block_B].start_bottom_states;
                        si!=m_blocks[index_block_B].start_non_bottom_states;
                      ++si)
      {
        // mCRL2complexity_gj(&m_states[*si], add_work(..., max_C), *this);
            // subsumed by the above call
        bool found=false;
        const outgoing_transitions_it end_it=((*si)+1>=m_states.size())?m_outgoing_transitions.end():m_states[(*si)+1].start_outgoing_transitions;
        for(outgoing_transitions_it tti=m_states[*si].start_outgoing_transitions;
                                     !found && tti!=end_it;
                                     ++tti)
        {
          // mCRL2complexity_gj(&m_transitions[tti->transition], add_work(..., max_C), *this);
          // subsumed by the above call
          const transition& t=m_aut.get_transitions()[m_BLC_transitions[tti->transition]];
          assert(t.from() == *si);
          if (m_aut.is_tau(m_aut.apply_hidden_label_map(t.label())) && m_blocks[m_states[t.to()].block].constellation==old_constellation)
          {
            found =true;
          }
          assert(tti <= tti->start_same_saC);
          tti = tti->start_same_saC;
        }
        if (!found)
        {
          // This state has no constellation-inert tau transition to the old constellation.
          m_U.add_todo(*si);
        }
        else
        {
          // The state *si has a tau transition to the old constellation that has just become constellation-non-inert.
          m_R.add_todo(*si);
          m_states[*si].counter=Rmarked;
          assert(si == m_states[*si].ref_states_in_blocks);
          assert(first_unmarked_bottom_state <= si);
          assert(si < m_blocks[index_block_B].start_non_bottom_states);
          swap_states_in_states_in_block(first_unmarked_bottom_state, si); // Move marked states to the front.
          first_unmarked_bottom_state++;
        }
      }
      if (first_unmarked_bottom_state != m_blocks[index_block_B].start_non_bottom_states)
      {
        // Splitting can commence.
        return true;
      }
      else
      {
        // Splitting is not possible. Reset the counter in m_states.
        m_U.clear();
        clear_state_counters(true);
        m_R.clear();
        return false;
      }
    } */


    // Select a block that is not the largest block in the constellation.
    // It is advantageous to select the smallest block.
    // The constellation ci is returned.
    block_index select_and_remove_a_block_in_a_non_trivial_constellation(constellation_index& ci)
    {
      // Do the minimal checking, i.e., only check two blocks in a constellation.
      ci=m_non_trivial_constellations.back();
      block_index index_block_B = m_states[*m_constellations[ci].start_const_states].block;          // The first block.
      block_index second_block_B= m_states[*std::prev(m_constellations[ci].end_const_states)].block; // The last block.

      if (number_of_states_in_block(index_block_B)<=number_of_states_in_block(second_block_B))
      {
        m_constellations[ci].start_const_states = m_blocks[index_block_B].end_states;
      }
      else
      {
        m_constellations[ci].end_const_states = m_blocks[second_block_B].start_bottom_states;
        index_block_B=second_block_B;
      }
      return index_block_B;
    }


    void refine_partition_until_it_becomes_stable()
    {
      // This represents the while loop in Algorithm 1 from line 1.6 to 1.25.

      // The instruction below has complexity O(|Act|);
      // calM will contain the m_BLC_transitions slices that need stabilization:
      std::vector<std::pair<BLC_list_iterator, BLC_list_iterator> > calM;
      // Algorithm 1, line 1.6.
      while (!m_non_trivial_constellations.empty())
      {
        /* static time_t last_log_time=time(nullptr)-1;
        time_t new_log_time = 0;
        if (time(&new_log_time)>last_log_time)
        {
          mCRL2log(log::verbose) << "Refining. There are " << m_blocks.size() << " blocks and " << m_constellations.size() << " constellations.\n";
          last_log_time=last_log_time = new_log_time;
        } */
        // print_data_structures("MAIN LOOP");
        assert(check_data_structures("MAIN LOOP"));
        assert(check_stability("MAIN LOOP"));

        // Algorithm 1, line 1.7.
        constellation_index ci=-1;
        block_index index_block_B=select_and_remove_a_block_in_a_non_trivial_constellation(ci);
// std::cerr << "REMOVE BLOCK " << index_block_B << " from constellation " << ci << "\n";

        // Algorithm 1, line 1.8.
        if (m_states[*m_constellations[ci].start_const_states].block == m_states[*std::prev(m_constellations[ci].end_const_states)].block)
        {
          // Constellation has become trivial.
          assert(m_non_trivial_constellations.back()==ci);
          m_non_trivial_constellations.pop_back();
        }
        m_constellations.emplace_back(m_blocks[index_block_B].start_bottom_states, m_blocks[index_block_B].end_states);
        const constellation_index old_constellation=m_blocks[index_block_B].constellation;
        assert(old_constellation == ci);
        const constellation_index new_constellation=m_constellations.size()-1;
        m_blocks[index_block_B].constellation=new_constellation;
        #ifdef CHECK_COMPLEXITY_GJ
          // m_constellations[new_constellation].work_counter = m_constellations[old_constellation].work_counter;
          const unsigned max_C = check_complexity::log_n - check_complexity::ilog2(number_of_states_in_constellation(new_constellation));
          mCRL2complexity_gj(&m_blocks[index_block_B], add_work(check_complexity::refine_partition_until_it_becomes_stable__find_splitter, max_C), *this);
        #endif
        // Here the variables block_to_constellation and the doubly linked list L_B->C in blocks must be still be updated.
        // This happens further below.

        // Algorithm 1, line 1.9.
        block_label_to_size_t_map block_label_to_cotransition;

        for(typename std::vector<state_index>::iterator i=m_blocks[index_block_B].start_bottom_states;
                                                        i!=m_blocks[index_block_B].end_states; ++i)
        {
          // mCRL2complexity_gj(m_states[*i], add_work(check_complexity::..., max_C), *this);
              // subsumed under the above counter
          // and visit the incoming transitions.
          const std::vector<transition>::iterator end_it=
                          ((*i)+1==m_states.size())?m_aut.get_transitions().end()
                                                   :m_states[(*i)+1].start_incoming_transitions;
          for(std::vector<transition>::iterator j=m_states[*i].start_incoming_transitions; j!=end_it; ++j)
          {
            const transition& t=*j;
            const transition_index t_index = std::distance(m_aut.get_transitions().begin(),j);
            // mCRL2complexity_gj(&m_transitions[t_index], add_work(check_complexity::..., max_C), *this);
                // subsumed under the above counter

            // Update the state-action-constellation (saC) references in m_outgoing_transitions.
            const outgoing_transitions_it old_pos = m_transitions[t_index].ref_outgoing_transitions;
            const outgoing_transitions_it end_same_saC = old_pos->start_same_saC < old_pos ? old_pos : old_pos->start_same_saC;
            const outgoing_transitions_it new_pos = end_same_saC->start_same_saC;
            assert(m_states[t.from()].start_outgoing_transitions <= new_pos);
            assert(new_pos <= old_pos);
            if (old_pos != new_pos)
            {
              std::swap(old_pos->transition,new_pos->transition);
              m_transitions[m_BLC_transitions[old_pos->transition]].ref_outgoing_transitions = old_pos;
              m_transitions[m_BLC_transitions[new_pos->transition]].ref_outgoing_transitions = new_pos;
            }
            if (new_pos < end_same_saC)
            {
              end_same_saC->start_same_saC = std::next(new_pos);
            }
            // correct start_same_saC provisionally: make them at least point at each other.
            // In the new saC-slice, all transitions point to the first one, except the first one: that shall point at the last one.
            new_pos->start_same_saC = new_pos;
            if (m_states[t.from()].start_outgoing_transitions < new_pos)
            {
              // Check if t is the first transition in the new saC slice:
              const transition_index prev_t_index = m_BLC_transitions[std::prev(new_pos)->transition];
              const transition& prev_t = m_aut.get_transitions()[prev_t_index];
              assert(prev_t.from() == t.from());
              if (m_states[prev_t.to()].block == index_block_B && label_or_divergence(prev_t) == label_or_divergence(t))
              {
                // prev_t also belongs to the new saC slice.
                new_pos->start_same_saC = std::prev(new_pos)->start_same_saC;
                assert(m_states[t.from()].start_outgoing_transitions <= new_pos->start_same_saC);
                assert(new_pos->start_same_saC < new_pos);
                assert(std::prev(new_pos) == new_pos->start_same_saC->start_same_saC);
                new_pos->start_same_saC->start_same_saC = new_pos;
              }
            }
          }
        }
        calM.clear();

        // Walk through all states in block B
        for(typename std::vector<state_index>::iterator i=m_blocks[index_block_B].start_bottom_states;
                                                        i!=m_blocks[index_block_B].end_states; ++i)
        {
          // mCRL2complexity_gj(m_states[*i], add_work(check_complexity::..., max_C), *this);
              // subsumed under the above counter

          // and visit the incoming transitions.
          const std::vector<transition>::iterator end_it=
                          ((*i)+1==m_states.size())?m_aut.get_transitions().end()
                                                   :m_states[(*i)+1].start_incoming_transitions;
          for(std::vector<transition>::iterator j=m_states[*i].start_incoming_transitions; j!=end_it; ++j)
          {
            const transition& t=*j;
            const transition_index t_index=std::distance(m_aut.get_transitions().begin(),j);
            // mCRL2complexity_gj(&m_transitions[t_index], add_work(check_complexity::..., max_C), *this);
                // subsumed under the above counter

            // Give the saC slice of this transition its final correction
            const outgoing_transitions_it out_pos = m_transitions[t_index].ref_outgoing_transitions;
            const outgoing_transitions_it start_new_saC = out_pos->start_same_saC;
            outgoing_transitions_it start_old_saC = start_new_saC;
            if (start_new_saC < out_pos)
            {
              if (out_pos < start_new_saC->start_same_saC)
              {
                out_pos->start_same_saC = start_new_saC->start_same_saC;
              }
              start_old_saC = start_new_saC->start_same_saC;
            }

            // mark one cotransition from the same source state (if it exists)
            // If the main transitions are constellation-inert (i.e. they go
            // from index_block_B to index_block_B), we shall not stabilize
            // under them, so we shall also not mark these transitions.
            // (This is also needed to ensure that the BLC set for the special
            // new-bottom-state split does not contain any spurious markings.)
            // If the main transitions start in the old constellation, the
            // co-transitions are constellation-inert and we do not need to
            // stabilize under them.

            // Perhaps we should only mark a co-transition if the source state is a bottom state.
            // Otherwise it might interfere with the handling of new bottom states.

            if (!m_blocks[m_states[t.from()].block].contains_new_bottom_states &&
                (!is_inert_during_init(t) || (m_states[t.from()].block != index_block_B &&
                                              m_blocks[m_states[t.from()].block].constellation != old_constellation)))
            {
              ++start_old_saC;
              if (start_old_saC < m_outgoing_transitions.end() &&
                  (t.from()+1 >= m_aut.num_states() || start_old_saC < m_states[t.from()+1].start_outgoing_transitions))
              {
                const transition& old_t = m_aut.get_transitions()[m_BLC_transitions[start_old_saC->transition]];
                assert(old_t.from() == t.from());
                assert(m_states[t.to()].block == index_block_B);
                if (label_or_divergence(t) == label_or_divergence(old_t) &&
                    old_constellation == m_blocks[m_states[old_t.to()].block].constellation)
                {
// std::cerr << "Marking " << m_transitions[m_BLC_transitions[start_old_saC->transition]].debug_id(*this) << " as a cotransition.\n";
                  mark_BLC_transition(start_old_saC);
                }
              }
            }

            // Update the block_label_to_cotransition map.
            if (block_label_to_cotransition.find(std::pair(m_states[t.from()].block,label_or_divergence(t))) == block_label_to_cotransition.end())
            {
              // Not found. Add a transition from the LBC_list to block_label_to_cotransition
              // that goes to C\B, or the null_transition if no such transition exists, which prevents searching
              // the list again. Except if t.from is in C\B and a=tau, because in that case it is a (former) constellation-inert transition.
              bool found=false;

              if (!is_inert_during_init(t) || m_blocks[m_states[t.from()].block].constellation!=ci)
              {
                const BLC_indicators& ind = *m_transitions[t_index].transitions_per_block_to_constellation;
                BLC_list_iterator transition_walker=ind.start_same_BLC;

                while (points_into_BLC_set(transition_walker, ind))
                {
                  const transition& tw=m_aut.get_transitions()[*transition_walker];
                  assert(m_states[tw.from()].block == m_states[t.from()].block);
                  assert(label_or_divergence(tw) == label_or_divergence(t));
                  if (m_blocks[m_states[tw.to()].block].constellation==ci)
                  {
                    found=true;
                    block_label_to_cotransition[std::pair(m_states[t.from()].block,label_or_divergence(t))] = *transition_walker;
                    break;
                  }
                  assert(m_blocks[m_states[tw.to()].block].constellation==new_constellation);
                  mCRL2complexity_gj(&m_transitions[*transition_walker], add_work(check_complexity::refine_partition_until_it_becomes_stable__find_cotransition, max_C), *this);
                  ++transition_walker;
                }
              }
              if (!found)
              {
                block_label_to_cotransition[std::pair(m_states[t.from()].block,label_or_divergence(t))] = null_transition;
              }
            }
            // Update the doubly linked list L_B->C in blocks as the constellation is split in B and C\B.
            if (update_the_doubly_linked_list_LBC_new_constellation(index_block_B, t, t_index) &&
                !(is_inert_during_init(t) &&
                  (assert(m_states[t.to()].block == index_block_B),
                   m_states[t.from()].block == index_block_B)))
            {
              // a new BLC set has been constructed, insert its start position into calM.
              BLC_list_iterator BLC_pos = m_BLC_transitions.begin() + m_transitions[t_index].ref_outgoing_transitions->transition;
              assert(t_index == *BLC_pos);
              calM.emplace_back(BLC_pos, BLC_pos);
              // The end-position (the second element in the pair) will need to be corrected later.
            }
          }
        }

        // mark all states in main splitters and correct the end-positions of calM entries
        for (std::pair<BLC_list_iterator, BLC_list_iterator>& calM_elt: calM)
        {
          linked_list <BLC_indicators>::iterator ind = m_transitions[*calM_elt.first].transitions_per_block_to_constellation;
          mCRL2complexity_gj(ind, add_work(check_complexity::refine_partition_until_it_becomes_stable__correct_end_of_calM, max_C), *this);
          assert(ind->start_same_BLC == calM_elt.first);
          assert(ind->start_same_BLC == calM_elt.second);
          assert(!has_marked_transitions(*ind));
          calM_elt.second = ind->start_marked_BLC;
          ind->start_marked_BLC = ind->start_same_BLC;
          // The mCRL2complexity_gj call above assigns work to every transition in ind, so we are allowed to mark all transitions at once.
        }

        // ---------------------------------------------------------------------------------------------
        // First carry out a co-split of B with respect to C\B and an action tau.
        if (m_branching)
        {
          linked_list<BLC_indicators>::iterator tau_co_splitter = find_inert_co_transition_for_block(index_block_B, old_constellation);

          // Algorithm 1, line 1.19.
          if (m_blocks[index_block_B].block_to_constellation.end() != tau_co_splitter)
          {
            tau_co_splitter->start_marked_BLC = tau_co_splitter->start_same_BLC;
            // We have to give credit for marking all transitions in the splitter at once:
            mCRL2complexity_gj(tau_co_splitter, add_work(check_complexity::refine_partition_until_it_becomes_stable__prepare_cosplit, max_C), *this);
            // The routine below has a side effect, as it sets m_R for all bottom states of block B.
            const BLC_list_iterator splitter_end_unmarked_BLC = tau_co_splitter->start_marked_BLC;
            std::vector<state_index>::iterator first_unmarked_bottom_state = not_all_bottom_states_are_touched(tau_co_splitter
                            #ifndef NDEBUG
                              , splitter_end_unmarked_BLC
                            #endif
                            );
            if (first_unmarked_bottom_state < m_blocks[index_block_B].start_non_bottom_states)
            {
            // Algorithm 1, line 1.10.

// std::cerr << "DO A TAU CO SPLIT " << old_constellation << "\n";
              splitB(tau_co_splitter, first_unmarked_bottom_state,
                        splitter_end_unmarked_BLC,
                        [&block_label_to_cotransition, ci, this]
                          (const block_index old_block,
                           const block_index new_block,
                           const transition_index moved_transition,
                           const transition_index alternative_transition)
                          {
                            maintain_block_label_to_cotransition(
                                    old_block,
                                    new_block,
                                    moved_transition,
                                    alternative_transition, block_label_to_cotransition,
                                    ci);
                          });
            }
            else
            {
              // not_all_bottom_states_are_touched() should have cleaned up
              assert(!has_marked_transitions(*tau_co_splitter));
            }
          }
        }
        // Algorithm 1, line 1.10.
        for (std::pair<BLC_list_iterator, BLC_list_iterator> calM_elt: calM)
        {
          // mCRL2complexity_gj(..., add_work(..., max_C), *this);
              // not needed as the inner loop is always executed at least once.
          //print_data_structures("Main loop");
          assert(check_data_structures("Main loop", false, false));
          check_stability("Main loop", &calM, &calM_elt, &block_label_to_cotransition);
          // Algorithm 1, line 1.11.
          assert(calM_elt.first < calM_elt.second);
          do
          {
            linked_list<BLC_indicators>::iterator splitter = m_transitions[*std::prev(calM_elt.second)].transitions_per_block_to_constellation;
            // mCRL2complexity_gj(splitter, add_work(..., max_C), *this);
                // not needed, as all transitions in calM are transitions into the small new constellation.
            assert(!points_into_BLC_set(calM_elt.second, *splitter));
            assert(splitter->start_same_BLC < calM_elt.second);
            assert(points_into_BLC_set(std::prev(calM_elt.second), *splitter));
            if (splitter->start_same_BLC != splitter->start_marked_BLC)
            {
              assert(splitter->start_same_BLC <= splitter->start_marked_BLC);
              assert(points_into_BLC_set(std::prev(splitter->start_marked_BLC), *splitter));
            }
            calM_elt.second = splitter->start_same_BLC;

            const transition& first_t = m_aut.get_transitions()[*splitter->start_same_BLC];
            const label_index a = label_or_divergence(first_t);
            assert(m_blocks[m_states[first_t.to()].block].constellation == new_constellation);
//std::cerr << "INVESTIGATE ACTION " << (m_aut.num_action_labels() == a ? "(tau-self-loops)" : m_aut.action_label(a)) << " source block " << bi << " target block " << index_block_B << "\n";
            block_index Bpp = m_states[first_t.from()].block;
            if (m_blocks[Bpp].contains_new_bottom_states ||
                (is_inert_during_init(first_t) && m_blocks[Bpp].constellation == new_constellation))
            {
              // The block Bpp contains new bottom states, and it is not necessary to spend any work on it now.
              // We will later stabilize it in stabilizeB().

              // Or the BLC slice only contains constellation-inert transitions,
              // and we do not need a main split. (The co-split has already been
              // executed earlier.)
              assert(!has_marked_transitions(*splitter));

              #ifndef NDEBUG
                // ensure that co-splitter transitions are unmarked:
                typename block_label_to_size_t_map::const_iterator bltc_it=block_label_to_cotransition.find(std::pair(Bpp,a));
                if (bltc_it!=block_label_to_cotransition.end() &&
                    bltc_it->second!=null_transition)
                {
                  const transition& co_t = m_aut.get_transitions()[bltc_it->second];
                  assert(m_states[co_t.from()].block == Bpp);
                  assert(m_blocks[m_states[co_t.to()].block].constellation == old_constellation);
                  assert(!(is_inert_during_init(co_t) && m_blocks[Bpp].constellation == old_constellation));
                  linked_list<BLC_indicators>::iterator co_splitter = m_transitions[bltc_it->second].transitions_per_block_to_constellation;
                  assert(m_states[m_aut.get_transitions()[*co_splitter->start_same_BLC].from()].block == Bpp);
                  assert(!has_marked_transitions(*co_splitter));
                }
              #endif
            }
            else
            {
              assert(has_marked_transitions(*splitter));
              assert(splitter->start_marked_BLC == splitter->start_same_BLC);
              // Check whether the bottom states of Bpp are not all included in Mleft.
              const BLC_list_iterator splitter_end_unmarked_BLC = splitter->start_marked_BLC;
              std::vector<state_index>::iterator first_unmarked_bottom_state = not_all_bottom_states_are_touched(splitter
                        #ifndef NDEBUG
                          , splitter_end_unmarked_BLC
                        #endif
                        );
              if (first_unmarked_bottom_state < m_blocks[Bpp].start_non_bottom_states)
              {
// std::cerr << "PERFORM A MAIN SPLIT \n";
              // Algorithm 1, line 1.12.
                // std::size_t dummy_number=0;
                Bpp = splitB(splitter,
                                        first_unmarked_bottom_state,
                                        splitter_end_unmarked_BLC,
                                        [&block_label_to_cotransition, ci, this]
                                          (const block_index old_block,
                                           const block_index new_block,
                                           const transition_index moved_transition,
                                           const transition_index alternative_transition)
                                          {
                                            maintain_block_label_to_cotransition(
                                                    old_block,
                                                    new_block,
                                                    moved_transition,
                                                    alternative_transition,
                                                    block_label_to_cotransition,
                                                    ci);
                                          });
                if (null_block == Bpp)
                {
                  // all bottom states in the R-subblock are new bottom states.
                  // Then no co-split is needed to stabilize the old bottom states in the R-block!
                  continue;
                }
                assert(0 <= Bpp); assert(Bpp < m_blocks.size());
                // Algorithm 1, line 1.13.
                // Algorithm 1, line 1.14 is implicitly done in the call of splitB above.
              }
              else
              {
                // not_all_bottom_states_are_marked() should have cleaned up
                assert(!has_marked_transitions(*splitter));
              }
              // Algorithm 1, line 1.17 and line 1.18.
// mCRL2log(log::debug) << "BLOCK THAT IS SPLITTER " << Bpp << "\n";
              typename block_label_to_size_t_map::const_iterator bltc_it=block_label_to_cotransition.find(std::pair(Bpp,a));
              // std::vector<state_index>::iterator first_unmarked_bottom_state;

              if (bltc_it!=block_label_to_cotransition.end() &&
                  bltc_it->second!=null_transition)
              {
// mCRL2log(log::debug) << "CO-TRANSITION  " << ptr(bltc_it->second) << "\n";
#ifndef NDEBUG
                const transition& co_t = m_aut.get_transitions()[bltc_it->second];
                assert(m_states[co_t.from()].block == Bpp);
                assert(m_blocks[m_states[co_t.to()].block].constellation == old_constellation);
                assert(!(is_inert_during_init(co_t) && m_blocks[Bpp].constellation == old_constellation));
#endif
                linked_list<BLC_indicators>::iterator co_splitter = m_transitions[bltc_it->second].transitions_per_block_to_constellation;
                assert(m_states[m_aut.get_transitions()[*co_splitter->start_same_BLC].from()].block == Bpp);
                // Algorithm 1, line 1.19.

                const BLC_list_iterator splitter_end_unmarked_BLC = co_splitter->start_marked_BLC;
                first_unmarked_bottom_state = not_all_bottom_states_are_touched(co_splitter
                            #ifndef NDEBUG
                              , splitter_end_unmarked_BLC
                            #endif
                            );
                if (first_unmarked_bottom_state < m_blocks[Bpp].start_non_bottom_states)
                {
// std::cerr << "PERFORM A MAIN CO-SPLIT \n";
                  splitB(co_splitter,
                        first_unmarked_bottom_state,
                        splitter_end_unmarked_BLC,
                        [&block_label_to_cotransition, ci, this]
                          (const block_index old_block,
                           const block_index new_block,
                           const transition_index moved_transition,
                           const transition_index alternative_transition)
                          {
                            maintain_block_label_to_cotransition(
                                    old_block,
                                    new_block,
                                    moved_transition,
                                    alternative_transition, block_label_to_cotransition,
                                    ci);
                          });
                }
                // Algorithm 1, line 1.20 and 1.21. P is updated implicitly when splitting Bpp.

              }
            }
          }
          while (calM_elt.first < calM_elt.second);
        }

// print_data_structures("Before stabilize");
        assert(check_data_structures("Before stabilize", false, false));
        check_stability("Before stabilize");
        stabilizeB();
      }
    }
};

/* ************************************************************************* */
/*                                                                           */
/*                             I N T E R F A C E                             */
/*                                                                           */
/* ************************************************************************* */


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
    // Line 1.2: Find tau-SCCs and contract each of them to a single state
mCRL2log(log::verbose) << "Start SCC\n";
    mcrl2::utilities::execution_timer timer;
    if (branching)
    {
        timer.start("scc_reduce");
        scc_reduce(l, preserve_divergence);
        timer.finish("scc_reduce");
    }

    // Now apply the branching bisimulation reduction algorithm.  If there
    // are no taus, this will automatically yield strong bisimulation.
mCRL2log(log::verbose) << "Start Partitioning\n";
    timer.start("reduction");
    bisim_partitioner_gj<LTS_TYPE> bisim_part(l, branching, preserve_divergence);
    timer.finish("reduction");

    // Assign the reduced LTS
    timer.start("quotient");
mCRL2log(log::verbose) << "Start finalizing\n";
    bisim_part.finalize_minimized_LTS();
    timer.finish("quotient");
    timer.report();
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
        mCRL2log(log::warning) << "The GJ24 branching bisimulation "
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
    bisim_partitioner_gj<LTS_TYPE> bisim_part(l1, branching, preserve_divergence);

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


} // end namespace detail
} // end namespace lts
} // end namespace mcrl2

#endif // ifndef LIBLTS_BISIM_GJ_H
