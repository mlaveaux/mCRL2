// Author(s): Muck van Weerdenburg, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts/detail/liblts_scc.h

#ifndef _LIBLTS_SCC_H
#define _LIBLTS_SCC_H

#include <unordered_set>

#include "mcrl2/lts/lts.h"
#include "mcrl2/utilities/logger.h"

#ifdef MCRL2_ENABLE_MODULES
    import utilities;
#endif

namespace mcrl2::lts
{

namespace detail
{
  // An indexed sorted vector below contains the outgoing or incoming tau transitions per state,
  // grouped per state. The input consists of an automaton with transitions. The incoming/outcoming
  // tau transitions are grouped by state in the m_states_with_outgoing_or_incoming_tau_transition. 
  // It is as long as the lts aut has tau transitions. The vector m_indices is as long as the number
  // of states. For each state it contains the place in the other vector where its tau transitions
  // start. So, the tau transitions reside at position indices[s] to indices[s+1]. These indices
  // can be acquired using the functions lowerbound and upperbound. 
  // This data structure is chosen due to its minimal memory and time footprint. 
  template <class LTS_TYPE>
  class indexed_sorted_vector_for_tau_transitions
  {
    protected:
      using state_type = std::size_t;
      using label_type = std::size_t;

      std::vector <state_type> m_states_with_outgoing_or_incoming_tau_transition;
      std::vector <size_t> m_indices;

    public:

      indexed_sorted_vector_for_tau_transitions(const LTS_TYPE& aut, bool outgoing)
       : m_indices(aut.num_states()+1,0)
      {
        // First count the number of outgoing transitions per state and put it in indices.
        for(const transition& t: aut.get_transitions())
        {
          if (aut.is_tau(aut.apply_hidden_label_map(t.label())))
          {
            m_indices[outgoing?t.from():t.to()]++;
          }
        }

        // Calculate the m_indices where the states with outgoing/incoming tau transition must be placed.
        // Put the starting index for state i at position i-1. When placing the transitions these indices
        // are decremented properly. 
        
        size_t sum=0;
        for(state_type& i: m_indices)  // The vector is changed. This must be a reference. 
        {
          sum=sum+i;
          i=sum;
        }

        // Now declare enough space for all transitions and store them in reverse order, while
        // at the same time decrementing the indices in m_indices. 
        m_states_with_outgoing_or_incoming_tau_transition.resize(sum);
        for(const transition& t: aut.get_transitions())
        {
          if (aut.is_tau(aut.apply_hidden_label_map(t.label())))
          {
            if (outgoing)
            {
              assert(t.from()<m_indices.size());
              assert(m_indices[t.from()]>0);
              m_indices[t.from()]--;
              assert(m_indices[t.from()] < m_states_with_outgoing_or_incoming_tau_transition.size());
              m_states_with_outgoing_or_incoming_tau_transition[m_indices[t.from()]]=t.to();
            }
            else
            {
              assert(t.to()<m_indices.size());
              assert(m_indices[t.to()]>0);
              m_indices[t.to()]--;
              assert(m_indices[t.to()] < m_states_with_outgoing_or_incoming_tau_transition.size());
              m_states_with_outgoing_or_incoming_tau_transition[m_indices[t.to()]]=t.from();
            }
          }
        }
        assert(m_indices.at(aut.num_states())==m_states_with_outgoing_or_incoming_tau_transition.size());
      }

      // Get the indexed transitions. 
      const std::vector<state_type>& get_transitions() const
      {
        return m_states_with_outgoing_or_incoming_tau_transition;
      }
    
      // Get the lowest index of incoming/outging transitions stored in m_states_with_outgoing_or_incoming_tau_transition.
      size_t lowerbound(const state_type s) const
      {
        assert(s+1<m_indices.size());
        return m_indices[s];
      }

      // Get 1 beyond the higest index of incoming/outging transitions stored in m_states_with_outgoing_or_incoming_tau_transition.
      size_t upperbound(const state_type s) const
      {
        assert(s+1<m_indices.size());
        return m_indices[s+1];
      }

      // Drastically clear the vectors by resetting its memory usage to minimal. 
      void clear()   
      {
        std::vector <state_type>().swap(m_states_with_outgoing_or_incoming_tau_transition);
        std::vector <size_t>().swap(m_indices);
        
      }
  };

/// \brief This class contains an scc partitioner removing inert tau loops.

template < class LTS_TYPE>
class scc_partitioner
{

  public:
    /** \brief Creates an scc partitioner for an LTS.
     *  \details This scc partitioner calculates a partition
     *  of the state space of the transition system l using
     *  a straightforward algorithm found in A.V. Aho, J.E. Hopcroft
     *  and J.D. Ullman, Data structures and algorithms. Addison Wesley,
     *  1987 on page 224. All states that reside on a loop of internal
     *  actions are put in the same equivalence class. The function l.is_tau
     *  is used to determine whether an action is internal. Partitioning is
     *  done immediately when an instance of this class is created.
     *  When applying the function \ref replace_transition_system the
     *  automaton l is replaced by (aka shrinked to) the automaton modulo the
     *  calculated partition.
     *  \param[in] l reference to an LTS. */
    scc_partitioner(LTS_TYPE& l);

    /** \brief Destroys this partitioner. */
    ~scc_partitioner()=default;

    /** \brief The lts for which this partioner is created is replaced by the lts modulo
     *        the calculated partition.
     *  \details The number of states of the new lts becomes the number of
     *        equivalence classes. In each transition the start and end state
     *        are replaced by its equivalence class. If a transition has a tau
     *        label (which is checked using the function l.is_tau) then it
     *        is preserved if either from and to state are different, or
     *        \e preserve_divergence_loops is true. All non tau transitions are
     *        always preserved. The label numbers for preserved transitions are
     *        not changed. 
     *
     * \param[in] preserve_divergence_loops If true preserve a tau loop on states that
     *     were part of a larger tau loop in the input transition system. Otherwise idle
     *     tau loops are removed. */
    void replace_transition_system(const bool preserve_divergence_loops);

    /** \brief Gives the number of bisimulation equivalence classes of the LTS.
     *  \return The number of bisimulation equivalence classes of the LTS.
     */
    std::size_t num_eq_classes() const;

    /** \brief Gives the equivalence class number of a state.
     *  The equivalence class numbers range from 0 upto (and excluding)
     *  \ref num_eq_classes().
     *  \param[in] s A state number.
     *  \return The number of the equivalence class to which \e s
     *    belongs. */
    std::size_t get_eq_class(const std::size_t s) const;

    /** \brief Returns whether two states are in the same bisimulation
     *     equivalence class.
     * \param[in] s A state number.
     * \param[in] t A state number.
     * \retval true if \e s and \e t are in the same bisimulation
     *    equivalence class;
     * \retval false otherwise. */
    bool in_same_class(const std::size_t s, const std::size_t t) const;

  private:

    using state_type = std::size_t;
    using label_type = std::size_t;

    LTS_TYPE& aut;

    std::vector < state_type > block_index_of_a_state;
    std::vector < state_type > dfsn2state;
    state_type equivalence_class_index = 0;

    void group_components(const state_type t,
                          const state_type equivalence_class_index,
                          const indexed_sorted_vector_for_tau_transitions<LTS_TYPE>& src_tgt_src,
                          std::vector < bool >& visited);
    void dfs_numbering(const state_type t,
                       const indexed_sorted_vector_for_tau_transitions<LTS_TYPE>& src_tgt,
                       std::vector < bool >& visited);

};


template < class LTS_TYPE>
scc_partitioner<LTS_TYPE>::scc_partitioner(LTS_TYPE& l)
  :aut(l),
    block_index_of_a_state(aut.num_states(),0)
{
  mCRL2log(log::debug) << "Tau loop (SCC) partitioner created for " << l.num_states() << " states and " <<
              l.num_transitions() << " transitions" << std::endl;

  dfsn2state.reserve(aut.num_states());

  // Initialise the data structures used in the recursive DFS procedure.
  std::vector<bool> visited(aut.num_states(),false); 
  indexed_sorted_vector_for_tau_transitions<LTS_TYPE> src_tgt(aut,true); // Group the tau transitions ordered per outgoing states. 

  // Number the states via a depth first search
  for (state_type i=0; i<aut.num_states(); ++i)
  {
    dfs_numbering(i,src_tgt,visited);
  }
  src_tgt.clear();

  indexed_sorted_vector_for_tau_transitions<LTS_TYPE> tgt_src(aut,false);
  for (std::vector < state_type >::reverse_iterator i=dfsn2state.rbegin();
       i!=dfsn2state.rend(); ++i)
  {
    if (visited[*i])  // Visited is used inversely here.
    {
      group_components(*i,equivalence_class_index,tgt_src,visited);
      equivalence_class_index++;
    }
  }
  mCRL2log(log::debug) << "Tau loop (SCC) partitioner reduces lts to " << equivalence_class_index << " states." << std::endl;

  dfsn2state.clear();
}


template < class LTS_TYPE>
void scc_partitioner<LTS_TYPE>::replace_transition_system(const bool preserve_divergence_loops)
{
  // Put all the non inert transitions in a set. Add the transitions that form a self
  // loop. Such transitions only exist in case divergence preserving branching bisimulation is
  // used. A set is used to remove double occurrences of transitions.
  /* std::unordered_set < transition > resulting_transitions;
  for (const transition& t: aut.get_transitions())
  {
    if (!aut.is_tau(aut.apply_hidden_label_map(t.label())) ||
        preserve_divergence_loops ||
        block_index_of_a_state[t.from()]!=block_index_of_a_state[t.to()])
    {
      resulting_transitions.insert(    // insert is faster than emplace. 
        transition(
          block_index_of_a_state[t.from()],
          aut.apply_hidden_label_map(t.label()),
          block_index_of_a_state[t.to()])); 
    }
  }

  aut.clear_transitions(resulting_transitions.size());
  // Copy the transitions from the set into the transition system.

  for (const transition& t: resulting_transitions)
  {
    aut.add_transition(t);
  } */
  
  for (transition& t: aut.get_transitions())
  { 
    t.set_from(block_index_of_a_state[t.from()]);
    t.set_label(aut.apply_hidden_label_map(t.label()));
    t.set_to(block_index_of_a_state[t.to()]);
  }
  
  sort_transitions(aut.get_transitions(),std::set<transition::size_type>(),tgt_lbl_src);
  
  // Compress the transitions while removing double occurrences, and if needed self-loops. 
  constexpr std::size_t non_existent=-1;
  std::size_t old_index=non_existent;
  for (const transition& t: aut.get_transitions())
  {
    if (!aut.is_tau(t.label()) ||
        preserve_divergence_loops ||
        t.from()!=t.to())
    {
      if (old_index==non_existent)
      {
        old_index=0;
        aut.get_transitions()[old_index]=t;
      }
      else if (t!=aut.get_transitions()[old_index])
      {
        old_index++;
        aut.get_transitions()[old_index]=t;
      }
    }
  }
  aut.set_num_transitions(old_index+1);

  // Merge the states, by setting the state labels of each state to the concatenation of the state labels of its
  // equivalence class. 
  if (aut.has_state_info())   /* If there are no state labels this step can be ignored */
  {
    /* Create a vector for the new labels */
    std::vector<typename LTS_TYPE::state_label_t> new_labels(num_eq_classes());

    for(std::size_t i=aut.num_states(); i>0; )
    {
      --i;
      const std::size_t new_index=block_index_of_a_state[i];
      new_labels[new_index]=new_labels[new_index]+aut.state_label(i);
    }

    for(std::size_t i=0; i<num_eq_classes(); ++i)
    {
      aut.set_state_label(i,new_labels[i]);
    }
  }
  
  aut.set_num_states(num_eq_classes()); 
  aut.set_initial_state(get_eq_class(aut.initial_state()));
}

template < class LTS_TYPE>
std::size_t scc_partitioner<LTS_TYPE>::num_eq_classes() const
{
  return equivalence_class_index;
}

template < class LTS_TYPE>
std::size_t scc_partitioner<LTS_TYPE>::get_eq_class(const std::size_t s) const
{
  return block_index_of_a_state[s];
}

template < class LTS_TYPE>
bool scc_partitioner<LTS_TYPE>::in_same_class(const std::size_t s, const std::size_t t) const
{
  return get_eq_class(s)==get_eq_class(t);
}

// Private methods of scc_partitioner

template < class LTS_TYPE>
void scc_partitioner<LTS_TYPE>::group_components(
  const state_type s,
  const state_type equivalence_class_index,
  const indexed_sorted_vector_for_tau_transitions<LTS_TYPE>& tgt_src,
  std::vector < bool >& visited)
{
  if (!visited[s])
  {
    return;
  }
  visited[s] = false;
  const size_t u=tgt_src.upperbound(s);  // only calculate the upperbound once. 
  for(state_type i=tgt_src.lowerbound(s); i<u; ++i)
  {
    group_components(tgt_src.get_transitions()[i],equivalence_class_index,tgt_src,visited);
  }
  block_index_of_a_state[s]=equivalence_class_index;
}

template < class LTS_TYPE>
void scc_partitioner<LTS_TYPE>::dfs_numbering(
  const state_type s,
  const indexed_sorted_vector_for_tau_transitions<LTS_TYPE>& src_tgt,
  std::vector < bool >& visited)
{
  if (visited[s])
  {
    return;
  }
  visited[s] = true;
  const size_t u=src_tgt.upperbound(s);  // only calculate the upperbound once. 
  for(state_type i=src_tgt.lowerbound(s); i<u; ++i)
  {
    dfs_numbering(src_tgt.get_transitions()[i],src_tgt,visited);
  }
  dfsn2state.push_back(s);
}

} // namespace detail

template < class LTS_TYPE>
void scc_reduce(LTS_TYPE& l,const bool preserve_divergence_loops = false)
{
  detail::scc_partitioner<LTS_TYPE> scc_part(l);
  scc_part.replace_transition_system(preserve_divergence_loops);
}

}

#endif // _LIBLTS_SCC_H
