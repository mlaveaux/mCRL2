// Author(s): Muck van Weerdenburg, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

// This file contains the merge algorithm that merges two lts's.
// Merges an LTS L with this LTS (say K) and stores the resulting LTS
// (say M) in this LTS datastructure, effectively replacing K.
// Conceptually, we just take the union of the sets of states and the
// sets of transitions of K and L:
//   States_M      = States_K + States_L
//   Transitions_M = Transitions_K + Transitions_L
// where + denotes set union.
// However, this assumes that States_K and States_L are disjoint,
// which is generally not the case. More specifically we have:
//   States_K = { 0, ..., N_K - 1 }   and
//   States_L = { 0, ..., N_L - 1 }
// for some N_K, N_L > 0.
// Therefore, state i of L will be numbered |N_K| + i in the resulting
// LTS M and state i of K will be numbered i in M. This yields:
//   States_M = { 0, ..., N_K + N_L - 1 }.


#ifndef MCRL2_LTS_LIBLTS_MERGE_H
#define MCRL2_LTS_LIBLTS_MERGE_H

#include "mcrl2/lts/lts.h"

#include <map>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::lts::detail
{

template <class LTS_TYPE>
void merge(LTS_TYPE& l1, const LTS_TYPE& l2)
{
  const std::size_t old_nstates=l1.num_states();
  l1.set_num_states(l1.num_states() + l2.num_states());

  // The resulting LTS will have state information only if BOTH LTSs
  // currently have state information.
  if (l1.has_state_info() && l2.has_state_info())
  {
    for (std::size_t i=0; i<l2.num_states(); ++i)
    {
      l1.add_state(l2.state_label(i));
    }
  }
  else
  {
    // remove state information from this LTS, if any
    l1.clear_state_labels();
  }

  // Before we can set the label data in a new transitions
  // array, we first have to collect the labels of both LTSs in a
  // map, of which the second element indicates the new index of each action label.

  using action_label_type = typename LTS_TYPE::action_label_t;
  using label_index = typename LTS_TYPE::labels_size_type;
  using insert_type = typename std::pair<typename std::map<action_label_type, label_index>::const_iterator, bool>;
  std::map < action_label_type,label_index > labs;

  // Put the labels of the LTS l1 in a map.
  for (std::size_t i = 0; i < l1.num_action_labels(); ++i)
  {
    labs.insert(std::pair <typename LTS_TYPE::action_label_t,typename LTS_TYPE::labels_size_type>
                (l1.action_label(i),i));
  }
  // Add the labels for the LTS l2, and put them there with a new index if it was
  // not added yet.
  // Furthermore, update the hidden_label_set. 
  // If label a1 is mapped to a2 in l2, then this must be the same
  // in l1. It may be that label a1 did not exist yet in which case it needs
  // to be added too.

  for (std::size_t i=0; i<l2.num_action_labels(); ++i)
  {
    typename LTS_TYPE::labels_size_type new_index;
    const insert_type it=labs.insert(std::pair < action_label_type,label_index >
                               (l2.action_label(i),l1.num_action_labels()));
    if (it.second) 
    {
      // New element has been inserted.
      new_index=l1.add_action(l2.action_label(i));
      if (l2.is_tau(l2.apply_hidden_label_map(i)))
      {
        l1.hidden_label_set().insert(new_index);
      }

    }
    else 
    {
      new_index=it.first->second; // Old index to which i is mapped.

      // If label i occurred in l1 and both were not mapped to the hidden label, raise an exception.
      if (l1.is_tau(l1.apply_hidden_label_map(new_index)) != l2.is_tau(l2.apply_hidden_label_map(i)))
      {
        throw mcrl2::runtime_error("The action " + pp(l2.action_label(i)) + " has incompatible hidden actions " +
                                       pp(l1.action_label(l1.apply_hidden_label_map(new_index))) + " and " + 
                                       pp(l2.action_label(l2.apply_hidden_label_map(i))) + ".");
      }
    }
    assert(new_index==it.first->second);
  }

  // Update the label numbers of all transitions of the LTS l1 to reflect
  // the new indices as given by labs.
  std::vector<transition>& trans1=l1.get_transitions();
  for (std::vector<transition>::iterator r=trans1.begin(); r!=trans1.end(); ++r)
  {
    r->set_label(labs[l1.action_label(r->label())]);
  }

  // Now add the transition labels of LTS l2
  // Now add the source and target states of the transitions of LTS l2.
  // The labels will be added below, depending on whether there is label
  // information in both LTSs.
  const std::vector<transition>& trans2=l2.get_transitions();
  for (std::vector<transition>::const_iterator r=trans2.begin(); r!=trans2.end(); ++r)
  {
    const transition transition_to_add=*r;
    l1.add_transition(transition(transition_to_add.from()+old_nstates,
                                 labs[l2.action_label(transition_to_add.label())],
                                 transition_to_add.to()+old_nstates));
  }

}
} // namespace mcrl2::lts::detail

#endif
