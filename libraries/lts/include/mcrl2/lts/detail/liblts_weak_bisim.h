// Author(s): Jan Friso Groote (upon request by Matthew Hennessy)
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts/detail/liblts_weak_bisim.h
/// \brief This file defines an algorithm for weak bisimulation, by
///        calculating the transitive tau closure and apply strong
///        bisimulation afterwards. In order to apply this algorithm
///        it is advisable to first apply a branching bisimulation reduction.

#ifndef MCRL2_LTS_WEAK_BISIM_H
#define MCRL2_LTS_WEAK_BISIM_H
#include "mcrl2/lts/detail/liblts_scc.h"
#include "mcrl2/lts/detail/liblts_tau_star_reduce.h"
#include "mcrl2/lts/detail/liblts_merge.h"
#include "mcrl2/lts/lts_aut.h"
#include "mcrl2/lts/lts_fsm.h"
#include "mcrl2/lts/lts_dot.h"

namespace mcrl2::lts::detail
{

/** \brief Reduce LTS l with respect to (divergence-preserving) weak bisimulation.
 * \param[in,out] l The transition system that is reduced.
 * \param[in] preserve_divergences Indicates whether loops of internal actions on states must be preserved. If false
 *            these are removed. If true these are preserved.  */
template < class LTS_TYPE>
void weak_bisimulation_reduce(
  LTS_TYPE& l,
  const bool preserve_divergences = false)
{
  if (1 < l.num_states())
  {
    bisimulation_reduce_dnj(l, true, preserve_divergences);   //< Apply branching bisimulation to l.
  }

  std::size_t divergence_label;
  if (preserve_divergences)
  {
    divergence_label=mark_explicit_divergence_transitions(l);
  }
  if (1 < l.num_states())
  {
    reflexive_transitive_tau_closure(l);                      // Apply transitive tau closure to l.
    bisimulation_reduce_dnj(l, false, false);                 // Apply strong bisimulation to l.
  }
  scc_reduce(l);                                              // Remove tau loops.
  remove_redundant_transitions(l);                            // Remove transitions s -a-> s' if also s-a->-tau->s' or s-tau->-a->s' is present.
                                                              // Note that this is correct, because l is reduced modulo strong bisimulation and
                                                              // does not contain tau loops.
  if (preserve_divergences)
  {
    unmark_explicit_divergence_transitions(l,divergence_label);
  }
}


/** \brief Checks whether the initial states of two LTSs are weakly bisimilar.
 * \details The LTSs l1 and l2 are not usable anymore after this call.
 *          The space consumption is O(n) and running time is dominated by the
 *          transitive closure (after branching bisimulation).
 * \param[in,out] l1 A first transition system.
 * \param[in,out] l2 A second transistion system.
 * \param[preserve_divergences] If true and branching is true, preserve tau loops on states.
 * \retval True iff the initial states of the current transition system and l2 are (divergence preserving) (branching) bisimilar */
template < class LTS_TYPE>
bool destructive_weak_bisimulation_compare(
  LTS_TYPE& l1,
  LTS_TYPE& l2,
  const bool preserve_divergences=false)
{
  weak_bisimulation_reduce(l1,preserve_divergences);
  weak_bisimulation_reduce(l2,preserve_divergences);
  return destructive_bisimulation_compare_dnj(l1,l2);
}


/** \brief Checks whether the initial states of two LTSs are weakly bisimilar.
 *  \details The LTSs l1 and l2 are first duplicated and subsequently
 *           reduced modulo bisimulation. If memory space is a concern, one could consider to
 *           use destructive_weak_bisimulation_compare.  The running time
 *           of this routine is dominated by the transitive closure
 *           (after branching bisimulation).  It uses O(m+n) memory
 *           in addition to the copies of l1 and l2, where n is the
 *           number of states and m is the number of transitions.
 * \param[in,out] l1 A first transition system.
 * \param[in,out] l2 A second transistion system.
 * \param[preserve_divergences] If true and branching is true, preserve tau loops on states.
 * \retval True iff the initial states of the current transition system and l2 are (divergence preserving) (branching) bisimilar */
template < class LTS_TYPE>
bool weak_bisimulation_compare(
  const LTS_TYPE& l1,
  const LTS_TYPE& l2,
  const bool preserve_divergences = false)
{
  LTS_TYPE l1_copy(l1);
  LTS_TYPE l2_copy(l2);
  return destructive_weak_bisimulation_compare(l1_copy, l2_copy,
                                                         preserve_divergences);
}

}

#endif  // MCRL2_LTS_WEAK_BISIM_H
