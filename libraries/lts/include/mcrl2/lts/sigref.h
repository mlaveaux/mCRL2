// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file sigref.h
/// \brief Implementation of LTS reductions using the signature refinement
///        approach of S. Blom and S. Orzan.
/// \author Jeroen Keiren

#ifndef MCRL2_LTS_SIGREF_H
#define MCRL2_LTS_SIGREF_H

#include "mcrl2/lts/lts_utilities.h"
#include "mcrl2/utilities/execution_timer.h"

namespace mcrl2
{
namespace lts
{

/** \brief A signature is a pair of an action label and a block */
typedef std::set<std::pair<std::size_t, std::size_t> > signature_t;

/** \brief Base class for signature computation */
template < class LTS_T >
class signature
{
protected:
  /** \brief The labelled transition system for which the signature is computed */
  const LTS_T& m_lts;

  /** \brief Signature stored per state */
  std::vector<signature_t> m_sig;

public:
  /** \brief Constructor
    */
  signature(const LTS_T& lts_)
    : m_lts(lts_), m_sig(m_lts.num_states(), signature_t())
  {}

  /** \brief Compute a new signature based on \a partition.
    * \param[in] partition The current partition
    */
  virtual void compute_signature(const std::vector<std::size_t>& partition) = 0;

  /** \brief Compute the transitions for the quotient according to \a partition.
    * \param[in] partition The partition that is used to compute the quotient
    * \param[out] transitions A set to which the transitions of the quotient are written
    */
  virtual void quotient_transitions(std::set<transition>& transitions, const std::vector<std::size_t>& partition)
  {
    for(std::vector<transition>::const_iterator i = m_lts.get_transitions().begin(); i != m_lts.get_transitions().end(); ++i)
    {
      transitions.insert(transition(partition[i->from()], i->label(), partition[i->to()]));
    }
  }

  /** \brief Return the signature for state \a i.
    * \param[in] i The state for which to return the signature.
    * \pre i < m_lts.num_states().
    */
  virtual const signature_t& get_signature(std::size_t i) const
  {
    return m_sig[i];
  }
};

/** \brief Class for computing the signature for strong bisimulation */
template < class LTS_T >
class signature_bisim: public signature<LTS_T>
{
protected:
  using signature<LTS_T>::m_lts;
  using signature<LTS_T>::m_sig;

public:
  /** \brief Constructor */
  signature_bisim(const LTS_T& lts_)
    : signature<LTS_T>(lts_)
  {
    mCRL2log(log::verbose) << "initialising signature computation for strong bisimulation" << std::endl;
  }

  /** \overload */
  virtual void
  compute_signature(const std::vector<std::size_t>& partition)
  {
    // compute signatures
    m_sig = std::vector<signature_t>(m_lts.num_states(), signature_t());
    for(std::vector<transition>::const_iterator i = m_lts.get_transitions().begin(); i != m_lts.get_transitions().end(); ++i)
    {
      m_sig[i->from()].insert(std::make_pair(m_lts.apply_hidden_label_map(i->label()), partition[i->to()]));
    }
  }

};

/** \brief Class for computing the signature for branching bisimulation */
template < class LTS_T >
class signature_branching_bisim: public signature<LTS_T>
{
protected:
  using signature<LTS_T>::m_lts;
  using signature<LTS_T>::m_sig;

  /** \brief Store the incoming transitions per state */
  outgoing_transitions_per_state_t m_prev_transitions;

  /** \brief Insert function
    * \param[in] partition The current partition
    * \param[in] t source state
    * \param[in] label_ transition label
    * \param[in] block target block
    *
    * Insert function as described in S. Blom, S. Orzan,
    * "Distributed Branching Bisimulation Reduction of State Spaces",
    * Proc. PDMC 2003.
    *
    * Inserts the pair (label_, block) in the signature of t, as well as
    * the signatures of all tau-predecessors of t within the same block.
    */
  void insert(const std::vector<std::size_t>& partition, const std::size_t t, const std::size_t label_, const std::size_t block)
  {
    if(m_sig[t].insert(std::make_pair(label_, block)).second)
    {
      // std::pair<outgoing_transitions_per_state_t::const_iterator, outgoing_transitions_per_state_t::const_iterator> pred_range
      //    = m_prev_transitions.equal_range(t);

      // for(outgoing_transitions_per_state_t::const_iterator i = pred_range.first; i != pred_range.second; ++i)
      // for(const outgoing_pair_t& p: m_prev_transitions[t])
      for (std::size_t i=m_prev_transitions.lowerbound(t); i<m_prev_transitions.upperbound(t); ++i)
      {
        const outgoing_pair_t& p = m_prev_transitions.get_transitions()[i];
        if(m_lts.is_tau(m_lts.apply_hidden_label_map(label(p))) && partition[t] == partition[to(p)])
        {
          insert(partition, to(p), label_, block);
        }
      }
    }
  }

public:
  /** \brief Constructor  */
  signature_branching_bisim(const LTS_T& lts_)
    : signature<LTS_T>(lts_),
      m_prev_transitions(lts_.get_transitions(),lts_.num_states(),false)  // transitions stored backward. 
  {
    mCRL2log(log::verbose) << "initialising signature computation for branching bisimulation" << std::endl;
  }

  /** \overload */
  virtual void compute_signature(const std::vector<std::size_t>& partition)
  {
    // compute signatures
    m_sig = std::vector<signature_t>(m_lts.num_states(), signature_t());
    for(std::vector<transition>::const_iterator i = m_lts.get_transitions().begin(); i != m_lts.get_transitions().end(); ++i)
    {
      if (!(m_lts.is_tau(m_lts.apply_hidden_label_map(i->label())) && (partition[i->from()] == partition[i->to()])))
      {
        insert(partition, i->from(), m_lts.apply_hidden_label_map(i->label()), partition[i->to()]);
      }
    }
  }

  /** \overload */
  virtual void quotient_transitions(std::set<transition>& transitions, const std::vector<std::size_t>& partition)
  {
    for(std::vector<transition>::const_iterator i = m_lts.get_transitions().begin(); i != m_lts.get_transitions().end(); ++i)
    {
      if(partition[i->from()] != partition[i->to()] || !m_lts.is_tau(m_lts.apply_hidden_label_map(i->label())))
      {
        transitions.insert(transition(partition[i->from()], m_lts.apply_hidden_label_map(i->label()), partition[i->to()]));
      }
    }
  }
};

/** \brief Class for computing the signature for divergence preserving branching bisimulation */
template < class LTS_T >
class signature_divergence_preserving_branching_bisim: public signature_branching_bisim<LTS_T>
{
protected:
  using signature_branching_bisim<LTS_T>::m_lts;
  using signature_branching_bisim<LTS_T>::m_sig;
  using signature_branching_bisim<LTS_T>::insert;

  /** \brief Record for each vertex whether it is in a tau-scc */
  std::vector<bool> m_divergent;

  /** \brief Iterative implementation of Tarjan's SCC algorithm.
   *
   * based on an earlier implementation by Sjoerd Cranen
   *
   * For the non-trivial tau-sccs (i.e. SCCs with more than one state, or
   * single-state SCCs with a tau-loop, we set \a m_divergent to true.
   */
  void compute_tau_sccs()
  {
    std::size_t unused = 1, lastscc = 1;
    std::vector<std::size_t> scc(m_lts.num_states(), 0);
    std::vector<std::size_t> low(m_lts.num_states(), 0);
    std::stack<std::size_t> stack;
    std::stack<std::size_t> sccstack;

    // Record forward transition relation sorted by state.
    outgoing_transitions_per_state_t m_lts_succ_transitions(m_lts.get_transitions(),m_lts.num_states(),true);

    for (std::size_t i = 0; i < m_lts.num_states(); ++i)
    {
      if (scc[i] == 0)
        stack.push(i);
      while (!stack.empty())
      {
        const std::size_t vi = stack.top();

        // Outgoing transitions of vi.
        // std::pair<outgoing_transitions_per_state_t::const_iterator, outgoing_transitions_per_state_t::const_iterator> succ_range
        //  = m_lts_succ_transitions.equal_range(vi);

        if (low[vi] == 0 && scc[vi] == 0)
        {
          scc[vi] = unused;
          low[vi] = unused++;
          sccstack.push(vi);

          // for (outgoing_transitions_per_state_t::const_iterator t = succ_range.first; t != succ_range.second; ++t)
          // for (const outgoing_pair_t& t: m_lts_succ_transitions[vi]) 
          // for (const outgoing_pair_t& t: m_lts_succ_transitions[vi]) 
          for (std::size_t i=m_lts_succ_transitions.lowerbound(vi); i<m_lts_succ_transitions.upperbound(vi); ++i)
          {
            const outgoing_pair_t& t=m_lts_succ_transitions.get_transitions()[i];
            if ((low[to(t)] == 0) && (scc[to(t)] == 0) && (m_lts.is_tau(m_lts.apply_hidden_label_map(label(t)))))
            {
              stack.push(to(t));
            }
          }
        }
        else
        {
          // for (outgoing_transitions_per_state_t::const_iterator t = succ_range.first; t != succ_range.second; ++t)
          for (std::size_t i=m_lts_succ_transitions.lowerbound(vi); i<m_lts_succ_transitions.upperbound(vi); ++i)
          {
            const outgoing_pair_t& t=m_lts_succ_transitions.get_transitions()[i];
            if ((low[to(t)] != 0) && (m_lts.is_tau(m_lts.apply_hidden_label_map(label(t)))))
              low[vi] = low[vi] < low[to(t)] ? low[vi] : low[to(t)];
          }
          if (low[vi] == scc[vi])
          {
            std::size_t tos, scc_id = lastscc++;
            std::vector<std::size_t> this_scc;
            do
            {
              tos = sccstack.top();
              low[tos] = 0;
              scc[tos] = scc_id;
              sccstack.pop();
              this_scc.push_back(tos);
            }
            while (tos != vi);

            // We only consider non-trivial sccs, hence
            // if the scc consists of a single schate, check whether it has a tau-loop
            if(this_scc.size() == 1)
            {
              // for(outgoing_transitions_per_state_t::const_iterator i = succ_range.first; i != succ_range.second; ++i)
              // for (const outgoing_pair_t& i: m_lts_succ_transitions[vi]) 
              for (std::size_t i_=m_lts_succ_transitions.lowerbound(vi); i_<m_lts_succ_transitions.upperbound(vi); ++i_)
              {
                const outgoing_pair_t& i=m_lts_succ_transitions.get_transitions()[i_];
                if(vi == to(i) && m_lts.is_tau(m_lts.apply_hidden_label_map(label(i))))
                {
                  m_divergent[tos] = true;
                  break;
                }
              }
            }
            else
            {
              m_divergent[tos] = true;
            }
          }
          stack.pop();
        }
      }
    }
  }

public:
  /** \brief Constructor
    *
    * This initialises \a m_divergent to record for each vertex whether it is
    * in a tau-scc.
    */
  signature_divergence_preserving_branching_bisim(const LTS_T& lts_)
    : signature_branching_bisim<LTS_T>(lts_),
      m_divergent(lts_.num_states(), false)
  {
    mCRL2log(log::verbose) << "initialising signature computation for divergence preserving branching bisimulation" << std::endl;
    compute_tau_sccs();
  }

  /** \overload
    *
    * Compute the signature as in branching bisimulation. In addition, add the
    * (tau, B) for edges s -tau-> t for which s,t in B and m_divergent[t]
    */
  virtual void compute_signature(const std::vector<std::size_t>& partition)
  {
    // compute signatures
    m_sig = std::vector<signature_t>(m_lts.num_states(), signature_t());
    for(std::vector<transition>::const_iterator i = m_lts.get_transitions().begin(); i != m_lts.get_transitions().end(); ++i)
    {
      if(!(partition[i->from()] == partition[i->to()] && m_lts.is_tau(m_lts.apply_hidden_label_map(i->label())))
         || m_divergent[i->to()])
      {
        insert(partition, i->from(), m_lts.apply_hidden_label_map(i->label()), partition[i->to()]);
      }
    }
  }

  /** \overload */
  virtual void quotient_transitions(std::set<transition>& transitions, const std::vector<std::size_t>& partition)
  {
    for(std::vector<transition>::const_iterator i = m_lts.get_transitions().begin(); i != m_lts.get_transitions().end(); ++i)
    {
      if(!(partition[i->from()] == partition[i->to()] && m_lts.is_tau(m_lts.apply_hidden_label_map(i->label())))
         || m_sig[i->from()].find(std::make_pair(m_lts.apply_hidden_label_map(i->label()), partition[i->to()])) != m_sig[i->from()].end())
      {
        transitions.insert(transition(partition[i->from()], m_lts.apply_hidden_label_map(i->label()), partition[i->to()]));
      }
    }
  }
};


/** \brief Signature based reductions for labelled transition systems.
  *
  * The implementation is based on the description in
  * S. Blom, S. Orzan. "Distributed Branching Bisimulation Reduction of State
  * Spaces", in Proc. PDMC 2003.
  *
  * The specific signature is a parameter of the algorithm.
  */
template < class LTS_T, typename Signature >
class sigref
{

protected:
  /** \brief Current partition; for each state (std::size_t) the block in which
             it resides is recorded. */
  std::vector<std::size_t> m_partition;

  /** \brief The number of blocks in the current partition */
  std::size_t m_count;

  /** \brief The LTS that we are reducing */
  LTS_T& m_lts;

  /** \brief Instance of a class performing the signature computation for the
             current equivalence */
  Signature m_signature;

  /** \brief Print a signature (for debugging purposes) */
  std::string print_sig(const signature_t& sig)
  {
    std::stringstream os;
    os << "{ ";
    for(signature_t::const_iterator i = sig.begin(); i != sig.end() ; ++i)
    {
      os << " (" << pp(m_lts.action_label(i->first)) << ", " << i->second << ") ";
    }
    os << " }";
    return os.str();
  }

  /** \brief Compute the partition. Repeatedly updates the signatures, and
             the partition, until the partition stabilises */
  void compute_partition()
  {
    std::size_t count_prev = m_count;
    std::size_t iterations = 0;

    sort_transitions(m_lts.get_transitions(), m_lts.hidden_label_set(), mcrl2::lts::lbl_tgt_src);

    do
    {
      mCRL2log(log::verbose) << "Iteration " << iterations
                                       << " currently have " << m_count << " blocks" << std::endl;

      m_signature.compute_signature(m_partition);

      count_prev = m_count;

      // Map signatures to block numbers
      std::map<signature_t, std::size_t> hashtable;
      m_count = 0;
      for(std::size_t i = 0; i < m_lts.num_states(); ++i)
      {
        if(hashtable.find(m_signature.get_signature(i)) == hashtable.end())
        {
          mCRL2log(log::debug) << "Adding block for signature " << print_sig(m_signature.get_signature(i)) << std::endl;
          hashtable[m_signature.get_signature(i)] = m_count++;
        }
      }

      // Map states to block numbers
      for(std::size_t i = 0; i < m_lts.num_states(); ++i)
      {
        m_partition[i] = hashtable[m_signature.get_signature(i)];
      }

      ++iterations;

    } while (count_prev != m_count);

    mCRL2log(log::verbose) << "Done after " << iterations << " iterations with " << m_count << " blocks" << std::endl;
  }

  /** \brief Perform the quotient with respect to the partition that has
             been computed */
  void quotient()
  {
    // Assign the reduced LTS
    m_lts.set_num_states(m_count);
    m_lts.set_initial_state(m_partition[m_lts.initial_state()]);

    // Compute quotient transitions
    // implemented in the signature class because it differs per equivalence.
    std::set<transition> transitions;
    m_signature.quotient_transitions(transitions, m_partition);

    // Set quotient transitions
    m_lts.clear_transitions();
    for(std::set<transition>::const_iterator i = transitions.begin(); i != transitions.end(); ++i)
    {
      m_lts.add_transition(*i);
    }
  }

public:
  /** \brief Constructor
    * \param[in] lts_ The LTS that is being reduced
    */
  sigref(LTS_T& lts_)
    : m_partition(std::vector<std::size_t>(lts_.num_states(), 0)),
      m_count(0),
      m_lts(lts_),
      m_signature(lts_)
  {}

  /** \brief Perform the reduction, modulo the equivalence for which the
    *        signature has been passed in as template parameter
    */
  void run()
  {
    // No need for state labels in the reduced LTS.
    m_lts.clear_state_labels();
    mcrl2::utilities::execution_timer timer;

    // Now apply the branching bisimulation reduction algorithm.  If there
    // are no taus, this will automatically yield strong bisimulation.
    timer.start("reduction");
    compute_partition();
    timer.finish("reduction");
    quotient();
  }
};

} // namespace lts
} // namespace mcrl2

#endif // MCRL2_LTS_SIGREF_H
