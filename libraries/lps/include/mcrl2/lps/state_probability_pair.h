// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

/** \file
 *
 * \brief This file contains a class with a state/probability pair. 
 * \details A state probability pair consists of a state and a probability.
 *          The probability should be larger than 0 and smaller or equal to 1.
 * \author Jan Friso Groote
 */


#ifndef MCRL2_LPS_STATE_PROBABILITY_PAIR_H
#define MCRL2_LPS_STATE_PROBABILITY_PAIR_H

#include <cassert>
#include <iostream>
#include "mcrl2/utilities/hash_utility.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm.cxx"
#else
  import atermpp;
#endif

namespace mcrl2::lps
{

template < class STATE, class PROBABILITY >
class state_probability_pair
{
  protected:
    // A state 
    STATE m_state;
    // A data expression of sort Real that indicates the probability of this state
    PROBABILITY m_probability;
  public:

    /// \brief constructor.
    /// \param state The state.
    /// \param probability The probability of the occurrence of this state. 
    state_probability_pair(const STATE& state, const PROBABILITY& probability)
     : m_state(state), 
       m_probability(probability)
    {
      assert(probability!=PROBABILITY::zero());
    } 

    /// \brief Copy constructor;
    /// \param p The state probability pair to be copied.
    state_probability_pair(const state_probability_pair& p)=default; 
    state_probability_pair(state_probability_pair&& p)=default; 

    /// \brief Standard assignment.
    /// \param p The state probability pair to be assigned.
    /// \return A reference to the assigned object.
    state_probability_pair& operator=(const state_probability_pair& p)=default; 
    state_probability_pair& operator=(state_probability_pair&& p)=default; 

    /** \brief Standard equality operator.
     *  \result Returns true iff the probabilistic states are equal.
     */
    bool operator==(const state_probability_pair& other) const
    {
      if constexpr (std::is_convertible_v<PROBABILITY, atermpp::aterm>)
      {
        // The probabilities are compared as aterms, and not based on their value, as comparing
        // probabilities using their value is expensive as it requires an application of the rewriter. 
        return m_state==other.m_state && 
               static_cast<atermpp::aterm>(m_probability)==static_cast<atermpp::aterm>(other.m_probability);
      }
      else
      {
        static_assert(!std::is_convertible_v<PROBABILITY, atermpp::aterm>);
        return m_state==other.m_state && m_probability==other.m_probability;
      }
    }

    /// \brief Get the state from a state probability pair.
    const STATE& state() const
    {
      return m_state;
    }

    /// \brief Get the state in a state probability pair.
    STATE& state() 
    {
      return m_state;
    }

    /// \brief get the probability from a state proability pair.
    const PROBABILITY& probability() const
    {
      return m_probability;
    }

    /// \brief Set the probability in a state probability pair.
    PROBABILITY& probability() 
    {
      return m_probability;
    }

};

} // namespace mcrl2::lps


namespace std
{

template < class STATE, class PROBABILITY >
struct hash<mcrl2::lps::state_probability_pair<STATE,PROBABILITY> >
{
  std::size_t operator()(const mcrl2::lps::state_probability_pair<STATE,PROBABILITY>& p) const
  {
    hash<STATE> state_hasher;
    hash<PROBABILITY> probability_hasher;
    return mcrl2::utilities::detail::hash_combine(state_hasher(p.state()), probability_hasher(p.probability()));
  }
};

} // end namespace std

#endif // MCRL2_LPS_STATE_PROBABILITY_PAIR_H
