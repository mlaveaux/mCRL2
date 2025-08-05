// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

/** \file probabilistic_lts.h
 *
 * \brief The file containing the core class for transition systems.
 * \details This is the LTS library's main header file. It declares all publicly
 * accessible data structures that form the interface of the library.
 * \author Jan Friso Groote
 */

#ifndef MCRL2_LTS_PROBABILISTIC_STATE_H
#define MCRL2_LTS_PROBABILISTIC_STATE_H

#include "mcrl2/lps/state_probability_pair.h"

#include <sstream>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::lts
{

/** \brief A class that contains a probabilistic state.
    \details A probabilistic state is essentially a sequence
             of pairs of a state and a probability. The probability
             indicates the likelyhood with which that particular
             state can be reached. The sum of all probabilities
             in a probabilistic state is one. 
             Internally, it is either stored as single state, 
             or as a vector of states. 
             When using a probabilistic state, a user should
             be aware which of the two internal representations are used.
             After calling "construct_internal_vector_representation"
             it is guaranteed that the state is stored as a vector,
             and iterators over begin/end, and cbegin/cend can be used.
*/

template < class STATE, class PROBABILITY >
class probabilistic_state
{
  public:

    friend std::hash<probabilistic_state>;

    using state_probability_pair = typename lps::state_probability_pair<STATE, PROBABILITY>;
    using state_t = STATE;
    using probability_t = PROBABILITY;
    using iterator = typename std::vector<state_probability_pair>::iterator;
    using const_iterator = typename std::vector<state_probability_pair>::const_iterator;
    using reverse_iterator = typename std::vector<state_probability_pair>::reverse_iterator;
    using const_reverse_iterator = typename std::vector<state_probability_pair>::const_reverse_iterator;

  protected:

    STATE m_single_state;                                      // If there is only one state, this is stored as
                                                               // a base state with implicit probability 1. 
    std::vector<state_probability_pair> m_probabilistic_state; // A vector with states and their associated probability.
                                                               // The sum of all probabilities in the vector must be one.
                                             
  public:

    /** \brief Default constructor
     **/
    probabilistic_state()
      : m_single_state(STATE(-1))
    {
    }

    /** \brief Constructor of a probabilistic state from a non probabilistic state.
     *  \param[in] s A state.
     *  \return The generated probabilistic state.
     **/
    explicit probabilistic_state(const STATE& s)
      : m_single_state(s)
    {
    }

    /** \brief Copy constructor **/
    probabilistic_state(const probabilistic_state& other)
      : m_single_state(other.m_single_state),
        m_probabilistic_state(other.m_probabilistic_state)
    {
      shrink_to_fit();
    }

    /** \brief Copy assignment constructor **/
    probabilistic_state& operator=(const probabilistic_state& other)
    {
      m_single_state= other.m_single_state;
      m_probabilistic_state = other.m_probabilistic_state;
      shrink_to_fit();
      return *this;
    }

    /** \brief Creates a probabilistic state on the basis of state_probability_pairs.
     * \param[in] begin Iterator to the first state_probability_pair.
     * \param[in] end Iterator to the last state_probability_pair.
     * \return Resulting probabilistic state.  **/
    template <class STATE_PROBABILITY_PAIR_ITERATOR>
    probabilistic_state(const STATE_PROBABILITY_PAIR_ITERATOR begin, const STATE_PROBABILITY_PAIR_ITERATOR end)
      : m_single_state(STATE(-1)),
        m_probabilistic_state(begin,end)
    {
      assert(begin!=end);
      shrink_to_fit();
    }

    /** \brief Standard equality operator.
     *  \result Returns true iff the probabilistic states are equal.
     */
    bool operator==(const probabilistic_state& other) const
    {
      if (m_single_state==STATE(-1))
      {
        if (other.m_single_state==STATE(-1))
        {
          assert(m_probabilistic_state.size()>0);
          assert(other.m_probabilistic_state.size()>0);
          return m_probabilistic_state==other.m_probabilistic_state;
        }
        else 
        {
          return m_probabilistic_state.size()==1 &&
                 other.m_single_state==m_probabilistic_state.front().state();
        }
      }
      else
      {
        if (other.m_single_state==STATE(-1))
        {
          return other.m_probabilistic_state.size()==1 &&
                 m_single_state==other.m_probabilistic_state.front().state();
        }
        else
        {
          assert(m_probabilistic_state.size()==0);
          assert(other.m_probabilistic_state.size()==0);
          return m_single_state==other.m_single_state;
        }
      }
    }

    /** \brief Standard equality operator.
     *  \result Returns true iff the probabilistic states are equal.
     */
    bool operator!=(const probabilistic_state& other) const
    {
      return !(operator==(other)); 
    }

    /** \brief Swap this probabilistic state.
     * \param[in] s A probabilistic state. */
    void swap(probabilistic_state& other) noexcept
    {
      std::swap(m_single_state,other.m_single_state);
      m_probabilistic_state.swap(other.m_probabilistic_state);
    }

    /** \brief Guarantee that this probabilistic state is internally stored as
     **        a vector, such that begin/end, cbegin/cend and size are properly 
     **        defined. */
    void construct_internal_vector_representation()
    {
      if (m_single_state!=STATE(-1))
      {
        assert(m_probabilistic_state.size()==0);
        m_probabilistic_state.emplace_back(m_single_state,PROBABILITY::one());
        m_single_state=STATE(-1);
      }
    }

    /** \brief Set this probabilistic state to a single state with probability one.
     * \details It is assumed that the given state probability pair does not have
               any element.
     * \param[in] s The state. */
    void set(const STATE& s)
    {
      assert(m_probabilistic_state.size()==0);
      m_single_state=s; 
    }

    /** \brief Get a probabilistic state if is is simple, i.e., consists of a single state. 
     * \details It is assumed that the given state probability pair does not have
               any element.
     * \param[in] s The state. */
    STATE get() const
    {
      if (m_probabilistic_state.size()>1)
      {
        throw mcrl2::runtime_error("Probabilistic state is not simple\n");
      }
      if (m_probabilistic_state.size()==1)
      {
        assert(m_probabilistic_state.front().state()!=STATE(-1));
        return m_probabilistic_state.front().state();
      }
      assert(m_single_state!=STATE(-1));
      return m_single_state; 
    }

    /** \brief Add a state with a probability to the probabilistic state
     * \param[in] s The state to be added.
     * \param[in] p The probability of this state.
     **/
    void add(const STATE& s, const PROBABILITY& p)
    {
      m_probabilistic_state.emplace_back(s,p);
    }

    /** \brief If a probabilistic state is ready, shrinking it to minimal size might be useful
     *         to reduce its memory usage. A requirement is that the sum of the probabilities must
     *         be one.
     */
    void shrink_to_fit()
    {
      m_probabilistic_state.shrink_to_fit();
      assert((m_probabilistic_state.size()>1 && m_single_state==std::size_t(-1)) ||
             (m_probabilistic_state.size()==0 && m_single_state!=std::size_t(-1)));
#ifndef NDEBUG
      PROBABILITY sum=PROBABILITY::zero();
      for(const state_probability_pair& p: m_probabilistic_state)
      {
        assert(p.probability()>PROBABILITY::zero());
        assert(p.probability()<=PROBABILITY::one());
        sum=sum+p.probability();
      }
      // SVC studio copies empty probabilistic states.
      assert(m_probabilistic_state.size()==0 || sum==PROBABILITY::one());
#endif
    }

    /** \brief Gets the number of probabilistic states in the vector representation of this state.
     **        If the state is stored as a simple number this size returns 0. So, 0 means that there
     **        is one state, which should then be obtained using get(). If the size is larger than 0
     **        the state is stored in a vector, and it must be accessed throught the iterators begin()
     **        and n. 
     * \return The number of probabilistic states, where 0 means there is one simple state. */
    std::size_t size() const
    {
      return m_probabilistic_state.size();
    }

    /** \brief Makes the probabilistic state empty.
     */
    void clear()
    {
      m_single_state=STATE(-1);
      m_probabilistic_state.clear();
    }

    /** \brief Gets an iterator over pairs of state and probability. This can only be used when the state is stored 
     **        internally as a vector.  
     *  \return The iterator pointing at the first state probability pair. */
    const_iterator begin() const
    {
      return m_probabilistic_state.begin();
    }

    /** \brief Gets the end iterator over pairs of state and probability.
     ** \return The iterator pointing beyond the last state probability pair in a probabilistic state. */
    const_iterator end() const
    {
      return m_probabilistic_state.end();
    }

    /** \brief Gets an iterator over pairs of state and probability. This can only be used if the state is internally stored as a vector. 
     ** \return The iterator pointing at the first state probability pair. */
    iterator begin()
    {
      return m_probabilistic_state.begin();
    }

    /** \brief Gets the end iterator over pairs of state and probability.
     ** \return The iterator pointing beyond the last state probability pair in a probabilistic state. */
    iterator end()
    {
      return m_probabilistic_state.end();
    }

    /** \brief Gets a reverse iterator over pairs of state and probability. This can only be used when the state is stored 
     **        internally as a vector.  
     *  \return The iterator pointing at the last state probability pair. */
    const_reverse_iterator rbegin() const
    {
      return m_probabilistic_state.rbegin();
    }

    /** \brief Gets the reverse end iterator over pairs of state and probability.
     ** \return The iterator pointing before the first state probability pair in a probabilistic state. */
    const_reverse_iterator rend() const
    {
      return m_probabilistic_state.rend();
    }

    /** \brief Gets a reverse iterator over pairs of state and probability. This can only be used if the state is internally stored as a vector. 
     ** \return The iterator pointing at the last state probability pair. */
    reverse_iterator rbegin()
    {
      return m_probabilistic_state.rbegin();
    }

    /** \brief Gets the reverse end iterator over pairs of state and probability.
     ** \return The iterator pointing before the first state probability pair in a probabilistic state. */
    reverse_iterator rend()
    {
      return m_probabilistic_state.rend();
    }
};


/* \brief A pretty print operator on action labels, returning it as a string.
 * */
template < class STATE, class PROBABILITY >
inline std::string pp(const probabilistic_state<STATE, PROBABILITY>& l)
{
  std::stringstream str;
  if (l.size()<=1)
  {
    str << "[ " << l.get() << ": 1.0 ]";
  }
  else
  {
    str << "[ ";
    bool first=true;
    for(const lps::state_probability_pair<STATE, PROBABILITY>& p: l)
    {
      if (first)
      {
        first=false;
      }
      else
      {
        str << "; ";
      }
      str << p.state() << ": " << p.probability();
    }
    str << " ]";
  }

  return str.str();
}

/** \brief Pretty print to an outstream.
 */
template < class STATE, class PROBABILITY >
inline
std::ostream& operator<<(std::ostream& out, const probabilistic_state<STATE, PROBABILITY>& l)
{
  return out << pp(l);
}



} // namespace mcrl2::lts


namespace std
{

/// \brief Specialization of the standard std::hash function.
/// \details It is essential that this hash function yields the same has for a singular
///          state in a distribution, stored as a number with implicit probability 1, or
///          as a vector of length 1. 
template < class STATE, class PROBABILITY >
struct hash< mcrl2::lts::probabilistic_state<STATE, PROBABILITY> >
{
  std::size_t operator()(const mcrl2::lts::probabilistic_state<STATE, PROBABILITY>& p) const
  {
    if (p.m_single_state!=STATE(-1))
    {
      assert(p.m_probabilistic_state.size()==0);
      hash<STATE> state_hasher;
      hash<PROBABILITY> probability_hasher;
      return mcrl2::utilities::detail::hash_combine(0, 
                   mcrl2::utilities::detail::hash_combine(state_hasher(p.m_single_state), 
                                                          probability_hasher(PROBABILITY::one())));
    }
    hash<vector<typename mcrl2::lps::state_probability_pair< STATE, PROBABILITY > > > hasher;
    return hasher(p.m_probabilistic_state);
  }
};

} // namespace std

#endif // MCRL2_LTS_PROBABILISTIC_STATE_H
