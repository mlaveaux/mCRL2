// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \details This file contains two functions to efficiently transform 
//           a plain_lts into a probabilistic lts and vice versa. The
//           lts from which the transformation takes place is destroyed.

#ifndef MCRL2_LTS_DETAIL_SWAP_TO_FROM_PROBABILISTIC_LTS_H
#define MCRL2_LTS_DETAIL_SWAP_TO_FROM_PROBABILISTIC_LTS_H

#include "mcrl2/lts/probabilistic_lts.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::lts::detail
{

template < class STATE_LABEL_T, class ACTION_LABEL_T, class PROBABILISTIC_STATE_T, class LTS_BASE >
void swap_to_non_probabilistic_lts(
        probabilistic_lts<STATE_LABEL_T,ACTION_LABEL_T, PROBABILISTIC_STATE_T, LTS_BASE>& l_probabilistic, 
        lts<STATE_LABEL_T,ACTION_LABEL_T, LTS_BASE>& l_plain)
{
  static_cast<lts<STATE_LABEL_T, ACTION_LABEL_T, LTS_BASE>& >(l_probabilistic).swap(l_plain);
  
  if (l_probabilistic.initial_probabilistic_state().size()<=1)
  { 
    l_plain.set_initial_state(l_probabilistic.initial_probabilistic_state().get());
  }
  else
  {
    throw mcrl2::runtime_error("Initial state is probabilistic and cannot be transformed into a non probabilistic state.");
  }
  
  // Adapt the probabilistic target states to non probabilistic target states.
  std::size_t transition_number=1;
  for(transition& t: l_plain.get_transitions())
  {
    std::size_t probabilistic_target_state_number=t.to();
    if (l_probabilistic.probabilistic_state(probabilistic_target_state_number).size()>1)
    {
      throw mcrl2::runtime_error("Transition " + std::to_string(transition_number) + " is probabilistic.");
    }
    else
    {
      t=transition(t.from(), t.label(), l_probabilistic.probabilistic_state(probabilistic_target_state_number).get());
    }
    transition_number++;
  }

  
}

template < class STATE_LABEL_T, class ACTION_LABEL_T, class PROBABILISTIC_STATE_T, class LTS_BASE >
void translate_to_probabilistic_lts(
        const lts<STATE_LABEL_T,ACTION_LABEL_T, LTS_BASE>& l_plain,
        probabilistic_lts<STATE_LABEL_T,ACTION_LABEL_T, PROBABILISTIC_STATE_T, LTS_BASE>& l_probabilistic) 
{
  static_cast<lts<STATE_LABEL_T,ACTION_LABEL_T, LTS_BASE>& >(l_probabilistic)=l_plain;
  l_probabilistic.set_initial_probabilistic_state(PROBABILISTIC_STATE_T(l_plain.initial_state()));
  for(std::size_t i=0; i<l_plain.num_states(); ++i)
  {
    l_probabilistic.add_probabilistic_state(PROBABILISTIC_STATE_T(i));
  }
}

} // detail
// lts
// mcrl2

#endif  // MCRL2_LTS_DETAIL_SWAP_TO_FROM_PROBABILISTIC_LTS_H
