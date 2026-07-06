// Author(s): Frank Stappers, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts2lps.h
/// \brief This file contains functions to translate lts's into
///        LPS format.


#ifndef MCRL2_LTS_LTS2LPS_H
#define MCRL2_LTS_LTS2LPS_H

#include "mcrl2/lps/io.h"
#include "mcrl2/lts/detail/lts_convert.h"

namespace mcrl2::lts
{

using namespace mcrl2::lps;
using namespace mcrl2::data;

/// \brief transform an lts in lts format into a linear process.
/// \param l A labelled transition system in lts format.
/// \return The function returns a linear process with the same behaviour as the lts.
inline lps::specification transform_lts2lps(const lts_lts_t& l)
{
  action_summand_vector action_summands;
  const variable process_parameter("x",mcrl2::data::sort_pos::pos());
  const variable_list process_parameters = { process_parameter };
  const std::set< data::variable> global_variables;
  // Add a single delta.
  const deadlock_summand_vector deadlock_summands(1,deadlock_summand(variable_list(), sort_bool::true_(), deadlock()));
  const process_initializer initial_process(data::data_expression_list{sort_pos::pos(l.initial_state()+1)}
  );

  const std::vector<transition>& trans=l.get_transitions();
  for (const auto& tran: trans)
  {
    const lps::multi_action actions = l.action_label(tran.label());

    assignment_list assignments;
    if (tran.from() != tran.to())
    {
      assignments = push_back(assignments, assignment(process_parameter, sort_pos::pos(tran.to() + 1)));
    }

    const action_summand summand(variable_list(),
      equal_to(process_parameter, sort_pos::pos(tran.from() + 1)),
      actions,
      assignments);
    action_summands.push_back(summand);
  }
  
  const linear_process lps1(process_parameters,deadlock_summands,action_summands);
  const specification spec(l.data(),l.action_label_declarations(),global_variables,lps1,initial_process);
  
  return spec;
}

/// \brief transform an lts in aut format into a linear process.
/// \param l1 A labelled transition system in aut format.
/// \param data A separate data specification.
/// \param action_labels A list containing the action labels used in the lts.
/// \param process_parameters The process parameters of the current process. 
/// \return The function returns a linear process with the same behaviour as the lts.
inline lps::specification transform_lts2lps(const lts_aut_t& l1,
    const data_specification& data,
    const process::action_label_list& action_labels,
    const variable_list& process_parameters)
{
  lts_lts_t l2;
  mcrl2::lts::detail::lts_convert(l1,l2,data,action_labels,process_parameters);  
  return transform_lts2lps(l2);
}

/// \brief transform an lts in fsm format into a linear process.
/// \param l1 A labelled transition system in fsm format.
/// \param data A separate data specification.
/// \param action_labels A list containing the action labels used in the lts.
/// \param process_parameters The process parameters of the current process. XXXXX Is this needed????
/// \return The function returns a linear process with the same behaviour as the lts.

inline lps::specification transform_lts2lps(const lts_fsm_t& l1,
    const data_specification& data,
    const process::action_label_list& action_labels,
    const variable_list& process_parameters)
{
  lts_lts_t l2;
  mcrl2::lts::detail::lts_convert(l1,l2,data,action_labels,process_parameters);  
  return transform_lts2lps(l2);
}

} // namespace mcrl2::lts

#endif // MCRL2_LTS_LTS2LPS_H
