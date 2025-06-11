// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/stochastic_linear_process.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_STOCHASTIC_LINEAR_PROCESS_H
#define MCRL2_LPS_STOCHASTIC_LINEAR_PROCESS_H

#include "mcrl2/lps/linear_process.h"
#include "mcrl2/lps/stochastic_action_summand.h"
#include "mcrl2/lps/stochastic_process_initializer.h"

namespace mcrl2 {

namespace lps {

namespace detail {

template <>
inline
stochastic_action_summand make_action_summand<stochastic_action_summand>(const data::variable_list& summation_variables,
                                                                         const data::data_expression& condition,
                                                                         const multi_action& a,
                                                                         const data::assignment_list& assignments,
                                                                         const stochastic_distribution& distribution
                                                                        )
{
  return stochastic_action_summand(summation_variables, condition, a, assignments, distribution);
}

inline
stochastic_action_summand_vector convert_action_summands(const action_summand_vector& action_summands)
{
  stochastic_action_summand_vector result;
  for (const action_summand& s: action_summands)
  {
    result.push_back(stochastic_action_summand(s));
  }
  return result;
}

} // namespace detail

/// \brief linear process.
class stochastic_linear_process: public linear_process_base<stochastic_action_summand>
{
  typedef linear_process_base<stochastic_action_summand> super;

  public:
    /// \brief Constructor.
    stochastic_linear_process()
    { }

    /// \brief Constructor.
    stochastic_linear_process(const atermpp::aterm& t, bool stochastic_distributions_allowed = true)
      : super(t, stochastic_distributions_allowed)
    { }

    /// \brief Constructor.
    stochastic_linear_process(const data::variable_list& process_parameters,
                              const deadlock_summand_vector& deadlock_summands,
                              const stochastic_action_summand_vector& action_summands
                             )
      : super(process_parameters, deadlock_summands, action_summands)
    { }

    /// \brief Constructor.
    explicit stochastic_linear_process(const linear_process& other)
      : super(other.process_parameters(), other.deadlock_summands(), detail::convert_action_summands(other.action_summands()))
    { }
};

//--- start generated class stochastic_linear_process ---//
// prototype declaration
std::string pp(const stochastic_linear_process& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const stochastic_linear_process& x)
{
  return out << lps::pp(x);
}
//--- end generated class stochastic_linear_process ---//

// template function overloads
std::set<data::variable> find_all_variables(const lps::stochastic_linear_process& x);
std::set<data::variable> find_free_variables(const lps::stochastic_linear_process& x);

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_STOCHASTIC_LINEAR_PROCESS_H
