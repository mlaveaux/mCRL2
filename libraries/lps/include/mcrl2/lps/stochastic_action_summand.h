// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/stochastic_action_summand.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_STOCHASTIC_ACTION_SUMMAND_H
#define MCRL2_LPS_STOCHASTIC_ACTION_SUMMAND_H

#include "mcrl2/lps/action_summand.h"

namespace mcrl2 {

namespace lps {

/// \brief LPS summand containing a multi-action.
class stochastic_action_summand: public action_summand
{
  protected:
    /// \brief The distribution of the summand
    stochastic_distribution m_distribution;

  public:
    /// \brief Constructor.
    stochastic_action_summand()
    {}

    /// \brief Constructor.
    stochastic_action_summand(const data::variable_list& summation_variables, const data::data_expression& condition, const lps::multi_action& action,
          const data::assignment_list& assignments, const stochastic_distribution& distribution)
      : action_summand(summation_variables, condition, action, assignments), m_distribution(distribution)
    {}

    /// \brief Constructor.
    stochastic_action_summand(const action_summand& s)
      : action_summand(s)
    {}

    /// Move semantics
    stochastic_action_summand(const stochastic_action_summand&) noexcept = default;
    stochastic_action_summand(stochastic_action_summand&&) noexcept = default;
    stochastic_action_summand& operator=(const stochastic_action_summand&) noexcept = default;
    stochastic_action_summand& operator=(stochastic_action_summand&&) noexcept = default;

    /// \brief Returns the distribution of this summand.
    const stochastic_distribution& distribution() const
    {
      return m_distribution;
    }

    /// \brief Returns the distribution of this summand.
    stochastic_distribution& distribution()
    {
      return m_distribution;
    }

    /// \brief Swaps the contents
    void swap(stochastic_action_summand& other)
    {
      action_summand::swap(other);
      using std::swap;
      swap(m_distribution, other.m_distribution);
    }
};

//--- start generated class stochastic_action_summand ---//
/// \\brief list of stochastic_action_summands
typedef atermpp::term_list<stochastic_action_summand> stochastic_action_summand_list;

/// \\brief vector of stochastic_action_summands
typedef std::vector<stochastic_action_summand>    stochastic_action_summand_vector;

// prototype declaration
std::string pp(const stochastic_action_summand& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const stochastic_action_summand& x)
{
  return out << lps::pp(x);
}

/// \\brief swap overload
inline void swap(stochastic_action_summand& t1, stochastic_action_summand& t2)
{
  t1.swap(t2);
}
//--- end generated class stochastic_action_summand ---//

/// \brief Comparison operator for action summands.
inline
bool operator<(const stochastic_action_summand& x, const stochastic_action_summand& y)
{
  if (x.distribution() == y.distribution())
  {
    return static_cast<const action_summand&>(x) < static_cast<const action_summand&>(y);
  }
  return x.distribution() < y.distribution();
}

/// \brief Equality operator of stochastic action summands
inline
bool operator==(const stochastic_action_summand& x, const stochastic_action_summand& y)
{
  return static_cast<const action_summand&>(x) == static_cast<const action_summand&>(y) && x.distribution() == y.distribution();
}

/// \brief Conversion to aterm.
inline
atermpp::aterm action_summand_to_aterm(const stochastic_action_summand& s)
{
  atermpp::aterm result = atermpp::aterm(core::detail::function_symbol_LinearProcessSummand(),
                       s.summation_variables(),
                       s.condition(),
                       s.multi_action(),
                       s.multi_action().time(),
                       s.assignments(),
                       s.distribution()
                     );
  return result;
}

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_STOCHASTIC_ACTION_SUMMAND_H
