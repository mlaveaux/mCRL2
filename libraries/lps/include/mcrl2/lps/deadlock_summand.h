// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/deadlock_summand.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_DEADLOCK_SUMMAND_H
#define MCRL2_LPS_DEADLOCK_SUMMAND_H

#include "mcrl2/lps/deadlock.h"
#include "mcrl2/lps/stochastic_distribution.h"
#include "mcrl2/lps/summand.h"

namespace mcrl2 {

namespace lps {

/// \brief LPS summand containing a deadlock.
class deadlock_summand: public summand_base
{
  protected:
    /// \brief The super class
    typedef summand_base super;

    /// \brief The deadlock of the summand
    lps::deadlock m_deadlock;

  public:
    /// \brief Constructor.
    // TODO: check if the default constructor results in a deadlock summand
    deadlock_summand()
    {}

    /// \brief Constructor.
    deadlock_summand(const data::variable_list& summation_variables, const data::data_expression& condition, const lps::deadlock& delta)
      : summand_base(summation_variables, condition),
        m_deadlock(delta)
    {}

    /// Move semantics
    deadlock_summand(const deadlock_summand&) noexcept = default;
    deadlock_summand(deadlock_summand&&) noexcept = default;
    deadlock_summand& operator=(const deadlock_summand&) noexcept = default;
    deadlock_summand& operator=(deadlock_summand&&) noexcept = default;

    /// \brief Returns the deadlock of this summand.
    const lps::deadlock& deadlock() const
    {
      return m_deadlock;
    }

    /// \brief Returns the deadlock of this summand.
    lps::deadlock& deadlock()
    {
      return m_deadlock;
    }

    /// \brief Returns true if time is available.
    /// \return True if time is available.
    bool has_time() const
    {
      return m_deadlock.has_time();
    }

    /// \brief Swaps the contents
    void swap(deadlock_summand& other)
    {
      summand_base::swap(other);
      using std::swap;
      swap(m_deadlock, other.m_deadlock);
    }
};

//--- start generated class deadlock_summand ---//
/// \\brief list of deadlock_summands
typedef atermpp::term_list<deadlock_summand> deadlock_summand_list;

/// \\brief vector of deadlock_summands
typedef std::vector<deadlock_summand>    deadlock_summand_vector;

// prototype declaration
std::string pp(const deadlock_summand& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const deadlock_summand& x)
{
  return out << lps::pp(x);
}

/// \\brief swap overload
inline void swap(deadlock_summand& t1, deadlock_summand& t2)
{
  t1.swap(t2);
}
//--- end generated class deadlock_summand ---//

/// \brief Conversion to atermappl.
/// \return The deadlock summand converted to aterm format.
inline
atermpp::aterm deadlock_summand_to_aterm(const deadlock_summand& s)
{
  atermpp::aterm result = atermpp::aterm(core::detail::function_symbol_LinearProcessSummand(),
                       s.summation_variables(),
                       s.condition(),
                       atermpp::aterm(core::detail::function_symbol_Delta()),
                       s.deadlock().time(),
                       data::assignment_list(),
                       stochastic_distribution()
                     );
  return result;
}

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_DEADLOCK_SUMMAND_H
