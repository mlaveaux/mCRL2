// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/deadlock.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_DEADLOCK_H
#define MCRL2_LPS_DEADLOCK_H

#include "mcrl2/data/print.h"

namespace mcrl2
{

namespace lps
{

/// \brief Represents a deadlock
/// \details A deadlock is 'delta' with an optional time tag.
class deadlock
{
    friend class deadlock_summand;

  protected:
    /// \brief The time of the deadlock. If <tt>m_time == data::undefined_real()</tt>
    /// the multi action has no time.
    data::data_expression m_time;

  public:
    /// \brief Constructor
    deadlock(data::data_expression time = data::undefined_real())
      : m_time(time)
    {}

    /// \brief Returns true if time is available.
    /// \return True if time is available.
    bool has_time() const
    {
      return m_time != data::undefined_real();
    }

    /// \brief Returns the time.
    /// \return The time.
    const data::data_expression& time() const
    {
      return m_time;
    }

    /// \brief Returns the time.
    /// \return The time.
    data::data_expression& time()
    {
      return m_time;
    }

    /// \brief Returns a string representation of the deadlock
    std::string to_string() const
    {
      return std::string("delta") + (has_time() ? (" @ " + data::pp(m_time)) : "");
    }

    /// \brief Comparison operator
    bool operator==(const deadlock& other) const
    {
      return m_time == other.m_time;
    }

    /// \brief Comparison operator
    bool operator!=(const deadlock& other) const
    {
      return !(*this == other);
    }

    /// \brief Swaps the contents
    void swap(deadlock& other)
    {
      using std::swap;
      swap(m_time, other.m_time);
    }
};

//--- start generated class deadlock ---//
/// \\brief list of deadlocks
typedef atermpp::term_list<deadlock> deadlock_list;

/// \\brief vector of deadlocks
typedef std::vector<deadlock>    deadlock_vector;

// prototype declaration
std::string pp(const deadlock& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const deadlock& x)
{
  return out << lps::pp(x);
}

/// \\brief swap overload
inline void swap(deadlock& t1, deadlock& t2)
{
  t1.swap(t2);
}
//--- end generated class deadlock ---//

// template function overloads
std::set<data::variable> find_all_variables(const lps::deadlock& x);
std::set<data::variable> find_free_variables(const lps::deadlock& x);

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_DEADLOCK_H
