// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/process_identifier.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_PROCESS_IDENTIFIER_H
#define MCRL2_PROCESS_PROCESS_IDENTIFIER_H

#include "mcrl2/process/action_name_multiset.h"

namespace mcrl2
{

namespace process
{

//--- start generated class process_identifier ---//
/// \\brief A process identifier
class process_identifier: public atermpp::aterm
{
  public:


    /// Move semantics
    process_identifier(const process_identifier&) noexcept = default;
    process_identifier(process_identifier&&) noexcept = default;
    process_identifier& operator=(const process_identifier&) noexcept = default;
    process_identifier& operator=(process_identifier&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[1]);
    }
//--- start user section process_identifier ---//
    /// \brief Default constructor.
    process_identifier()
      : atermpp::aterm(core::detail::default_values::ProcVarId)
    {}

    /// \brief Constructor.
    /// \param term A term
    explicit process_identifier(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_ProcVarId(*this));
    }

    /// \brief Constructor.
    process_identifier(const core::identifier_string& name, const data::variable_list& variables)
    {
      atermpp::make_term_appl(*this, core::detail::function_symbol_ProcVarId(), name, variables);
    }

    /// \brief Constructor.
    process_identifier(const std::string& name, const data::variable_list& variables)
    {
      atermpp::make_term_appl(*this, core::detail::function_symbol_ProcVarId(), core::identifier_string(name), variables);
    }
//--- end user section process_identifier ---//
};

/// \\brief Make_process_identifier constructs a new term into a given address.
/// \\ \param t The reference into which the new process_identifier is constructed. 
template <class... ARGUMENTS>
inline void make_process_identifier(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ProcVarId(), args...);
}

/// \\brief list of process_identifiers
typedef atermpp::term_list<process_identifier> process_identifier_list;

/// \\brief vector of process_identifiers
typedef std::vector<process_identifier>    process_identifier_vector;

/// \\brief Test for a process_identifier expression
/// \\param x A term
/// \\return True if \\a x is a process_identifier expression
inline
bool is_process_identifier(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ProcVarId;
}

// prototype declaration
std::string pp(const process_identifier& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const process_identifier& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(process_identifier& t1, process_identifier& t2)
{
  t1.swap(t2);
}
//--- end generated class process_identifier ---//

// template function overloads
std::string pp(const process_identifier_list& x);
std::string pp(const process_identifier_vector& x);
void normalize_sorts(process_identifier_vector& x, const data::sort_specification& sortspec);

} // namespace process

} // namespace mcrl2

#endif // MCRL2_PROCESS_PROCESS_IDENTIFIER_H
