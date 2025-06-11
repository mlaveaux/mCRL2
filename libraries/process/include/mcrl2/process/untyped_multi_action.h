// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/untyped_multi_action.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_UNTYPED_MULTI_ACTION_H
#define MCRL2_PROCESS_UNTYPED_MULTI_ACTION_H

#include "mcrl2/data/untyped_data_parameter.h"

namespace mcrl2 {

namespace process {

//--- start generated class untyped_multi_action ---//
/// \\brief An untyped multi action or data application
class untyped_multi_action: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    untyped_multi_action()
      : atermpp::aterm(core::detail::default_values::UntypedMultiAction)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit untyped_multi_action(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_UntypedMultiAction(*this));
    }

    /// \\brief Constructor Z12.
    explicit untyped_multi_action(const data::untyped_data_parameter_list& actions)
      : atermpp::aterm(core::detail::function_symbol_UntypedMultiAction(), actions)
    {}

    /// Move semantics
    untyped_multi_action(const untyped_multi_action&) noexcept = default;
    untyped_multi_action(untyped_multi_action&&) noexcept = default;
    untyped_multi_action& operator=(const untyped_multi_action&) noexcept = default;
    untyped_multi_action& operator=(untyped_multi_action&&) noexcept = default;

    const data::untyped_data_parameter_list& actions() const
    {
      return atermpp::down_cast<data::untyped_data_parameter_list>((*this)[0]);
    }
};

/// \\brief Make_untyped_multi_action constructs a new term into a given address.
/// \\ \param t The reference into which the new untyped_multi_action is constructed. 
template <class... ARGUMENTS>
inline void make_untyped_multi_action(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_UntypedMultiAction(), args...);
}

/// \\brief list of untyped_multi_actions
typedef atermpp::term_list<untyped_multi_action> untyped_multi_action_list;

/// \\brief vector of untyped_multi_actions
typedef std::vector<untyped_multi_action>    untyped_multi_action_vector;

/// \\brief Test for a untyped_multi_action expression
/// \\param x A term
/// \\return True if \\a x is a untyped_multi_action expression
inline
bool is_untyped_multi_action(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::UntypedMultiAction;
}

// prototype declaration
std::string pp(const untyped_multi_action& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const untyped_multi_action& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(untyped_multi_action& t1, untyped_multi_action& t2)
{
  t1.swap(t2);
}
//--- end generated class untyped_multi_action ---//

} // namespace process

} // namespace mcrl2

#endif // MCRL2_PROCESS_UNTYPED_MULTI_ACTION_H
