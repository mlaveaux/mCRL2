// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/action_label.h
/// \brief The class action_label.

#ifndef MCRL2_PROCESS_ACTION_LABEL_H
#define MCRL2_PROCESS_ACTION_LABEL_H

#include "mcrl2/data/data_specification.h"

namespace mcrl2::process
{

//--- start generated class action_label ---//
/// \\brief An action label
class action_label: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    action_label()
      : atermpp::aterm(core::detail::default_values::ActId)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit action_label(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_ActId(*this));
    }

    /// \\brief Constructor Z12.
    action_label(const core::identifier_string& name, const data::sort_expression_list& sorts)
      : atermpp::aterm(core::detail::function_symbol_ActId(), name, sorts)
    {}

    /// \\brief Constructor Z1.
    action_label(const std::string& name, const data::sort_expression_list& sorts)
      : atermpp::aterm(core::detail::function_symbol_ActId(), core::identifier_string(name), sorts)
    {}

    /// Move semantics
    action_label(const action_label&) noexcept = default;
    action_label(action_label&&) noexcept = default;
    action_label& operator=(const action_label&) noexcept = default;
    action_label& operator=(action_label&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const data::sort_expression_list& sorts() const
    {
      return atermpp::down_cast<data::sort_expression_list>((*this)[1]);
    }
};

/// \\brief The function make_action_label constructs a new term into a given address.
/// \\ \param t The reference into which the new action_label is constructed.
template <class... ARGUMENTS>
inline void make_action_label(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActId(), args...);
}

/// \\brief list of action_labels
using action_label_list = atermpp::term_list<action_label>;

/// \\brief vector of action_labels
using action_label_vector = std::vector<action_label>;

/// \\brief Test for a action_label expression
/// \\param x A term
/// \\return True if \\a x is a action_label expression
inline
bool is_action_label(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActId;
}

// prototype declaration
std::string pp(const action_label& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const action_label& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(action_label& t1, action_label& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated class action_label ---//

// template function overloads
std::string pp(const action_label_list& x, bool precedence_aware = true);
std::string pp(const action_label_vector& x, bool precedence_aware = true);
action_label_list normalize_sorts(const action_label_list& x, const data::sort_specification& sortspec);
std::set<data::sort_expression> find_sort_expressions(const process::action_label_list& x);

/// Compare two actions names.
/// This is implemented using string comparison, such that the comparison is
/// stable across different runs.
struct action_name_compare
{
  bool operator()(const core::identifier_string& s1, const core::identifier_string& s2) const
  {
    return std::string(s1) < std::string(s2);
  }
};

/// Total ordering on action labels: name first (string order), then sorts.
inline bool operator<(const action_label& a1, const action_label& a2)
{
  const std::string& s1 = a1.name();
  const std::string& s2 = a2.name();
  return s1 < s2 || (s1 == s2 && a1.sorts() < a2.sorts());
}

} // namespace mcrl2::process



#endif // MCRL2_PROCESS_ACTION_LABEL_H
