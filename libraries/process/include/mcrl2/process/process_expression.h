// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/process_expression.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_PROCESS_EXPRESSION_H
#define MCRL2_PROCESS_PROCESS_EXPRESSION_H

#include "mcrl2/data/assignment.h"
#include "mcrl2/data/untyped_data_parameter.h"
#include "mcrl2/process/action_label.h"
#include "mcrl2/process/communication_expression.h"
#include "mcrl2/process/process_identifier.h"
#include "mcrl2/process/rename_expression.h"

namespace mcrl2::process
{

//--- start generated classes ---//
/// \\brief A process expression
class process_expression: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    process_expression()
      : atermpp::aterm(core::detail::default_values::ProcExpr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit process_expression(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_ProcExpr(*this));
    }

    /// \\brief Constructor Z6.
    process_expression(const data::untyped_data_parameter& x)
      : atermpp::aterm(x)
    {}

    /// Move semantics
    process_expression(const process_expression&) noexcept = default;
    process_expression(process_expression&&) noexcept = default;
    process_expression& operator=(const process_expression&) noexcept = default;
    process_expression& operator=(process_expression&&) noexcept = default;
};

/// \\brief list of process_expressions
using process_expression_list = atermpp::term_list<process_expression>;

/// \\brief vector of process_expressions
using process_expression_vector = std::vector<process_expression>;

// prototypes
inline bool is_action(const atermpp::aterm& x);
inline bool is_process_instance(const atermpp::aterm& x);
inline bool is_process_instance_assignment(const atermpp::aterm& x);
inline bool is_delta(const atermpp::aterm& x);
inline bool is_tau(const atermpp::aterm& x);
inline bool is_sum(const atermpp::aterm& x);
inline bool is_block(const atermpp::aterm& x);
inline bool is_hide(const atermpp::aterm& x);
inline bool is_rename(const atermpp::aterm& x);
inline bool is_comm(const atermpp::aterm& x);
inline bool is_allow(const atermpp::aterm& x);
inline bool is_sync(const atermpp::aterm& x);
inline bool is_at(const atermpp::aterm& x);
inline bool is_seq(const atermpp::aterm& x);
inline bool is_if_then(const atermpp::aterm& x);
inline bool is_if_then_else(const atermpp::aterm& x);
inline bool is_bounded_init(const atermpp::aterm& x);
inline bool is_merge(const atermpp::aterm& x);
inline bool is_left_merge(const atermpp::aterm& x);
inline bool is_choice(const atermpp::aterm& x);
inline bool is_stochastic_operator(const atermpp::aterm& x);
inline bool is_untyped_process_assignment(const atermpp::aterm& x);

/// \\brief Test for a process_expression expression
/// \\param x A term
/// \\return True if \\a x is a process_expression expression
inline
bool is_process_expression(const atermpp::aterm& x)
{
  return data::is_untyped_data_parameter(x) ||
         process::is_action(x) ||
         process::is_process_instance(x) ||
         process::is_process_instance_assignment(x) ||
         process::is_delta(x) ||
         process::is_tau(x) ||
         process::is_sum(x) ||
         process::is_block(x) ||
         process::is_hide(x) ||
         process::is_rename(x) ||
         process::is_comm(x) ||
         process::is_allow(x) ||
         process::is_sync(x) ||
         process::is_at(x) ||
         process::is_seq(x) ||
         process::is_if_then(x) ||
         process::is_if_then_else(x) ||
         process::is_bounded_init(x) ||
         process::is_merge(x) ||
         process::is_left_merge(x) ||
         process::is_choice(x) ||
         process::is_stochastic_operator(x) ||
         process::is_untyped_process_assignment(x);
}

// prototype declaration
std::string pp(const process_expression& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const process_expression& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(process_expression& t1, process_expression& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief An action
class action: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    action()
      : process_expression(core::detail::default_values::Action)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit action(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Action(*this));
    }

    /// \\brief Constructor Z14.
    action(const action_label& label, const data::data_expression_list& arguments)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Action(), label, arguments))
    {}

    /// Move semantics
    action(const action&) noexcept = default;
    action(action&&) noexcept = default;
    action& operator=(const action&) noexcept = default;
    action& operator=(action&&) noexcept = default;

    const action_label& label() const
    {
      return atermpp::down_cast<action_label>((*this)[0]);
    }

    const data::data_expression_list& arguments() const
    {
      return atermpp::down_cast<data::data_expression_list>((*this)[1]);
    }
};

/// \\brief Make_action constructs a new term into a given address.
/// \\ \param t The reference into which the new action is constructed. 
template <class... ARGUMENTS>
inline void make_action(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Action(), args...);
}

/// \\brief list of actions
using action_list = atermpp::term_list<action>;

/// \\brief vector of actions
using action_vector = std::vector<action>;

/// \\brief Test for a action expression
/// \\param x A term
/// \\return True if \\a x is a action expression
inline
bool is_action(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Action;
}

// prototype declaration
std::string pp(const action& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const action& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(action& t1, action& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief A process
class process_instance: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    process_instance()
      : process_expression(core::detail::default_values::Process)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit process_instance(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Process(*this));
    }

    /// \\brief Constructor Z14.
    process_instance(const process_identifier& identifier, const data::data_expression_list& actual_parameters)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Process(), identifier, actual_parameters))
    {}

    /// Move semantics
    process_instance(const process_instance&) noexcept = default;
    process_instance(process_instance&&) noexcept = default;
    process_instance& operator=(const process_instance&) noexcept = default;
    process_instance& operator=(process_instance&&) noexcept = default;

    const process_identifier& identifier() const
    {
      return atermpp::down_cast<process_identifier>((*this)[0]);
    }

    const data::data_expression_list& actual_parameters() const
    {
      return atermpp::down_cast<data::data_expression_list>((*this)[1]);
    }
};

/// \\brief Make_process_instance constructs a new term into a given address.
/// \\ \param t The reference into which the new process_instance is constructed. 
template <class... ARGUMENTS>
inline void make_process_instance(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Process(), args...);
}

/// \\brief Test for a process_instance expression
/// \\param x A term
/// \\return True if \\a x is a process_instance expression
inline
bool is_process_instance(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Process;
}

// prototype declaration
std::string pp(const process_instance& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const process_instance& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(process_instance& t1, process_instance& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief A process assignment
class process_instance_assignment: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    process_instance_assignment()
      : process_expression(core::detail::default_values::ProcessAssignment)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit process_instance_assignment(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_ProcessAssignment(*this));
    }

    /// \\brief Constructor Z14.
    process_instance_assignment(const process_identifier& identifier, const data::assignment_list& assignments)
      : process_expression(atermpp::aterm(core::detail::function_symbol_ProcessAssignment(), identifier, assignments))
    {}

    /// Move semantics
    process_instance_assignment(const process_instance_assignment&) noexcept = default;
    process_instance_assignment(process_instance_assignment&&) noexcept = default;
    process_instance_assignment& operator=(const process_instance_assignment&) noexcept = default;
    process_instance_assignment& operator=(process_instance_assignment&&) noexcept = default;

    const process_identifier& identifier() const
    {
      return atermpp::down_cast<process_identifier>((*this)[0]);
    }

    const data::assignment_list& assignments() const
    {
      return atermpp::down_cast<data::assignment_list>((*this)[1]);
    }
};

/// \\brief Make_process_instance_assignment constructs a new term into a given address.
/// \\ \param t The reference into which the new process_instance_assignment is constructed. 
template <class... ARGUMENTS>
inline void make_process_instance_assignment(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ProcessAssignment(), args...);
}

/// \\brief Test for a process_instance_assignment expression
/// \\param x A term
/// \\return True if \\a x is a process_instance_assignment expression
inline
bool is_process_instance_assignment(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ProcessAssignment;
}

// prototype declaration
std::string pp(const process_instance_assignment& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const process_instance_assignment& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(process_instance_assignment& t1, process_instance_assignment& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The value delta
class delta: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    delta()
      : process_expression(core::detail::default_values::Delta)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit delta(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Delta(*this));
    }

    /// Move semantics
    delta(const delta&) noexcept = default;
    delta(delta&&) noexcept = default;
    delta& operator=(const delta&) noexcept = default;
    delta& operator=(delta&&) noexcept = default;
};

/// \\brief Test for a delta expression
/// \\param x A term
/// \\return True if \\a x is a delta expression
inline
bool is_delta(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Delta;
}

// prototype declaration
std::string pp(const delta& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const delta& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(delta& t1, delta& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The value tau
class tau: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    tau()
      : process_expression(core::detail::default_values::Tau)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit tau(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Tau(*this));
    }

    /// Move semantics
    tau(const tau&) noexcept = default;
    tau(tau&&) noexcept = default;
    tau& operator=(const tau&) noexcept = default;
    tau& operator=(tau&&) noexcept = default;
};

/// \\brief Test for a tau expression
/// \\param x A term
/// \\return True if \\a x is a tau expression
inline
bool is_tau(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Tau;
}

// prototype declaration
std::string pp(const tau& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const tau& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(tau& t1, tau& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The sum operator
class sum: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    sum()
      : process_expression(core::detail::default_values::Sum)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit sum(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Sum(*this));
    }

    /// \\brief Constructor Z14.
    sum(const data::variable_list& variables, const process_expression& operand)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Sum(), variables, operand))
    {}

    /// Move semantics
    sum(const sum&) noexcept = default;
    sum(sum&&) noexcept = default;
    sum& operator=(const sum&) noexcept = default;
    sum& operator=(sum&&) noexcept = default;

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[0]);
    }

    const process_expression& operand() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_sum constructs a new term into a given address.
/// \\ \param t The reference into which the new sum is constructed. 
template <class... ARGUMENTS>
inline void make_sum(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Sum(), args...);
}

/// \\brief Test for a sum expression
/// \\param x A term
/// \\return True if \\a x is a sum expression
inline
bool is_sum(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Sum;
}

// prototype declaration
std::string pp(const sum& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const sum& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(sum& t1, sum& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The block operator
class block: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    block()
      : process_expression(core::detail::default_values::Block)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit block(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Block(*this));
    }

    /// \\brief Constructor Z14.
    block(const core::identifier_string_list& block_set, const process_expression& operand)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Block(), block_set, operand))
    {}

    /// Move semantics
    block(const block&) noexcept = default;
    block(block&&) noexcept = default;
    block& operator=(const block&) noexcept = default;
    block& operator=(block&&) noexcept = default;

    const core::identifier_string_list& block_set() const
    {
      return atermpp::down_cast<core::identifier_string_list>((*this)[0]);
    }

    const process_expression& operand() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_block constructs a new term into a given address.
/// \\ \param t The reference into which the new block is constructed. 
template <class... ARGUMENTS>
inline void make_block(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Block(), args...);
}

/// \\brief Test for a block expression
/// \\param x A term
/// \\return True if \\a x is a block expression
inline
bool is_block(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Block;
}

// prototype declaration
std::string pp(const block& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const block& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(block& t1, block& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The hide operator
class hide: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    hide()
      : process_expression(core::detail::default_values::Hide)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit hide(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Hide(*this));
    }

    /// \\brief Constructor Z14.
    hide(const core::identifier_string_list& hide_set, const process_expression& operand)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Hide(), hide_set, operand))
    {}

    /// Move semantics
    hide(const hide&) noexcept = default;
    hide(hide&&) noexcept = default;
    hide& operator=(const hide&) noexcept = default;
    hide& operator=(hide&&) noexcept = default;

    const core::identifier_string_list& hide_set() const
    {
      return atermpp::down_cast<core::identifier_string_list>((*this)[0]);
    }

    const process_expression& operand() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_hide constructs a new term into a given address.
/// \\ \param t The reference into which the new hide is constructed. 
template <class... ARGUMENTS>
inline void make_hide(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Hide(), args...);
}

/// \\brief Test for a hide expression
/// \\param x A term
/// \\return True if \\a x is a hide expression
inline
bool is_hide(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Hide;
}

// prototype declaration
std::string pp(const hide& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const hide& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(hide& t1, hide& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The rename operator
class rename: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    rename()
      : process_expression(core::detail::default_values::Rename)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit rename(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Rename(*this));
    }

    /// \\brief Constructor Z14.
    rename(const rename_expression_list& rename_set, const process_expression& operand)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Rename(), rename_set, operand))
    {}

    /// Move semantics
    rename(const rename&) noexcept = default;
    rename(rename&&) noexcept = default;
    rename& operator=(const rename&) noexcept = default;
    rename& operator=(rename&&) noexcept = default;

    const rename_expression_list& rename_set() const
    {
      return atermpp::down_cast<rename_expression_list>((*this)[0]);
    }

    const process_expression& operand() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_rename constructs a new term into a given address.
/// \\ \param t The reference into which the new rename is constructed. 
template <class... ARGUMENTS>
inline void make_rename(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Rename(), args...);
}

/// \\brief Test for a rename expression
/// \\param x A term
/// \\return True if \\a x is a rename expression
inline
bool is_rename(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Rename;
}

// prototype declaration
std::string pp(const rename& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const rename& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(rename& t1, rename& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The communication operator
class comm: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    comm()
      : process_expression(core::detail::default_values::Comm)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit comm(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Comm(*this));
    }

    /// \\brief Constructor Z14.
    comm(const communication_expression_list& comm_set, const process_expression& operand)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Comm(), comm_set, operand))
    {}

    /// Move semantics
    comm(const comm&) noexcept = default;
    comm(comm&&) noexcept = default;
    comm& operator=(const comm&) noexcept = default;
    comm& operator=(comm&&) noexcept = default;

    const communication_expression_list& comm_set() const
    {
      return atermpp::down_cast<communication_expression_list>((*this)[0]);
    }

    const process_expression& operand() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_comm constructs a new term into a given address.
/// \\ \param t The reference into which the new comm is constructed. 
template <class... ARGUMENTS>
inline void make_comm(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Comm(), args...);
}

/// \\brief Test for a comm expression
/// \\param x A term
/// \\return True if \\a x is a comm expression
inline
bool is_comm(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Comm;
}

// prototype declaration
std::string pp(const comm& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const comm& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(comm& t1, comm& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The allow operator
class allow: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    allow()
      : process_expression(core::detail::default_values::Allow)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit allow(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Allow(*this));
    }

    /// \\brief Constructor Z14.
    allow(const action_name_multiset_list& allow_set, const process_expression& operand)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Allow(), allow_set, operand))
    {}

    /// Move semantics
    allow(const allow&) noexcept = default;
    allow(allow&&) noexcept = default;
    allow& operator=(const allow&) noexcept = default;
    allow& operator=(allow&&) noexcept = default;

    const action_name_multiset_list& allow_set() const
    {
      return atermpp::down_cast<action_name_multiset_list>((*this)[0]);
    }

    const process_expression& operand() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_allow constructs a new term into a given address.
/// \\ \param t The reference into which the new allow is constructed. 
template <class... ARGUMENTS>
inline void make_allow(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Allow(), args...);
}

/// \\brief Test for a allow expression
/// \\param x A term
/// \\return True if \\a x is a allow expression
inline
bool is_allow(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Allow;
}

// prototype declaration
std::string pp(const allow& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const allow& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(allow& t1, allow& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The synchronization operator
class sync: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    sync()
      : process_expression(core::detail::default_values::Sync)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit sync(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Sync(*this));
    }

    /// \\brief Constructor Z14.
    sync(const process_expression& left, const process_expression& right)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Sync(), left, right))
    {}

    /// Move semantics
    sync(const sync&) noexcept = default;
    sync(sync&&) noexcept = default;
    sync& operator=(const sync&) noexcept = default;
    sync& operator=(sync&&) noexcept = default;

    const process_expression& left() const
    {
      return atermpp::down_cast<process_expression>((*this)[0]);
    }

    const process_expression& right() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_sync constructs a new term into a given address.
/// \\ \param t The reference into which the new sync is constructed. 
template <class... ARGUMENTS>
inline void make_sync(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Sync(), args...);
}

/// \\brief Test for a sync expression
/// \\param x A term
/// \\return True if \\a x is a sync expression
inline
bool is_sync(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Sync;
}

// prototype declaration
std::string pp(const sync& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const sync& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(sync& t1, sync& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The at operator
class at: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    at()
      : process_expression(core::detail::default_values::AtTime)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit at(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_AtTime(*this));
    }

    /// \\brief Constructor Z14.
    at(const process_expression& operand, const data::data_expression& time_stamp)
      : process_expression(atermpp::aterm(core::detail::function_symbol_AtTime(), operand, time_stamp))
    {}

    /// Move semantics
    at(const at&) noexcept = default;
    at(at&&) noexcept = default;
    at& operator=(const at&) noexcept = default;
    at& operator=(at&&) noexcept = default;

    const process_expression& operand() const
    {
      return atermpp::down_cast<process_expression>((*this)[0]);
    }

    const data::data_expression& time_stamp() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[1]);
    }
};

/// \\brief Make_at constructs a new term into a given address.
/// \\ \param t The reference into which the new at is constructed. 
template <class... ARGUMENTS>
inline void make_at(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_AtTime(), args...);
}

/// \\brief Test for a at expression
/// \\param x A term
/// \\return True if \\a x is a at expression
inline
bool is_at(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::AtTime;
}

// prototype declaration
std::string pp(const at& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const at& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(at& t1, at& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The sequential composition
class seq: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    seq()
      : process_expression(core::detail::default_values::Seq)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit seq(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Seq(*this));
    }

    /// \\brief Constructor Z14.
    seq(const process_expression& left, const process_expression& right)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Seq(), left, right))
    {}

    /// Move semantics
    seq(const seq&) noexcept = default;
    seq(seq&&) noexcept = default;
    seq& operator=(const seq&) noexcept = default;
    seq& operator=(seq&&) noexcept = default;

    const process_expression& left() const
    {
      return atermpp::down_cast<process_expression>((*this)[0]);
    }

    const process_expression& right() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_seq constructs a new term into a given address.
/// \\ \param t The reference into which the new seq is constructed. 
template <class... ARGUMENTS>
inline void make_seq(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Seq(), args...);
}

/// \\brief Test for a seq expression
/// \\param x A term
/// \\return True if \\a x is a seq expression
inline
bool is_seq(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Seq;
}

// prototype declaration
std::string pp(const seq& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const seq& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(seq& t1, seq& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The if-then operator
class if_then: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    if_then()
      : process_expression(core::detail::default_values::IfThen)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit if_then(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_IfThen(*this));
    }

    /// \\brief Constructor Z14.
    if_then(const data::data_expression& condition, const process_expression& then_case)
      : process_expression(atermpp::aterm(core::detail::function_symbol_IfThen(), condition, then_case))
    {}

    /// Move semantics
    if_then(const if_then&) noexcept = default;
    if_then(if_then&&) noexcept = default;
    if_then& operator=(const if_then&) noexcept = default;
    if_then& operator=(if_then&&) noexcept = default;

    const data::data_expression& condition() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[0]);
    }

    const process_expression& then_case() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_if_then constructs a new term into a given address.
/// \\ \param t The reference into which the new if_then is constructed. 
template <class... ARGUMENTS>
inline void make_if_then(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_IfThen(), args...);
}

/// \\brief Test for a if_then expression
/// \\param x A term
/// \\return True if \\a x is a if_then expression
inline
bool is_if_then(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::IfThen;
}

// prototype declaration
std::string pp(const if_then& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const if_then& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(if_then& t1, if_then& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The if-then-else operator
class if_then_else: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    if_then_else()
      : process_expression(core::detail::default_values::IfThenElse)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit if_then_else(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_IfThenElse(*this));
    }

    /// \\brief Constructor Z14.
    if_then_else(const data::data_expression& condition, const process_expression& then_case, const process_expression& else_case)
      : process_expression(atermpp::aterm(core::detail::function_symbol_IfThenElse(), condition, then_case, else_case))
    {}

    /// Move semantics
    if_then_else(const if_then_else&) noexcept = default;
    if_then_else(if_then_else&&) noexcept = default;
    if_then_else& operator=(const if_then_else&) noexcept = default;
    if_then_else& operator=(if_then_else&&) noexcept = default;

    const data::data_expression& condition() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[0]);
    }

    const process_expression& then_case() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }

    const process_expression& else_case() const
    {
      return atermpp::down_cast<process_expression>((*this)[2]);
    }
};

/// \\brief Make_if_then_else constructs a new term into a given address.
/// \\ \param t The reference into which the new if_then_else is constructed. 
template <class... ARGUMENTS>
inline void make_if_then_else(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_IfThenElse(), args...);
}

/// \\brief Test for a if_then_else expression
/// \\param x A term
/// \\return True if \\a x is a if_then_else expression
inline
bool is_if_then_else(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::IfThenElse;
}

// prototype declaration
std::string pp(const if_then_else& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const if_then_else& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(if_then_else& t1, if_then_else& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The bounded initialization
class bounded_init: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    bounded_init()
      : process_expression(core::detail::default_values::BInit)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit bounded_init(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_BInit(*this));
    }

    /// \\brief Constructor Z14.
    bounded_init(const process_expression& left, const process_expression& right)
      : process_expression(atermpp::aterm(core::detail::function_symbol_BInit(), left, right))
    {}

    /// Move semantics
    bounded_init(const bounded_init&) noexcept = default;
    bounded_init(bounded_init&&) noexcept = default;
    bounded_init& operator=(const bounded_init&) noexcept = default;
    bounded_init& operator=(bounded_init&&) noexcept = default;

    const process_expression& left() const
    {
      return atermpp::down_cast<process_expression>((*this)[0]);
    }

    const process_expression& right() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_bounded_init constructs a new term into a given address.
/// \\ \param t The reference into which the new bounded_init is constructed. 
template <class... ARGUMENTS>
inline void make_bounded_init(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_BInit(), args...);
}

/// \\brief Test for a bounded_init expression
/// \\param x A term
/// \\return True if \\a x is a bounded_init expression
inline
bool is_bounded_init(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::BInit;
}

// prototype declaration
std::string pp(const bounded_init& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const bounded_init& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(bounded_init& t1, bounded_init& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The merge operator
class merge: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    merge()
      : process_expression(core::detail::default_values::Merge)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit merge(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Merge(*this));
    }

    /// \\brief Constructor Z14.
    merge(const process_expression& left, const process_expression& right)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Merge(), left, right))
    {}

    /// Move semantics
    merge(const merge&) noexcept = default;
    merge(merge&&) noexcept = default;
    merge& operator=(const merge&) noexcept = default;
    merge& operator=(merge&&) noexcept = default;

    const process_expression& left() const
    {
      return atermpp::down_cast<process_expression>((*this)[0]);
    }

    const process_expression& right() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_merge constructs a new term into a given address.
/// \\ \param t The reference into which the new merge is constructed. 
template <class... ARGUMENTS>
inline void make_merge(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Merge(), args...);
}

/// \\brief Test for a merge expression
/// \\param x A term
/// \\return True if \\a x is a merge expression
inline
bool is_merge(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Merge;
}

// prototype declaration
std::string pp(const merge& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const merge& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(merge& t1, merge& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The left merge operator
class left_merge: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    left_merge()
      : process_expression(core::detail::default_values::LMerge)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit left_merge(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_LMerge(*this));
    }

    /// \\brief Constructor Z14.
    left_merge(const process_expression& left, const process_expression& right)
      : process_expression(atermpp::aterm(core::detail::function_symbol_LMerge(), left, right))
    {}

    /// Move semantics
    left_merge(const left_merge&) noexcept = default;
    left_merge(left_merge&&) noexcept = default;
    left_merge& operator=(const left_merge&) noexcept = default;
    left_merge& operator=(left_merge&&) noexcept = default;

    const process_expression& left() const
    {
      return atermpp::down_cast<process_expression>((*this)[0]);
    }

    const process_expression& right() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_left_merge constructs a new term into a given address.
/// \\ \param t The reference into which the new left_merge is constructed. 
template <class... ARGUMENTS>
inline void make_left_merge(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_LMerge(), args...);
}

/// \\brief Test for a left_merge expression
/// \\param x A term
/// \\return True if \\a x is a left_merge expression
inline
bool is_left_merge(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::LMerge;
}

// prototype declaration
std::string pp(const left_merge& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const left_merge& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(left_merge& t1, left_merge& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The choice operator
class choice: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    choice()
      : process_expression(core::detail::default_values::Choice)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit choice(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_Choice(*this));
    }

    /// \\brief Constructor Z14.
    choice(const process_expression& left, const process_expression& right)
      : process_expression(atermpp::aterm(core::detail::function_symbol_Choice(), left, right))
    {}

    /// Move semantics
    choice(const choice&) noexcept = default;
    choice(choice&&) noexcept = default;
    choice& operator=(const choice&) noexcept = default;
    choice& operator=(choice&&) noexcept = default;

    const process_expression& left() const
    {
      return atermpp::down_cast<process_expression>((*this)[0]);
    }

    const process_expression& right() const
    {
      return atermpp::down_cast<process_expression>((*this)[1]);
    }
};

/// \\brief Make_choice constructs a new term into a given address.
/// \\ \param t The reference into which the new choice is constructed. 
template <class... ARGUMENTS>
inline void make_choice(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Choice(), args...);
}

/// \\brief Test for a choice expression
/// \\param x A term
/// \\return True if \\a x is a choice expression
inline
bool is_choice(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Choice;
}

// prototype declaration
std::string pp(const choice& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const choice& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(choice& t1, choice& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The distribution operator
class stochastic_operator: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    stochastic_operator()
      : process_expression(core::detail::default_values::StochasticOperator)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit stochastic_operator(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_StochasticOperator(*this));
    }

    /// \\brief Constructor Z14.
    stochastic_operator(const data::variable_list& variables, const data::data_expression& distribution, const process_expression& operand)
      : process_expression(atermpp::aterm(core::detail::function_symbol_StochasticOperator(), variables, distribution, operand))
    {}

    /// Move semantics
    stochastic_operator(const stochastic_operator&) noexcept = default;
    stochastic_operator(stochastic_operator&&) noexcept = default;
    stochastic_operator& operator=(const stochastic_operator&) noexcept = default;
    stochastic_operator& operator=(stochastic_operator&&) noexcept = default;

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[0]);
    }

    const data::data_expression& distribution() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[1]);
    }

    const process_expression& operand() const
    {
      return atermpp::down_cast<process_expression>((*this)[2]);
    }
};

/// \\brief Make_stochastic_operator constructs a new term into a given address.
/// \\ \param t The reference into which the new stochastic_operator is constructed. 
template <class... ARGUMENTS>
inline void make_stochastic_operator(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StochasticOperator(), args...);
}

/// \\brief Test for a stochastic_operator expression
/// \\param x A term
/// \\return True if \\a x is a stochastic_operator expression
inline
bool is_stochastic_operator(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StochasticOperator;
}

// prototype declaration
std::string pp(const stochastic_operator& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const stochastic_operator& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(stochastic_operator& t1, stochastic_operator& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief An untyped process assginment
class untyped_process_assignment: public process_expression
{
  public:
    /// \\brief Default constructor X3.
    untyped_process_assignment()
      : process_expression(core::detail::default_values::UntypedProcessAssignment)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit untyped_process_assignment(const atermpp::aterm& term)
      : process_expression(term)
    {
      assert(core::detail::check_term_UntypedProcessAssignment(*this));
    }

    /// \\brief Constructor Z14.
    untyped_process_assignment(const core::identifier_string& name, const data::untyped_identifier_assignment_list& assignments)
      : process_expression(atermpp::aterm(core::detail::function_symbol_UntypedProcessAssignment(), name, assignments))
    {}

    /// \\brief Constructor Z2.
    untyped_process_assignment(const std::string& name, const data::untyped_identifier_assignment_list& assignments)
      : process_expression(atermpp::aterm(core::detail::function_symbol_UntypedProcessAssignment(), core::identifier_string(name), assignments))
    {}

    /// Move semantics
    untyped_process_assignment(const untyped_process_assignment&) noexcept = default;
    untyped_process_assignment(untyped_process_assignment&&) noexcept = default;
    untyped_process_assignment& operator=(const untyped_process_assignment&) noexcept = default;
    untyped_process_assignment& operator=(untyped_process_assignment&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const data::untyped_identifier_assignment_list& assignments() const
    {
      return atermpp::down_cast<data::untyped_identifier_assignment_list>((*this)[1]);
    }
};

/// \\brief Make_untyped_process_assignment constructs a new term into a given address.
/// \\ \param t The reference into which the new untyped_process_assignment is constructed. 
template <class... ARGUMENTS>
inline void make_untyped_process_assignment(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_UntypedProcessAssignment(), args...);
}

/// \\brief Test for a untyped_process_assignment expression
/// \\param x A term
/// \\return True if \\a x is a untyped_process_assignment expression
inline
bool is_untyped_process_assignment(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::UntypedProcessAssignment;
}

// prototype declaration
std::string pp(const untyped_process_assignment& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const untyped_process_assignment& x)
{
  return out << process::pp(x);
}

/// \\brief swap overload
inline void swap(untyped_process_assignment& t1, untyped_process_assignment& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated classes ---//

// template function overloads
std::string pp(const process_expression_list& x, bool precedence_aware = true);
std::string pp(const process_expression_vector& x, bool precedence_aware = true);
std::set<data::sort_expression> find_sort_expressions(const process::process_expression& x);
std::string pp(const action_list& x, bool precedence_aware = true);
std::string pp(const action_vector& x, bool precedence_aware = true);
action normalize_sorts(const action& x, const data::sort_specification& sortspec);
action translate_user_notation(const action& x);
process::process_expression translate_user_notation(const process::process_expression& x);
std::set<data::variable> find_all_variables(const action& x);
std::set<data::variable> find_free_variables(const action& x);

/// \brief Compares the signatures of two actions
/// \param a An action
/// \param b An action
/// \return Returns true if the actions a and b have the same label, and
/// the sorts of the arguments of a and b are equal.
inline
bool equal_signatures(const action& a, const action& b)
{
  if (a.label() != b.label())
  {
    return false;
  }

  const data::data_expression_list& a_args = a.arguments();
  const data::data_expression_list& b_args = b.arguments();

  if (a_args.size() != b_args.size())
  {
    return false;
  }

  return std::equal(a_args.begin(), a_args.end(), b_args.begin(), [](const data::data_expression& x, const data::data_expression& y) { return x.sort() == y.sort(); });
}

} // namespace mcrl2::process

namespace std
{

/// \brief Standard has function for actions.
template <>
struct hash<mcrl2::process::action>
{
  std::size_t operator()(const mcrl2::process::action& t) const
  {
    return std::hash<atermpp::aterm>()(t);
  }

};


} // namespace std

#endif // MCRL2_PROCESS_PROCESS_EXPRESSION_H
