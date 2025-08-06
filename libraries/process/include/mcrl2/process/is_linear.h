// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/is_linear.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_IS_LINEAR_H
#define MCRL2_PROCESS_IS_LINEAR_H

#include <iostream>

#include "mcrl2/process/traverser.h"

namespace mcrl2::process
{

namespace detail
{

/// \brief Returns true if the process instance assignment a matches with the process equation eq.
inline
bool check_process_instance_assignment(const process_equation& eq, const process_instance_assignment& inst)
{
  if (inst.identifier() != eq.identifier())
  {
    return false;
  }
  const data::variable_list& v = eq.formal_parameters();

  // check if the left hand sides of the assignments exist
  for (const auto& a: inst.assignments())
  {
    if (std::find(v.begin(), v.end(), a.lhs()) == v.end())
    {
      return false;
    }
  }
  return true;
}

/// \brief Returns true if the process instance a matches with the process equation eq.
inline
bool check_process_instance(const process_equation& eq, const process_instance& init)
{
  if (eq.identifier() != init.identifier())
  {
    return false;
  }
  const data::variable_list& v = eq.formal_parameters();
  const data::data_expression_list& e = init.actual_parameters();
  data::variable_list::const_iterator i = v.begin();
  data::data_expression_list::const_iterator j = e.begin();
  for (; i != v.end(); ++i, ++j)
  {
    if (i->sort() != j->sort())
    {
      return false;
    }
  }
  return true;
}

/// \brief Returns true if the argument is a process instance
inline
bool is_process(const process_expression& x)
{
  return is_process_instance(x)
         || is_process_instance_assignment(x)
         ;
}

/// \brief Returns true if the argument is a process instance, optionally wrapped in a stochastic
/// distribution.
inline
bool is_stochastic_process(const process_expression& x)
{
  if (is_process(x))
  {
    return true;
  }
  if (is_stochastic_operator(x))
  {
    return is_process(atermpp::down_cast<stochastic_operator>(x).operand());
  }
  return false;
}

/// \brief Returns true if the argument is a deadlock
inline
bool is_timed_deadlock(const process_expression& x)
{
  return is_delta(x)
         || is_at(x)
         ;
}

/// \brief Returns true if the argument is a multi-action
inline
bool is_multiaction(const process_expression& x)
{
  return is_tau(x)
         || is_sync(x)
         || is_action(x)
         ;
}

/// \brief Returns true if the argument is a multi-action
inline
bool is_timed_multiaction(const process_expression& x)
{
  return is_at(x)
         || is_multiaction(x);
}

/// \brief Returns true if the argument is an action prefix
inline
bool is_action_prefix(const process_expression& x)
{
  return is_seq(x)
         || is_timed_multiaction(x);
}

/// \brief Returns true if the argument is a conditional deadlock
inline
bool is_conditional_deadlock(const process_expression& x)
{
  return is_if_then(x)
         || is_timed_deadlock(x);
}

/// \brief Returns true if the argument is a conditional action prefix.
inline
bool is_conditional_action_prefix(const process_expression& x)
{
  return is_if_then(x)
         || is_action_prefix(x);
}

/// \brief Returns true if the argument is an alternative composition
inline
bool is_alternative(const process_expression& x)
{
  return is_sum(x)
         || is_conditional_action_prefix(x)
         || is_conditional_deadlock(x)
         ;
}

/// \brief Returns true if the argument is a linear process
inline
bool is_linear_process_term(const process_expression& x)
{
  return is_choice(x)
         || is_alternative(x)
         ;
}

/// \brief Checks if a process equation is linear.
/// Use the is_linear() member function for this.
struct linear_process_expression_traverser: public process_expression_traverser<linear_process_expression_traverser>
{
  using super = process_expression_traverser<linear_process_expression_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  /// \brief The process equation that is checked.
  process_equation eqn;

  /// \brief Exception that is thrown by linear_process_expression_traverser.
  struct non_linear_process_error: public mcrl2::runtime_error
  {
    explicit non_linear_process_error(const std::string& msg)
      : mcrl2::runtime_error(msg)
    {}
  };

  linear_process_expression_traverser(const process_equation& eqn_ = process_equation())
   : eqn(eqn_)
  {}

  // TODO: check if this function is needed
  void enter(const process::process_instance& x)
  {
    if (!detail::check_process_instance(eqn, x))
    {
      throw non_linear_process_error(process::pp(x) + " is not a valid process instance");
    }
  }

  void enter(const process::process_instance_assignment& x)
  {
    if (!detail::check_process_instance_assignment(eqn, x))
    {
      throw non_linear_process_error(process::pp(x) + " is not a valid process instance assignment");
    }
  }

  void enter(const process::sum& x)
  {
    if (!is_alternative(x.operand()))
    {
      throw non_linear_process_error(process::pp(x.operand()) + " is not an alternative expression");
    }
  }

  void enter(const process::block& x)
  {
    throw non_linear_process_error("block expression " + process::pp(x) + " encountered");
  }

  void enter(const process::hide& x)
  {
    throw non_linear_process_error("hide expression " + process::pp(x) + " encountered");
  }

  void enter(const process::rename& x)
  {
    throw non_linear_process_error("rename expression " + process::pp(x) + " encountered");
  }

  void enter(const process::comm& x)
  {
    throw non_linear_process_error("comm expression " + process::pp(x) + " encountered");
  }

  void enter(const process::allow& x)
  {
    throw non_linear_process_error("allow expression " + process::pp(x) + " encountered");
  }

  void enter(const process::sync& x)
  {
    if (!is_multiaction(x.left()) || !is_multiaction(x.right()))
    {
      if (!is_multiaction(x.left()))
      {
        throw non_linear_process_error(process::pp(x.left()) + " is not a multi action");
      }
      else
      {
        throw non_linear_process_error(process::pp(x.right()) + " is not a multi action");
      }
    }
  }

  void enter(const process::at& x)
  {
    if (!is_multiaction(x.operand()) && !is_delta(x.operand()))
    {
      throw non_linear_process_error(process::pp(x.operand()) + " is not a multi action and not a deadlock");
    }
  }

  void enter(const process::seq& x)
  {
    if (!is_timed_multiaction(x.left()) || !is_stochastic_process(x.right()))
    {
      throw non_linear_process_error(process::pp(x.left()) + " is not a timed multi action and not a process");
    }
    process_expression right = x.right();
    if (is_stochastic_operator(right))
    {
      right = atermpp::down_cast<stochastic_operator>(right).operand();
    }
    if (is_process_instance(right))
    {
      const auto& right_ = atermpp::down_cast<process_instance>(right);
      if (right_.identifier() != eqn.identifier())
      {
        throw non_linear_process_error(process::pp(right_) + " has an unexpected identifier");
      }
    }
    else if (is_process_instance_assignment(right))
    {
      const auto& right_ = atermpp::down_cast<process_instance_assignment>(right);
      if (right_.identifier() != eqn.identifier())
      {
        throw non_linear_process_error(process::pp(right_) + " has an unexpected identifier");
      }
    }
    else
    {
      std::cerr << "seq right hand side: " << process::pp(right) << std::endl;
      throw std::runtime_error("unexpected error in visit_seq");
    }
  }

  void enter(const process::if_then& x)
  {
    if (!is_action_prefix(x.then_case()) && !is_timed_deadlock(x.then_case()))
    {
      throw non_linear_process_error(process::pp(x) + " is not an action prefix and not a timed deadlock");
    }
  }

  void enter(const process::if_then_else& x)
  {
    throw non_linear_process_error("if then else expression " + process::pp(x) + " encountered");
  }

  void enter(const process::bounded_init& x)
  {
    throw non_linear_process_error("bounded init expression " + process::pp(x) + " encountered");
  }

  void enter(const process::merge& x)
  {
    throw non_linear_process_error("merge expression " + process::pp(x) + " encountered");
  }

  void enter(const process::left_merge& x)
  {
    throw non_linear_process_error("left merge expression " + process::pp(x) + " encountered");
  }


  /// \brief Returns true if the process equation e is linear.
  bool is_linear(const process_expression& x, bool verbose = false)
  {
    try
    {
      this->apply(x);
    }
    catch (non_linear_process_error& p)
    {
      if (verbose)
      {
        std::clog << p.what() << std::endl;
      }
      return false;
    }
    return true;
  }
};

} // namespace detail

/// \brief Returns true if the process specification is linear.
inline
bool is_linear(const process_specification& p, bool verbose = false)
{
  if (p.equations().size() != 1)
  {
    if (verbose)
    {
      std::clog << "warning: the number of equations is not equal to 1" << std::endl;
    }
    return false;
  }
  const process_equation& eqn = p.equations().front();
  detail::linear_process_expression_traverser visitor(eqn);
  {
    if (!visitor.is_linear(eqn.expression(), verbose))
    {
      if (verbose)
      {
        std::clog << "warning: the first equation is not linear" << std::endl;
      }
      return false;
    }
    if (!detail::is_stochastic_process(p.init()))
    {
      if (verbose)
      {
        std::clog << "warning: the initial process " << process::pp(p.init()) << " is not a process instance or a process instance assignment" << std::endl;
      }
      return false;
    }
  }
  return true;
}

/// \brief Returns true if the process equation is linear.
inline
bool is_linear(const process_equation& eqn)
{
  detail::linear_process_expression_traverser f(eqn);
  return f.is_linear(eqn.expression());
}

/// \brief Returns true if the process expression is linear.
/// \param x A process expression.
/// \param eqn The linear equation belonging to the indicated process.
inline
bool is_linear(const process_expression& x, const process_equation& eqn)
{
  detail::linear_process_expression_traverser f(eqn);
  return f.is_linear(x);
}

} // namespace mcrl2::process

#endif // MCRL2_PROCESS_IS_LINEAR_H
