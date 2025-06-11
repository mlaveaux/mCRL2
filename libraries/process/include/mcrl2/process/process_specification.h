// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/process/process_specification.h
/// \brief add your file description here.

#ifndef MCRL2_PROCESS_PROCESS_SPECIFICATION_H
#define MCRL2_PROCESS_PROCESS_SPECIFICATION_H

#include "mcrl2/core/load_aterm.h"
#include "mcrl2/core/print.h"
#include "mcrl2/data/data_io.h"
#include "mcrl2/data/detail/io.h"
#include "mcrl2/process/process_equation.h"

namespace mcrl2
{

/// \brief The main namespace for the Process library.
namespace process
{

class process_specification;
atermpp::aterm process_specification_to_aterm(const process_specification& spec);
void complete_data_specification(process_specification&);

// template function overloads
std::string pp(const process_specification& x, bool precedence_aware);
void normalize_sorts(process_specification& x, const data::sort_specification& sortspec);
void translate_user_notation(process::process_specification& x);
std::set<data::sort_expression> find_sort_expressions(const process::process_specification& x);
std::set<core::identifier_string> find_identifiers(const process::process_specification& x);
std::set<data::variable> find_free_variables(const process::process_specification& x);

/// \brief Test for a process specification expression
/// \param x A term
/// \return True if \a x is a process specification expression
inline
bool is_process_specification(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ProcSpec;
}

/// \brief Process specification consisting of a data specification, action labels, a sequence of process equations and a process initialization.
//<ProcSpec>     ::= ProcSpec(<DataSpec>, <ActSpec>, <GlobVarSpec>, <ProcEqnSpec>, <ProcInit>)
class process_specification
{
  protected:
    /// \brief The data specification of the specification
    data::data_specification m_data;

    /// \brief The action specification of the specification
    process::action_label_list m_action_labels;

    /// \brief The set of global variables
    std::set<data::variable> m_global_variables;

    /// \brief The equations of the specification
    std::vector<process_equation> m_equations;

    /// \brief The initial state of the specification
    process_expression m_initial_process;

    /// \brief Initializes the specification with an aterm.
    /// \param t A term
    void construct_from_aterm(const atermpp::aterm& t)
    {
      atermpp::aterm::iterator i = t.begin();
      m_data            = data::data_specification(*i++);
      m_action_labels   = atermpp::down_cast<process::action_label_list>((*i++)[0]);
      data::variable_list global_variables = atermpp::down_cast<data::variable_list>((*i++)[0]);
      m_global_variables = std::set<data::variable>(global_variables.begin(),global_variables.end());
      process_equation_list l = atermpp::down_cast<process_equation_list>((*i++)[0]);
      atermpp::aterm init = (*i);
      m_initial_process = process_expression(init[0]);
      m_equations       = std::vector<process_equation>(l.begin(), l.end());
    }

  public:
    /// \brief Constructor.
    process_specification()
    {}

    /// \brief Constructor.
    /// \param t A term containing an aterm representation of a process specification.
    process_specification(atermpp::aterm t)
    {
      assert(core::detail::check_term_ProcSpec(t));
      construct_from_aterm(t);
      complete_data_specification(*this);
    }

    /// \brief Constructor that sets the global variables to empty;
    process_specification(data::data_specification data, process::action_label_list action_labels, process_equation_list equations, process_expression init)
      : m_data(data),
        m_action_labels(action_labels),
        m_equations(equations.begin(), equations.end()),
        m_initial_process(init)
    {}

    /// \brief Constructor of a process specification.
    process_specification(
      data::data_specification data,
      process::action_label_list action_labels,
      data::variable_list global_variables,
      process_equation_list equations,
      process_expression init)
      : m_data(data),
        m_action_labels(action_labels),
        m_global_variables(global_variables.begin(),global_variables.end()),
        m_equations(equations.begin(), equations.end()),
        m_initial_process(init)
    {}

    /// \brief Returns the data specification
    /// \return The data specification
    const data::data_specification& data() const
    {
      return m_data;
    }

    /// \brief Returns the data specification
    /// \return The data specification
    data::data_specification& data()
    {
      return m_data;
    }

    /// \brief Returns the action label specification
    /// \return The action label specification
    const process::action_label_list& action_labels() const
    {
      return m_action_labels;
    }

    /// \brief Returns the action label specification
    /// \return The action label specification
    process::action_label_list& action_labels()
    {
      return m_action_labels;
    }

    /// \brief Returns the declared free variables of the process specification.
    /// \return The declared free variables of the process specification.
    const std::set<data::variable>& global_variables() const
    {
      return m_global_variables;
    }

    /// \brief Returns the declared free variables of the process specification.
    /// \return The declared free variables of the process specification.
    std::set<data::variable>& global_variables()
    {
      return m_global_variables;
    }

    /// \brief Returns the equations of the process specification
    /// \return The equations of the process specification
    const std::vector<process_equation>& equations() const
    {
      return m_equations;
    }

    /// \brief Returns the equations of the process specification
    /// \return The equations of the process specification
    std::vector<process_equation>& equations()
    {
      return m_equations;
    }

    /// \brief Returns the initialization of the process specification
    /// \return The initialization of the process specification
    const process_expression& init() const
    {
      return m_initial_process;
    }

    /// \brief Returns the initialization of the process specification
    /// \return The initialization of the process specification
    process_expression& init()
    {
      return m_initial_process;
    }
};

//--- start generated class process_specification ---//
// prototype declaration
std::string pp(const process_specification& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const process_specification& x)
{
  return out << process::pp(x);
}
//--- end generated class process_specification ---//

/// \brief Adds all sorts that appear in the process specification spec
///  to the data specification of spec.
/// \param spec A process specification
inline
void complete_data_specification(process_specification& spec)
{
  std::set<data::sort_expression> s = process::find_sort_expressions(spec);
  spec.data().add_context_sorts(s);
}

/// \brief Conversion to aterm.
/// \return The specification converted to aterm format.
/// \param spec A process specification
inline
atermpp::aterm process_specification_to_aterm(const process_specification& spec)
{
  return atermpp::aterm(core::detail::function_symbol_ProcSpec(),
           data::detail::data_specification_to_aterm(spec.data()),
           atermpp::aterm(core::detail::function_symbol_ActSpec(), spec.action_labels()),
           atermpp::aterm(core::detail::function_symbol_GlobVarSpec(), data::variable_list(spec.global_variables().begin(),spec.global_variables().end())),
           atermpp::aterm(core::detail::function_symbol_ProcEqnSpec(), process_equation_list(spec.equations().begin(), spec.equations().end())),
           atermpp::aterm(core::detail::function_symbol_ProcessInit(), spec.init())
         );
}

/// \brief Equality operator
inline
bool operator==(const process_specification& spec1, const process_specification& spec2)
{
  return process_specification_to_aterm(spec1) == process_specification_to_aterm(spec2);
}

/// \brief Inequality operator
inline
bool operator!=(const process_specification& spec1, const process_specification& spec2)
{
  return !(spec1 == spec2);
}

} // namespace process

} // namespace mcrl2

#endif // MCRL2_PROCESS_PROCESS_SPECIFICATION_H


