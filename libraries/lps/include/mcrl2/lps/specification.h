// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/specification.h
/// \brief The class specification.

#ifndef MCRL2_LPS_SPECIFICATION_H
#define MCRL2_LPS_SPECIFICATION_H

#include "mcrl2/data/data_io.h"
#include "mcrl2/data/detail/io.h"
#include "mcrl2/lps/process_initializer.h"
#include "mcrl2/lps/linear_process.h"

#include <cerrno>

namespace mcrl2
{

/// \brief The main namespace for the LPS library.
namespace lps
{

// prototype definitions
template <typename Object> bool check_well_typedness(const Object& o);
template <typename LinearProcess, typename InitialProcessExpression> class specification_base;
template <typename LinearProcess, typename InitialProcessExpression> atermpp::aterm specification_to_aterm(const specification_base<LinearProcess, InitialProcessExpression>& spec);
class specification;
void complete_data_specification(specification& spec);

// template function overloads
bool check_well_typedness(const specification& x);
void normalize_sorts(lps::specification& x, const data::sort_specification& sortspec);

/// \brief Test for a specification expression
/// \param x A term
/// \return True if \a x is a specification expression
inline
bool is_specification(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::LinProcSpec;
}

template <typename LinearProcess, typename InitialProcessExpression>
class specification_base
{
  protected:
    /// \brief The data specification of the specification
    data::data_specification m_data;

    /// \brief The action specification of the specification
    process::action_label_list m_action_labels;

    /// \brief The set of global variables
    std::set<data::variable> m_global_variables;

    /// \brief The linear process of the specification
    LinearProcess m_process;

    /// \brief The initial state of the specification
    InitialProcessExpression m_initial_process;

  public:
    /// \brief The process type
    typedef LinearProcess process_type;

    /// \brief The initial process type
    typedef InitialProcessExpression initial_process_type;

    /// \brief Constructor.
    specification_base()
    { }

    /// \brief Constructor.
    /// \param data A data specification
    /// \param action_labels A sequence of action labels
    /// \param global_variables A set of global variables
    /// \param lps A linear process
    /// \param initial_process A process initializer
    specification_base(const data::data_specification& data,
                       const process::action_label_list& action_labels,
                       const std::set<data::variable>& global_variables,
                       const LinearProcess& lps,
                       const InitialProcessExpression& initial_process)
       :
      m_data(data),
      m_action_labels(action_labels),
      m_global_variables(global_variables),
      m_process(lps),
      m_initial_process(initial_process)
    {
      assert(lps.process_parameters().size()==initial_process.expressions().size());
    }

    /// \brief Returns the linear process of the specification.
    /// \return The linear process of the specification.
    const LinearProcess& process() const
    {
      return m_process;
    }

    /// \brief Returns a reference to the linear process of the specification.
    /// \return The linear process of the specification.
    LinearProcess& process()
    {
      return m_process;
    }

    /// \brief Returns the data specification.
    /// \return The data specification.
    const data::data_specification& data() const
    {
      return m_data;
    }

    /// \brief Returns a reference to the data specification.
    /// \return The data specification.
    data::data_specification& data()
    {
      return m_data;
    }

    /// \brief Returns a sequence of action labels.
    /// This sequence contains all action labels occurring in the specification (but it can have more).
    /// \return A sequence of action labels.
    const process::action_label_list& action_labels() const
    {
      return m_action_labels;
    }

    /// \brief Returns a sequence of action labels.
    /// This sequence contains all action labels occurring in the specification (but it can have more).
    /// \return A sequence of action labels.
    process::action_label_list& action_labels()
    {
      return m_action_labels;
    }

    /// \brief Returns the declared free variables of the LPS.
    /// \return The declared free variables of the LPS.
    const std::set<data::variable>& global_variables() const
    {
      return m_global_variables;
    }

    /// \brief Returns the declared free variables of the LPS.
    /// \return The declared free variables of the LPS.
    std::set<data::variable>& global_variables()
    {
      return m_global_variables;
    }

    /// \brief Returns the initial process.
    /// \return The initial process.
    const InitialProcessExpression& initial_process() const
    {
      return m_initial_process;
    }

    /// \brief Returns a reference to the initial process.
    /// \return The initial process.
    InitialProcessExpression& initial_process()
    {
      return m_initial_process;
    }
};

/// \brief Linear process specification.
class specification: public specification_base<linear_process, process_initializer>
{
  protected:
    typedef specification_base<linear_process, process_initializer> super;

  public:
    /// \brief Constructor.
    specification() = default;

    /// \brief Constructor.
    /// \param data A data specification
    /// \param action_labels A sequence of action labels
    /// \param global_variables A set of global variables
    /// \param lps A linear process
    /// \param initial_process A process initializer
    specification(const data::data_specification& data,
                  const process::action_label_list& action_labels,
                  const std::set<data::variable>& global_variables,
                  const linear_process& lps,
                  const process_initializer& initial_process)
      : super(data, action_labels, global_variables, lps, initial_process)
    {
      complete_data_specification(*this);
    }
};

//--- start generated class specification ---//
// prototype declaration
std::string pp(const specification& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const specification& x)
{
  return out << lps::pp(x);
}
//--- end generated class specification ---//

std::string pp_with_summand_numbers(const specification& x);

// template function overloads
std::set<data::sort_expression> find_sort_expressions(const lps::specification& x);
std::set<data::variable> find_all_variables(const lps::specification& x);
std::set<data::variable> find_free_variables(const lps::specification& x);
std::set<data::function_symbol> find_function_symbols(const lps::specification& x);
std::set<core::identifier_string> find_identifiers(const lps::specification& x);
std::set<process::action_label> find_action_labels(const lps::specification& x);

/// \brief Conversion to aterm.
/// \return The specification converted to aterm format.
template <typename LinearProcess, typename InitialProcessExpression>
atermpp::aterm specification_to_aterm(const specification_base<LinearProcess, InitialProcessExpression>& spec)
{
  return atermpp::aterm(core::detail::function_symbol_LinProcSpec(),
           data::detail::data_specification_to_aterm(spec.data()),
           atermpp::aterm(core::detail::function_symbol_ActSpec(), spec.action_labels()),
           atermpp::aterm(core::detail::function_symbol_GlobVarSpec(), data::variable_list(spec.global_variables().begin(),spec.global_variables().end())),
           linear_process_to_aterm(spec.process()),
           spec.initial_process()
         );
}

/// \brief Equality operator
inline
bool operator==(const specification& spec1, const specification& spec2)
{
  return specification_to_aterm(spec1) == specification_to_aterm(spec2);
}

/// \brief Inequality operator
inline
bool operator!=(const specification& spec1, const specification& spec2)
{
  return !(spec1 == spec2);
}

/// \brief Adds all sorts that appear in the process of l to the data specification of l.
/// \param spec A linear process specification
inline
void complete_data_specification(specification& spec)
{
  std::set<data::sort_expression> s = lps::find_sort_expressions(spec);
  spec.data().add_context_sorts(s);
}

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_SPECIFICATION_H
