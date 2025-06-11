// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_LPS_STOCHASTIC_SPECIFICATION_H
#define MCRL2_LPS_STOCHASTIC_SPECIFICATION_H

#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/stochastic_linear_process.h"

namespace mcrl2 {

namespace lps {

class stochastic_specification;
void complete_data_specification(stochastic_specification& spec);

// template function overloads
std::set<data::sort_expression> find_sort_expressions(const lps::stochastic_specification& x);
std::set<data::variable> find_all_variables(const lps::stochastic_specification& x);
std::set<data::variable> find_free_variables(const lps::stochastic_specification& x);
std::set<data::function_symbol> find_function_symbols(const lps::stochastic_specification& x);
std::set<core::identifier_string> find_identifiers(const lps::stochastic_specification& x);

// template function overloads
bool check_well_typedness(const stochastic_specification& x);
void normalize_sorts(stochastic_specification& x, const data::sort_specification& sortspec);

/// \brief Linear process specification.
class stochastic_specification: public specification_base<stochastic_linear_process, stochastic_process_initializer>
{
  protected:
    typedef specification_base<stochastic_linear_process, stochastic_process_initializer> super;

  public:
    /// \brief Constructor.
    stochastic_specification()
    { }

    /// \brief Constructor.
    /// \param data A data specification
    /// \param action_labels A sequence of action labels
    /// \param global_variables A set of global variables
    /// \param lps A linear process
    /// \param initial_process A process initializer
    stochastic_specification(const data::data_specification& data,
                             const process::action_label_list& action_labels,
                             const std::set<data::variable>& global_variables,
                             const stochastic_linear_process& lps,
                             const stochastic_process_initializer& initial_process)
      : super(data, action_labels, global_variables, lps, initial_process)
    {
      complete_data_specification(*this);
    }

    /// \brief Constructor. This constructor is explicit as implicit conversions of this kind is a source of bugs.
    explicit stochastic_specification(const specification& other)
      : super(other.data(), 
        other.action_labels(), 
        other.global_variables(), 
        stochastic_linear_process(other.process()),
        stochastic_process_initializer(other.initial_process().expressions(),stochastic_distribution()))
    { } 
};

//--- start generated class stochastic_specification ---//
// prototype declaration
std::string pp(const stochastic_specification& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const stochastic_specification& x)
{
  return out << lps::pp(x);
}
//--- end generated class stochastic_specification ---//

inline
bool operator==(const stochastic_specification& spec1, const stochastic_specification& spec2)
{
  return specification_to_aterm(spec1) == specification_to_aterm(spec2);
}

/// \brief Inequality operator
inline
bool operator!=(const stochastic_specification& spec1, const stochastic_specification& spec2)
{
  return !(spec1 == spec2);
}

std::string pp_with_summand_numbers(const stochastic_specification& x);

/// \brief Adds all sorts that appear in the process of l to the data specification of l.
/// \param spec A linear process specification
inline
void complete_data_specification(stochastic_specification& spec)
{
  std::set<data::sort_expression> s = lps::find_sort_expressions(spec);
  spec.data().add_context_sorts(s);
}

/// \brief Converts a stochastic specification to a specification. Throws an exception if
/// non-empty distributions are encountered.
inline
specification remove_stochastic_operators(const stochastic_specification& spec)
{
  specification result;
  result.data() = spec.data();
  result.action_labels() = spec.action_labels();
  result.global_variables() = spec.global_variables();

  if (spec.initial_process().distribution().is_defined())
  {
    throw mcrl2::runtime_error("The initial state has a non-empty stochastic distribution " + pp(spec.initial_process().distribution()) + ".\n" + 
                               "Transformation of this stochastic lps to a plain lps fails.");
  }
  result.initial_process() = process_initializer(spec.initial_process().expressions());

  linear_process& proc = result.process();
  proc.process_parameters() = spec.process().process_parameters();
  proc.deadlock_summands() = spec.process().deadlock_summands();

  action_summand_vector v;

  for (const stochastic_action_summand& s: spec.process().action_summands())
  {
    if (s.distribution().is_defined())
    {
      throw mcrl2::runtime_error("There is an action summand that has a non-empty stochastic distribution " + pp(s.distribution()) + ".\n" + 
                                 "Transformation of this stochastic lps to a plain lps fails.");;
    }
    v.push_back(s);
  }
  proc.action_summands() = v;
  return result;
}

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_STOCHASTIC_SPECIFICATION_H
