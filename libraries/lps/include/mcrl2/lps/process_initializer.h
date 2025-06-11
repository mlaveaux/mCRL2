// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/process_initializer.h
/// \brief The class process_initializer.

#ifndef MCRL2_LPS_PROCESS_INITIALIZER_H
#define MCRL2_LPS_PROCESS_INITIALIZER_H

#include "mcrl2/data/print.h"
#include "mcrl2/data/replace.h"
#include "mcrl2/data/substitutions/assignment_sequence_substitution.h"
#include "mcrl2/lps/stochastic_distribution.h"
#include "mcrl2/process/action_label.h"

namespace mcrl2
{

namespace lps
{

/// \brief A process initializer
class process_initializer: public atermpp::aterm
{
  public:
    /// \brief Default constructor.
    process_initializer()
      : atermpp::aterm(core::detail::default_values::LinearProcessInit)
    {}

    /// \brief Constructor.
    /// \param term A term.
    /// \param check_distribution Check whether the initial state is plain or a state distribution.
    explicit process_initializer(const atermpp::aterm& term, bool check_distribution = true)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_LinearProcessInit(*this));
      const lps::stochastic_distribution& dist = atermpp::down_cast<lps::stochastic_distribution>(term[1]);
      if (check_distribution && dist.is_defined())
      {
        throw mcrl2::runtime_error("initial state with non-empty stochastic distribution encountered");
      }
    }

    /// \brief Constructor.
    explicit process_initializer(const data::data_expression_list& expressions)
      : atermpp::aterm(core::detail::function_symbol_LinearProcessInit(), expressions, stochastic_distribution())
    {}

    /// Move semantics
    process_initializer(const process_initializer&) noexcept = default;
    process_initializer(process_initializer&&) noexcept = default;
    process_initializer& operator=(const process_initializer&) noexcept = default;
    process_initializer& operator=(process_initializer&&) noexcept = default;

    data::data_expression_list expressions() const
    {
      return atermpp::down_cast<data::data_expression_list>((*this)[0]);
    }
};


template <class EXPRESSION_LIST>
inline void make_process_initializer(atermpp::aterm& t, EXPRESSION_LIST args)
{
  make_term_appl(t, core::detail::function_symbol_LinearProcessInit(), args, stochastic_distribution());
}


//--- start generated class process_initializer ---//
/// \\brief list of process_initializers
typedef atermpp::term_list<process_initializer> process_initializer_list;

/// \\brief vector of process_initializers
typedef std::vector<process_initializer>    process_initializer_vector;

/// \\brief Test for a process_initializer expression
/// \\param x A term
/// \\return True if \\a x is a process_initializer expression
inline
bool is_process_initializer(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::LinearProcessInit;
}

// prototype declaration
std::string pp(const process_initializer& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const process_initializer& x)
{
  return out << lps::pp(x);
}

/// \\brief swap overload
inline void swap(process_initializer& t1, process_initializer& t2)
{
  t1.swap(t2);
}
//--- end generated class process_initializer ---//

// template function overloads
std::set<data::variable> find_free_variables(const lps::process_initializer& x);
std::set<process::action_label> find_action_labels(const lps::process_initializer& x);

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_PROCESS_INITIALIZER_H
