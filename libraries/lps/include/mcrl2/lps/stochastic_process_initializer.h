// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/stochastic_process_initializer.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_STOCHASTIC_PROCESS_INITIALIZER_H
#define MCRL2_LPS_STOCHASTIC_PROCESS_INITIALIZER_H

#include "mcrl2/lps/process_initializer.h"

namespace mcrl2 {

namespace lps {

/// \brief A stochastic process initializer
class stochastic_process_initializer: public process_initializer
{
  typedef process_initializer super;

  public:
    /// \brief Default constructor.
    stochastic_process_initializer()
      : super(core::detail::default_values::LinearProcessInit)
    {}

    /// \brief Constructor.
    /// \param term A term
    explicit stochastic_process_initializer(const atermpp::aterm& term)
      : super(term, false)
    {
      assert(core::detail::check_term_LinearProcessInit(*this));
    }

    /// \brief Constructor.
    stochastic_process_initializer(const data::data_expression_list& expressions, const stochastic_distribution& distribution)
      : super(atermpp::aterm(core::detail::function_symbol_LinearProcessInit(), expressions, distribution), false)
    {}

    const stochastic_distribution& distribution() const
    {
      return atermpp::down_cast<stochastic_distribution>((*this)[1]);
    }
};

template <class... ARGUMENTS>
inline void make_stochastic_process_initializer(atermpp::aterm& t, ARGUMENTS... args)
{
  make_term_appl(t, core::detail::function_symbol_LinearProcessInit(), args...);
}


//--- start generated class stochastic_process_initializer ---//
/// \\brief list of stochastic_process_initializers
typedef atermpp::term_list<stochastic_process_initializer> stochastic_process_initializer_list;

/// \\brief vector of stochastic_process_initializers
typedef std::vector<stochastic_process_initializer>    stochastic_process_initializer_vector;

/// \\brief Test for a stochastic_process_initializer expression
/// \\param x A term
/// \\return True if \\a x is a stochastic_process_initializer expression
inline
bool is_stochastic_process_initializer(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::LinearProcessInit;
}

// prototype declaration
std::string pp(const stochastic_process_initializer& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const stochastic_process_initializer& x)
{
  return out << lps::pp(x);
}

/// \\brief swap overload
inline void swap(stochastic_process_initializer& t1, stochastic_process_initializer& t2)
{
  t1.swap(t2);
}
//--- end generated class stochastic_process_initializer ---//

// template function overloads
std::set<data::variable> find_free_variables(const lps::stochastic_process_initializer& x);

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_STOCHASTIC_PROCESS_INITIALIZER_H
