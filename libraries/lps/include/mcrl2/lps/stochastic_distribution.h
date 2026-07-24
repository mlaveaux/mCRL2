// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/stochastic_distribution.h
/// \brief The class stochastic_distribution.

#ifndef MCRL2_LPS_STOCHASTIC_DISTRIBUTION_H
#define MCRL2_LPS_STOCHASTIC_DISTRIBUTION_H

#include "mcrl2/data/variable.h"

namespace mcrl2::lps {

//--- start generated class stochastic_distribution ---//
/// \\brief A stochastic distribution
class stochastic_distribution: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    stochastic_distribution()
      : atermpp::aterm(core::detail::default_values::Distribution)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit stochastic_distribution(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_Distribution(*this));
    }

    /// \\brief Constructor Z12.
    stochastic_distribution(const data::variable_list& variables, const data::data_expression& distribution)
      : atermpp::aterm(core::detail::function_symbol_Distribution(), variables, distribution)
    {}

    /// Move semantics
    stochastic_distribution(const stochastic_distribution&) noexcept = default;
    stochastic_distribution(stochastic_distribution&&) noexcept = default;
    stochastic_distribution& operator=(const stochastic_distribution&) noexcept = default;
    stochastic_distribution& operator=(stochastic_distribution&&) noexcept = default;

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[0]);
    }

    const data::data_expression& distribution() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[1]);
    }
//--- start user section stochastic_distribution ---//
    /// \brief Returns true if the distribution is defined, i.e. it contains a valid distribution.
    /// This is encoded by a non-empty list of variables.
    bool is_defined() const
    {
      return !variables().empty();
    }
//--- end user section stochastic_distribution ---//
};

/// \\brief The function make_stochastic_distribution constructs a new term into a given address.
/// \\ \param t The reference into which the new stochastic_distribution is constructed.
template <class... ARGUMENTS>
inline void make_stochastic_distribution(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_Distribution(), args...);
}

/// \\brief list of stochastic_distributions
using stochastic_distribution_list = atermpp::term_list<stochastic_distribution>;

/// \\brief vector of stochastic_distributions
using stochastic_distribution_vector = std::vector<stochastic_distribution>;

/// \\brief Test for a stochastic_distribution expression
/// \\param x A term
/// \\return True if \\a x is a stochastic_distribution expression
inline
bool is_stochastic_distribution(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::Distribution;
}

// prototype declaration
std::string pp(const stochastic_distribution& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const stochastic_distribution& x)
{
  return out << lps::pp(x);
}

/// \\brief swap overload
inline void swap(stochastic_distribution& t1, stochastic_distribution& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated class stochastic_distribution ---//

} // namespace mcrl2::lps



#endif // MCRL2_LPS_STOCHASTIC_DISTRIBUTION_H
