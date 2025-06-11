// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/abstraction.h
/// \brief The class abstraction.

#ifndef MCRL2_DATA_ABSTRACTION_H
#define MCRL2_DATA_ABSTRACTION_H

#include "mcrl2/data/binder_type.h"
#include "mcrl2/data/data_expression.h"

namespace mcrl2
{

namespace data
{

/// \brief An abstraction expression.
class abstraction: public data_expression
{
  public:
    /// \brief Default constructor.
    abstraction()
      : data_expression(core::detail::default_values::Binder)
    {}

    /// \brief Constructor.
    /// \param term A term
    explicit abstraction(const atermpp::aterm& term)
      : data_expression(term)
    {
      assert(core::detail::check_term_Binder(*this));
    }

    /// \brief Constructor.
    abstraction(const binder_type& binding_operator, const variable_list& variables, const data_expression& body)
      : data_expression(atermpp::aterm(core::detail::function_symbol_Binder(), binding_operator, variables, body))
    {}

    /// \brief Constructor.
    template <typename Container>
    abstraction(const binder_type& binding_operator, const Container& variables, const data_expression& body, typename atermpp::enable_if_container<Container, variable>::type* = nullptr)
      : data_expression(atermpp::aterm(core::detail::function_symbol_Binder(), binding_operator, variable_list(variables.begin(), variables.end()), body))
    {}

    /// Move semantics
    abstraction(const abstraction&) noexcept = default;
    abstraction(abstraction&&) noexcept = default;
    abstraction& operator=(const abstraction&) noexcept = default;
    abstraction& operator=(abstraction&&) noexcept = default;

    const binder_type& binding_operator() const
    {
      return atermpp::down_cast<const binder_type>((*this)[0]);
    }

    const variable_list& variables() const
    {
      return atermpp::down_cast<const variable_list>((*this)[1]);
    }

    const data_expression& body() const
    {
      return atermpp::down_cast<const data_expression>((*this)[2]);
    }
};

template <class... ARGUMENTS>
void make_abstraction(atermpp::aterm& result, ARGUMENTS... arguments)
{
  make_term_appl(result, core::detail::function_symbol_Binder(), arguments...);
}

/* inline void make_abstraction(atermpp::aterm& result, const binder_type& binding_operator, const variable_list& variables, const data_expression& body)
{
  make_term_appl(result, core::detail::function_symbol_Binder(), binding_operator, variables, body);
} */

//--- start generated class abstraction ---//
// prototype declaration
std::string pp(const abstraction& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const abstraction& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(abstraction& t1, abstraction& t2)
{
  t1.swap(t2);
}
//--- end generated class abstraction ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_ABSTRACTION_H

