// Author(s): Jeroen Keiren, Jeroen van der Wulp
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/structured_sort_constructor_argument.h
/// \brief The class structured_sort_constructor_arguments.

#ifndef MCRL2_DATA_STRUCTURED_SORT_CONSTRUCTOR_ARGUMENT_H
#define MCRL2_DATA_STRUCTURED_SORT_CONSTRUCTOR_ARGUMENT_H

#include "mcrl2/data/function_symbol.h"

namespace mcrl2
{

namespace data
{

//--- start generated class structured_sort_constructor_argument ---//
/// \\brief An argument of a constructor of a structured sort
class structured_sort_constructor_argument: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    structured_sort_constructor_argument()
      : atermpp::aterm(core::detail::default_values::StructProj)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit structured_sort_constructor_argument(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_StructProj(*this));
    }

    /// \\brief Constructor Z12.
    structured_sort_constructor_argument(const core::identifier_string& name, const sort_expression& sort)
      : atermpp::aterm(core::detail::function_symbol_StructProj(), name, sort)
    {}

    /// \\brief Constructor Z1.
    structured_sort_constructor_argument(const std::string& name, const sort_expression& sort)
      : atermpp::aterm(core::detail::function_symbol_StructProj(), core::identifier_string(name), sort)
    {}

    /// Move semantics
    structured_sort_constructor_argument(const structured_sort_constructor_argument&) noexcept = default;
    structured_sort_constructor_argument(structured_sort_constructor_argument&&) noexcept = default;
    structured_sort_constructor_argument& operator=(const structured_sort_constructor_argument&) noexcept = default;
    structured_sort_constructor_argument& operator=(structured_sort_constructor_argument&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const sort_expression& sort() const
    {
      return atermpp::down_cast<sort_expression>((*this)[1]);
    }
//--- start user section structured_sort_constructor_argument ---//
    /// \brief Constructor
    ///
    /// \param[in] sort The sort of the argument.
    structured_sort_constructor_argument(const sort_expression& sort)
      : atermpp::aterm(core::detail::function_symbol_StructProj(), core::empty_identifier_string(), sort)
    {}

    /// \brief Constructor.
    ///
    /// \overload to work around problem that MSVC reinterprets char* or char[] as core::identifier_string
    template < std::size_t S >
    structured_sort_constructor_argument(const char(&name)[S], const sort_expression& sort)
      : atermpp::aterm(core::detail::function_symbol_StructProj(), core::identifier_string(name), sort)
    {}
//--- end user section structured_sort_constructor_argument ---//
};

/// \\brief Make_structured_sort_constructor_argument constructs a new term into a given address.
/// \\ \param t The reference into which the new structured_sort_constructor_argument is constructed. 
template <class... ARGUMENTS>
inline void make_structured_sort_constructor_argument(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StructProj(), args...);
}

/// \\brief list of structured_sort_constructor_arguments
typedef atermpp::term_list<structured_sort_constructor_argument> structured_sort_constructor_argument_list;

/// \\brief vector of structured_sort_constructor_arguments
typedef std::vector<structured_sort_constructor_argument>    structured_sort_constructor_argument_vector;

/// \\brief Test for a structured_sort_constructor_argument expression
/// \\param x A term
/// \\return True if \\a x is a structured_sort_constructor_argument expression
inline
bool is_structured_sort_constructor_argument(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StructProj;
}

// prototype declaration
std::string pp(const structured_sort_constructor_argument& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const structured_sort_constructor_argument& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(structured_sort_constructor_argument& t1, structured_sort_constructor_argument& t2)
{
  t1.swap(t2);
}
//--- end generated class structured_sort_constructor_argument ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_STRUCTURED_SORT_CONSTRUCTOR_ARGUMENT_H

