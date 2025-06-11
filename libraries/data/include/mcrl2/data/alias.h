// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/alias.h
/// \brief The class alias.

#ifndef MCRL2_DATA_ALIAS_H
#define MCRL2_DATA_ALIAS_H

#include "mcrl2/data/basic_sort.h"

namespace mcrl2
{

namespace data
{

//--- start generated class alias ---//
/// \\brief A sort alias
class alias: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    alias()
      : atermpp::aterm(core::detail::default_values::SortRef)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit alias(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_SortRef(*this));
    }

    /// \\brief Constructor Z12.
    alias(const basic_sort& name, const sort_expression& reference)
      : atermpp::aterm(core::detail::function_symbol_SortRef(), name, reference)
    {}

    /// Move semantics
    alias(const alias&) noexcept = default;
    alias(alias&&) noexcept = default;
    alias& operator=(const alias&) noexcept = default;
    alias& operator=(alias&&) noexcept = default;

    const basic_sort& name() const
    {
      return atermpp::down_cast<basic_sort>((*this)[0]);
    }

    const sort_expression& reference() const
    {
      return atermpp::down_cast<sort_expression>((*this)[1]);
    }
};

/// \\brief Make_alias constructs a new term into a given address.
/// \\ \param t The reference into which the new alias is constructed. 
template <class... ARGUMENTS>
inline void make_alias(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_SortRef(), args...);
}

/// \\brief list of aliass
typedef atermpp::term_list<alias> alias_list;

/// \\brief vector of aliass
typedef std::vector<alias>    alias_vector;

/// \\brief Test for a alias expression
/// \\param x A term
/// \\return True if \\a x is a alias expression
inline
bool is_alias(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortRef;
}

// prototype declaration
std::string pp(const alias& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const alias& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(alias& t1, alias& t2)
{
  t1.swap(t2);
}
//--- end generated class alias ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_ALIAS_H

