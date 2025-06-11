// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/sort_expression.h
/// \brief The class sort_expression.

#ifndef MCRL2_DATA_SORT_EXPRESSION_H
#define MCRL2_DATA_SORT_EXPRESSION_H

#include "mcrl2/core/detail/default_values.h"
#include "mcrl2/core/detail/soundness_checks.h"

namespace mcrl2
{

namespace data
{

/// \brief Returns true if the term t is a basic sort
inline bool is_basic_sort(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortId;
}

/// \brief Returns true if the term t is a function sort
inline bool is_function_sort(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortArrow;
}

/// \brief Returns true if the term t is a container sort
inline bool is_container_sort(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortCons;
}

/// \brief Returns true if the term t is a structured sort
inline bool is_structured_sort(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::SortStruct;
}

/// \brief Returns true if the term t is the unknown sort
inline bool is_untyped_sort(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::UntypedSortUnknown;
}

/// \brief Returns true if the term t is an expression for multiple possible sorts
inline bool is_untyped_possible_sorts(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::UntypedSortsPossible;
}

//--- start generated class sort_expression ---//
/// \\brief A sort expression
class sort_expression: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    sort_expression()
      : atermpp::aterm(core::detail::default_values::SortExpr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit sort_expression(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_SortExpr(*this));
    }

    /// Move semantics
    sort_expression(const sort_expression&) noexcept = default;
    sort_expression(sort_expression&&) noexcept = default;
    sort_expression& operator=(const sort_expression&) noexcept = default;
    sort_expression& operator=(sort_expression&&) noexcept = default;
//--- start user section sort_expression ---//
    /// \brief Returns the target sort of this expression.
    /// \return For a function sort D->E, return the target sort of E. Otherwise return this sort.
    inline
    const sort_expression& target_sort() const
    {
      if (is_function_sort(*this))
      {
        return atermpp::down_cast<const sort_expression>((*this)[1]);
      }
      else
      {
        return *this;
      }
    }
//--- end user section sort_expression ---//
};

/// \\brief list of sort_expressions
typedef atermpp::term_list<sort_expression> sort_expression_list;

/// \\brief vector of sort_expressions
typedef std::vector<sort_expression>    sort_expression_vector;

// prototype declaration
std::string pp(const sort_expression& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const sort_expression& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(sort_expression& t1, sort_expression& t2)
{
  t1.swap(t2);
}
//--- end generated class sort_expression ---//

/// \brief Test for a sort_expression expression
/// \param x A term
/// \return True if it is a sort_expression expression
inline
bool is_sort_expression(const atermpp::aterm& x)
{
  return is_basic_sort(x)             ||
         is_function_sort(x)          ||
         is_container_sort(x)         ||
         is_structured_sort(x)        ||
         is_untyped_sort(x)           ||
         is_untyped_possible_sorts(x);
}

// template function overloads
std::string pp(const sort_expression_list& x, bool precedence_aware = true);
std::string pp(const sort_expression_vector& x, bool precedence_aware = true);
std::set<data::sort_expression> find_sort_expressions(const data::sort_expression& x);

} // namespace data

} // namespace mcrl2

namespace std
{

template <>
struct hash<mcrl2::data::sort_expression>
{
  std::size_t operator()(const mcrl2::data::sort_expression& x) const
  {
    return hash<atermpp::aterm>()(x);
  }
};

} // namespace std

#endif // MCRL2_DATA_SORT_EXPRESSION_H

