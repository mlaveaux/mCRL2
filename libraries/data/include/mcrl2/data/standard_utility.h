// Author(s): Jeroen van der Wulp
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/standard_utility.h
/// \brief Provides utilities for working with data expressions of standard sorts

#ifndef MCRL2_DATA_STANDARD_UTILITY_H
#define MCRL2_DATA_STANDARD_UTILITY_H

#include <queue>

#include "mcrl2/data/real.h"
#include "mcrl2/utilities/detail/join.h"



namespace mcrl2::data
{

namespace sort_bool
{
/// \brief Constructs expression of type Bool from an integral type
/// \param b A Boolean
inline data_expression bool_(bool b)
{
  return (b) ? sort_bool::true_() : sort_bool::false_();
}

/// \brief Determines whether b is a Boolean constant
/// \param b A data expression
inline bool is_boolean_constant(data_expression const& b)
{
  return sort_bool::is_true_function_symbol(b) ||
         sort_bool::is_false_function_symbol(b);
}
} // namespace sort_bool

/// \brief Returns true iff the expression represents a standard sort.
/// \param[in] s a sort expression.
inline
bool
is_system_defined(const sort_expression& s)
{
  return sort_bool::is_bool(s) || sort_real::is_real(s)
         || sort_int::is_int(s) || sort_nat::is_nat(s) || sort_pos::is_pos(s)
         || is_container_sort(s) || is_structured_sort(s);
}

/** \brief A collection of utilities for lazy expression construction
 *
 * The basic idea is to keep expressions that result from application of
 * any of the container operations by applying the usual rules of logic.
 *
 * For example and(true, x) as in @c and applied to @c true and @c x yields x.
 **/
namespace lazy
{
/// \brief Returns an expression equivalent to not p
/// \param p A data expression
/// \return The value <tt>!p</tt>
inline data_expression not_(data_expression const& p)
{
  if (p == sort_bool::true_())
  {
    return sort_bool::false_();
  }
  else if (p == sort_bool::false_())
  {
    return sort_bool::true_();
  }

  return sort_bool::not_(p);
}

/// \brief Returns an expression equivalent to p and q
/// \param p A data expression
/// \param q A data expression
/// \return The value <tt>p && q</tt>
inline data_expression or_(data_expression const& p, data_expression const& q)
{
  if ((p == sort_bool::true_()) || (q == sort_bool::true_()))
  {
    return sort_bool::true_();
  }
  else if ((p == q) || (p == sort_bool::false_()))
  {
    return q;
  }
  else if (q == sort_bool::false_())
  {
    return p;
  }

  return sort_bool::or_(p, q);
}

/// \brief Returns an expression equivalent to p or q
/// \param p A data expression
/// \param q A data expression
/// \return The value p || q
inline data_expression and_(data_expression const& p, data_expression const& q)
{
  if ((p == sort_bool::false_()) || (q == sort_bool::false_()))
  {
    return sort_bool::false_();
  }
  else if ((p == q) || (p == sort_bool::true_()))
  {
    return q;
  }
  else if (q == sort_bool::true_())
  {
    return p;
  }

  return sort_bool::and_(p, q);
}

/// \brief Returns an expression equivalent to p implies q
/// \param p A data expression
/// \param q A data expression
/// \return The value p || q
inline data_expression implies(data_expression const& p, data_expression const& q)
{
  if ((p == sort_bool::false_()) || (q == sort_bool::true_()) || (p == q))
  {
    return sort_bool::true_();
  }
  else if (p == sort_bool::true_())
  {
    return q;
  }
  else if (q == sort_bool::false_())
  {
    return sort_bool::not_(p);
  }

  return sort_bool::implies(p, q);
}

/// \brief Returns an expression equivalent to p == q
/// \param p A data expression
/// \param q A data expression
/// \return The value p == q
inline data_expression equal_to(data_expression const& p, data_expression const& q)
{
  if (p == q)
  {
    return sort_bool::true_();
  }

  return data::equal_to(p, q);
}

/// \brief Returns an expression equivalent to p == q
/// \param p A data expression
/// \param q A data expression
/// \return The value ! p == q
inline data_expression not_equal_to(data_expression const& p, data_expression const& q)
{
  if (p == q)
  {
    return sort_bool::false_();
  }

  return data::not_equal_to(p, q);
}

/// @brief Returns an expression equivalent to if(cond,then,else_)
/// @return The value if(cond,then,else_)
inline data_expression if_(const data_expression& cond, const data_expression& then, const data_expression& else_)
{
  if (cond == sort_bool::true_() || then == else_)
  {
    return then;
  }
  else if (cond == sort_bool::false_())
  {
    return else_;
  }

  return data::if_(cond, then, else_);
}

/// \brief Returns or applied to the sequence of data expressions [first, last)
/// \param first Start of a sequence of data expressions
/// \param last End of a sequence of data expressions
/// \return Or applied to the sequence of data expressions [first, last)
template < typename ForwardTraversalIterator >
inline data_expression join_or(ForwardTraversalIterator first, ForwardTraversalIterator last)
{
  return utilities::detail::join(first, last, lazy::or_, static_cast< data_expression const& >(sort_bool::false_()));
}

/// \brief Returns and applied to the sequence of data expressions [first, last)
/// \param first Start of a sequence of data expressions
/// \param last End of a sequence of data expressions
/// \return And applied to the sequence of data expressions [first, last)
template < typename ForwardTraversalIterator >
inline data_expression join_and(ForwardTraversalIterator first, ForwardTraversalIterator last)
{
  return utilities::detail::join(first, last, lazy::and_, static_cast< data_expression const& >(sort_bool::true_()));
}
} // namespace lazy

} // namespace mcrl2::data



#endif

