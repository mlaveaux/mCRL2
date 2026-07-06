// Author(s): Jeroen  van der Wulp
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/detail/construction_utility.h
/// \brief construction utility

#ifndef MCRL2_CORE_DETAIL_CONSTRUCTION_UTILITY_H
#define MCRL2_CORE_DETAIL_CONSTRUCTION_UTILITY_H

#include "mcrl2/core/identifier_string.h"



namespace mcrl2::core::detail
{

// Component that helps applying the Singleton design pattern
template < typename Derived, typename Expression = atermpp::aterm >
class singleton_expression 
{
  public:
    static const Expression& instance()
    {
      static Expression single_instance = Expression(Derived::initialise());
      return single_instance;
    }

    singleton_expression(const singleton_expression&) = delete;
    singleton_expression& operator=(const singleton_expression &) = delete;
    singleton_expression(singleton_expression &&) = delete;
    singleton_expression & operator=(singleton_expression &&) = delete;

  protected:
    singleton_expression() = default;
};

template < typename Derived >
struct singleton_identifier : public singleton_expression< Derived, core::identifier_string >
  {};

} // namespace mcrl2::core::detail




#endif

