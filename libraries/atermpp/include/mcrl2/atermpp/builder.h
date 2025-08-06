// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/builder.h
/// \brief add your file description here.

#ifndef MCRL2_ATERMPP_BUILDER_H
#define MCRL2_ATERMPP_BUILDER_H

#include "mcrl2/atermpp/aterm_int.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm_list.cxx"
#else
  import atermpp;
#endif

namespace atermpp {

template <typename Derived>
struct builder
{
  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }

  // Enter object
  template <typename T>
  void enter(const T&)
  {}

  // Leave object
  template <typename T>
  void leave(const T&)
  {}

  template <class T>
  void apply(T& result, const aterm_int& x)
  {
    derived().enter(x);
    derived().leave(x);
    result=x;
  }

  template <class T>
  void apply(T& result, const aterm_list& x)
  {
    derived().enter(x);
    make_term_list(static_cast<aterm_list&>(result), x.begin(), x.end(), [&](aterm& r, const aterm& v) { return derived().apply(r, v); } ) ;
    derived().leave(x);

    // return mcrl2::workaround::return_std_move(result);
  }

  /* template <class T>
  void apply(T& result, const aterm& x)
  {
    derived().enter(x);
    make_term_appl(result, x.function() , x.begin(), x.end(), [&](aterm& r, const aterm& v) { return derived().apply(r, v); } );
    derived().leave(x);

    // return mcrl2::workaround::return_std_move(result);
  } */

  template <class T>
  void apply(T& result, const aterm& x)
  {
    derived().enter(x);
    if (x.type_is_list())
    {
      derived().apply(result, atermpp::down_cast<aterm_list>(x));
    }
    else if (x.type_is_int())
    {
      derived().apply(result, atermpp::down_cast<aterm_int>(x));
    }
    else // It is an application, other than a list or an int. 
    {
      make_term_appl(result, x.function() , x.begin(), x.end(), [&](aterm& r, const aterm& v) { return derived().apply(r, v); } );
    }
    derived().leave(x);
  }
};

} // namespace atermpp

#endif // MCRL2_ATERMPP_BUILDER_H
