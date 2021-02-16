// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_IMPLEMENTATION_H
#define MCRL2_ATERMPP_ATERM_IMPLEMENTATION_H
#pragma once

#include "mcrl2/atermpp/detail/global_aterm_pool.h"

namespace atermpp
{

inline aterm::aterm() noexcept
{
#ifndef MCRL2_ATERMPP_REFERENCE_COUNTED
  detail::g_thread_term_pool().register_variable(this);
#endif
}

inline aterm::~aterm() noexcept
{
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  decrement_reference_count();
#else
  detail::g_thread_term_pool().remove_variable(this);
#endif
}

inline aterm::aterm(const detail::_aterm *t) noexcept
{
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  t->increment_reference_count();
#else
  detail::g_thread_term_pool().register_variable(this);
#endif
  m_term = t;
}

inline aterm::aterm(const aterm& other) noexcept
 : unprotected_aterm(other.m_term)
{
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  increment_reference_count();
#else
  detail::g_thread_term_pool().register_variable(this);
#endif
}

inline aterm::aterm(aterm&& other) noexcept
 : unprotected_aterm(other.m_term)
{
#ifndef MCRL2_ATERMPP_REFERENCE_COUNTED
  detail::g_thread_term_pool().register_variable(this);
#endif
  other.m_term=nullptr;
}

inline aterm_container::aterm_container()
{
#ifndef MCRL2_ATERMPP_REFERENCE_COUNTED
  detail::g_thread_term_pool().register_container(this);
#endif
}

inline aterm_container::~aterm_container()
{
#ifndef MCRL2_ATERMPP_REFERENCE_COUNTED
  detail::g_thread_term_pool().remove_container(this);
#endif
}

inline void add_creation_hook(const function_symbol& function, term_callback callback)
{
  detail::g_term_pool().add_creation_hook(function, callback);
}

inline void add_deletion_hook(const function_symbol& function, term_callback callback)
{
  detail::g_term_pool().add_deletion_hook(function, callback);
}

} // namespace atermpp

#endif // MCRL2_ATERMPP_TERM_IMPLEMENTATION_H
