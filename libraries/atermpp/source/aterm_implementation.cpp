// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/atermpp/aterm_io.h"

using namespace atermpp;
using namespace atermpp::detail;

/// \brief Check for reasonably sized aterm (32 bits, 4 bytes)
///        This check might break on perfectly valid architectures
///        that have char == 2 bytes, and sizeof(header_type) == 2
static_assert(sizeof(std::size_t) == sizeof(_aterm*), "The size of an aterm pointer is not equal to the size of type std::size_t. Cannot compile the MCRL2 toolset for this platform.");
static_assert(sizeof(std::size_t) >= 4,"The size of std::size_t should at least be four bytes. Cannot compile the toolset for this platform.");

aterm::aterm() noexcept
{
#ifndef MCRL2_ATERMPP_REFERENCE_COUNTED
  g_thread_term_pool().register_variable(this);
#endif
}

aterm::~aterm() noexcept
{
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  decrement_reference_count();
#else
  g_thread_term_pool().remove_variable(this);
#endif
}

aterm::aterm(const detail::_aterm *t) noexcept
{
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  t->increment_reference_count();
#else
  g_thread_term_pool().register_variable(this);
#endif
  m_term = t;
}

aterm::aterm(const aterm& other) noexcept
 : unprotected_aterm(other.m_term)
{
#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  increment_reference_count();
#else
  g_thread_term_pool().register_variable(this);
#endif
}

aterm::aterm(aterm&& other) noexcept
 : unprotected_aterm(other.m_term)
{
#ifndef MCRL2_ATERMPP_REFERENCE_COUNTED
  g_thread_term_pool().register_variable(this);
#endif
  other.m_term=nullptr;
}

void atermpp::add_creation_hook(const function_symbol& function, term_callback callback)
{
  g_term_pool().add_creation_hook(function, callback);
}

void atermpp::add_deletion_hook(const function_symbol& function, term_callback callback)
{  
  g_term_pool().add_deletion_hook(function, callback);
}

aterm_stream::~aterm_stream() {}

aterm_istream::~aterm_istream() {}

aterm_ostream::~aterm_ostream() {}

typename std::aligned_storage<sizeof(aterm_pool), alignof(aterm_pool)>::type atermpp::detail::g_aterm_pool_storage = {};
