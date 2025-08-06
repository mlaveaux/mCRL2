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

#include "mcrl2/atermpp/detail/aterm_container.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/global_aterm_pool.cxx"
#else
  import atermpp;
#endif

namespace atermpp
{

namespace detail
{
template <typename T>
const reference_aterm<T, typename std::enable_if_t<std::is_base_of_v<aterm_core, T>>>&
reference_aterm<T, typename std::enable_if_t<std::is_base_of_v<aterm_core, T>>>::operator=(
    const unprotected_aterm_core& other) noexcept
{
  mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
  m_term = address(other);
  return *this;
  }

  template <typename T>
  const reference_aterm<T, typename std::enable_if_t<std::is_base_of_v<aterm_core, T>>>&
  reference_aterm<T, typename std::enable_if_t<std::is_base_of_v<aterm_core, T>>>::operator=(
      unprotected_aterm_core&& other) noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_term = address(other);
    return *this;
  }

  template<typename T, typename Allocator>
  void aterm_allocator<T,Allocator>::deallocate(T* p, size_type n)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_allocator.deallocate(p, n);
  }

  inline
  aterm_container::aterm_container(std::function<void(term_mark_stack&)> mark_func, std::function<std::size_t()> size_func)
    : mark_func(mark_func),
      size_func(size_func)
  {
    g_thread_term_pool().register_container(this);
  }

  inline
  aterm_container::~aterm_container()
  {
    g_thread_term_pool().deregister_container(this);
  }
}

inline aterm_core::aterm_core() noexcept
{
  detail::g_thread_term_pool().register_variable(this);
}

inline aterm_core::~aterm_core() noexcept
{
  detail::g_thread_term_pool().deregister_variable(this);
}

inline aterm_core::aterm_core(const detail::_aterm *t) noexcept
{
  detail::g_thread_term_pool().register_variable(this);
  m_term = t;
}

inline aterm_core::aterm_core(const aterm_core& other) noexcept
 : unprotected_aterm_core(other.m_term)
{
  detail::g_thread_term_pool().register_variable(this);
}

inline aterm_core::aterm_core(aterm_core&& other) noexcept
 : unprotected_aterm_core(other.m_term)
{
  detail::g_thread_term_pool().register_variable(this);
}

inline aterm_core& aterm_core::operator=(const aterm_core& other) noexcept
{
  mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
  m_term = other.m_term;
  return *this;
}

inline aterm_core& aterm_core::assign(const aterm_core& other,
  detail::thread_aterm_pool& pool) noexcept
{
  mcrl2::utilities::shared_guard guard = pool.lock_shared();
  m_term = other.m_term;
  return *this;
}


template <bool CHECK_BUSY_FLAG /* =true*/>
inline aterm_core& aterm_core::unprotected_assign(const aterm_core& other) noexcept
{
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe && CHECK_BUSY_FLAG)
  {
    assert(detail::g_thread_term_pool().is_shared_locked());
  }
  m_term = other.m_term;
  return *this;
}


inline aterm_core& aterm_core::operator=(aterm_core&& other) noexcept
{
  mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
  // Using hash set protection it is cheaper just to move the value to the new term.
  m_term = other.m_term;    
  return *this;
}


} // namespace atermpp

#endif // MCRL2_ATERMPP_TERM_IMPLEMENTATION_H
