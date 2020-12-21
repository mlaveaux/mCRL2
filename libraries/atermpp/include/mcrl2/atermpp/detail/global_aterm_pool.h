// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_DETAIL_GLOBAL_ATERM_POOL_H_
#define MCRL2_ATERMPP_DETAIL_GLOBAL_ATERM_POOL_H_

#include "mcrl2/atermpp/detail/aterm_pool.h"
#include "mcrl2/atermpp/detail/thread_aterm_pool.h"

namespace atermpp
{
namespace detail
{

/// \brief Storage for a global term pool that is not initialized.
extern typename std::aligned_storage<sizeof(aterm_pool), alignof(aterm_pool)>::type g_aterm_pool_storage;

/// \brief A reference to the global term pool storage
static aterm_pool& g_aterm_pool_instance = *static_cast<aterm_pool*>(static_cast<void*>(&g_aterm_pool_storage));

inline thread_aterm_pool& g_thread_term_pool()
{
  static thread_local thread_aterm_pool per_thread_pool(g_aterm_pool_instance);
  return per_thread_pool;
}

/// \brief obtain a reference to the global aterm pool.
/// \details provides lazy initialization which should be used when instantiating
///          global terms and function symbols.
template<bool lazy = false>
inline aterm_pool& g_term_pool()
{
  if (lazy)
  {
    static bool initialized = false;
    if (!initialized)
    {
      new (&g_aterm_pool_instance) aterm_pool();
      initialized = true;
    }
  }

  return g_aterm_pool_instance;
}

} // namespace detail
} // namespace atermpp

#endif // MCRL2_ATERMPP_DETAIL_GLOBAL_ATERM_POOL_H_
