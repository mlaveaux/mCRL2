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
MCRL2_MODULE;

#include <cstddef>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/aterm_pool.cxx"
  #include "mcrl2/atermpp/detail/thread_aterm_pool.cxx"
#else
  export module atermpp:detail.global_aterm_pool;

  import :detail.aterm_pool;
  import :detail.thread_aterm_pool;
#endif

namespace atermpp::detail
{

/// \brief Storage for a global term pool that is not initialized.
alignas(aterm_pool)
extern std::byte g_aterm_pool_storage[sizeof(aterm_pool)];

/// \brief A reference to the global term pool storage
static aterm_pool& g_aterm_pool_instance = *reinterpret_cast<aterm_pool*>(&g_aterm_pool_storage);

/// \brief obtain a reference to the global aterm pool.
/// \details provides lazy initialization which should be used when instantiating
///          global terms and function symbols.
MCRL2_MODULE_EXPORT
template<bool lazy = false>
inline aterm_pool& g_term_pool()
{
  if constexpr (lazy)
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

} // namespace atermpp::detail

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/aterm_implementation.cxx"
#endif 

#endif // MCRL2_ATERMPP_DETAIL_GLOBAL_ATERM_POOL_H_
