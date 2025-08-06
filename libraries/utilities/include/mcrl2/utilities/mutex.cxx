// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_MUTEX_H_
#define MCRL2_UTILITIES_MUTEX_H_
MCRL2_MODULE;

#include <mutex>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/configuration.cxx"
#else
  export module utilities:mutex;

  import :configuration;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities
{

/// \brief This is simply an exclusive lock based on the standard library with the ability to perform no locks when 
///        thread safety is not desired.
class mutex 
{
public:
  void lock() {
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe) {
        m_mutex.lock();
    }
  }

  void unlock()
  {        
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe) {
        m_mutex.unlock();
    }
  }

private:
    std::mutex m_mutex;
};

} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_MUTEX_H_