// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_COPYABLE_MUTEX_H
#define MCRL2_UTILITIES_COPYABLE_MUTEX_H

#include <mutex>

namespace mcrl2::utilities
{

/// \brief The default mutex is not copyable or movable, because copying classes
///        concurrently requires additional attention. However, this also means
///        that for all classes that use mutexes need explicit copy/move constructors
///        and assignments. For classes that are only copied during initialisation
///        in the main thread we can use the copyable_mutex.
class copyable_mutex : public std::mutex
{  
public:
  copyable_mutex() {};
  copyable_mutex(const copyable_mutex&) {};
  copyable_mutex& operator=(const copyable_mutex&) { return *this; };

  copyable_mutex(copyable_mutex&&) {};
  copyable_mutex& operator=(copyable_mutex&&) { return *this; };
};

}

#endif // MCRL2_UTILITIES_COPYABLE_MUTEX_H
