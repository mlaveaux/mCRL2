// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#define BOOST_TEST_MODULE shared_mutex_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/utilities/configuration.h"
#include "mcrl2/utilities/shared_mutex.h"

#include <mutex>
#include <shared_mutex>
#include <type_traits>

using namespace mcrl2::utilities;

// The mutex must be usable with the standard library lock guards, so it must satisfy
// the corresponding named requirements.
static_assert(detail::IsLockable<shared_mutex>);
static_assert(detail::IsSharedLockable<shared_mutex>);

// The guards must be the standard library guards.
static_assert(std::is_same_v<lock_guard, std::unique_lock<shared_mutex>>);
static_assert(std::is_same_v<shared_guard, std::shared_lock<shared_mutex>>);

BOOST_AUTO_TEST_CASE(exclusive_lock_unlock)
{
  shared_mutex mutex;

  // The standard library interface: lock() returns void and unlock() releases.
  mutex.lock();
  mutex.unlock();

  // try_lock returns true on success, as std::shared_mutex::try_lock does.
  BOOST_CHECK(mutex.try_lock());
  mutex.unlock();
}

BOOST_AUTO_TEST_CASE(shared_lock_unlock_recursive)
{
  shared_mutex mutex;

  BOOST_CHECK(!mutex.is_shared_locked());

  // Recursive shared locking is an explicitly documented extension.
  mutex.lock_shared();
  BOOST_CHECK(mutex.is_shared_locked());

  mutex.lock_shared();
  BOOST_CHECK(mutex.is_shared_locked());

  mutex.unlock_shared();
  BOOST_CHECK(mutex.is_shared_locked());

  mutex.unlock_shared();
  BOOST_CHECK(!mutex.is_shared_locked());

  BOOST_CHECK(mutex.try_lock_shared());
  BOOST_CHECK(mutex.is_shared_locked());
  mutex.unlock_shared();
  BOOST_CHECK(!mutex.is_shared_locked());
}

BOOST_AUTO_TEST_CASE(standard_guards)
{
  shared_mutex mutex;

  {
    // std::lock_guard works on the exclusive interface.
    std::lock_guard<shared_mutex> guard(mutex);
  }

  {
    // std::unique_lock works on the exclusive interface.
    std::unique_lock<shared_mutex> guard(mutex);
    BOOST_CHECK(guard.owns_lock());

    guard.unlock();
    BOOST_CHECK(!guard.owns_lock());

    guard.lock();
    BOOST_CHECK(guard.owns_lock());
  }

  {
    // std::unique_lock with try_to_lock uses try_lock.
    std::unique_lock<shared_mutex> guard(mutex, std::try_to_lock);
    BOOST_CHECK(guard.owns_lock());
  }

  {
    // std::shared_lock works on the shared interface.
    std::shared_lock<shared_mutex> guard(mutex);
    BOOST_CHECK(guard.owns_lock());
    BOOST_CHECK(mutex.is_shared_locked());

    guard.unlock();
    BOOST_CHECK(!guard.owns_lock());
    BOOST_CHECK(!mutex.is_shared_locked());

    guard.lock();
    BOOST_CHECK(guard.owns_lock());
    BOOST_CHECK(mutex.is_shared_locked());
  }

  BOOST_CHECK(!mutex.is_shared_locked());

  {
    // The aliases used throughout the code base are exactly these guards.
    shared_guard shared(mutex);
    BOOST_CHECK(mutex.is_shared_locked());
  }
  {
    lock_guard exclusive(mutex);
  }

  {
    // Guards are movable, as in the standard library.
    shared_guard outer = [&]()
    {
      return shared_guard(mutex);
    }();
    BOOST_CHECK(outer.owns_lock());
    BOOST_CHECK(mutex.is_shared_locked());
  }
  BOOST_CHECK(!mutex.is_shared_locked());
}

BOOST_AUTO_TEST_CASE(exclusive_blocks_shared)
{
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
  {
    // Two mutexes sharing the same group; holding one exclusively forbids
    // shared locking of the others, which try_lock_shared reports without blocking.
    shared_mutex first;
    shared_mutex second(first);

    first.lock();
    BOOST_CHECK(!second.try_lock_shared());
    first.unlock();

    BOOST_CHECK(second.try_lock_shared());
    second.unlock_shared();
  }
}
