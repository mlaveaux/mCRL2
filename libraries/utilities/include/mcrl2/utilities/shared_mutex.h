// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_DETAIL_SHARED_MUTEX_H
#define MCRL2_UTILITIES_DETAIL_SHARED_MUTEX_H

#include <algorithm>
#include <atomic>
#include <cassert>
#include <concepts>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <vector>

#include "mcrl2/utilities/configuration.h"
#include "mcrl2/utilities/hardware_interference_size.h"

namespace mcrl2::utilities
{

// Forward declaration.
class shared_mutex;

struct shared_mutex_data 
{
  /// \brief The list of other mutexes.
  std::vector<shared_mutex*> other;

  /// \brief Mutex for adding/removing shared_guards.
  std::mutex mutex;  
  
  /// Adds a shared mutex to the data.
  inline
  void register_mutex(shared_mutex* shared_mutex)
  {
    std::lock_guard guard(mutex);
    other.emplace_back(shared_mutex);
  }
  
  // Removes a shared mutex from the data
  inline void unregister_mutex(shared_mutex* shared_mutex) 
  {
    std::lock_guard guard(mutex);
    auto it = std::find(other.begin(), other.end(), shared_mutex);
    assert(it != other.end());

    other.erase(it);
  }
};

/// An implementation of a shared mutex (also called readers-write lock in the literature) based on
/// the notion of busy and forbidden flags.
class alignas(hardware_destructive_interference_size) shared_mutex
{
public:
  shared_mutex()
    : m_shared(std::make_shared<shared_mutex_data>())
  {
    m_shared->register_mutex(this);
  }

  ~shared_mutex() 
  { 
    m_shared->unregister_mutex(this);
  }

  /// The copy/move constructor/assignment should not be called while any lock_guard or shared_guard is alive.
  shared_mutex(const shared_mutex& other)
    : m_shared(other.m_shared)
  {
    m_shared->register_mutex(this);
  }

  shared_mutex(shared_mutex&& other) noexcept
      : m_shared(other.m_shared)
  {
    m_shared->register_mutex(this);
  }

  shared_mutex& operator=(const shared_mutex& other)
  {
    if (this != &other)
    {
      // Remove ourselves, and register into the other shared.
      m_shared->unregister_mutex(this);

      m_shared = other.m_shared;
      m_shared->register_mutex(this);
    }
    return *this;
  }

  shared_mutex& operator=(shared_mutex&& other) noexcept
  {
    if (this != &other)
    {
      m_shared->unregister_mutex(this);

      m_shared = other.m_shared;
      m_shared->register_mutex(this);
      m_shared->unregister_mutex(&other);
    }
    return *this;
  }

  /// \brief Obtain exclusive access, and stop all other threads that use this mutex.
  /// \details Equivalent to std::shared_mutex::lock. Blocks until exclusive access is acquired.
  /// \pre The calling thread does not hold a (shared) lock on this mutex.
  /// \post All threads sharing this mutex are suspended outside of their shared sections.
  inline
  void lock()
  {
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
    {
      // Only one thread can halt everything.
      m_shared->mutex.lock();
      set_forbidden_flags();
    }
  }

  /// \brief Try to obtain exclusive access without blocking on other exclusive locks.
  /// \details Equivalent to std::shared_mutex::try_lock.
  /// \pre The calling thread does not hold a (shared) lock on this mutex.
  /// \returns True iff exclusive access was acquired.
  [[nodiscard]]
  inline
  bool try_lock()
  {
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
    {
      // Try to obtain the lock. If the surrounding mutex cannot be locked, this fails.
      if (!m_shared->mutex.try_lock())
      {
        return false;
      }

      set_forbidden_flags();
    }
    return true;
  }

  /// \brief Release exclusive access.
  /// \details Equivalent to std::shared_mutex::unlock.
  /// \pre The calling thread holds the exclusive lock.
  inline
  void unlock()
  {
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
    {
      for (shared_mutex* mutex : m_shared->other)
      {
        mutex->set_forbidden(false);
      }

      m_shared->mutex.unlock();
    }
  }

  /// \brief Acquires a shared lock on this mutex.
  /// \details Equivalent to std::shared_mutex::lock_shared, except that recursive shared
  ///          locking by the owning thread is explicitly allowed: every lock_shared() must
  ///          be balanced by exactly one unlock_shared().
  /// \post No thread can acquire the exclusive lock until the last unlock_shared().
  inline
  void lock_shared()
  {
    if (mcrl2::utilities::detail::GlobalThreadSafe && m_lock_depth == 0)
    {
      assert(!m_busy_flag);
      m_busy_flag.store(true);

      // Wait for the forbidden flag to become false.
      while (m_forbidden_flag.load())
      {
        m_busy_flag = false;

        // Wait for the global lock.
        m_shared->mutex.lock();
        m_shared->mutex.unlock();

        m_busy_flag = true;
      }
    }

    ++m_lock_depth;
  }

  /// \brief Tries to acquire a shared lock on this mutex without blocking.
  /// \details Equivalent to std::shared_mutex::try_lock_shared.
  /// \returns True iff the shared lock was acquired.
  [[nodiscard]]
  inline
  bool try_lock_shared()
  {
    if (mcrl2::utilities::detail::GlobalThreadSafe && m_lock_depth == 0)
    {
      assert(!m_busy_flag);
      m_busy_flag.store(true);

      if (m_forbidden_flag.load())
      {
        // An exclusive lock is being held or requested.
        m_busy_flag.store(false);
        return false;
      }
    }

    ++m_lock_depth;
    return true;
  }

  /// \brief Releases a shared lock on this mutex.
  /// \details Equivalent to std::shared_mutex::unlock_shared.
  /// \pre The calling thread holds a shared lock on this mutex.
  inline
  void unlock_shared()
  {
    assert(!mcrl2::utilities::detail::GlobalThreadSafe || m_lock_depth > 0);

    --m_lock_depth;
    if (mcrl2::utilities::detail::GlobalThreadSafe && m_lock_depth == 0)
    {
      assert(m_busy_flag);
      m_busy_flag.store(false, std::memory_order_release);
    }
  }

  /// \returns True iff the shared mutex is in the shared section
  bool is_shared_locked() const
  {
    return m_lock_depth != 0;
  }

private:
  inline
  void set_forbidden_flags()
  {
    // Shared and exclusive sections MUST be disjoint.
    assert(!m_busy_flag);

    assert(std::find(m_shared->other.begin(), m_shared->other.end(), this) != m_shared->other.end());

    // Indicate that threads must wait.
    for (shared_mutex* mutex : m_shared->other)
    {
      if (mutex != this)
      {
        mutex->set_forbidden(true);
      }
    }

    // Wait for all pools to indicate that they are not busy.
    for (const shared_mutex* mutex : m_shared->other)
    {
      if (mutex != this)
      {
        mutex->wait_for_busy();
      }
    }
  }

  /// \returns True iff the shared mutex has its busy flag set. 
  inline
  bool is_busy() const
  {
    return m_busy_flag.load();
  }

  /// \brief Waits for the busy flag to become false.
  inline
  void wait_for_busy() const
  {
    while (m_busy_flag.load()) { /* wait */ };
  }

  inline
  void set_forbidden(bool value)
  {
    m_forbidden_flag.store(value);
  }

  /// \brief A boolean flag indicating whether this thread is working inside the global aterm pool.
  std::atomic<bool> m_busy_flag = false;
  std::atomic<bool> m_forbidden_flag = false;

  /// \brief It can happen that un/lock_shared calls are nested, so keep track of the nesting depth and only
  ///        actually perform un/locking at the root.
  std::size_t m_lock_depth = 0;

  std::shared_ptr<shared_mutex_data> m_shared;
};

namespace detail
{

/// \brief The Cpp17Lockable named requirement; the standard library provides no concept for it.
template<typename Mutex>
concept IsLockable = requires(Mutex m) {
  { m.lock() } -> std::same_as<void>;
  { m.unlock() } -> std::same_as<void>;
  { m.try_lock() } -> std::same_as<bool>;
};

/// \brief The Cpp17SharedLockable named requirement; the standard library provides no concept for it.
template<typename Mutex>
concept IsSharedLockable = requires(Mutex m) {
  { m.lock_shared() } -> std::same_as<void>;
  { m.try_lock_shared() } -> std::same_as<bool>;
  { m.unlock_shared() } -> std::same_as<void>;
};

} // namespace detail

// The interface must behave as std::shared_mutex, such that the standard library lock guards apply.
static_assert(detail::IsLockable<shared_mutex>,
    "shared_mutex must satisfy Cpp17Lockable, like std::shared_mutex, for use with std::unique_lock");
static_assert(detail::IsSharedLockable<shared_mutex>,
    "shared_mutex must satisfy Cpp17SharedLockable, like std::shared_mutex, for use with std::shared_lock");

/// \brief An exclusive lock guard for the shared_mutex, as in the standard library.
using lock_guard = std::unique_lock<shared_mutex>;

/// \brief A shared lock guard for the shared_mutex, as in the standard library.
using shared_guard = std::shared_lock<shared_mutex>;

} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_DETAIL_SHARED_MUTEX_H
