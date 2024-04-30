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

#include <assert.h>
#include <atomic>
#include <memory>
#include <shared_mutex>

#include "mcrl2/utilities/noncopyable.h"
#include "mcrl2/utilities/configuration.h"


namespace mcrl2::utilities
{

// Forward declaration.
class shared_mutex;

/// A shared lock guard for the shared_mutex.
class shared_guard : private mcrl2::utilities::noncopyable
{
public:
  
  /// Locks the guard again explicitly.
  inline
  void lock_shared();

  /// Unlocks the acquired shared guard explicitly. Otherwise, performed in destructor.
  inline
  void unlock_shared();

  ~shared_guard()
  {
    if (is_locked)
    {
      unlock_shared();
    }
  }

private:
  friend class shared_mutex;

  shared_guard(shared_mutex& mutex)
    : m_mutex(mutex)
  {}

  shared_mutex& m_mutex;
  bool is_locked = true;
};

/// An exclusive lock guard for the shared_mutex.
class lock_guard : private mcrl2::utilities::noncopyable
{
public:
  /// Unlocks the acquired shared guard explicitly. Otherwise, performed in destructor.
  void unlock();

  ~lock_guard()
  {
    if (is_locked)
    {
      unlock();
    }
  }

private:
  friend class shared_mutex;

  lock_guard(shared_mutex& mutex)
    : m_mutex(mutex)
  {}

  shared_mutex& m_mutex;
  bool is_locked = true;
};

const std::size_t NUM_SHARED_MUTEXES = 8;

namespace 
{
  struct padding
  {
    alignas(64)
    std::shared_mutex mutex;
  };

  static_assert(sizeof(padding) == 64);
}

struct shared_mutex_data 
{
  /// \brief The list of other mutexes.
  std::array<padding, NUM_SHARED_MUTEXES> mutexes;
  std::atomic<std::size_t> m_next = 0;
  
  /// Returns the next index in the queue
  inline
  std::size_t register_mutex()
  {
    m_next = (m_next + 1) % NUM_SHARED_MUTEXES;
    return m_next;
  }
};

/// An implementation of a shared mutex (also called readers-write lock in the literature) based on
/// the notion of busy and forbidden flags.
class shared_mutex
{
public:
  shared_mutex()
    : m_shared(std::make_shared<shared_mutex_data>())
  {
    m_index = m_shared->register_mutex();
  }

  shared_mutex(const shared_mutex& other)
    : m_shared(other.m_shared)
  {
    m_index = m_shared->register_mutex();
  }
  
  // Obtain exclusive access to the busy-forbidden lock.
  inline
  lock_guard lock()
  {
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
    {
      // Shared and exclusive sections MUST be disjoint.
      assert(m_lock_depth == 0);

      m_lock_depth = 1;

      for (auto& padding : m_shared->mutexes)
      {
        padding.mutex.lock();
      }
    }

    return lock_guard(*this);
  }

  // Release exclusive access to the busy-forbidden lock.
  inline
  void unlock()
  {
    unlock_impl();
  }

  /// Acquires a shared lock on this instance, returns a shared guard that keeps the lock until it is destroyed.
  /// Or alternative, unlock_shared is called explicitly.
  inline
  shared_guard lock_shared()
  {
    lock_shared_impl();
    return shared_guard(*this);
  }

  /// \returns True iff the shared mutex is in the shared section
  bool is_shared_locked() const
  {
    return m_lock_depth != 0;
  }

private:
  friend class lock_guard;
  friend class shared_guard;
  friend class shared_mutex_pool;
  
  void unlock_impl()
  {
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
    {
      // Unlock in the reverse order.
      for (auto it = m_shared->mutexes.rbegin(); it != m_shared->mutexes.rend(); ++it)
      {
        it->mutex.unlock();
      }

      assert(m_lock_depth == 1);
      m_lock_depth = 0;
    }
  }
    
  inline
  void lock_shared_impl()
  {
    if (mcrl2::utilities::detail::GlobalThreadSafe && m_lock_depth == 0)
    {
      m_shared->mutexes[m_index].mutex.lock_shared();
    }

    ++m_lock_depth;
  }

  // Release shared access to the busy-forbidden lock.
  inline
  void unlock_shared()
  {
    assert(!mcrl2::utilities::detail::GlobalThreadSafe || m_lock_depth > 0);

    --m_lock_depth;
    if (mcrl2::utilities::detail::GlobalThreadSafe && m_lock_depth == 0)
    {
      m_shared->mutexes[m_index].mutex.unlock_shared();
    }
  }

  /// \brief A boolean flag indicating whether this thread is working inside the global aterm pool.
  std::size_t m_index = 0;

  /// \brief It can happen that un/lock_shared calls are nested, so keep track of the nesting depth and only
  ///        actually perform un/locking at the root.
  std::size_t m_lock_depth = 0;

  std::shared_ptr<shared_mutex_data> m_shared;
};

inline
void shared_guard::lock_shared()   
{    
  // Uses the internal implementation since we don't need a shared_guard.
  m_mutex.lock_shared_impl();
  is_locked = true;
}

inline
void shared_guard::unlock_shared()
{    
  m_mutex.unlock_shared();
  is_locked = false;
}

inline
void lock_guard::unlock()
{
  m_mutex.unlock();
  is_locked = false;
}

} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_DETAIL_SHARED_MUTEX_H
