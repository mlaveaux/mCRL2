// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file utilities/detail/indexed_set.cpp
/// \brief This file contains some constants and functions shared
///        between indexed_sets and tables.

#ifndef MCRL2_UTILITIES_DETAIL_INDEXED_SET_H
#define MCRL2_UTILITIES_DETAIL_INDEXED_SET_H

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/indexed_set.cxx"    // necessary for header test.
#endif

namespace mcrl2::utilities
{
namespace detail
{

static const std::size_t STEP = 1; ///< The position on which the next hash entry is searched.

/// in the hashtable we use the following constant to indicate free positions.
static constexpr std::size_t EMPTY(std::numeric_limits<std::size_t>::max());

static constexpr std::size_t RESERVED(std::numeric_limits<std::size_t>::max()-1);

static constexpr float max_load_factor = 0.6f; ///< The load factor before the hash table is resized.

static constexpr std::size_t PRIME_NUMBER = 999953;

#ifndef NDEBUG  // Numbers are small in debug mode for more intensive checks. 
  static constexpr std::size_t minimal_hashtable_size = 16; 
#else
  static constexpr std::size_t minimal_hashtable_size = 2048;
#endif
  static constexpr std::size_t RESERVATION_FRACTION = 8;        // If the reserved keys entries are exploited, 1/RESERVATION_FRACTION new
                                                                // keys are reserved. This is an expensive operation, as it is executed
                                                                // using a lock_exclusive. 

static_assert(minimal_hashtable_size/RESERVATION_FRACTION != 0);
static_assert(minimal_hashtable_size>=8);       ///< With a max_load of 0.75 the minimal size of the hashtable must be 8.

} // namespace detail

#define INDEXED_SET_TEMPLATE template <class Key, bool ThreadSafe, typename Hash, typename Equals, typename Allocator, typename KeyTable>
#define INDEXED_SET indexed_set<Key, ThreadSafe, Hash, Equals, Allocator, KeyTable>

INDEXED_SET_TEMPLATE
inline void INDEXED_SET::reserve_indices(const std::size_t thread_index)
{
  lock_guard guard = m_shared_mutexes[thread_index].lock();

  if (m_next_index + m_shared_mutexes.size() >= m_keys.size())   // otherwise another process already reserved entries, and nothing needs to be done. 
  {
    assert(m_next_index <= m_keys.size());
    m_keys.resize(m_keys.size() + std::max(m_keys.size() / detail::RESERVATION_FRACTION, m_shared_mutexes.size()));  // Increase with at least the number of threads. 

    while ((detail::max_load_factor * m_hashtable.size()) < m_keys.size())
    {
       resize_hashtable();
    }
  }
}

INDEXED_SET_TEMPLATE
inline typename INDEXED_SET::size_type INDEXED_SET::put_in_hashtable(
                  const key_type& key, 
                  std::size_t value, 
                  std::size_t& new_position)
{
  // Find a place to insert key and find whether key already exists.
  assert(m_hashtable.size() > 0);

  new_position = ((m_hasher(key) * detail::PRIME_NUMBER) >> 2) % m_hashtable.size();
  [[maybe_unused]] // Not used in release mode
  std::size_t start = new_position;

  while (true)
  {
    std::size_t index = m_hashtable[new_position];
    assert(index == detail::EMPTY || index == detail::RESERVED || index < m_keys.size());

    if (index == detail::EMPTY)
    {
      // Found an empty spot, insert a new index belonging to key,
      std::size_t pos=detail::EMPTY;
      if (reinterpret_cast<std::atomic<std::size_t>*>(&m_hashtable[new_position])->compare_exchange_strong(pos,value))
      {
        return value;
      }
      index=pos;             // Insertion failed, but another process put an alternative value "pos"
                             // at this position. 
    }

    // If the index is RESERVED, we go into a busy loop, as another process 
    // will shortly change the RESERVED value into a sensible index. 
    if (index != detail::RESERVED) 
    {
      assert(index!=detail::EMPTY);
      if (m_equals(m_keys[index], key))
      {
        // key is already in the set, return position of key.
        assert(index<m_next_index && m_next_index<=m_keys.size());
        return index;
      }
      assert(m_hashtable.size()>0);
      new_position = (new_position + detail::STEP) % m_hashtable.size();
      assert(new_position != start); // In this case the hashtable is full, which should never happen.
    }
  }

  // not reached. 
  std::abort();
  return detail::EMPTY;
}

INDEXED_SET_TEMPLATE
inline void INDEXED_SET::resize_hashtable()
{
  m_hashtable.assign(m_hashtable.size() * 2, detail::EMPTY);
  size_t index = 0;
  for (const Key& k: m_keys)
  {
    if (index<m_next_index)
    {
      std::size_t new_position;  // The resulting new_position is not used here. 
      put_in_hashtable(k, index, new_position);
    }
    else 
    {
      break;
    }
    ++index;
  }
}

INDEXED_SET_TEMPLATE
inline INDEXED_SET::indexed_set()
  : indexed_set(1, detail::minimal_hashtable_size)   // Run with one main thread. 
{} 

INDEXED_SET_TEMPLATE
inline INDEXED_SET::indexed_set(std::size_t number_of_threads)
  : indexed_set(number_of_threads, detail::minimal_hashtable_size)
{
  assert(number_of_threads != 0);
}

INDEXED_SET_TEMPLATE
inline INDEXED_SET::indexed_set(
           std::size_t number_of_threads,
           std::size_t initial_size,
           const hasher& hasher,
           const key_equal& equals)
  : m_hashtable(std::max(initial_size, detail::minimal_hashtable_size), detail::EMPTY), 
    m_mutex(new std::mutex()),
    m_hasher(hasher),
    m_equals(equals)
{
  assert(number_of_threads != 0);

  // Insert the main mutex.
  m_shared_mutexes.emplace_back();

  for (std::size_t i = 1; i < ((number_of_threads == 1) ? 1 : number_of_threads + 1); ++i)
  {
    // Copy the mutex n times for all the other threads.
    m_shared_mutexes.emplace_back(m_shared_mutexes[0]);
  }
}

INDEXED_SET_TEMPLATE
inline typename INDEXED_SET::size_type INDEXED_SET::index(const key_type& key, const std::size_t thread_index) const
{
  shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
  assert(m_hashtable.size() > 0);

  std::size_t start = ((m_hasher(key) * detail::PRIME_NUMBER) >> 2) % m_hashtable.size();
  std::size_t position = start;
  do
  {
    std::size_t index = m_hashtable[position];
    if (index == detail::EMPTY)
    {
      return npos; // Not found.
    }
    // If the index is RESERVED, go into a busy loop. Another thread will 
    // change this RESERVED index shortly into a sensible index. 
    if (index != detail::RESERVED)
    {
      assert(index < m_keys.size());
      if (m_equals(key, m_keys[index]))
      {
        assert(index<m_next_index && m_next_index <= m_keys.size());
        return index;
      }

      assert(m_hashtable.size() > 0);
      position = (position + detail::STEP) % m_hashtable.size();
      assert(position != start); // The hashtable is full. This should never happen.
    }
  }
  while (true);

  std::abort();
  return npos; // Dummy return.
}

INDEXED_SET_TEMPLATE
inline typename INDEXED_SET::const_iterator INDEXED_SET::find(const key_type& key, const std::size_t thread_index) const
{
  const std::size_t idx = index(key, thread_index);
  if (idx < m_keys.size())
  {
    return begin(thread_index) + idx;
  }

  return end(thread_index);
}


INDEXED_SET_TEMPLATE
inline const Key& INDEXED_SET::at(std::size_t index) const
{
  if (index >= m_next_index)
  {
    throw std::out_of_range("indexed_set: index too large: " + std::to_string(index) + " > " + std::to_string(m_next_index) + ".");
  }

  return m_keys[index];
}

INDEXED_SET_TEMPLATE
inline const Key& INDEXED_SET::operator[](std::size_t index) const
{
  assert(index<m_keys.size());
  const Key& key = m_keys[index];
  return key;
}

INDEXED_SET_TEMPLATE
inline void INDEXED_SET::clear(const std::size_t thread_index)
{
  lock_guard guard = m_shared_mutexes[thread_index].lock();
  m_hashtable.assign(m_hashtable.size(), detail::EMPTY);

  m_keys.clear();
  m_next_index.store(0);
}


INDEXED_SET_TEMPLATE
inline std::pair<typename INDEXED_SET::size_type, bool> INDEXED_SET::insert(const Key& key, const std::size_t thread_index)
{
  shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
  assert(m_next_index <= m_keys.size());
  if (m_next_index + m_shared_mutexes.size() >= m_keys.size())
  {
    guard.unlock_shared();
    reserve_indices(thread_index);
    guard.lock_shared();
  }
  std::size_t new_position;
  const std::size_t index = put_in_hashtable(key, detail::RESERVED, new_position);
  
  if (index != detail::RESERVED) // Key already exists.
  {
    assert(index < m_next_index && m_next_index <= m_keys.size());
    return std::make_pair(index, false);
  }

  const std::size_t new_index = m_next_index.fetch_add(1);
  assert(new_index < m_keys.size());
  m_keys[new_index] = key; 

  std::atomic_thread_fence(std::memory_order_seq_cst);   // Necessary for ARM. std::memory_order_acquire and 
                                                            // std::memory_order_release appear to work, too.
  m_hashtable[new_position] = new_index;


  assert(new_index < m_next_index && m_next_index <= m_keys.size());
  return std::make_pair(new_index, true);
}

#undef INDEXED_SET_TEMPLATE 
#undef INDEXED_SET 


} // namespace mcrl2::utilities



#endif // MCRL2_UTILITIES_DETAIL_INDEXED_SET_H
