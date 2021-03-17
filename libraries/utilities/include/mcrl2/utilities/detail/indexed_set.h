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
#pragma once

#include "mcrl2/utilities/indexed_set.h"    // necessary for header test. 

namespace mcrl2
{
namespace utilities
{
namespace detail
{

static const std::size_t STEP = 1; ///< The position on which the next hash entry is searched.

/// in the hashtable we use the following constant to indicate free positions.
static constexpr std::size_t EMPTY(std::numeric_limits<std::size_t>::max());

static constexpr float max_load_factor = 0.75f; ///< The load factor before the hash table is resized.

static constexpr size_t minimal_hashtable_size = 8;  ///< With a max_load of 0.75 the minimal size of the hashtable must be 8.

static constexpr std::size_t PRIME_NUMBER = 999953;

} // namespace detail

template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline
std::size_t indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::put_in_hashtable(const key_type& key, std::size_t value)
{
  // Find a place to insert key and find whether key already exists.
  std::size_t start = (m_hasher(key) * detail::PRIME_NUMBER) % m_hashtable.size();
  std::size_t position = start;

  while (true)
  {
    std::size_t index = m_hashtable[position];
    assert(index == detail::EMPTY || index < m_keys.size());

    if (index == detail::EMPTY)
    {
      // Found an empty spot, insert a new index belonging to key,
      m_hashtable[position] = value;
      return value;
    }

    if (m_equals(m_keys[index], key))
    {
      // key is already in the set, return position of key.
      return index;
    }
    position = (position + detail::STEP) % m_hashtable.size();

    assert(position != start); // In this case the hashtable is full, which should never happen.
  }

  return position;
}


template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline void indexed_set<Key, Hash, Equals, Allocator, ThreadSafe>::resize_hashtable()
{
  if constexpr (ThreadSafe) { m_mutex.lock(); }
  m_hashtable = std::vector<std::size_t>(m_hashtable.size() * 2, detail::EMPTY);

  size_t index = 0;
  for (const Key& k : m_keys)
  {
    put_in_hashtable(k, index);
    ++index;
  }
  if constexpr (ThreadSafe) { m_mutex.unlock(); }
}

template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::indexed_set()
  : indexed_set(128)
{
} 

template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::indexed_set(std::size_t initial_size,
  const hasher& hasher,
  const key_equal& equals)
      : m_hashtable(std::max(initial_size, detail::minimal_hashtable_size), detail::EMPTY),
        m_hasher(hasher),
        m_equals(equals)
{}

template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline typename indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::size_type indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::index(const key_type& key) const
{
  std::size_t start = (m_hasher(key) * detail::PRIME_NUMBER) % m_hashtable.size();
  std::size_t position = start;

  if constexpr (ThreadSafe) { m_mutex.lock(); }
  do
  {
    std::size_t index = m_hashtable[position];
    if (index == detail::EMPTY)
    {
      if constexpr (ThreadSafe) { m_mutex.unlock(); }
      return npos; // Not found.
    }
    assert(index < m_keys.size());
    if (m_equals(key, m_keys[index]))
    {
      if constexpr (ThreadSafe) { m_mutex.unlock(); }
      return index;
    }

    position = (position + detail::STEP) % m_hashtable.size();
    assert(position!=start); // The hashtable is full. This should never happen.
  }
  while (true);

  if constexpr (ThreadSafe) { m_mutex.unlock(); }
  return npos; // Not found.
}

template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline typename indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::const_iterator indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::find(const key_type& key) const
{
  const std::size_t idx = index(key);
  if (idx < m_keys.size())
  {
    return begin() + idx;
  }

  return end();
}


template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline const Key& indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::at(std::size_t index) const
{
  if (index >= m_keys.size())
  {
    throw std::out_of_range("indexed_set: index too large: " + std::to_string(index) + " > " + std::to_string(m_keys.size()) + ".");
  }

  if constexpr (ThreadSafe) { m_mutex.lock(); }
  const Key& key = m_keys[index];
  if constexpr (ThreadSafe) { m_mutex.unlock(); }

  return key;
}

template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline const Key& indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::operator[](std::size_t index) const
{
  assert(index<m_keys.size());

  if constexpr (ThreadSafe) { m_mutex.lock(); }
  const Key& key = m_keys[index];
  if constexpr (ThreadSafe) { m_mutex.unlock(); }

  return key;
}

template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline void indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::clear()
{
  m_hashtable.assign(m_hashtable.size(), detail::EMPTY);
  m_keys.clear();
}


template <class Key, typename Hash, typename Equals, typename Allocator, bool ThreadSafe>
inline std::pair<std::size_t, bool> indexed_set<Key,Hash,Equals,Allocator,ThreadSafe>::insert(const Key& key)
{
  if constexpr (ThreadSafe) { m_mutex.lock(); }

  const std::size_t new_index = m_keys.size();
  const std::size_t index = put_in_hashtable(key, new_index);

  if (index != new_index) // Key already exists.
  {
    if constexpr (ThreadSafe) { m_mutex.unlock(); }
    return std::make_pair(index, false);
  }

  m_keys.push_back(key);

  if ((detail::max_load_factor * m_hashtable.size()) < m_keys.size())
  {
    if constexpr (ThreadSafe) { m_mutex.unlock(); }
    resize_hashtable();
  }

  if constexpr (ThreadSafe) { m_mutex.unlock(); }
  return std::make_pair(index, true);
}


} // namespace utilities

} // namespace mcrl2

#endif // MCRL2_UTILITIES_DETAIL_INDEXED_SET_H
