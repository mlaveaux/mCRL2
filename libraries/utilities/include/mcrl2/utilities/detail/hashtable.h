// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_DETAIL_HASHTABLE_H
#define MCRL2_UTILITIES_DETAIL_HASHTABLE_H
#pragma once

#include "mcrl2/utilities/hashtable.h"    // necessary for header test.
#include "mcrl2/utilities/indexed_set.h"    // necessary for header test.

namespace mcrl2
{
namespace utilities
{

template <class Key, typename Hash, typename Equals, typename Allocator>
inline void hashtable<Key, Hash, Equals, Allocator>::rehash(std::size_t size)
{
  // Copy the old hashtable.
  std::vector<Key> old = std::move(m_hashtable);

  m_hashtable = std::vector<Key>(size, nullptr);
  m_buckets_mask = m_hashtable.size() - 1;

  for (const Key& key : old)
  {
    auto it = begin() + get_index(key);
    // Find the first empty position.
    while (*it != nullptr)
    {
      ++it;
      if (it == end())
      {
        it = begin();
      }

      assert(it != begin() + get_index(key));
    }

    // Found an empty spot, insert a new index belonging to key,
    *it = key;
  }
}

template <class Key, typename Hash, typename Equals, typename Allocator>
inline hashtable<Key,Hash,Equals,Allocator>::hashtable()
  : hashtable(128)
{} 

template <class Key, typename Hash, typename Equals, typename Allocator>
inline hashtable<Key,Hash,Equals,Allocator>::hashtable(std::size_t initial_size,
  const hasher& hasher,
  const key_equal& equals)
      : m_hashtable(std::max(initial_size, detail::minimal_hashtable_size), nullptr),
        m_hasher(hasher),
        m_equals(equals)
{
  m_buckets_mask = m_hashtable.size() - 1;
}

template <class Key, typename Hash, typename Equals, typename Allocator>
inline void hashtable<Key,Hash,Equals,Allocator>::clear()
{
  m_hashtable.clear();
}

template <class Key, typename Hash, typename Equals, typename Allocator>
inline std::pair<typename hashtable<Key,Hash,Equals,Allocator>::iterator, bool> hashtable<Key,Hash,Equals,Allocator>::insert(const Key& key)
{
  ++m_number_of_elements;
  // Resize hashtable if necessary.
  if (2 * m_number_of_elements >= m_hashtable.size())
  {
    rehash(2 * m_hashtable.size());
  }

  auto it = begin() + get_index(key);

  // Find the first empty position.
  while (*it != nullptr)
  {
    ++it;
    if (it == end())
    {
      it = begin();
    }

    assert(it != begin() + get_index(key));
  }

  // Found an empty spot, insert a new index belonging to key,
  *it = key;
  return std::make_pair(it, true);
}


template <class Key, typename Hash, typename Equals, typename Allocator>
inline typename hashtable<Key,Hash,Equals,Allocator>::iterator hashtable<Key,Hash,Equals,Allocator>::erase(const Key& key)
{
  auto it = begin() + get_index(key);

  // Find the key.
  while (!m_equals(*it, key))
  {
    ++it;
    if (it == end())
    {
      it = begin();
    }

    if (it == begin() + get_index(key))
    {
      // This is a workaround for the thread_local hashsets that are deleted before the global aterms.
      return it;
    }
  }

  *it = nullptr;
  --m_number_of_elements;
  return it;
}

// PRIVATE FUNCTIONS

template <class Key, typename Hash, typename Equals, typename Allocator>
inline std::size_t hashtable<Key,Hash,Equals,Allocator>::get_index(const Key& key)
{
  return m_hasher(key) & m_buckets_mask;
}


} // namespace utilities

} // namespace mcrl2

#endif // MCRL2_UTILITIES_DETAIL_INDEXED_SET_H
