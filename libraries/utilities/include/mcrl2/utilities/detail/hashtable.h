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
inline void hashtable<Key, Hash, Equals, Allocator>::resize_if_needed()
{
  // Resize hashtable if necessary.
  if (2 * m_number_of_elements >= m_hashtable.size())
  {
    // Copy the old hashtable.
    std::vector<Key> old = std::move(m_hashtable);

    m_hashtable = std::vector<Key>(old.size() * 2, nullptr);
    m_buckets_mask = m_hashtable.size() - 1;

    for (const Key& key : old)
    {
      // Find a place to insert key and find whether key already exists.
      std::size_t start = get_index(key);
      std::size_t position = start;

      do
      {
        const Key& element = m_hashtable[position];

        if (element == nullptr)
        {
          // Found an empty spot, insert a new index belonging to key,
          m_hashtable[position] = key;
          break;
        }

        position = next_index(position);
      }      
      while (start != position);
    }
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
  resize_if_needed();

  // Find a place to insert key and find whether key already exists.
  std::size_t start = get_index(key);
  std::size_t position = start;

  do
  {
    const Key& element = m_hashtable[position];

    if (element == nullptr)
    {
      // Found an empty spot, insert a new index belonging to key,
      m_hashtable[position] = key;
      ++m_number_of_elements;
      return std::make_pair(m_hashtable.begin() + position, true);
    }

    position = next_index(position);
  }
  while (position != start);

  return std::make_pair(m_hashtable.end(), false);
}


template <class Key, typename Hash, typename Equals, typename Allocator>
inline typename hashtable<Key,Hash,Equals,Allocator>::iterator hashtable<Key,Hash,Equals,Allocator>::erase(const Key& key)
{
  // Find the key.
  std::size_t start = get_index(key);
  std::size_t position = start;

  do
  {
    const Key& element = m_hashtable[position];

    if (element == key)
    {
      m_hashtable[position] = nullptr;
      --m_number_of_elements;

      // key is already in the set, return position of key.
      return m_hashtable.begin() + position;
    }

    position = next_index(position);
  }
  while (position != start);

  return m_hashtable.end();
}

// PRIVATE FUNCTIONS

template <class Key, typename Hash, typename Equals, typename Allocator>
inline std::size_t hashtable<Key,Hash,Equals,Allocator>::get_index(const Key& key)
{
  return std::hash<Key>()(key) * detail::PRIME_NUMBER & m_buckets_mask;
}

template <class Key, typename Hash, typename Equals, typename Allocator>
inline std::size_t hashtable<Key,Hash,Equals,Allocator>::next_index(std::size_t index)
{
  return (index + 1) & m_buckets_mask;
}

} // namespace utilities

} // namespace mcrl2

#endif // MCRL2_UTILITIES_DETAIL_INDEXED_SET_H
