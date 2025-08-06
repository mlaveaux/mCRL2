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

#include "mcrl2/utilities/hashtable.h"    // necessary for header test.

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/power_of_two.cxx" 
  #include "mcrl2/utilities/indexed_set.cxx"    // necessary for header test.
#else
  import utilities;
#endif

namespace mcrl2::utilities
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
    const std::size_t key_index = get_index(key);
    auto it = begin() + key_index;
    // Find the first empty position.
    while (*it != nullptr)
    {
      ++it;
      if (it == end())
      {
        it = begin();
      }

      assert(it != begin() + key_index);
    }

    // Found an empty spot, insert a new index belonging to key,
    *it = key;
  }
}

template <class Key, typename Hash, typename Equals, typename Allocator>
inline hashtable<Key,Hash,Equals,Allocator>::hashtable()
  : hashtable(128)
{
} 

template <class Key, typename Hash, typename Equals, typename Allocator>
inline hashtable<Key,Hash,Equals,Allocator>::hashtable(std::size_t initial_size,
  const hasher& hasher,
  const key_equal& equals)
      : m_hashtable(std::max(initial_size, detail::minimal_hashtable_size), nullptr),
        m_hasher(hasher),
        m_equals(equals)
{
  assert(utilities::is_power_of_two(initial_size));
  m_buckets_mask = m_hashtable.size() - 1;
}

template <class Key, typename Hash, typename Equals, typename Allocator>
inline void hashtable<Key,Hash,Equals,Allocator>::clear()
{
  m_hashtable.clear();
}

template <class Key, typename Hash, typename Equals, typename Allocator>
bool hashtable<Key,Hash,Equals,Allocator>::must_resize()
{
  return (2 * m_number_of_elements >= m_hashtable.size());
}

template <class Key, typename Hash, typename Equals, typename Allocator>
void hashtable<Key,Hash,Equals,Allocator>::resize()
{
  rehash(2 * m_hashtable.size());
}

template <class Key, typename Hash, typename Equals, typename Allocator>
inline std::pair<typename hashtable<Key,Hash,Equals,Allocator>::iterator, bool> hashtable<Key,Hash,Equals,Allocator>::insert(const Key& key)
{
  // Resize hashtable must be done explicitly.
  assert(!must_resize());
  ++m_number_of_elements;

  const std::size_t key_index = get_index(key);
  auto it = begin() + key_index;

  // Find the first empty position.
  while (*it != nullptr)
  {
    ++it;
    if (it == end())
    {
      it = begin();
    }

    assert(it != begin() + key_index);
  }

  // Found an empty spot, insert a new index belonging to key,
  *it = key;
  return std::make_pair(it, true);
}


template <class Key, typename Hash, typename Equals, typename Allocator>
inline typename hashtable<Key,Hash,Equals,Allocator>::iterator hashtable<Key,Hash,Equals,Allocator>::erase(const Key& key)
{
  const std::size_t key_index = get_index(key);
  auto it = begin() + key_index;

  // Find the key.
  while (!m_equals(*it, key))
  {
    ++it;
    if (it == end())
    {
      it = begin();
    }

    if (it == begin() + key_index)
    {
      // An element not in the hashset is begin removed.
      // When optimizing, the gcc compiler tends to generate
      // destructions of non generated aterms. If this is
      // repaired, this safety escape can be removed. 
      return it; 
    }
    assert(it != begin() + key_index);
  }

  *it = nullptr;
  --m_number_of_elements;
  return it;
}


template <class Key, typename Hash, typename Equals, typename Allocator>
inline typename hashtable<Key,Hash,Equals,Allocator>::iterator hashtable<Key,Hash,Equals,Allocator>::find(const Key& key)
{
  for (auto it = begin(); it != end(); ++it)
  {
    if (*it == key)
    {
      return it;
    }
  }

  return end();
}

// PRIVATE FUNCTIONS

template <class Key, typename Hash, typename Equals, typename Allocator>
inline std::size_t hashtable<Key,Hash,Equals,Allocator>::get_index(const Key& key)
{
  return m_hasher(key) & m_buckets_mask;
}


} // namespace mcrl2::utilities



#endif // MCRL2_UTILITIES_DETAIL_INDEXED_SET_H
