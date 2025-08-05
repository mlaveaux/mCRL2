// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef MCRL2_UTILITIES_UNORDERED_MAP_IMPLEMENTATION_H
#define MCRL2_UTILITIES_UNORDERED_MAP_IMPLEMENTATION_H
#pragma once

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/unordered_map.cxx" // necessary for header test. 
#else
  // Implementation unit
  module utilities:unordered_map;
#endif

#define MCRL2_UNORDERED_MAP_TEMPLATES template<typename Key, typename T, typename Hash, typename Equals, typename Allocator, bool ThreadSafe, bool Resize>
#define MCRL2_UNORDERED_MAP_CLASS unordered_map<Key, T, Hash, Equals, Allocator, ThreadSafe, Resize>

namespace mcrl2::utilities
{

MCRL2_UNORDERED_MAP_TEMPLATES
template<typename ...Args>
auto MCRL2_UNORDERED_MAP_CLASS::try_emplace(const key_type& key, Args&&... args) -> std::pair<iterator, bool>
{
  m_set.rehash_if_needed();

  std::size_t bucket = m_set.find_bucket_index(key);
  // If the key can be found, then return it and otherwise add it in the same bucket.
  iterator it = m_set.find_impl(bucket, key);
  if (it != end())
  {
    return std::make_pair(it, false);
  }
  else
  {
    auto [x, y] = m_set.emplace_impl(bucket, key, std::forward<Args>(args)...);
    return std::make_pair(iterator(x), y);
  }
}

MCRL2_UNORDERED_MAP_TEMPLATES
T& MCRL2_UNORDERED_MAP_CLASS::operator[](const key_type& key)
{
  // Insert a new object and return a reference to it;
  // auto pair = m_set.emplace(std::make_pair<const Key, T>(key, mapped_type()));
  auto pair = m_set.emplace(std::make_pair(key, mapped_type()));

  // The second element of the pair is mutable.
  return const_cast<mapped_type&>((*pair.first).second);
}

MCRL2_UNORDERED_MAP_TEMPLATES
const T& MCRL2_UNORDERED_MAP_CLASS::at(const key_type& key) const
{
  auto it = m_set.find(key);
  assert(it != m_set.end());
  return (*it).second;
}

#undef MCRL2_UNORDERED_MAP_TEMPLATES
#undef MCRL2_UNORDERED_MAP_CLASS

} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_UNORDERED_MAP_IMPLEMENTATION_H
