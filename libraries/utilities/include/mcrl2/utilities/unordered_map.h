// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_UNORDERED_MAP_H
#define MCRL2_UTILITIES_UNORDERED_MAP_H

#include "mcrl2/utilities/unordered_set.h"

#include <utility>
#include <functional>

namespace mcrl2::utilities
{

/// \brief A class for a map of keys to values in T based using the simple hash table set implementation.
template<typename Key,
         typename T,
         typename Hash,
         typename KeyEqual,
         typename Allocator,
         bool ThreadSafe>
class unordered_map
{
public:
  using key_type = Key;
  using mapped_type = T;
  using value_type = std::pair<const Key, T>;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  using hasher = Hash;
  using key_equal = KeyEqual;
  using allocator_type = typename std::allocator_traits<Allocator>::template rebind_alloc<value_type>;

  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = typename std::allocator_traits<Allocator>::pointer;
  using const_pointer = typename std::allocator_traits<Allocator>::const_pointer;

private:

  // Hashes only the keys of each pair.
  struct PairHash
  {
    using is_transparent = void;

    PairHash() = default;

    explicit PairHash(const hasher& hash)
      : hash(hash)
    {}

    hasher hash;

    std::size_t operator()(const value_type& pair) const
    {
      return hash(pair.first);
    }

    std::size_t operator()(const key_type& key, const mapped_type&) const
    {
      return hash(key);
    }

    std::size_t operator()(const key_type& key) const
    {
      return hash(key);
    }
  };

  // Compares only the keys of each pair.
  struct PairEquals
  {
    using is_transparent = void;

    PairEquals() = default;

    explicit PairEquals(const key_equal& equals)
      : equals(equals)
    {}

    key_equal equals;

    bool operator()(const value_type& first, const value_type& second) const
    {
      return equals(first.first, second.first);
    }

    bool operator()(const value_type& first, const key_type& key, const mapped_type&) const
    {
      return equals(first.first, key);
    }

    bool operator()(const value_type& first, const key_type& key) const
    {
      return equals(first.first, key);
    }
  };

  using Set = unordered_set<value_type, PairHash, PairEquals, allocator_type, ThreadSafe>;
  using bucket_type = typename Set::bucket_type;

  Set m_set; ///< The underlying set storing <key, value> pairs.

public:
  using iterator = typename Set::template unordered_set_iterator<bucket_type, false>;
  using const_iterator = typename Set::const_iterator;
  using local_iterator = typename bucket_type::iterator;
  using const_local_iterator = typename Set::const_local_iterator;

  unordered_map()
    : m_set(0, PairHash(hasher()), PairEquals(key_equal()))
  {}

  /// \brief Constructs an unordered_map that can store initial_size number of elements before resizing.
  unordered_map(std::size_t initial_size,
    const hasher& hash = hasher(),
    const key_equal& equals = key_equal())
    : m_set(initial_size, PairHash(hash), PairEquals(equals))
  {}

  /// \returns A reference to the local node allocator.
  const allocator_type& get_allocator() const noexcept { return m_set.get_allocator(); }
  allocator_type& get_allocator() noexcept { return m_set.get_allocator(); }

  iterator begin() { return m_set.begin(); }
  iterator end() { return m_set.end(); }

  const_iterator begin() const { return m_set.begin(); }
  const_iterator end() const { return m_set.end(); }

  const_iterator cbegin() const { return m_set.begin(); }
  const_iterator cend() const { return m_set.end(); }

  /// \returns True iff the set is empty.
  bool empty() const noexcept { return m_set.empty(); }

  /// \returns The number of elements.
  size_type size() const { return m_set.size(); }
  size_type max_size() const noexcept { return m_set.max_size(); }

  /// \brief Clears the content.
  void clear() { m_set.clear(); }

  /// \brief Inserts elements.
  std::pair<iterator, bool> insert(const value_type& pair) { auto[x, y] = m_set.emplace(pair); return std::make_pair(iterator(x), y); }

  template<typename ...Args>
  std::pair<iterator, bool> emplace(Args&&... args) { auto[x, y] = m_set.emplace(std::forward<Args>(args)...); return std::make_pair(iterator(x), y); }

  template<typename ...Args>
  std::pair<iterator, bool> try_emplace(const key_type& key, Args&&... args);

  /// \brief Erases elements.
  void erase(const key_type& key) { const_iterator it = m_set.find(key); m_set.erase(it); }
  iterator erase(const_iterator it) { return m_set.erase(it); }

  /// \brief Provides access to the value associated with the given key.
  const T& at(const key_type& key) const;

  /// \brief Provides access to the value associated with the given key, constructs a default
  ///        value whenever the key was undefined.
  mapped_type& operator[](const key_type& key);

  /// \returns The number of elements matching the specified key.
  size_type count(const key_type& key) const { return m_set.count(key); }

  /// \returns Element with the specified key.
  template<typename ...Args>
  iterator find(const Args&... args) { return m_set.find(args...); }

  template<typename ...Args>
  const_iterator find(const Args&... args) const { return m_set.find(args...); }

  /// \returns An iterator to the elements in the given bucket with index n.
  local_iterator begin(size_type n) { return m_set.begin(n); }
  local_iterator end(size_type n) { return m_set.end(n); }

  const_local_iterator begin(size_type n) const { return m_set.begin(n); }
  const_local_iterator end(size_type n) const { return m_set.end(n); }

  const_local_iterator cbegin(size_type n) const { return m_set.begin(n); }
  const_local_iterator cend(size_type n) const { return m_set.end(n); }

  /// \returns The number of buckets.
  size_type bucket_count() const noexcept { return m_set.bucket_count(); }
  size_type max_bucket_count() const noexcept { return m_set.max_bucket_count(); }

  size_type bucket_size(size_type n) const noexcept { return m_set.bucket_size(n); }
  size_type bucket(const key_type& key) const noexcept { return m_set.bucket(key); }

  float load_factor() const { return static_cast<float>(size()) / bucket_count(); }
  float max_load_factor() const { return m_set.max_load_factor(); }
  void max_load_factor(float factor) { m_set.max_load_factor(factor); }

  /// \brief Resize the number buckets to at least number_of_buckets.
  void rehash(size_type number_of_buckets);

  /// \brief Resizes the set to the given number of elements.
  void reserve(size_type count) { rehash(std::ceil(static_cast<float>(count) / max_load_factor())); }

  hasher hash_function() const { return m_set.hash_function(); }
  key_equal key_eq() const { return m_set.key_eq(); }

  /// \brief Number of elements that can be stored before rehash.
  /// \details Nonstandard.
  size_type capacity() { return m_set.capacity(); }
};

/// \brief A specialization for large unordered maps that uses the block_allocator internally by default.
template<typename Key,
         typename T,
         typename Hash = std::hash<Key>,
         typename Equals = std::equal_to<Key>,
         typename Allocator = mcrl2::utilities::block_allocator<Key>,
         bool ThreadSafe = false>
using unordered_map_large = unordered_map<Key, T, Hash, Equals, Allocator, ThreadSafe>;

} // namespace mcrl2::utilities

#include "mcrl2/utilities/detail/unordered_map_implementation.h"

#endif // MCRL2_UTILITIES_UNORDERED_MAP_H
