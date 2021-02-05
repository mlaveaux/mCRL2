// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_UNORDERED_SET_H
#define MCRL2_UTILITIES_UNORDERED_SET_H

#include "mcrl2/utilities/block_allocator.h"
#include "mcrl2/utilities/logger.h"
#include "mcrl2/utilities/detail/bucket_list.h"

#include <cmath>
#include <mutex>

namespace mcrl2::utilities
{

/// \brief Enables lockfree implementation of emplace.
constexpr static bool EnableLockfreeInsertion = true;

/// \brief Number of buckets per mutex.
constexpr static long BucketsPerMutex = 256;

/// \brief Prints various information for unordered_set like data structures.
template<typename T>
void print_performance_statistics(const T& unordered_set);

// Forward declaration of an unordered_map.
template<typename Key,
         typename T,
         typename Hash = std::hash<Key>,
         typename KeyEqual = std::equal_to<Key>,
         typename Allocator = std::allocator<Key>,
         bool ThreadSafe = false>
class unordered_map;

/// \brief A unordered_set with a subset of the interface of std::unordered_set that only stores a single pointer for each element.
/// \details Only supports input iterators (not bidirectional) compared to std::unordered_set. Furthermore, iterating over all elements
///          in the set is O(n + m), where n is the number of elements in the set and m the number of empty buckets. Also incrementing the
///          iterator is O(m) as opposed to O(1) as the standard mandates.
///
///          Additionally, the unordered_set supports allocators that have a specialized allocate_args(args...) to vary the allocation size based
///          on the arguments used. This is required to store _aterm_appl classes with the function symbol arity determined at runtime.
///
///          Threadsafe enables concurrent emplace calls, and Resize enables automatically resizing if needed.
///
/// \todo Does not implement std::unordered_map equal_range and swap.
template<typename Key,
         typename Hash = std::hash<Key>,
         typename Equals = std::equal_to<Key>,
         typename Allocator = std::allocator<Key>,
         bool ThreadSafe = false,
         bool Resize = true>
class unordered_set
{
  static_assert (!(ThreadSafe && Resize), "ThreadSafe cannot be enabled together with automatic resizing.");

private:
  /// \brief Combine the bucket list and a lock that locks modifications to the bucket list.
  using bucket_type = detail::bucket_list<Key, Allocator>;
  using bucket_iterator = typename std::vector<bucket_type>::iterator;
  using const_bucket_iterator = typename std::vector<bucket_type>::const_iterator;
  using mutex_type = std::mutex;

  template<typename Key_, typename T, typename Hash_, typename KeyEqual, typename Allocator_, bool ThreadSafe_>
  friend class unordered_map;

public:
  /// \brief An iterator over all elements in the unordered set.
  template<typename Bucket, bool Constant>
  class unordered_set_iterator : std::iterator_traits<Key>
  {
  private:
    friend class unordered_set;

    template<typename Key_, typename T, typename Hash_, typename KeyEqual, typename Allocator_, bool ThreadSafe_>
    friend class unordered_map;

    using bucket_it = typename std::vector<Bucket>::const_iterator;
    using key_it_type = typename Bucket::const_iterator;

  public:
    using value_type = Key;
    using reference = typename std::conditional<Constant, const Key&, Key&>::type;
    using pointer = typename std::conditional<Constant, const Key*, Key*>::type;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    unordered_set_iterator() = default;

    operator unordered_set_iterator<Bucket, true>() const
    {
      return unordered_set_iterator<Bucket, true>(m_bucket_it, m_bucket_end, m_key_before_it, m_key_it);
    }

    unordered_set_iterator& operator++()
    {
      ++m_key_before_it;
      ++m_key_it;
      goto_next_bucket();
      return *this;
    }

    unordered_set_iterator operator++(int)
    {
      unordered_set_iterator copy(*this);
      ++(*this);
      return *copy;
    }

    reference operator*() const
    {
      return const_cast<reference>(*m_key_it);
    }

    pointer operator->() const
    {
      return const_cast<pointer>(&(*m_key_it));
    }

    bool operator!=(const unordered_set_iterator& other) const
    {
      return m_key_it != other.m_key_it || m_bucket_it != other.m_bucket_it;
    }

    bool operator==(const unordered_set_iterator& other) const
    {
      return !(*this != other);
    }

  private:
    /// \brief Construct an iterator over all keys passed in this bucket and all remaining buckets.
    unordered_set_iterator(bucket_it it, bucket_it end, key_it_type before_it, key_it_type key) :
      m_bucket_it(it), m_bucket_end(end), m_key_before_it(before_it), m_key_it(key)
    {}

    /// \brief Construct the begin iterator (over all elements).
    unordered_set_iterator(bucket_it it, bucket_it end) :
      m_bucket_it(it), m_bucket_end(end), m_key_before_it((*it).before_begin()), m_key_it((*it).begin())
    {
      goto_next_bucket();
    }

    /// \brief Construct the end iterator
    explicit unordered_set_iterator(bucket_it it) :
      m_bucket_it(it)
    {}

    operator unordered_set_iterator<Bucket, false>() const
    {
      return unordered_set_iterator<Bucket, false>(m_bucket_it, m_bucket_end, m_key_before_it, m_key_it);
    }

    /// \returns A reference to the before key iterator.
    key_it_type& key_before_it() { return m_key_before_it; }

    /// \returns A reference to the key iterator.
    key_it_type& key_it() { return m_key_it; }

    /// \returns A reference to the bucket iterator.
    bucket_it& get_bucket_it() { return m_bucket_it; }

    /// \brief Iterate to the next non-empty bucket.
    void goto_next_bucket()
    {
      // Find the first bucket that is not empty.
      while(!(m_key_it != detail::EndIterator))
      {
        // Take the next bucket and reset the key iterator.
        ++m_bucket_it;

        if (m_bucket_it != m_bucket_end)
        {
          m_key_it = (*m_bucket_it).begin();
          m_key_before_it = (*m_bucket_it).before_begin();
        }
        else
        {
          // Reached the end of the buckets.
          break;
        }
      }

      // The current bucket contains elements or we are at the end.
      assert(m_bucket_it == m_bucket_end || m_key_it != detail::EndIterator);
    }

    bucket_it m_bucket_it;
    bucket_it m_bucket_end;
    key_it_type m_key_before_it;
    key_it_type m_key_it; // Invariant: m_key_it != EndIterator.
  };

  using key_type = Key;
  using value_type = Key;
  using hasher = Hash;
  using key_equal = Equals;
  using allocator_type = typename bucket_type::NodeAllocator;

  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = typename std::allocator_traits<Allocator>::pointer;
  using const_pointer = typename std::allocator_traits<Allocator>::const_pointer;

  using const_local_iterator = typename bucket_type::const_iterator;
  using local_iterator = typename bucket_type::const_iterator;
  using const_iterator = unordered_set_iterator<bucket_type, true>;
  using iterator = const_iterator;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  unordered_set() { rehash(0); }

  /// \brief Constructs an unordered_set that contains bucket_count number of buckets.
  explicit unordered_set(size_type bucket_count,
    const hasher& hash = hasher(),
    const key_equal& equals = key_equal())
    : m_hash(hash),
      m_equals(equals)
  {
    rehash(bucket_count);
  }

  // Copy operators.
  unordered_set(const unordered_set& set);
  unordered_set& operator=(const unordered_set& set);

  // Default move operators.
  unordered_set(unordered_set&& other) = default;
  unordered_set& operator=(unordered_set&& other) = default;

  ~unordered_set();

  /// \returns A reference to the local node allocator.
  const allocator_type& get_allocator() const noexcept { return m_allocator; }
  allocator_type& get_allocator() noexcept { return m_allocator; }

  /// \returns An iterator over all keys.
  iterator begin() { return iterator(m_buckets.begin(), m_buckets.end()); }
  iterator end() { return iterator(m_buckets.end()); }

  /// \returns A const iterator over all keys.
  const_iterator begin() const { return const_iterator(m_buckets.begin(), m_buckets.end()); }
  const_iterator end() const { return const_iterator(m_buckets.end()); }

  /// \returns A const iterator over all keys.
  const_iterator cbegin() const { return const_iterator(m_buckets.begin(), m_buckets.end()); }
  const_iterator cend() const { return const_iterator(m_buckets.end()); }

  /// \returns True iff the set is empty.
  bool empty() const noexcept { return size() == 0; }

  /// \returns The amount of elements stored in this set.
  size_type size() const noexcept { return m_number_of_elements; }
  size_type max_size() const noexcept { return m_buckets.max_size(); }

  /// \brief Removes all elements from the set.
  /// \details Does not free the vector of buckets itself.
  void clear();

  /// \brief Inserts an element Key(args...) into the set if it did not already exist.
  /// \returns A pair of the iterator pointing to the element and a boolean that is true iff
  ///         a new element was inserted (as opposed to it already existing in the set).
  /// \threadsafe
  template<typename ...Args>
  std::pair<iterator, bool> emplace(Args&&... args);

  /// \brief Erases the given key_type(args...) from the unordered set.
  template<typename...Args>
  void erase(const Args&... args);

  /// \brief Erases the element pointed to by the iterator.
  /// \returns An iterator to the next key.
  iterator erase(const_iterator it);

  /// \brief Counts the number of occurrences of the given key (1 when it exists and 0 otherwise).
  template<typename ...Args>
  size_type count(const Args&... args) const;

  /// \brief Searches whether an object key_type(args...) occurs in the set.
  /// \returns An iterator to the matching element or the end when this object does not exist.
  template<typename...Args>
  const_iterator find(const Args&... args) const;

  template<typename...Args>
  iterator find(const Args&... args);

  /// \returns An iterator to the elements in the given bucket with index n.
  local_iterator begin(size_type n) { return m_buckets[n].begin(); }
  local_iterator end(size_type n) { return m_buckets[n].end(); }

  const_local_iterator begin(size_type n) const { return m_buckets[n].begin(); }
  const_local_iterator end(size_type n) const { return m_buckets[n].end(); }

  const_local_iterator cbegin(size_type n) const { return m_buckets[n].begin(); }
  const_local_iterator cend(size_type n) const { return m_buckets[n].end(); }

  /// \returns The number of buckets.
  size_type bucket_count() const noexcept { return m_buckets.size(); }
  size_type max_bucket_count() const noexcept { return m_buckets.max_size(); }

  size_type bucket_size(size_type n) const noexcept { return std::distance(m_buckets[n].begin(), m_buckets[n].end()); }
  size_type bucket(const key_type& key) const noexcept { return std::distance(m_buckets.begin(), find_bucket_index(key)); }

  float load_factor() const { return static_cast<float>(size()) / bucket_count(); }
  float max_load_factor() const { return m_max_load_factor; }
  void max_load_factor(float factor) { m_max_load_factor = factor; }

  /// \brief Resize the number buckets to at least number_of_buckets.
  void rehash(size_type number_of_buckets);

  /// \brief Resizes the set to the given number of elements.
  void reserve(size_type count) { rehash(static_cast<std::size_t>(std::ceil(static_cast<float>(count) / max_load_factor()))); }

  hasher hash_function() const { return m_hash; }
  key_equal key_eq() const { return m_equals; }

  /// \returns The number of elements that can be present in the set before resizing.
  /// \details Not standard.
  size_type capacity() const noexcept { return m_buckets.size(); }

  /// \brief Resizes the hash table if necessary.
  /// \details Not standard.
  void rehash_if_needed();

private:
  template<typename Key_, typename T, typename Hash_, typename KeyEqual, typename Allocator_, bool ThreadSafe_>
  friend class unordered_map;

  // Check for the existence of the is_transparent type.
  template <typename... >
  using void_t = void;

  template <typename X, typename = void>
  struct is_transparent : std::false_type { };

  template <typename X>
  struct is_transparent<X, void_t<typename X::is_transparent>>
  : std::true_type { };

  /// \brief Inserts T(args...) into the given bucket, assumes that it did not exists before.
  /// \threadsafe
  template<typename ...Args>
  std::pair<iterator, bool> emplace_impl(size_type bucket_index, Args&&... args);

  /// \brief Removes T(args...) from the set.
  template<typename ...Args>
  void erase_impl(const Args&... args);

  /// \returns The index of the bucket that might contain the element constructed by the given arguments.
  template<typename ...Args>
  size_type find_bucket_index(const Args&... args) const;

  /// \brief Searches for the element in the given bucket.
  template<typename ...Args>
  const_iterator find_impl(size_type bucket_index, const Args&... args) const;

  /// \brief True iff the hash and equals functions allow transparent lookup,
  static constexpr bool allow_transparent = is_transparent<Hash>() && is_transparent<Equals>();

  /// \brief The number of elements stored in this set.
  std::conditional_t<ThreadSafe, std::atomic<size_type>, size_type> m_number_of_elements = 0;

  /// \brief Always equal to m_buckets.size() - 1.
  size_type m_buckets_mask;

  std::vector<bucket_type> m_buckets;
  std::vector<std::mutex> m_bucket_mutexes;

  float m_max_load_factor = 1.0f;

  hasher m_hash = hasher();
  key_equal m_equals = key_equal();
  allocator_type m_allocator;
};

/// \brief A specialization for large unordered sets that uses the block_allocator internally by default.
template<typename Key,
         typename Hash = std::hash<Key>,
         typename Equals = std::equal_to<Key>,
         typename Allocator = mcrl2::utilities::block_allocator<Key>,
         bool ThreadSafe = false>
using unordered_set_large = unordered_set<Key, Hash, Equals, Allocator, ThreadSafe>;

} // namespace mcrl2::utilities

#include "mcrl2/utilities/detail/unordered_set_implementation.h"

#endif // MCRL2_UTILITIES_UNORDERED_SET_H
