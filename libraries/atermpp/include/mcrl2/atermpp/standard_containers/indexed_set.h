// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_ATERMPP_INDEXED_SET_H
#define MCRL2_ATERMPP_INDEXED_SET_H

#include "mcrl2/atermpp/detail/thread_aterm_pool.h"
#include "mcrl2/atermpp/standard_containers/deque.h"
#include "mcrl2/utilities/detail/container_utility.h"
#include "mcrl2/utilities/indexed_set.h"
#include "mcrl2/utilities/shared_mutex.h"


namespace atermpp
{

/// \brief A set that assigns each element an unique index, and protects its internal terms en masse.
/// \details The API of a mcrl2::utilities::indexed_set is exposed by explicit forwarding
///          (composition) instead of inheritance, such that every operation that touches
///          aterms is guaranteed to acquire the shared lock that keeps the garbage
///          collector out of the underlying key table while it is being accessed.
template<typename Key,
         bool ThreadSafe = false,
         typename Hash = std::hash<Key>,
         typename Equals = std::equal_to<Key>,
         typename Allocator = std::allocator<Key>,
         typename KeyTable = atermpp::deque<Key > >
class indexed_set
{
protected:
  using super = mcrl2::utilities::indexed_set<Key, ThreadSafe, Hash, Equals, Allocator, KeyTable>;

  super m_set;

public:
  using key_type = typename super::key_type;
  using size_type = typename super::size_type;
  using value_type = typename super::value_type;
  using hasher = typename super::hasher;
  using key_equal = typename super::key_equal;
  using iterator = typename super::iterator;
  using const_iterator = typename super::const_iterator;
  using reverse_iterator = typename super::reverse_iterator;
  using const_reverse_iterator = typename super::const_reverse_iterator;

  /// \brief Value returned when an element does not exist in the set.
  static constexpr size_type npos = super::npos;

  /// \brief Constructor of an empty indexed set. Starts with a hashtable of size 128.
  indexed_set() = default;

  /// \brief Constructor of an empty indexed set. Starts with a hashtable of size 128.
  indexed_set(std::size_t number_of_threads)
    : m_set(number_of_threads)
  {}

  /// \brief Constructor of an empty index set. Starts with a hashtable of the indicated size. 
  /// \param initial_hashtable_size The initial size of the hashtable.
  /// \param hash The hash function.
  /// \param equals The comparison function for its elements.
  indexed_set(std::size_t number_of_threads,
              std::size_t initial_hashtable_size,
              const typename super::hasher& hash = typename super::hasher(),
              const typename super::key_equal& equals = typename super::key_equal()) 
    : m_set(number_of_threads, initial_hashtable_size, hash, equals)
  {}
  
  /// \brief Clears the indexed set.
  void clear(std::size_t thread_index=0)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_set.clear(thread_index);
  }

  /// \brief Insert a key in the indexed set and return its index.
  /// \returns The index of the key and whether the element was newly inserted.
  [[nodiscard]] std::pair<size_type, bool> insert(const Key& key, std::size_t thread_index=0)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_set.insert(key, thread_index);
  }

  /// \returns The index of the given key, or npos when the key is not present.
  [[nodiscard]] size_type index(const key_type& key, std::size_t thread_index = 0) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_set.index(key, thread_index);
  }

  /// \returns Whether the given key is present.
  /// \details Performs the lookup under a single guard such that a concurrent insert
  ///          cannot invalidate the answer halfway through the check.
  [[nodiscard]] bool contains(const key_type& key, std::size_t thread_index = 0) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_set.index(key, thread_index) != npos;
  }

  /// \returns An iterator to the given key, or end() when the key is not present.
  [[nodiscard]] const_iterator find(const key_type& key, std::size_t thread_index = 0) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_set.find(key, thread_index);
  }

  /// \returns The key at the given index; throws std::out_of_range when absent.
  [[nodiscard]] const key_type& at(size_type index) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_set.at(index);
  }

  /// \returns The key at the given index.
  [[nodiscard]] const key_type& operator[](size_type index) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_set[index];
  }

  /// \returns The number of elements in the indexed set.
  [[nodiscard]] size_type size(std::size_t thread_index = 0) const
  {
    return m_set.size(thread_index);
  }

  iterator begin(std::size_t thread_index = 0) { return m_set.begin(thread_index); }
  iterator end(std::size_t thread_index = 0) { return m_set.end(thread_index); }

  const_iterator begin(std::size_t thread_index = 0) const { return m_set.begin(thread_index); }
  const_iterator end(std::size_t thread_index = 0) const { return m_set.end(thread_index); }

  const_iterator cbegin(std::size_t thread_index = 0) const { return m_set.cbegin(thread_index); }
  const_iterator cend(std::size_t thread_index = 0) const { return m_set.cend(thread_index); }

  reverse_iterator rbegin(std::size_t thread_index = 0) { return m_set.rbegin(thread_index); }
  reverse_iterator rend(std::size_t thread_index = 0) { return m_set.rend(thread_index); }

  const_reverse_iterator crbegin(std::size_t thread_index = 0) const { return m_set.crbegin(thread_index); }
  const_reverse_iterator crend(std::size_t thread_index = 0) const { return m_set.crend(thread_index); }
};

} // end namespace atermpp

namespace mcrl2::utilities::detail
{

// Specialization of a function defined in mcrl2/utilities/detail/container_utility.h.
// In utilities, atermpp is not known. 
template<typename Key,
         bool ThreadSafe,
         typename Hash,
         typename Equals,
         typename Allocator,
         typename KeyTable>
bool contains(const atermpp::indexed_set<Key, ThreadSafe, Hash, Equals, Allocator, KeyTable>& c, 
              const typename atermpp::indexed_set<Key, ThreadSafe, Hash, Equals, Allocator, KeyTable>::key_type& v,
              const std::size_t thread_index=0)
{
  return c.contains(v, thread_index);
}

} // namespace mcrl2::utilities::detail

#endif // MCRL2_ATERMPP_INDEXED_SET_H
