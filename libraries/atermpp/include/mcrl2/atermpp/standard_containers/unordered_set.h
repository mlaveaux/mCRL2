// Author(s): Jan Friso Groote, Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//
/// \file mcrl2/atermpp/standard_containers/unordered_set.h
/// \brief This file contains an unordered_set class that behaves
///        exactly as a standard unordered_set. It can only be used
///        to store class instances that derive from aterms.
///        The stored aterms are protected as a whole, i.e.,
///        time and memory is saved as individual protection
///        per element is unnecessary.

#ifndef MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_SET_H
#define MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_SET_H
#pragma once

#include <concepts>
#include <unordered_set>
#include "mcrl2/atermpp/detail/aterm_container.h"
#include "mcrl2/utilities/shared_mutex.h"

/// \brief The main namespace for the aterm++ library.
namespace atermpp
{

/// \brief An unordered_set class in which aterms can be stored.
/// \details The API of a std::unordered_set is exposed by explicit forwarding (composition)
///          instead of inheritance, such that every operation is guaranteed to acquire the
///          shared lock that keeps the garbage collector out of the container while it is used.
template < class Key,
           class Hash = std::hash<detail::reference_aterm<Key>>,
           class Pred = std::equal_to<detail::reference_aterm<Key>>,
           class Alloc = std::allocator<detail::reference_aterm<Key>> >
class unordered_set
{
protected:
  using super = std::unordered_set<detail::reference_aterm<Key>, Hash, Pred, Alloc>;

  super m_container;
  detail::generic_aterm_container<super> container_wrapper;

public:

  /// Standard typedefs.
  using key_type = typename super::key_type;
  using allocator_type = typename super::allocator_type;
  using value_type = typename super::value_type;
  using size_type = typename super::size_type;
  using node_type = typename super::node_type;
  using reference = typename super::reference;
  using const_reference = typename super::const_reference;
  using iterator = typename super::iterator;
  using const_iterator = typename super::const_iterator;
  using insert_return_type = typename super::insert_return_type;
  using hasher = typename super::hasher;
  using key_equal = typename super::key_equal;

  /// \brief Default constructor.
  unordered_set()
    : m_container(),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  explicit unordered_set(const allocator_type& alloc)
    : m_container(alloc),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  explicit unordered_set(size_type n, const allocator_type& alloc = allocator_type())
    : m_container(n, alloc),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  template <class InputIterator>
  unordered_set(InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type())
    : m_container(first, last, alloc),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  unordered_set(const unordered_set& x)
    : m_container(x.m_container),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  unordered_set(const unordered_set& x, const allocator_type& alloc)
    : m_container(x.m_container, alloc),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  unordered_set(unordered_set&& x) noexcept
    : m_container(std::move(x.m_container)),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  unordered_set(unordered_set&& x, const allocator_type& alloc)
    : m_container(std::move(x.m_container), alloc),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  unordered_set(std::initializer_list<value_type> il, const allocator_type& alloc = allocator_type())
    : m_container(il, alloc),
      container_wrapper(m_container)
  {}

  /// \brief Standard assignment.
  unordered_set& operator=(const unordered_set& other);

  /// \brief Standard move assignment.
  unordered_set& operator=(unordered_set&& other) noexcept;

  /// \brief Standard destructor.
  ~unordered_set() = default;

  /// \returns The allocator of the underlying container.
  [[nodiscard]] allocator_type get_allocator() const noexcept
  {
    return m_container.get_allocator();
  }

  void clear() noexcept;

  std::pair<iterator, bool> insert(const value_type& value);

  template<class P>
    requires std::constructible_from<detail::reference_aterm<Key>, P>
  std::pair<iterator, bool> insert(P&& value);

  iterator insert(const_iterator hint, const value_type& value);

  template<class P>
    requires std::constructible_from<detail::reference_aterm<Key>, P>
  iterator insert(const_iterator hint, P&& value);

  template<class InputIt>
  void insert(InputIt first, InputIt last);

  void insert(std::initializer_list<value_type> ilist);

  insert_return_type insert(node_type&& nh);

  iterator insert(const_iterator hint, node_type&& nh);

  template<class... Args>
  std::pair<iterator, bool> emplace(Args&&... args);

  template<class... Args>
  iterator emplace_hint(const_iterator hint, Args&&... args);

  iterator erase(iterator pos);

  // iterator erase(const_iterator pos); TODO: Enable this by making the iterator and const_iterator of aterms different. 
  //                                           A const_iterator now behaves as an iterator. 

  iterator erase(const_iterator first, const_iterator last);

  size_type erase(const Key& key);

  void swap(unordered_set& other) noexcept;

  node_type extract(const_iterator position);

  node_type extract(const Key& k);

  /// \returns Element with the specified key, or end() otherwise.
  [[nodiscard]] iterator find(const Key& key);
  [[nodiscard]] const_iterator find(const Key& key) const;

  /// \returns The number of elements matching the specified key.
  [[nodiscard]] size_type count(const Key& key) const;

  /// \returns Whether an element with the specified key is stored.
  [[nodiscard]] bool contains(const Key& key) const;

  [[nodiscard]] std::size_t size() const
  {
    // Concurrent read/write on the size.
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.size();
  }

  /// \returns An iterator over all keys.
  iterator begin();
  iterator end();

  /// \returns A const iterator over all keys.
  const_iterator begin() const;
  const_iterator end() const;

  /// \returns A const iterator over all keys.
  const_iterator cbegin() const;
  const_iterator cend() const;

  /// \returns True iff the set is empty.
  [[nodiscard]] bool empty() const noexcept;

  /// \returns The amount of elements stored in this set.
  [[nodiscard]] size_type max_size() const noexcept;

  void rehash(size_type count);

  void reserve(size_type count);

  [[nodiscard]] hasher hash_function() const { return m_container.hash_function(); }
  [[nodiscard]] key_equal key_eq() const { return m_container.key_eq(); }

  /// \brief Marks all stored terms during garbage collection; used when this
  ///        container is an element of another protected container.
  void mark(term_mark_stack& todo) const
  {
    for (const detail::reference_aterm<Key>& element : m_container)
    {
      element.mark(todo);
    }
  }

  bool operator==(const unordered_set& other) const { return m_container == other.m_container; }
};

namespace utilities
{

/// \brief An unordered_set class in which aterms can be stored.
/// \details The API of a mcrl2::utilities::unordered_set is exposed by explicit forwarding
///          (composition) instead of inheritance, such that every operation is guaranteed to
///          acquire the necessary locks.
template < class Key,
           class Hash = std::hash<detail::reference_aterm<Key>>,
           class Pred = std::equal_to<detail::reference_aterm<Key>>,
           class Alloc = std::allocator<detail::reference_aterm<Key>>,
           bool ThreadSafe = false >
class unordered_set
{
protected:
  using super = mcrl2::utilities::unordered_set<detail::reference_aterm<Key>, Hash, Pred, Alloc, ThreadSafe, false>;

  super m_container;
  detail::generic_aterm_container<super> container_wrapper;

public:

  /// Standard typedefs.
  using allocator_type = typename super::allocator_type;
  using value_type = typename super::value_type;
  using size_type = typename super::size_type;
  using reference = typename super::reference;
  using const_reference = typename super::const_reference;
  using iterator = typename super::iterator;
  using const_iterator = typename super::const_iterator;

  /// \brief Default constructor.
  unordered_set()
    : m_container(),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  explicit unordered_set(size_type n)
    : m_container(n),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  unordered_set(const unordered_set& x)
    : m_container(x.m_container),
      container_wrapper(m_container)
  {}

  /// \brief Constructor.
  unordered_set(unordered_set&& x) noexcept
    : m_container(std::move(x.m_container)),
      container_wrapper(m_container)
  {}

  /// \brief Standard assignment.
  unordered_set& operator=(const unordered_set& other);

  /// \brief Standard move assignment.
  unordered_set& operator=(unordered_set&& other) noexcept;

  /// \brief Standard destructor.
  ~unordered_set() = default;

  void clear() noexcept;

  /// \brief Standard find function in a set.
  /// \returns Element with the specified key.
  template<typename... Args>
  [[nodiscard]] iterator find(const Args&... args);

  /// \brief Standard find function in a set.
  /// \returns Element with the specified key.
  template<typename... Args>
  [[nodiscard]] const_iterator find(const Args&... args) const;

  /// \returns The number of elements matching the specified key.
  template<typename... Args>
  [[nodiscard]] size_type count(const Args&... args) const;

  /// \returns Whether an element with the specified key is stored.
  template<typename... Args>
  [[nodiscard]] bool contains(const Args&... args) const;

  template<typename... Args>
  std::pair<iterator, bool> emplace(Args&&... args);

  template<typename... Args>
  void erase(const Args&... args);

  iterator erase(const_iterator it);

  void swap(unordered_set& other) noexcept;

  [[nodiscard]] std::size_t size() const
  {
    // Concurrent read/write on the size.
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.size();
  }

  /// \returns An iterator over all keys.
  iterator begin();
  iterator end();

  /// \returns A const iterator over all keys.
  const_iterator begin() const;
  const_iterator end() const;

  /// \returns A const iterator over all keys.
  const_iterator cbegin() const;
  const_iterator cend() const;

  /// \returns True iff the set is empty.
  [[nodiscard]] bool empty() const noexcept;

  /// \returns The amount of elements stored in this set.
  [[nodiscard]] size_type max_size() const noexcept;

  /// \brief Marks all stored terms during garbage collection; used when this
  ///        container is an element of another protected container.
  void mark(term_mark_stack& todo) const
  {
    for (const detail::reference_aterm<Key>& element : m_container)
    {
      element.mark(todo);
    }
  }

protected:

  void rehash(std::size_t new_size);

  /// \brief Rehashes the container when the load factor is exceeded, re-checking the
  ///        load factor under the exclusive lock when ThreadSafe to avoid redundant rehashes.
  void rehash_if_needed();

  bool rehash_is_needed() const;
};

} // namespace utilities

} // namespace atermpp

#include "mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h"

#endif // MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_SET_H
