// Author(s): Jan Friso Groote, Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//
/// \file mcrl2/atermpp/standard_containers/unordered_map.h
/// \brief This file contains a unordered_map class that behaves 
///        exactly as a standard unordered_map. It can only be used
///        to store class instances that derive from aterms.
///        The stored aterms are protected as a whole, i.e.,
///        time and memory is saved as individual protection
///        per element is unnecessary. 

#ifndef MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_MAP_H
#define MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_MAP_H
#pragma once

#include <concepts>
#include <unordered_map>
#include "mcrl2/atermpp/detail/aterm_container.h"
#include "mcrl2/utilities/shared_mutex.h"

/// \brief The main namespace for the aterm++ library.
namespace atermpp
{

/// \brief A unordered_map class in which aterms can be stored.
/// \details The API of a std::unordered_map is exposed by explicit forwarding (composition)
///          instead of inheritance, such that every operation is guaranteed to acquire the
///          shared lock that keeps the garbage collector out of the container while it is used.
template < class Key, 
           class T,
           class Hash = std::hash<detail::markable_aterm<Key> >,
           class Pred = std::equal_to<detail::markable_aterm<Key> >,
           class Alloc = std::allocator< std::pair<const detail::markable_aterm<Key>, detail::markable_aterm<T> > > >

class unordered_map
{
protected:
  using super = std::unordered_map<detail::markable_aterm<Key>, detail::markable_aterm<T>, Hash, Pred, Alloc>;

  super m_container;
  detail::generic_aterm_container<super> container_wrapper;

public:
  
  /// Standard typedefs.
  using key_type = typename super::key_type;
  using mapped_type = typename super::mapped_type;
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
  unordered_map()
   : m_container(),
     container_wrapper(m_container)     
  {}

  /// \brief Constructor.
  explicit unordered_map(const allocator_type& alloc)
   : m_container(alloc),
     container_wrapper(m_container)     
  {}

  /// \brief Constructor.
  explicit unordered_map(size_type n, const allocator_type& alloc = allocator_type())
   : m_container(n, alloc),
     container_wrapper(m_container)     
  {}

  /// \brief Constructor.
  template <class InputIterator>
  unordered_map(InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type())
   : m_container(first, last, alloc),
     container_wrapper(m_container)     
  {}
    
  /// \brief Constructor.
  unordered_map(const unordered_map& x)
   : m_container(x.m_container),
     container_wrapper(m_container)     
  {}

  /// \brief Constructor.
  unordered_map(const unordered_map& x, const allocator_type& alloc)
   : m_container(x.m_container, alloc),
     container_wrapper(m_container)     
  {}
  
  /// \brief Constructor.
  unordered_map(unordered_map&& x) noexcept
      : m_container(std::move(x.m_container)),
        container_wrapper(m_container)
  {}


  /// \brief Constructor.
  unordered_map(unordered_map&& x, const allocator_type& alloc)
   : m_container(std::move(x.m_container), alloc),
     container_wrapper(m_container)
  {}

  /// \brief Constructor.
  unordered_map(std::initializer_list<value_type> il, const allocator_type& alloc = allocator_type())
    : m_container(il, alloc),
      container_wrapper(m_container)      
  {}

  /// \brief Standard assignment.
  unordered_map& operator=(const unordered_map& other);

  /// \brief Standard move assignment.
  unordered_map& operator=(unordered_map&& other) noexcept;

  /// \brief Standard destructor.
  ~unordered_map()=default;

  /// \returns The allocator of the underlying container.
  [[nodiscard]] allocator_type get_allocator() const noexcept
  {
    return m_container.get_allocator();
  }

  void clear() noexcept;

  /// \brief Inserts an element referring to a default value in the map. 
  std::pair<iterator,bool> insert( const value_type& value );

  template< class P >
    requires std::constructible_from<typename std::unordered_map<detail::markable_aterm<Key>, detail::markable_aterm<T>, Hash, Pred, Alloc>::value_type, P>
  std::pair<iterator,bool> insert( P&& value );

  iterator insert( const_iterator hint, const value_type& value );

  template< class P >
    requires std::constructible_from<typename std::unordered_map<detail::markable_aterm<Key>, detail::markable_aterm<T>, Hash, Pred, Alloc>::value_type, P>
  iterator insert( const_iterator hint, P&& value );

  template< class InputIt >
  void insert( InputIt first, InputIt last );

  void insert( std::initializer_list<value_type> ilist );

  insert_return_type insert( node_type&& nh );

  iterator insert( const_iterator hint, node_type&& nh );

  template <class M>
  std::pair<iterator, bool> insert_or_assign( const Key& k, M&& obj );

  template <class M>
  std::pair<iterator, bool> insert_or_assign( Key&& k, M&& obj );

  template <class M>
  iterator insert_or_assign( const_iterator hint, const Key& k, M&& obj );

  template <class M>
  iterator insert_or_assign( const_iterator hint, Key&& k, M&& obj );

  template< class... Args >
  std::pair<iterator,bool> emplace( Args&&... args );

  template <class... Args>
  iterator emplace_hint( const_iterator hint, Args&&... args );

  template< class... Args >
  std::pair<iterator, bool> try_emplace( const Key& k, Args&&... args );

  template< class... Args >
  std::pair<iterator, bool> try_emplace( Key&& k, Args&&... args );

  template< class... Args >
  iterator try_emplace( const_iterator hint, const Key& k, Args&&... args );

  template< class... Args >
  iterator try_emplace( const_iterator hint, Key&& k, Args&&... args );

  iterator erase( iterator pos );

  iterator erase( const_iterator pos );
  
  iterator erase( const_iterator first, const_iterator last );

  size_type erase( const Key& key );

  void swap(unordered_map& other) noexcept;

  node_type extract( const_iterator position );
  
  node_type extract( const Key& k );

  /// \brief Provides access to the value associated with the given key,
  ///        inserting a default value when the key is not yet present.
  mapped_type& operator[]( const Key& key );

  /// \brief Provides access to the value associated with the given key.
  [[nodiscard]] mapped_type& at( const Key& key );
  [[nodiscard]] const mapped_type& at( const Key& key ) const;

  /// \returns Element with the specified key, or end() otherwise.
  [[nodiscard]] iterator find( const Key& key );
  [[nodiscard]] const_iterator find( const Key& key ) const;

  /// \returns The number of elements matching the specified key.
  [[nodiscard]] size_type count( const Key& key ) const;

  /// \returns Whether an element with the specified key is stored.
  [[nodiscard]] bool contains( const Key& key ) const;

  [[nodiscard]] std::size_t size() const
  {
    // Concurrent read/write on the size.
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.size();
  }

  /// \returns An iterator over all keys.
  inline
  iterator begin();

  inline
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

  void rehash( size_type count );

  void reserve( size_type count );

  [[nodiscard]] hasher hash_function() const { return m_container.hash_function(); }
  [[nodiscard]] key_equal key_eq() const { return m_container.key_eq(); }

  /// \brief Marks all stored terms during garbage collection; used when this
  ///        container is an element of another protected container.
  void mark(term_mark_stack& todo) const
  {
    for (const auto& element : m_container)
    {
      element.first.mark(todo);
      element.second.mark(todo);
    }
  }

  bool operator==(const unordered_map& other) const { return m_container == other.m_container; }
};

namespace utilities
{
/// \brief A unordered_map class in which aterms can be stored.
/// \details The API of a mcrl2::utilities::unordered_map is exposed by explicit forwarding
///          (composition) instead of inheritance, such that every operation is guaranteed to
///          acquire the necessary locks.
template < class Key,
  class T,
  class Hash = std::hash<detail::markable_aterm<Key> >,
  class Pred = std::equal_to<detail::markable_aterm<Key> >,
  class Alloc = std::allocator< std::pair<const detail::markable_aterm<Key>, detail::markable_aterm<T> > >,
  bool ThreadSafe = false >

class unordered_map
{
  protected:
    using super = mcrl2::utilities::
        unordered_map<detail::markable_aterm<Key>, detail::markable_aterm<T>, Hash, Pred, Alloc, ThreadSafe, false>;

    super m_container;
    detail::generic_aterm_container<super> container_wrapper;

  public:

    /// Standard typedefs.
    using key_type = typename super::key_type;
    using mapped_type = typename super::mapped_type;
    using allocator_type = typename super::allocator_type;
    using value_type = typename super::value_type;
    using size_type = typename super::size_type;
    using reference = typename super::reference;
    using const_reference = typename super::const_reference;
    using iterator = typename super::iterator;
    using const_iterator = typename super::const_iterator;

    /// \brief Default constructor.
    unordered_map()
      : m_container(),
      container_wrapper(m_container)
    {}

    /// \brief Constructor.
    explicit unordered_map(size_type n)
      : m_container(n),
      container_wrapper(m_container)
    {}

    /// \brief Constructor.
    unordered_map(const unordered_map& x)
      : m_container(x.m_container),
      container_wrapper(m_container)
    {}

    /// \brief Constructor.
    unordered_map(unordered_map&& x) noexcept
        : m_container(std::move(x.m_container)),
          container_wrapper(m_container)
    {}

    /// \brief Standard assignment.
    unordered_map& operator=(const unordered_map& other);

    /// \brief Standard move assignment.
    unordered_map& operator=(unordered_map&& other) noexcept;

    /// \brief Standard destructor.
    ~unordered_map() = default;

    void clear() noexcept;

    /// \brief Standard find function in a map.
    /// \returns Element with the specified key.
    template<typename ...Args>
    iterator find(const Args&... args);

    /// \brief Standard find function in a map.
    /// \returns Element with the specified key.
    template<typename ...Args>
    const_iterator find(const Args&... args) const;

    std::pair<iterator, bool> insert(const value_type& value);

    template< class P >
    std::pair<iterator, bool> insert(P&& value);

    iterator insert(const_iterator hint, value_type&& value);

    template< class P >
    iterator insert(const_iterator hint, P&& value);

    template< class InputIt >
    void insert(InputIt first, InputIt last);

    void insert(std::initializer_list<value_type> ilist);

    /*
    insert_return_type insert(node_type&& nh);

    iterator insert(const_iterator hint, node_type&& nh);
    */
    
    template <class M>
    std::pair<iterator, bool> insert_or_assign(const Key& k, M&& obj);

    template <class M>
    std::pair<iterator, bool> insert_or_assign(Key&& k, M&& obj);

    template <class M>
    iterator insert_or_assign(const_iterator hint, const Key& k, M&& obj);

    template <class M>
    iterator insert_or_assign(const_iterator hint, Key&& k, M&& obj);

    template< class... Args >
    std::pair<iterator, bool> emplace(Args&&... args);

    template <class... Args>
    iterator emplace_hint(const_iterator hint, Args&&... args);

    template< class... Args >
    std::pair<iterator, bool> try_emplace(const Key& k, Args&&... args);

    template< class... Args >
    std::pair<iterator, bool> try_emplace(Key&& k, Args&&... args);

    template< class... Args >
    iterator try_emplace(const_iterator hint, const Key& k, Args&&... args);

    template< class... Args >
    iterator try_emplace(const_iterator hint, Key&& k, Args&&... args);

    iterator erase(iterator pos);

    iterator erase(const_iterator pos);

    iterator erase(const_iterator first, const_iterator last);

    size_type erase(const Key& key);

    void swap(unordered_map& other) noexcept;

    /// \brief Provides access to the value associated with the given key,
    ///        inserting a default value when the key is not yet present.
    mapped_type& operator[](const key_type& key);

    /// \brief Provides access to the value associated with the given key.
    [[nodiscard]] const mapped_type& at(const key_type& key) const;

    /// \returns The number of elements matching the specified key.
    [[nodiscard]] size_type count(const key_type& key) const;

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
      for (const auto& element : m_container)
      {
        element.first.mark(todo);
        element.second.mark(todo);
      }
    }

  protected:

    void rehash(std::size_t new_size);

    /// \brief Rehashes the container when the load factor is exceeded, re-checking the
    ///        load factor under the exclusive lock when ThreadSafe to avoid redundant rehashes.
    void rehash_if_needed();

    bool rehash_is_needed() const;
};
}

} // namespace atermpp
#endif // MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_MAP_H


