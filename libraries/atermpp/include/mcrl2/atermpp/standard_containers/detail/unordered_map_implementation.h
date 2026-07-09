// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//
/// \file mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h
/// \brief This file contains an implementation of the hash function 
///        to break circular header dependencies. 

#ifndef MCRL2_ATERMPP_STANDARD_CONTAINER_DETAIL_UNORDERED_MAP_IMPLEMENTATION_H
#define MCRL2_ATERMPP_STANDARD_CONTAINER_DETAIL_UNORDERED_MAP_IMPLEMENTATION_H

#include "mcrl2/atermpp/standard_containers/unordered_map.h"
#include "mcrl2/atermpp/detail/thread_aterm_pool.h"
#include "mcrl2/utilities/shared_mutex.h"

namespace atermpp
{

  template< class Key, class T, class Hash, class Pred, class Alloc >
  unordered_map<Key,T,Hash,Pred,Alloc>& unordered_map<Key,T,Hash,Pred,Alloc>::operator=(const unordered_map& other)
  {
    if (this != &other)
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      m_container = other.m_container;
    }
    return *this;
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  unordered_map<Key,T,Hash,Pred,Alloc>& unordered_map<Key,T,Hash,Pred,Alloc>::operator=(unordered_map&& other) noexcept
  {
    if (this != &other)
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      m_container = std::move(other.m_container);
    }
    return *this;
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >  
  void unordered_map<Key,T,Hash,Pred,Alloc>::clear() noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.clear();
  }

  /// \brief Inserts an element referring to a default value in the map. 
  template< class Key, class T, class Hash, class Pred, class Alloc >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator,bool> unordered_map<Key,T,Hash,Pred,Alloc>::insert( const value_type& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(value);
  }
    
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template< class P >
    requires std::constructible_from<typename std::unordered_map<detail::markable_aterm<Key>, detail::markable_aterm<T>, Hash, Pred, Alloc>::value_type, P>
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator,bool> unordered_map<Key,T,Hash,Pred,Alloc>::insert( P&& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(std::forward<P>(value));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::insert( const_iterator hint, const value_type& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(hint, value);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template< class P >
    requires std::constructible_from<typename std::unordered_map<detail::markable_aterm<Key>, detail::markable_aterm<T>, Hash, Pred, Alloc>::value_type, P>
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::insert( const_iterator hint, P&& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(hint, std::forward<P>(value));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template< class InputIt >
  void unordered_map<Key,T,Hash,Pred,Alloc>::insert( InputIt first, InputIt last )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.insert(first, last);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  void unordered_map<Key,T,Hash,Pred,Alloc>::insert( std::initializer_list<value_type> ilist )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.insert(ilist);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::insert_return_type unordered_map<Key,T,Hash,Pred,Alloc>::insert( node_type&& nh )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(std::move(nh));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::insert( const_iterator hint, node_type&& nh )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(hint, std::move(nh));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template <class M>
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc>::insert_or_assign( const Key& k, M&& obj )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert_or_assign(k, std::forward<M>(obj));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template <class M>
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc>::insert_or_assign( Key&& k, M&& obj )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert_or_assign(std::move(k), std::forward<M>(obj));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template <class M>
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::insert_or_assign( const_iterator hint, const Key& k, M&& obj )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert_or_assign(hint, k, std::forward<M>(obj)); 
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template <class M>
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::insert_or_assign( const_iterator hint, Key&& k, M&& obj )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert_or_assign(hint, std::move(k), std::forward<M>(obj));   
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template< class... Args >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator,bool> unordered_map<Key,T,Hash,Pred,Alloc>::emplace( Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.emplace(std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template <class... Args>
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::emplace_hint( const_iterator hint, Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.emplace_hint(hint, std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template< class... Args >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc>::try_emplace( const Key& k, Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.try_emplace(detail::markable_aterm<Key>(k), std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template< class... Args >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc>::try_emplace( Key&& k, Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.try_emplace(detail::markable_aterm<Key>(std::move(k)), std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template< class... Args >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::try_emplace( const_iterator hint, const Key& k, Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.try_emplace(hint, detail::markable_aterm<Key>(k), std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  template< class... Args >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::try_emplace( const_iterator hint, Key&& k, Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.try_emplace(hint, detail::markable_aterm<Key>(std::move(k)), std::forward<Args>(args)...);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::erase( iterator pos )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(pos);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::erase( const_iterator pos )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(pos);   
  }	
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::erase( const_iterator first, const_iterator last )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(first, last);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::size_type unordered_map<Key,T,Hash,Pred,Alloc>::erase( const Key& key )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(detail::markable_aterm<Key>(key));
  }

  template<class Key, class T, class Hash, class Pred, class Alloc>
  void unordered_map<Key, T, Hash, Pred, Alloc>::swap(unordered_map& other) noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.swap(other.m_container);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::node_type unordered_map<Key,T,Hash,Pred,Alloc>::extract( const_iterator position )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.extract(position);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::node_type unordered_map<Key,T,Hash,Pred,Alloc>::extract( const Key& k )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.extract(detail::markable_aterm<Key>(k));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::mapped_type& unordered_map<Key,T,Hash,Pred,Alloc>::operator[]( const Key& key )
  {
    // The element access may insert a default element and can trigger a rehash.
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container[detail::markable_aterm<Key>(key)];
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::mapped_type& unordered_map<Key,T,Hash,Pred,Alloc>::at( const Key& key )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.at(detail::markable_aterm<Key>(key));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  const typename unordered_map<Key,T,Hash,Pred,Alloc>::mapped_type& unordered_map<Key,T,Hash,Pred,Alloc>::at( const Key& key ) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.at(detail::markable_aterm<Key>(key));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::find( const Key& key )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.find(detail::markable_aterm<Key>(key));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc>::find( const Key& key ) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.find(detail::markable_aterm<Key>(key));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::size_type unordered_map<Key,T,Hash,Pred,Alloc>::count( const Key& key ) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.count(detail::markable_aterm<Key>(key));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  bool unordered_map<Key,T,Hash,Pred,Alloc>::contains( const Key& key ) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.contains(detail::markable_aterm<Key>(key));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  void unordered_map<Key,T,Hash,Pred,Alloc>::rehash( size_type count )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.rehash(count);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  void unordered_map<Key,T,Hash,Pred,Alloc>::reserve( size_type count )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.reserve(count);
  }
      
  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::begin() 
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.begin();
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::iterator unordered_map<Key,T,Hash,Pred,Alloc>::end()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.end();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc>::begin() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.begin();
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc>::end() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.end();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc>::cbegin() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.cbegin();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc>::cend() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.cend();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  bool unordered_map<Key,T,Hash,Pred,Alloc>::empty() const noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.empty();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc >
  typename unordered_map<Key,T,Hash,Pred,Alloc>::size_type unordered_map<Key,T,Hash,Pred,Alloc>::max_size() const noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.max_size();
  }
} // namespace atermpp

namespace atermpp::utilities {  

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>&
  unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::operator=(const unordered_map& other)
  {
    if (this != &other)
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      m_container = other.m_container;
    }
    return *this;
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>&
  unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::operator=(unordered_map&& other) noexcept
  {
    if (this != &other)
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      m_container = std::move(other.m_container);
    }
    return *this;
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  inline void unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::rehash(std::size_t new_size )
  { 
    if constexpr (ThreadSafe)
    {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      m_container.rehash(new_size);
    }
    else
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      m_container.rehash(new_size);
    }
  }

  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  inline bool unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::rehash_is_needed() const
  {
    return m_container.load_factor() >= m_container.max_load_factor();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  inline void unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::rehash_if_needed()
  {
    if constexpr (ThreadSafe)
    {
      if (rehash_is_needed())
      {
        mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
        // Re-check under the exclusive lock: another thread may have rehashed already.
        if (rehash_is_needed())
        {
          m_container.rehash(m_container.bucket_count() * 2);
        }
      }
    }
    else
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      if (rehash_is_needed())
      {
        m_container.rehash(m_container.bucket_count() * 2);
      }
    }
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  void unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::clear() noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.clear();
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template<typename ...Args>
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::find(const Args&... args)
  { 
    if constexpr (ThreadSafe)
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      return m_container.find(args...); 
    }
    return m_container.find(args...); 
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template<typename ...Args>
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::find(const Args&... args) const 
  {
    if constexpr (ThreadSafe)
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      return m_container.find(args...);
    }
    return m_container.find(args...);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::mapped_type&
  unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::operator[](const key_type& key)
  {
    // The element access may insert a default element and can trigger a rehash.
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container[key];
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  const typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::mapped_type&
  unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::at(const key_type& key) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.at(key);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::size_type
  unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::count(const key_type& key) const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.count(key);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert(const value_type& value)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(value);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template< class P >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert(P&& value)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(std::forward<P>(value));
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert(const_iterator hint, value_type&& value)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(hint, std::move(value));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template< class P >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert(const_iterator hint, P&& value)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(hint, std::forward<P>(value));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template< class InputIt >
  void unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert(InputIt first, InputIt last)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.insert(first, last);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  void unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert(std::initializer_list<value_type> ilist)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.insert(ilist);
  }

    
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template <class M>
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert_or_assign(const Key& k, M&& obj)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert_or_assign(k, std::forward<M>(obj));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template <class M>
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert_or_assign(Key&& k, M&& obj)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert_or_assign(std::move(k), std::forward<M>(obj));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template <class M>
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert_or_assign(const_iterator hint, const Key& k, M&& obj)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert_or_assign(hint, k, std::forward<M>(obj));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template <class M>
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::insert_or_assign(const_iterator hint, Key&& k, M&& obj)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert_or_assign(hint, std::move(k), std::forward<M>(obj));
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template< class... Args >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::emplace(Args&&... args)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.emplace(std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template <class... Args>
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::emplace_hint(const_iterator hint, Args&&... args)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.emplace_hint(hint, std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template< class... Args >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::try_emplace(const Key& k, Args&&... args)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.try_emplace(k, std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template< class... Args >
  std::pair<typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator, bool> unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::try_emplace(Key&& k, Args&&... args)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.try_emplace(std::move(k), std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template< class... Args >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::try_emplace(const_iterator hint, const Key& k, Args&&... args)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.try_emplace(hint, k, std::forward<Args>(args)...);
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  template< class... Args >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::try_emplace(const_iterator hint, Key&& k, Args&&... args)
  {
    rehash_if_needed();
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.try_emplace(hint, std::move(k), std::forward<Args>(args)...);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::erase(iterator pos)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(pos);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::erase(const_iterator pos)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(pos);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::erase(const_iterator first, const_iterator last)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(first, last);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::size_type unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::erase(const Key& key)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();

    // The underlying erase(key) returns void and requires the key to be present,
    // so find first and report the number of erased elements as std::unordered_map does.
    auto it = m_container.find(key);
    if (it == m_container.end())
    {
      return 0;
    }
    m_container.erase(typename super::const_iterator(it));
    return 1;
  }

  template<class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe>
  void unordered_map<Key, T, Hash, Pred, Alloc, ThreadSafe>::swap(unordered_map& other) noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.swap(other.m_container);
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::begin() 
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.begin();
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::end()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.end();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::begin() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.begin();
  }
  
  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::end() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.end();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::cbegin() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.cbegin();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::const_iterator unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::cend() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.cend();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  bool unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::empty() const noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.empty();
  }

  template< class Key, class T, class Hash, class Pred, class Alloc, bool ThreadSafe >
  typename unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::size_type unordered_map<Key,T,Hash,Pred,Alloc,ThreadSafe>::max_size() const noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.max_size();
  }
} // namespace atermpp::utilities

#endif // MCRL2_ATERMPP_STANDARD_CONTAINER_DETAIL_UNORDERED_MAP_IMPLEMENTATION_H
