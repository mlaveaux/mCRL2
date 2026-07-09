// Author(s): Jan Friso Groote, Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//
/// \file mcrl2/atermpp/standard_containers/deque.h
/// \brief This file contains a deque class that behaves 
///        exactly as a standard deque. It can only be used
///        to store class instances that derive from aterms.
///        The stored aterms are protected as a whole, i.e.,
///        time and memory is saved as individual protection
///        per element is unnecessary. 

#ifndef MCRL2_ATERMPP_STANDARD_CONTAINER_DEQUE_H
#define MCRL2_ATERMPP_STANDARD_CONTAINER_DEQUE_H

#include <compare>
#include <deque>
#include "mcrl2/atermpp/detail/aterm_container.h"
#include "mcrl2/atermpp/detail/global_aterm_pool.h"
#include "mcrl2/atermpp/detail/thread_aterm_pool.h"
#include "mcrl2/utilities/shared_mutex.h"

/// \brief The main namespace for the aterm++ library.
namespace atermpp
{

/// \brief A deque class in which aterms can be stored.
/// \details The API of a std::deque is exposed by explicit forwarding (composition) instead of
///          inheritance, such that every structural modification is guaranteed to acquire the
///          shared lock that keeps the garbage collector out of the container while it changes.
template < class T, class Alloc = std::allocator<detail::markable_aterm<T> > > 
class deque
{
protected:
  using super = std::deque<detail::markable_aterm<T>, Alloc>;

  super m_container;
  detail::generic_aterm_container<super> container_wrapper;

public:
  
  /// Standard typedefs.
  using allocator_type = typename super::allocator_type;
  using value_type = typename super::value_type;
  using size_type = typename super::size_type;
  using difference_type = typename super::difference_type;
  using reference = typename super::reference;
  using const_reference = typename super::const_reference;
  using pointer = typename super::pointer;
  using const_pointer = typename super::const_pointer;
  using iterator = typename super::iterator;
  using const_iterator = typename super::const_iterator;
  using reverse_iterator = typename super::reverse_iterator;
  using const_reverse_iterator = typename super::const_reverse_iterator;

  /// \brief Default constructor.
  deque()
   : m_container(),
     container_wrapper(m_container)     
  {}

  /// \brief Constructor.
  explicit deque(const allocator_type& alloc)
   : m_container(alloc),
     container_wrapper(m_container)     
  {}

  /// \brief Constructor.
  explicit deque(size_type n, const allocator_type& alloc = allocator_type())
   : m_container(n, alloc),
     container_wrapper(m_container)
  {}

  /// \brief Constructor.
  deque(size_type n, const value_type& val, const allocator_type& alloc = allocator_type())
   : m_container(n, detail::markable_aterm<T>(val), alloc),
     container_wrapper(m_container)    
  {}

  /// \brief Constructor.
  template <class InputIterator>
  deque(InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type())
   : m_container(first, last, alloc),
     container_wrapper(m_container)     
  {}
    
  /// \brief Constructor.
  deque(const deque& x)
   : m_container(x.m_container),
     container_wrapper(m_container)     
  {}

  /// \brief Constructor.
  deque(const deque& x, const allocator_type& alloc)
   : m_container(x.m_container, alloc),
     container_wrapper(m_container)     
  {}
  
  /// \brief Constructor.
  deque(deque&& x) noexcept
      : m_container(std::move(x.m_container)),
        container_wrapper(m_container)
  {}

  /// \brief Constructor.
  deque(deque&& x, const allocator_type& alloc)
   : m_container(std::move(x.m_container), alloc),
     container_wrapper(m_container)     
  {}

  /// \brief Constructor. 
  deque(std::initializer_list<value_type> il, const allocator_type& alloc = allocator_type())
    : m_container(il.begin(), il.end(), alloc),
      container_wrapper(m_container)      
  {}

  /// \brief Copy assignment operator.
  deque& operator=(const deque& other)
  {
    if (this != &other)
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      m_container = other.m_container;
    }
    return *this;
  }

  /// \brief Move assignment operator.
  deque& operator=(deque&& other) noexcept
  {
    if (this != &other)
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      m_container = std::move(other.m_container);
    }
    return *this;
  }

  /// \brief Standard destructor.
  ~deque() = default;

  /// \returns The allocator of the underlying container.
  [[nodiscard]] allocator_type get_allocator() const noexcept
  {
    return m_container.get_allocator();
  }

  void assign(size_type count, const value_type& value)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.assign(count, value);
  }

  template<class InputIt>
  void assign(InputIt first, InputIt last)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.assign(first, last);
  }

  void assign(std::initializer_list<value_type> ilist)
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.assign(ilist);
  }

  [[nodiscard]] reference at(size_type pos) { return m_container.at(pos); }
  [[nodiscard]] const_reference at(size_type pos) const { return m_container.at(pos); }

  [[nodiscard]] reference operator[](size_type pos) { return m_container[pos]; }
  [[nodiscard]] const_reference operator[](size_type pos) const { return m_container[pos]; }

  [[nodiscard]] reference front() { return m_container.front(); }
  [[nodiscard]] const_reference front() const { return m_container.front(); }

  [[nodiscard]] reference back() { return m_container.back(); }
  [[nodiscard]] const_reference back() const { return m_container.back(); }

  iterator begin() noexcept { return m_container.begin(); }
  const_iterator begin() const noexcept { return m_container.begin(); }
  const_iterator cbegin() const noexcept { return m_container.cbegin(); }

  iterator end() noexcept { return m_container.end(); }
  const_iterator end() const noexcept { return m_container.end(); }
  const_iterator cend() const noexcept { return m_container.cend(); }

  reverse_iterator rbegin() noexcept { return m_container.rbegin(); }
  const_reverse_iterator rbegin() const noexcept { return m_container.rbegin(); }
  const_reverse_iterator crbegin() const noexcept { return m_container.crbegin(); }

  reverse_iterator rend() noexcept { return m_container.rend(); }
  const_reverse_iterator rend() const noexcept { return m_container.rend(); }
  const_reverse_iterator crend() const noexcept { return m_container.crend(); }

  [[nodiscard]] bool empty() const noexcept
  {
    // Concurrent read/write on the size.
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.empty();
  }

  [[nodiscard]] size_type max_size() const noexcept
  {
    return m_container.max_size();
  }

  void shrink_to_fit()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.shrink_to_fit();
  }

  void clear() noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.clear();
  }

  iterator insert( const_iterator pos, const T& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(pos, value);
  }

  iterator insert( const_iterator pos, T&& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(pos, std::move(value));
  }
  
  iterator insert( const_iterator pos, size_type count, const T& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(pos, count, value);
  }
    
  template< class InputIt >
  iterator insert( const_iterator pos,
                  InputIt first, InputIt last )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(pos, first, last);  
  }
    
  iterator insert( const_iterator pos, std::initializer_list<T> ilist )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.insert(pos, ilist.begin(), ilist.end());
  }
  
  template< class... Args >
  iterator emplace( const_iterator pos, Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.emplace(pos, std::forward<Args>(args)...);
  }

  iterator erase( const_iterator pos )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(pos);
  }

  iterator erase( const_iterator first, const_iterator last )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.erase(first, last);    
  }

  void push_back( const T& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.push_back(value);  
  }
  
  void push_back( T&& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.push_back(std::move(value));
  }

  template< class... Args >
  reference emplace_back( Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.emplace_back(std::forward<Args>(args)...);
  }

  void pop_back()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.pop_back();
  }

  void push_front( const T& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.push_front(value);
  }

  void push_front( T&& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.push_front(std::move(value));
  }

  template< class... Args >
  reference emplace_front( Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.emplace_front(std::forward<Args>(args)...);
  }

  void pop_front()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.pop_front();
  }

  void resize( size_type count )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.resize(count);
  }
  
  void resize( size_type count, const value_type& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.resize(count, value);
  }

  [[nodiscard]] std::size_t size() const
  {
    // Concurrent read/write on the size.
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return m_container.size();
  }

  void swap( deque& other ) noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    m_container.swap(other.m_container); // Invalidates end() so must be protected.
  }

  /// \brief Marks all stored terms during garbage collection; used when this
  ///        container is an element of another protected container.
  void mark(term_mark_stack& todo) const
  {
    for (const detail::markable_aterm<T>& element : m_container)
    {
      element.mark(todo);
    }
  }

  bool operator==(const deque& other) const { return m_container == other.m_container; }

  auto operator<=>(const deque& other) const { return m_container <=> other.m_container; }
};

} // namespace atermpp
#endif // MCRL2_ATERMPP_STANDARD_CONTAINER_DEQUE_H
