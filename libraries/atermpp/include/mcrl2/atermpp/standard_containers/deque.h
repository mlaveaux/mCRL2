// Author(s): Jan Friso Groote, Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//
/// \file mcrl2/data/standard_containers/deque.h
/// \brief This file contains a deque class that behaves 
///        exactly as a standard deque. It can only be used
///        to store class instances that derive from aterms.
///        The stored aterms are protected as a whole, i.e.,
///        time and memory is saved as individual protection
///        per element is unnecessary. 

#ifndef MCRL2_ATERMPP_STANDARD_CONTAINER_DEQUE_H
#define MCRL2_ATERMPP_STANDARD_CONTAINER_DEQUE_H

#include <deque>
#include "mcrl2/atermpp/detail/aterm_container.h"
#include "mcrl2/atermpp/detail/global_aterm_pool.h"
#include "mcrl2/atermpp/detail/thread_aterm_pool.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/shared_mutex.cxx"
#else
  import utilities;
#endif

/// \brief The main namespace for the aterm++ library.
namespace atermpp
{

/// \brief A deque class in which aterms can be stored. 
template < class T, class Alloc = std::allocator<detail::reference_aterm<T> > > 
class deque : public std::deque< detail::reference_aterm<T>, Alloc >              
{
protected:
  using super = std::deque<detail::reference_aterm<T>, Alloc>;

  detail::generic_aterm_container<std::deque<detail::reference_aterm<T>, Alloc>> container_wrapper;

public:
  
  /// Standard typedefs.
  using allocator_type = typename super::allocator_type;
  using value_type = typename super::value_type;
  using size_type = typename super::size_type;
  using reference = typename super::reference;
  using iterator = typename super::iterator;
  using const_iterator = typename super::const_iterator;

  /// \brief Default constructor.
  deque()
   : super(),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  explicit deque (const allocator_type& alloc)
   : super::deque(alloc),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  explicit deque (size_type n, const allocator_type& alloc = allocator_type())
   : super::deque(n, alloc),
     container_wrapper(*this)
  {}

  /// \brief Constructor.
  deque(size_type n, const value_type& val, const allocator_type& alloc = allocator_type())
   : super::deque(n, detail::reference_aterm(val), alloc),
     container_wrapper(*this)    
  {}

  /// \brief Constructor.
  template <class InputIterator>
  deque(InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type())
   : super::deque(first, last, alloc),
     container_wrapper(*this)     
  {}
    
  /// \brief Constructor.
  deque(const deque& x)
   : super::deque(x),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  deque(const deque& x, const allocator_type& alloc)
   : super::deque(x, alloc),
     container_wrapper(*this)     
  {}
  
  /// \brief Constructor.
  deque(deque&& x) noexcept
      : super::deque(std::move(x)),
        container_wrapper(*this)
  {}

  /// \brief Constructor.
  deque(deque&& x, const allocator_type& alloc)
   : super::deque(std::move(x), alloc),
     container_wrapper(*this)     
  {}

  /// \brief Constructor. 
  deque(std::initializer_list<value_type> il, const allocator_type& alloc = allocator_type())
    : super::deque(il.begin(), il.end(), alloc),
      container_wrapper(*this)      
  {}

  /// \brief Copy assignment operator.
  deque& operator=(const deque& other) = default;

  /// \brief Move assignment operator.
  deque& operator=(deque&& other) = default;

  /// \brief Standard destructor.
  ~deque() = default;

  void shrink_to_fit()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::shrink_to_fit();
  }

  void clear() noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::clear();
  }

  iterator insert( const_iterator pos, const T& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, value);
  }

  iterator insert( const_iterator pos, T&& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, value);
  }
  
  iterator insert( const_iterator pos, size_type count, const T& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, count, value);
  }
    
  template< class InputIt >
  iterator insert( const_iterator pos,
                  InputIt first, InputIt last )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, first, last);  
  }
    
  iterator insert( const_iterator pos, std::initializer_list<T> ilist )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, ilist);
  }
  
  template< class... Args >
  iterator emplace( const_iterator pos, Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::emplace(pos, args...);
  }

  iterator erase( const_iterator pos )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::erase(pos);
  }

  iterator erase( const_iterator first, const_iterator last )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::erase(first, last);    
  }

  void push_back( const T& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::push_back(value);  
  }
  
  void push_back( T&& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::push_back(value);
  }

  template< class... Args >
  reference emplace_back( Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::emplace_back(args...);
  }

  void pop_back()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::pop_back();
  }

  void push_front( const T& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::push_front(value);
  }
		
  void push_front( T&& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::push_front(value);
  }

  template< class... Args >
  reference emplace_front( Args&&... args )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::emplace_front(args...);
  }

  void resize( size_type count )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::resize(count);
  }
  
  void resize( size_type count, const value_type& value )
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::resize(count, value);
  }

  std::size_t size() const
  {
    return super::size();
  }

  void swap( deque& other ) noexcept
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::swap(other); // Invalidates end() so must be protected.
  }
};

} // namespace atermpp
#endif // MCRL2_ATERMPP_STANDARD_CONTAINER_DEQUE_H
