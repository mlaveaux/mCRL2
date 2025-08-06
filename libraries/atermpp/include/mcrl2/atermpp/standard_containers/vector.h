// Author(s): Jan Friso Groote, Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//
/// \file mcrl2/data/standard_containers/vector.h
/// \brief This file contains a vector class that behaves 
///        exactly as a standard vector. It can only be used
///        to store class instances that derive from aterms.
///        The stored aterms are protected as a whole, i.e.,
///        time and memory is saved as individual protection
///        per element is unnecessary. 

#ifndef MCRL2_ATERMPP_STANDARD_CONTAINER_VECTOR_H
#define MCRL2_ATERMPP_STANDARD_CONTAINER_VECTOR_H

#include <vector>


#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/aterm_container.cxx"
  #include "mcrl2/atermpp/detail/thread_aterm_pool.cxx"
  #include "mcrl2/utilities/shared_mutex.cxx"
#else
  import atermpp;
  import utilities;
#endif

/// \brief The main namespace for the aterm++ library.
namespace atermpp
{

/// \brief A vector class in which aterms can be stored. 
template < class T, class Alloc = std::allocator<detail::reference_aterm<T> >, bool ThreadSafe = false > 
class vector : public std::vector< detail::reference_aterm<T>, Alloc >
{
protected:
  using super = std::vector<detail::reference_aterm<T>, Alloc>;

  detail::generic_aterm_container<std::vector<detail::reference_aterm<T>, Alloc>> container_wrapper;

public: 
  /// Standard typedefs.
  using allocator_type = typename super::allocator_type;
  using value_type = typename super::value_type;
  using size_type = typename super::size_type;
  using reference = typename super::reference;
  using iterator = typename super::iterator;
  using const_iterator = typename super::const_iterator;

  /// \brief Default constructor.
  vector()
   : super(),
     container_wrapper(*this)
  {}

  /// \brief Constructor.
  explicit vector (const allocator_type& alloc)
   : super::vector(alloc),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  explicit vector (size_type n, const allocator_type& alloc = allocator_type())
   : super::vector(n, alloc),
     container_wrapper(*this)     
  {}

  vector (size_type n, const value_type& val, const allocator_type& alloc = allocator_type())
   : super::vector(n, detail::reference_aterm(val), alloc),
     container_wrapper(*this)
  {}

  /// \brief Constructor.
  template <class InputIterator>
  vector (InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type())
   : super::vector(first, last, alloc),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  vector (const vector& x)
   : super::vector(x),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  vector (const vector& x, const allocator_type& alloc)
   : super::vector(x, alloc),
     container_wrapper(*this)     
  {}
  
  /// \brief Constructor.
  vector(vector&& x) noexcept
      : super::vector(std::move(x)),
        container_wrapper(*this)
  {}

  /// \brief Constructor.
  vector (vector&& x, const allocator_type& alloc)
   : super::vector(std::move(x), alloc),
     container_wrapper(*this)
  {}

  /// \brief Constructor. To be done later....
  vector (std::initializer_list<value_type> il, const allocator_type& alloc = allocator_type())
    : super::vector(il, alloc),
      container_wrapper(*this)      
  {}

  /// \brief Assignment operator.
  vector& operator=(const vector& x) = default;
  
  /// \brief Move assignment operator
  vector& operator=(vector&& x) = default;

  /// \brief Standard destructor.
  ~vector() = default;

  void shrink_to_fit()
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::shrink_to_fit();
    } else {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::shrink_to_fit();
    }
  }

  void clear() noexcept
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::clear();
    } else {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::clear();
    }
  }

  iterator insert( const_iterator pos, const T& value )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();

      // This is not thread safe otherwise since the length or end iterator is updated during this producedure.
      return super::insert(pos, value);
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, value);
  }

  iterator insert( const_iterator pos, T&& value )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      return super::insert(pos, value);
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, value);
  }
  
  iterator insert( const_iterator pos, size_type count, const T& value )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      return super::insert(pos, count, value);
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, count, value);
  }
    
  template< class InputIt >
  iterator insert( const_iterator pos,
                  InputIt first, InputIt last )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      return super::insert(pos, first, last); 
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, first, last);  
  }
    
  iterator insert( const_iterator pos, std::initializer_list<T> ilist )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      return super::insert(pos, ilist);
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::insert(pos, ilist);
  }
  
  template< class... Args >
  iterator emplace( const_iterator pos, Args&&... args )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      return super::emplace(pos, args...);   
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::emplace(pos, args...);   
  }

  iterator erase( const_iterator pos )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      return super::erase(pos);
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::erase(pos);
  }

  iterator erase( const_iterator first, const_iterator last )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      return super::erase(first, last);
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::erase(first, last);
  }

  void push_back( const T& value )
  {
    if constexpr (ThreadSafe) 
    {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::push_back(value);
    } 
    else 
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::push_back(value);
    }
  }

  void push_back( T&& value )
  {
    if constexpr (ThreadSafe) 
    {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::push_back(value);
    } 
    else 
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::push_back(value);
    }
  }

  template< class... Args >
  reference emplace_back( Args&&... args )
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      return super::emplace_back(args...);  
    }

    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::emplace_back(args...);  
  }

  void pop_back()
  {
    if constexpr (ThreadSafe) 
    {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::pop_back();
    } 
    else 
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::pop_back();
    }
  }

  void resize( size_type count )
  {
    if constexpr (ThreadSafe) 
    {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::resize(count);
    } 
    else 
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::resize(count);
    }
  }

  void resize( size_type count, const value_type& value )
  {
    if constexpr (ThreadSafe) 
    {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::resize(count, value);
    } 
    else 
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::resize(count, value);
    }
  }

  void swap( vector& other ) noexcept
  {
    if constexpr (ThreadSafe) 
    {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::swap(other); // Invalidates end() so must be protected.
    } 
    else 
    {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::swap(other); // Invalidates end() so must be protected.
    }
  }

  std::size_t size() const
  {
    // Concurrent read/write on the size.
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::size();
  }
};

} // namespace atermpp
#endif // MCRL2_ATERMPP_STANDARD_CONTAINER_VECTOR_H
