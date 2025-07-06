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

#ifndef MCRL2_ATERMPP_STANDARD_CONTAINER_STACK_H
#define MCRL2_ATERMPP_STANDARD_CONTAINER_STACK_H

#include "mcrl2/atermpp/unprotected_aterm.h"
#include "mcrl2/atermpp/detail/aterm_container.h"
#include "mcrl2/atermpp/standard_containers/deque.h"

#include <mcrl3/mutex.h>

/// \brief The main namespace for the aterm++ library.
namespace atermpp
{

/// \brief A deque class in which aterms can be stored. 
template < class T, 
           class Container = atermpp::deque< T > > 
class stack
{
private:
  Container m_container;

public:
  
  /// Standard typedefs.
  using allocator_type = typename Container::allocator_type;
  using value_type = typename Container::value_type;
  using size_type = typename Container::size_type;
  using reference = typename Container::reference;
  using const_reference = typename Container::const_reference;
  using iterator = typename Container::iterator;
  using const_iterator = typename Container::const_iterator;
  
  /// \brief Constructor.
  explicit stack(const Container& cont = Container())
   : m_container(cont)     
  {}

  /// \brief Constructor.
  explicit stack( Container&& cont )
   : m_container(cont)
  {}

  stack(const stack& other) = default;
  stack(stack&& other) = default;

  /// \brief Constructor.
  template< class InputIt >
  stack( InputIt first, InputIt last )
    : m_container(first, last)
  {}

  /// \brief Constructor.
  template <class InputIterator>
  stack(InputIterator first, InputIterator last,
        const allocator_type& alloc = allocator_type())
   : m_container(first, last, alloc)  
  {}

  template< class Alloc >
  explicit stack( const Alloc& alloc )
    : m_container(alloc)
  {}

  template< class Alloc >
  stack( const Container& cont, const Alloc& alloc )
    : m_container(cont, alloc)
  {}

  template< class Alloc >
  stack( Container&& cont, const Alloc& alloc )
    : m_container(cont, alloc)
  {}

  template< class Alloc >
  stack( const stack& other, const Alloc& alloc)
    : m_container(other.m_container, alloc)
  {}

  template< class Alloc >
  stack( stack&& other, const Alloc& alloc)
    : m_container(std::move(other.m_container), alloc)
  {}

  template< class InputIt, class Alloc >
  stack( InputIt first, InputIt last, const Alloc& alloc )
    : m_container(first, last, alloc)
  {}
    
  /// \brief Copy assignment operator.
  stack& operator=(const stack& other) = default;

  /// \brief Move assignment operator.
  stack& operator=(stack&& other) = default;

  /// \brief Standard destructor.
  ~stack()=default;

  reference top()
  {
    return m_container.back();
  }

  const_reference top() const
  {
    return m_container.back();
  }

  bool empty() const
  {
    return m_container.empty();
  }

  size_type size() const
  {
    return m_container.size();
  }

  void push( const value_type& value )
  {
    m_container.push_back(value);
  }

  void push( value_type&& value )
  {
    m_container.push_back(std::move(value));
  }

  template< class... Args >
  void emplace( Args&&... args )
  {
    m_container.emplace_back(std::forward<Args>(args)...);
  }

  void pop()
  {
    m_container.pop_back();
  }

  void swap( stack& other ) noexcept
  {
    using std::swap; swap(m_container, other.m_container);
  }

  void mark(std::stack<unprotected_aterm>& todo) const
  {
    m_container.mark(todo);
  }
};

} // namespace atermpp
#endif // MCRL2_ATERMPP_STANDARD_CONTAINER_STACK_H
