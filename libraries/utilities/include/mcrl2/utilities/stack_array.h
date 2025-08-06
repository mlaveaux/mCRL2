// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_STACK_VECTOR_H
#define MCRL2_UTILITIES_STACK_VECTOR_H

#include "mcrl2/utilities/detail/memory_utility.h"

#include <cstddef>
#include <iterator>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/noncopyable.cxx"
#else
  import utilities;
#endif

namespace mcrl2::utilities
{

/// \brief Provides (a subset of) the interface of std::array<T> for a portion of preallocated memory. Can be used to
///        interface with a portion of memory allocated on the stack, \see MCRL2_DECLARE_STACK_ARRAY. The advantage over
///        MCRL2_SPECIFIC_STACK_ALLOCATOR is that the lifetime of the underlying objects is bounded by the lifetime of the
///        stack_array.
template<typename T>
class stack_array : public noncopyable
{
public:
  using iterator = T*;
  using const_iterator = const T*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  /// \brief The given pointer should be able to hold N element of sizeof(T) bytes.
  stack_array(T* reserved_memory, std::size_t N)
    : m_reserved_memory(reserved_memory),
      m_size(N)
  {
    for (T& element : *this)
    {
      new (&element) T();
    }
  }

  ~stack_array()
  {
    for (T& element : *this)
    {
      element.~T();
    }
  }

  // The remaining interface of std::array

  iterator begin() { return data(); }
  const_iterator begin() const { return data(); }

  iterator end() { return data() + size();  }
  const_iterator end() const { return data() + size(); }

  T* data() { return m_reserved_memory; }
  const T* data() const { return m_reserved_memory; }

  bool empty() const { return size() != 0; }

  reverse_iterator rbegin() { return reverse_iterator(data() + size()); }
  const_reverse_iterator rbegin() const { return const_reverse_iterator(data() + size()); }

  reverse_iterator rend() { return reverse_iterator(data()); }
  const_reverse_iterator rend() const { return const_reverse_iterator(data()); }

  std::size_t size() const { return m_size; }

  std::size_t max_size() const { return m_size; }

  T& operator[](std::size_t index)
  {
    assert(index < size());
    return m_reserved_memory[index];
  }

private:
  T* m_reserved_memory;
  std::size_t m_size;
};

} // namespace mcrl2::utilities


/// \brief Define a (hopefully) unique name for the underlying reserved stack memory.
#define MCRL2_STACK_ARRAY_NAME(NAME) \
  NAME ## _reserved_stack_memory

/// \brief Declares a stack_array<TYPE> with the specified NAME that stores SIZE elements type TYPE.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MCRL2_DECLARE_STACK_ARRAY(NAME, TYPE, SIZE) \
  TYPE* MCRL2_STACK_ARRAY_NAME(NAME) = MCRL2_SPECIFIC_STACK_ALLOCATOR(TYPE, SIZE); \
  mcrl2::utilities::stack_array<TYPE> NAME (MCRL2_STACK_ARRAY_NAME(NAME), SIZE) 

#endif // MCRL2_UTILITIES_STACK_VECTOR_H
