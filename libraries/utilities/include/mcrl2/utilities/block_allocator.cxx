// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_BLOCK_ALLOCATOR_H_
#define MCRL2_UTILITIES_BLOCK_ALLOCATOR_H_
MCRL2_MODULE;

#include <memory>
#include <stddef.h>

#ifndef MCRL2_ENABLE_MODULES
  #include "memory_pool.cxx"
#else
  export module utilities:block_allocator;

  import :memory_pool;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities
{

/// \brief The block allocator provides the allocator interface for the memory pool class.
///        As such is can be used as an allocator for allocator aware storages.
/// \details Also provides several non-standard allocate functions specifically for the term pool.
template <class T,
          std::size_t ElementsPerBlock = 1024,
          bool ThreadSafe = false>
class block_allocator : public memory_pool<T, ElementsPerBlock, ThreadSafe>
{
private:
  using super = memory_pool<T, ElementsPerBlock, ThreadSafe>;

public:
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  template <class U>
  struct rebind
  {
    using other = block_allocator<U, ElementsPerBlock, ThreadSafe>;
  };

  block_allocator() = default;

  /// \details The unused parameter is to make the interface equivalent
  ///          to the allocator.
  T* allocate(size_type n, const void* hint = nullptr)
  {
    if (n != 1 || hint)
    {
      throw std::bad_alloc();
    }

    return super::allocate();
  }

  /// \details The unused parameter is to make the interface equivalent
  ///          to the allocator.
  void deallocate(T* p, size_type)
  {
    super::deallocate(p);
  }

  // Move assignment and construction is possible.
  block_allocator(block_allocator&&) = default;
  block_allocator& operator=(block_allocator&&) = default;
};

} // namespace mcrl2::utilities


#endif
