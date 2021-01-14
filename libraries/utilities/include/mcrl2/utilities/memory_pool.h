// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_MEMORY_POOL_H_
#define MCRL2_UTILITIES_MEMORY_POOL_H_

#include "mcrl2/utilities/detail/free_list.h"
#include "mcrl2/utilities/noncopyable.h"

#include <array>
#include <cstdint>
#include <forward_list>
#include <type_traits>
#include <mutex>

namespace mcrl2
{
namespace utilities
{

/// \brief The memory pool allocates elements of size T from blocks.
/// \details When ThreadSafe is true then the thread-safe guarantees will be satisfied.
template <class T, 
          std::size_t ElementsPerBlock = 1024, 
          bool ThreadSafe = false>
class memory_pool : private mcrl2::utilities::noncopyable
{
public:
  memory_pool() = default;

  /// \brief Triggers the (possibly non-trivial) destructor of all elements stored
  ///        in the pool.
  ~memory_pool()
  {
    m_freelist.destructive_mark();

    /// For all actual elements stored in the pool trigger the destructor.
    bool first_block = true;
    for (auto& block : m_blocks)
    {
      if (first_block)
      {
        // This is the first block, for that one only m_current_index elements were inserted.
        for (std::size_t index = 0; index < m_current_index; ++index)
        {
          auto& slot = block[index];
          if (!slot.is_marked())
          {
            reinterpret_cast<T*>(&slot)->~T();
          }
        }
        first_block = false;
      }
      else
      {
        for (auto& slot : block)
        {
          if (!slot.is_marked())
          {
            reinterpret_cast<T*>(&slot)->~T();
          }
        }
      }
    }

    assert(m_freelist.empty());
  }

  /// \brief Reuses memory from block and allocates a new block when
  ///        no slots are free.
  /// \returns A pointer to a block of memory that can store an object of type T.
  /// \threadsafe
  T* allocate()
  {
    // Only allow one thread to allocate at the time.
    if constexpr (ThreadSafe) { m_block_mutex.lock(); }

    if (!m_freelist.empty())
    {
      T& element = m_freelist.pop_front();
      if constexpr (ThreadSafe) { m_block_mutex.unlock(); }
      return &element;
    }

    if (m_current_index >= ElementsPerBlock)
    {
      // The whole buffer was used so allocate a new one.
      m_blocks.emplace_front();
      ++m_number_of_blocks;
      m_current_index = 0;
    }

    // The object was last written as this slot is not part of the freelist.
    T& slot = (m_blocks.front()[m_current_index++]).element();

    if constexpr (ThreadSafe) { m_block_mutex.unlock(); }
    assert(contains(&slot));
    return &slot;
  }

  /// \brief Free the memory used by the given pointer that has been allocated by this pool.
  void deallocate(T* pointer)
  {    
    if constexpr (ThreadSafe) { m_block_mutex.lock(); }

    assert(contains(pointer));
    m_freelist.push_front(reinterpret_cast<Slot&>(*pointer));

    if constexpr (ThreadSafe) { m_block_mutex.unlock(); }
  }

  /// \brief Frees blocks that are no longer storing elements of T.
  /// \returns The number of blocks that have been removed.
  std::size_t consolidate()
  {
    m_freelist.destructive_mark();
    std::size_t old_number_of_blocks = m_number_of_blocks;

    // Keep track of the last block to be able to erase the current block.
    auto block_before_it = m_blocks.before_begin();
    for (auto it = m_blocks.begin(); it != m_blocks.end();)
    {
      Block& block = *it;

      // Keep track of the head of the freelist.
      FreelistIt old_freelist_head = m_freelist.begin();

      // Indicate that current block only contains slots in the freelist.
      bool block_only_freelist = true;
      for (Slot& slot : block)
      {
        if (slot.is_marked())
        {
          m_freelist.push_front(slot);
        }
        else
        {
          // There is one slot that was not part of the freelist.
          block_only_freelist = false;
        }
      }

      // Erase blocks that only have elements in the freelist.
      if (block_only_freelist)
      {
        // Remove the slots in the freelist that point into this block.
        m_freelist.erase_after(m_freelist.before_begin(), old_freelist_head);

        // The current block only has elements in the freelist, so erase it.
        it = m_blocks.erase_after(block_before_it);
        --m_number_of_blocks;
      }
      else
      {
        block_before_it = it;
        ++it;
      }
    }

    return old_number_of_blocks - m_number_of_blocks;
  }

  /// \returns True when thi memory pool has space for at least one more element without allocating
  ///          new memory.
  bool has_free_slots() const noexcept
  {
    // The freelist is not full, i.e. points to some slot or the index has not yet reached the end.
    return m_current_index < ElementsPerBlock || m_freelist.empty();
  }

  /// \returns The total number of elements that could be stored in this memory pool.
  std::size_t capacity() const noexcept
  {
    return m_number_of_blocks * ElementsPerBlock;
  }

  // Move constructor and assignment are possible.
  memory_pool(memory_pool&& other) = default;
  memory_pool& operator=(memory_pool&& other) = default;

private:
  using Freelist = typename detail::free_list<T>;
  using FreelistIt = typename Freelist::iterator;
  using Slot = typename Freelist::slot;

  /// \brief An array that stores ElementsPerBlock number of objects of type T.
  using Block = std::array<Slot, ElementsPerBlock>;

  /// \brief The last slot in the first block that has never been returned by allocate.
  using SizeType = typename std::conditional<ThreadSafe, std::atomic<std::size_t>, std::size_t>::type;
  SizeType m_current_index = ElementsPerBlock;

  /// \brief Equal to the size of the blocks array to prevent iterating over the block list.
  SizeType m_number_of_blocks = 0;

  /// \brief The list of blocks allocated by this pool.
  std::forward_list<Block> m_blocks;

  /// \brief Ensures that the block list is only modified by a single thread.
  std::mutex m_block_mutex = {};

  /// \brief Indicates the head of the freelist.
  Freelist m_freelist;

  /// \returns Check whether the pointer is contained in this memory pool.
  bool contains(T* p)
  {
    assert(p != nullptr);

    std::uintptr_t pointer = reinterpret_cast<std::uintptr_t>(p);
    for (auto& block : m_blocks)
    {
      std::uintptr_t firstSlot = reinterpret_cast<std::uintptr_t>(block.data());
      if (firstSlot <= pointer && pointer < firstSlot + sizeof(Slot) * ElementsPerBlock)
      {
        assert((pointer - firstSlot) % sizeof(Slot) == 0);
        return true;
      }
    }

    return false;
  }

};

} // namespace utilities
} // namespace mcrl2

#endif // MCRL2_UTILITIES_MEMORY_POOL_H_
