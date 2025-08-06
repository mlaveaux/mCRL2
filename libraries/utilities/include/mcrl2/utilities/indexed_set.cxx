// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_INDEXED_SET_H
#define MCRL2_UTILITIES_INDEXED_SET_H
MCRL2_MODULE;

#include <cassert>
#include <deque>
#include <mutex>
#include <vector>
#include <utility>
#include <memory>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/unordered_map.cxx"
  #include "mcrl2/utilities/detail/atomic_wrapper.cxx"
  #include "mcrl2/utilities/shared_mutex.cxx"  
#else
  export module utilities:indexed_set;

  import :unordered_map;
  import :shared_mutex;
  import :detail.atomic_wrapper;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities
{

/// \brief A set that assigns each element an unique index.
template<typename Key,
         bool ThreadSafe = false,
         typename Hash = std::hash<Key>,
         typename Equals = std::equal_to<Key>,
         typename Allocator = std::allocator<Key>,
         typename KeyTable = std::deque< Key, Allocator > > 
class indexed_set
{
private:
  std::vector<detail::atomic_wrapper<std::size_t>> m_hashtable;
  KeyTable m_keys;

  /// \brief Mutex for the m_hashtable and m_keys data structures.
  mutable std::shared_ptr<std::mutex> m_mutex;
  mutable std::vector<shared_mutex> m_shared_mutexes;

  /// m_next_index indicates the next index that 
  //  has not yet been used. This allows to increase m_keys in 
  //  large steps, avoiding exclusive access too often.  
  detail::atomic_wrapper<size_t> m_next_index;

  Hash m_hasher;
  Equals m_equals;

  /// \brief Reserve indices that can be used. Doing this 
  ///        infrequently prevents obtaining an exclusive lock for the
  ///        indexed set too often. This operation requires a
  ///        resize of m_keys. 
  void reserve_indices(std::size_t thread_index);

  /// \brief Inserts the given (key, n) pair into the indexed set.
  std::size_t put_in_hashtable(const Key& key, std::size_t value, std::size_t& new_position);

  /// \brief Resizes the hash table to twice its current size.
  void resize_hashtable();

public:
  using key_type = Key;
  using size_type = std::size_t;
  using value_type = std::pair<const key_type, size_type>;
  using key_equal = Equals;
  using hasher = Hash;

  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;

  using iterator = typename KeyTable::iterator;
  using const_iterator = typename KeyTable::const_iterator;

  using reverse_iterator = typename KeyTable::reverse_iterator;
  using const_reverse_iterator = typename KeyTable::const_reverse_iterator;

  using difference_type = std::ptrdiff_t;

  /// \brief Value returned when an element does not exist in the set.
  /// \return Value indicating non existing element, equal to std::numeric_limits<std::size_t>::max(). 
  static constexpr size_type npos = std::numeric_limits<std::size_t>::max();

  /// \brief Constructor of an empty indexed set. Starts with a hashtable of size 128 and assumes one single thread.
  indexed_set();

  /// \brief Constructor of an empty indexed set. 
  /// \details With a single thread it delivers contiguous values for states.
  ///          With multiple threads some indices may be skipped. Each thread
  ///          reserves numbers, which it hands out. If a thread does not have
  ///          the opportunity to hand out all numbers, holes in the contiguous
  ///          numbering can occur. The holes are always of limited size. 
  /// \param number_of_threads The number of threads that use this index set. If the number is 1, it is treated
  ///        as a sequential set. If this number is larger than 1, the threads must be numbered
  ///        from 1 up and including number_of_threads. The number 0 cannot be used in that case. 
  indexed_set(std::size_t number_of_threads);

  /// \brief Constructor of an empty index set. Starts with a hashtable of the indicated size. 
  /// \details With one thread the numbering is contiguous. With multiple threads limited size holes
  ///          can occur in the numbering. 
  /// \param number_of_threads The number of threads that use this index set. This number is either 1, and then the implementation 
  ///        assumes that the thread has number 0, or it is larger than 1, and it is assumed that threads are numbered from 1 upwards. 
  /// \param initial_hashtable_size The initial size of the hashtable.
  /// \param hash The hash function.
  /// \param equals The comparison function for its elements.
  indexed_set(
    std::size_t number_of_threads, 
    std::size_t initial_hashtable_size, 
    const hasher& hash = hasher(),
    const key_equal& equals = key_equal());

  /// \brief Returns a reference to the mapped value.
  /// \details Returns an invalid value, larger or equal than the size of the indexed set, if there is no element with the given key.
  size_type index(const key_type& key, std::size_t thread_index = 0) const;

  /// \brief Returns a reference to the mapped value.
  /// \details Returns an out_of_range exception if there is no element with the given key.
  /// \param index The position in the indexed set.
  /// \return The value at position index.
  const key_type& at(const size_type index) const;

  /// \brief Operator that provides a const reference at the position indicated by index.
  /// \param index The position in the indexed set.
  /// \return The value at position index.
  /// \details threadsafe
  const key_type& operator[](const size_type index) const;

  /// \brief Forward iterator which runs through the elements from the lowest to the largest number.
  /// \details Complexity is constant per operation.
  iterator begin(std::size_t thread_index = 0) 
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    iterator i = m_keys.begin();
    return i;
  }

  /// \brief End of the forward iterator.
  iterator end(std::size_t thread_index = 0)
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    iterator i = m_keys.begin()+m_next_index;
    return i;
  }

  /// \brief Forward iterator which runs through the elements from the lowest to the largest number.
  /// \details Complexity is constant per operation.
  const_iterator begin(std::size_t thread_index = 0) const
  {
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    const_iterator i = m_keys.begin();
    return i;
  }

  /// \brief End of the forward iterator.
  const_iterator end(std::size_t thread_index = 0) const
  {
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    const_iterator i = m_keys.begin()+m_next_index;
    return i;
  }

  /// \brief const_iterator going through the elements in the set numbered from zero upwards. 
  const_iterator cbegin(std::size_t thread_index = 0) const
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    const_iterator i = m_keys.begin();
    return i;
  }

  /// \brief End of the forward const_iterator. 
  const_iterator cend(std::size_t thread_index = 0) const 
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    const_iterator i = m_keys.cbegin() + m_next_index;
    return i;
  }

  /// \brief Reverse iterator going through the elements in the set from the largest to the smallest index. 
  reverse_iterator rbegin(std::size_t thread_index = 0) 
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    reverse_iterator i = m_keys.rend() - m_next_index;
    return i;
  }

  /// \brief End of the reverse iterator. 
  reverse_iterator rend(std::size_t thread_index = 0)
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    reverse_iterator i = m_keys.rend();
    return i;
  }

  /// \brief Reverse const_iterator going through the elements from the highest to the lowest numbered element. 
  const_reverse_iterator crbegin(std::size_t thread_index = 0) const
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    const_reverse_iterator i=m_keys.crend() - m_next_index;
    return i;
  }

  /// \brief End of the reverse const_iterator. 
  const_reverse_iterator crend(std::size_t thread_index = 0) const 
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    reverse_iterator i = m_keys.crend();
    return i;
  }

  /// \brief Clears the indexed set by removing all its elements. It is not guaranteed that the memory is released too. 
  void clear(std::size_t thread_index=0);

  /// \brief Insert a key in the indexed set and return its index. 
  /// \details If the element was already in the set, the resulting bool is true, and the existing index is returned.
  ///         Otherwise, the key is inserted in the set, and the next available index is assigned to it. 
  /// \param  key The key to be inserted in the set.
  /// \return The index of the key and a boolean indicating whether the element was actually inserted.
  /// \details threadsafe
  std::pair<size_type, bool> insert(const key_type& key, std::size_t thread_index = 0);

  /// \brief Provides an iterator to the stored key in the indexed set.
  /// \param key The key that is sought.
  /// \return An iterator to the key, otherwise end().
  const_iterator find(const key_type& key, std::size_t thread_index = 0) const;

  /// \brief The number of elements in the indexed set.
  /// \return The number of elements in the indexed set. 
  /// \details threadsafe
  size_type size(std::size_t thread_index = 0) const
  { 
    shared_guard guard = m_shared_mutexes[thread_index].lock_shared();
    size_type result=m_next_index;
    return result;
  }
};

} // end namespacee mcrl2::utilities

#include "mcrl2/utilities/detail/indexed_set_implementation.cxx"

#endif // MCRL2_UTILITIES_INDEXED_SET_H
