// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_INDEXED_SET_H
#define MCRL2_UTILITIES_INDEXED_SET_H

#include <deque>

#include "mcrl2/utilities/unordered_map.h"
#include "mcrl2/utilities/copyable_mutex.h"

namespace mcrl2
{
namespace utilities
{

/// \brief A set that assigns each element an unique index.
template<typename Key,
         typename Hash = std::hash<Key>,
         typename Equals = std::equal_to<Key>,
         typename Allocator = std::allocator<Key>,
         bool ThreadSafe = false>
class indexed_set
{
private:
  std::vector<std::size_t> m_hashtable;
  std::deque<Key, Allocator> m_keys;

  Hash m_hasher;
  Equals m_equals;

  /// \brief Mutex for the m_hashtable and m_keys data structures.
  mutable mcrl2::utilities::copyable_mutex m_mutex;

  /// \brief Inserts the given (key, n) pair into the indexed set.
  std::size_t put_in_hashtable(const Key& key, std::size_t n);

  /// \brief Resizes the hash table to twice its current size.
  inline void resize_hashtable();


public:
  typedef Key key_type;
  typedef std::size_t size_type;
  typedef std::pair<const key_type, size_type> value_type;
  typedef Equals key_equal;
  typedef Hash hasher;

  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;

  typedef typename std::deque < Key >::iterator iterator;
  typedef typename std::deque < Key >::const_iterator const_iterator;

  typedef typename std::deque < Key >::reverse_iterator reverse_iterator;
  typedef typename std::deque < Key >::const_reverse_iterator const_reverse_iterator;

  typedef std::ptrdiff_t difference_type;
  
  /// \brief Value returned when an element does not exist in the set.
  /// \return Value indicating non existing element, equal to std::numeric_limits<std::size_t>::max(). 
  static constexpr size_type npos = std::numeric_limits<std::size_t>::max();

  /// \brief Constructor of an empty indexed set. Starts with a hashtable of size 128.
  indexed_set();

  /// \brief Constructor of an empty index set. Starts with a hashtable of the indicated size. 
  /// \param initial_hashtable_size The initial size of the hashtable.
  /// \param hash The hash function.
  /// \param equals The comparison function for its elements.
  indexed_set(std::size_t initial_hashtable_size,
    const hasher& hash = hasher(),
    const key_equal& equals = key_equal());

  /// \brief Returns a reference to the mapped value.
  /// \details Returns an invalid value, larger or equal than the size of the indexed set, if there is no element with the given key.
  /// \threadsafe
  size_type index(const key_type& key) const;

  /// \brief Returns a reference to the mapped value.
  /// \details Returns an out_of_range exception if there is no element with the given key.
  /// \param index The position in the indexed set.
  /// \return The value at position index.
  const key_type& at(const size_type index) const;

  /// \brief Operator that provides a const reference at the position indicated by index.
  /// \param index The position in the indexed set.
  /// \return The value at position index.
  const key_type& operator[](const size_type index) const;

  /// \brief Forward iterator which runs through the elements from the lowest to the largest number.
  /// \details Complexity is constant per operation.
  iterator begin() 
  { 
    return m_keys.begin(); 
  }

  /// \brief End of the forward iterator.
  iterator end()
  { 
    return m_keys.end(); 
  }

  /// \brief Forward iterator which runs through the elements from the lowest to the largest number.
  /// \details Complexity is constant per operation.
  const_iterator begin() const
  {
    return m_keys.begin();
  }

  /// \brief End of the forward iterator.
  const_iterator end() const
  {
    return m_keys.end();
  }

  /// \brief const_iterator going through the elements in the set numbered from zero upwards. 
  const_iterator cbegin() const
  { 
    return m_keys.cbegin(); 
  }

  /// \brief End of the forward const_iterator. 
  const_iterator cend() const 
  { 
    return m_keys.cend(); 
  }

  /// \brief Reverse iterator going through the elements in the set from the largest to the smallest index. 
  iterator rbegin() 
  { 
    return m_keys.rbegin(); 
  }

  /// \brief End of the reverse iterator. 
  iterator rend()
  { 
    return m_keys.rend(); 
  }

  /// \brief Reverse const_iterator going through the elements from the highest to the lowest numbered element. 
  const_iterator crbegin() const
  { 
    return m_keys.crbegin(); 
  }

  /// \brief End of the reverse const_iterator. 
  const_iterator crend() const 
  { 
    return m_keys.crend(); 
  }

  /// \brief Clears the indexed set by removing all its elements. It is not guaranteed that the memory is released too. 
  void clear();

  /// \brief Insert a key in the indexed set and return its index. 
  /// \details If the element was already in the set, the resulting bool is true, and the existing index is returned.
  ///         Otherwise, the key is inserted in the set, and the next available index is assigned to it. 
  /// \param  key The key to be inserted in the set.
  /// \return The index of the key and a boolean indicating whether the element was actually inserted.
  /// \threadsafe
  std::pair<size_type, bool> insert(const key_type& key);

  /// \brief Provides an iterator to the stored key in the indexed set.
  /// \param key The key that is sought.
  /// \return An iterator to the key, otherwise end().
  /// \threadsafe
  const_iterator find(const key_type& key) const;

  /// \brief The number of elements in the indexed set.
  /// \return The number of elements in the indexed set. 
  size_type size() const
  { 
    return m_keys.size();
  }
};

} // end namespace utilities
} // end namespace mcrl2

#include "mcrl2/utilities/detail/indexed_set.h"


#endif // MCRL2_UTILITIES_INDEXED_SET_H
