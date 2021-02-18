// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_HASHTABLE_H
#define MCRL2_UTILITIES_HASHTABLE_H

#include <vector>

#include "mcrl2/utilities/unordered_map.h"

namespace mcrl2
{
namespace utilities
{

/// \brief A set that assigns each element an unique index.
template<typename Key,
         typename Hash = std::hash<Key>,
         typename Equals = std::equal_to<Key>,
         typename Allocator = std::allocator<Key>>
class hashtable
{
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

  typedef typename std::vector < Key >::iterator iterator;
  typedef typename std::vector < Key >::const_iterator const_iterator;

  typedef typename std::vector < Key >::reverse_iterator reverse_iterator;
  typedef typename std::vector < Key >::const_reverse_iterator const_reverse_iterator;

  typedef std::ptrdiff_t difference_type;

  /// \brief Constructor of an empty indexed set. Starts with a hashtable of size 128.
  hashtable();

  /// \brief Constructor of an empty index set. Starts with a hashtable of the indicated size. 
  /// \param initial_hashtable_size The initial size of the hashtable.
  /// \param hash The hash function.
  /// \param equals The comparison function for its elements.
  hashtable(std::size_t initial_hashtable_size,
    const hasher& hash = hasher(),
    const key_equal& equals = key_equal());

  /// \brief Forward iterator which runs through the elements from the lowest to the largest number.
  /// \details Complexity is constant per operation.
  iterator begin() 
  { 
    return m_hashtable.begin();
  }

  /// \brief End of the forward iterator.
  iterator end()
  { 
    return m_hashtable.end();
  }

  /// \brief Forward iterator which runs through the elements from the lowest to the largest number.
  /// \details Complexity is constant per operation.
  iterator begin() const
  {
    return m_hashtable.begin();
  }

  /// \brief End of the forward iterator.
  iterator end() const
  {
    return m_hashtable.end();
  }

  /// \brief const_iterator going through the elements in the set numbered from zero upwards. 
  const_iterator cbegin() const
  { 
    return m_hashtable.cbegin();
  }

  /// \brief End of the forward const_iterator. 
  const_iterator cend() const 
  { 
    return m_hashtable.cend();
  }

  /// \brief Reverse iterator going through the elements in the set from the largest to the smallest index. 
  iterator rbegin() 
  { 
    return m_hashtable.rbegin();
  }

  /// \brief End of the reverse iterator. 
  iterator rend()
  { 
    return m_hashtable.rend();
  }

  /// \brief Reverse const_iterator going through the elements from the highest to the lowest numbered element. 
  const_iterator crbegin() const
  { 
    return m_hashtable.crbegin();
  }

  /// \brief End of the reverse const_iterator. 
  const_iterator crend() const 
  { 
    return m_hashtable.crend();
  }

  /// \brief Clears the indexed set by removing all its elements. It is not guaranteed that the memory is released too. 
  void clear();

  /// \brief Insert a key in the indexed set and return its index. 
  /// \details If the element was already in the set, the resulting bool is true, and the existing index is returned.
  ///         Otherwise, the key is inserted in the set, and the next available index is assigned to it. 
  /// \param  key The key to be inserted in the set.
  /// \return The index of the key and a boolean indicating whether the element was actually inserted. 
  std::pair<iterator, bool> insert(const key_type& key);

  /// \brief Erases the key assuming that this key is present in the hashtable..
  /// \returns An iterator to the next key.
  iterator erase(const key_type& key);

  /// \brief The number of elements in the indexed set.
  /// \return The number of elements in the indexed set. 
  size_type size() const
  { 
    return m_number_of_elements;
  }

private:
  /// \return The index where this key should be stored.
  std::size_t get_index(const key_type& key);

  std::vector<Key> m_hashtable;

  Hash m_hasher;
  Equals m_equals;

  std::size_t m_number_of_elements = 0;

  /// \brief Always equal to m_hashtable.size() - 1.
  size_type m_buckets_mask;

  /// \brief Resizes the hash table to the given power of two size.
  inline void rehash(std::size_t size);
};

} // end namespace utilities
} // end namespace mcrl2

#include "mcrl2/utilities/detail/hashtable.h"


#endif // MCRL2_UTILITIES_HASHTABLE_H
