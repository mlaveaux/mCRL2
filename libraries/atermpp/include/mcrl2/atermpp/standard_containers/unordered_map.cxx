// Author(s): Jan Friso Groote, Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//
/// \file mcrl2/data/standard_containers/unordered_map.h
/// \brief This file contains a unordered_map class that behaves 
///        exactly as a standard unordered_map. It can only be used
///        to store class instances that derive from aterms.
///        The stored aterms are protected as a whole, i.e.,
///        time and memory is saved as individual protection
///        per element is unnecessary. 

#ifndef MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_MAP_H
#define MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_MAP_H
MCRL2_MODULE;

#include <cstdint>
#include <unordered_map>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/aterm_container.cxx"
  #include "mcrl2/utilities/shared_mutex.cxx"
#else
  export module atermpp:standard_containers.unordered_map;

  import :detail.aterm_configuration;
  import :detail.aterm_container;

  import utilities;
#endif

/// \brief The main namespace for the aterm++ library.
MCRL2_MODULE_EXPORT namespace atermpp
{

/// \brief A unordered_map class in which aterms can be stored. 
template < class Key, 
           class T,
           class Hash = std::hash<detail::reference_aterm<Key> >,
           class Pred = std::equal_to<detail::reference_aterm<Key> >,
           class Alloc = std::allocator< std::pair<const detail::reference_aterm<Key>, detail::reference_aterm<T> > > >

class unordered_map : public std::unordered_map< detail::reference_aterm<Key>, detail::reference_aterm<T>, Hash, Pred, Alloc >
{
protected:
  using super = std::unordered_map<detail::reference_aterm<Key>, detail::reference_aterm<T>, Hash, Pred, Alloc>;

  detail::generic_aterm_container<std::unordered_map< detail::reference_aterm<Key>, detail::reference_aterm<T>, Hash, Pred, Alloc > > container_wrapper;

public:
  
  /// Standard typedefs.
  using allocator_type = typename super::allocator_type;
  using value_type = typename super::value_type;
  using size_type = typename super::size_type;
  using node_type = typename super::node_type;
  using reference = typename super::reference;
  using iterator = typename super::iterator;
  using const_iterator = typename super::const_iterator;
  using insert_return_type = typename super::insert_return_type;

  /// \brief Default constructor.
  unordered_map()
   : super(),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  explicit unordered_map (const allocator_type& alloc)
   : super::unordered_map(alloc),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  explicit unordered_map (size_type n, const allocator_type& alloc = allocator_type())
   : super::unordered_map(n, alloc),
     container_wrapper(*this)     
  {}

  unordered_map (size_type n, const value_type& val, const allocator_type& alloc = allocator_type())
   : super::unordered_map(n, detail::reference_aterm(val), alloc),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  template <class InputIterator>
  unordered_map (InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type())
   : super::unordered_map(first, last, alloc),
     container_wrapper(*this)     
  {}
    
  /// \brief Constructor.
  unordered_map (const unordered_map& x)
   : super::unordered_map(x),
     container_wrapper(*this)     
  {}

  /// \brief Constructor.
  unordered_map (const unordered_map& x, const allocator_type& alloc)
   : super::unordered_map(x, alloc),
     container_wrapper(*this)     
  {}
  
  /// \brief Constructor.
  unordered_map(unordered_map&& x) noexcept
      : super::unordered_map(std::move(x)),
        container_wrapper(*this)
  {}


  /// \brief Constructor.
  unordered_map (unordered_map&& x, const allocator_type& alloc)
   : super::unordered_map(std::move(x), alloc),
     container_wrapper(*this)
  {}

  /// \brief Constructor. To be done later....
  unordered_map (std::initializer_list<value_type> il, const allocator_type& alloc = allocator_type())
    : super::unordered_map(il, alloc),
      container_wrapper(*this)      
  {}

  /// \brief Standard assignment.
  unordered_map& operator=(const unordered_map& other)=default;

  /// \brief Standard move assignment.
  unordered_map& operator=(unordered_map&& other)=default;

  /// \brief Standard destructor.
  ~unordered_map()=default;

  void clear() noexcept;

  /// \brief Inserts an element referring to a default value in the map. 
  std::pair<iterator,bool> insert( const value_type& value );

  template< class P >
  std::pair<iterator,bool> insert( P&& value );

  iterator insert( const_iterator hint, const value_type& value );

  template< class P >
  iterator insert( const_iterator hint, P&& value );

  template< class InputIt >
  void insert( InputIt first, InputIt last );

  void insert( std::initializer_list<value_type> ilist );

  insert_return_type insert( node_type&& nh );

  iterator insert( const_iterator hint, node_type&& nh );

  template <class M>
  std::pair<iterator, bool> insert_or_assign( const Key& k, M&& obj );

  template <class M>
  std::pair<iterator, bool> insert_or_assign( Key&& k, M&& obj );

  template <class M>
  iterator insert_or_assign( const_iterator hint, const Key& k, M&& obj );

  template <class M>
  iterator insert_or_assign( const_iterator hint, Key&& k, M&& obj );

  template< class... Args >
  std::pair<iterator,bool> emplace( Args&&... args );

  template <class... Args>
  iterator emplace_hint( const_iterator hint, Args&&... args );

  template< class... Args >
  std::pair<iterator, bool> try_emplace( const Key& k, Args&&... args );

  template< class... Args >
  std::pair<iterator, bool> try_emplace( Key&& k, Args&&... args );

  template< class... Args >
  iterator try_emplace( const_iterator hint, const Key& k, Args&&... args );

  template< class... Args >
  iterator try_emplace( const_iterator hint, Key&& k, Args&&... args );

  iterator erase( iterator pos );

  iterator erase( const_iterator pos );
  
  iterator erase( const_iterator first, const_iterator last );

  size_type erase( const Key& key );

  void swap(unordered_map& other) noexcept;

  node_type extract( const_iterator position );
  
  node_type extract( const Key& k );

  template<class H2, class P2>
  void merge( std::unordered_map<Key, T, H2, P2, allocator_type>& source );

  template<class H2, class P2>
  void merge( std::unordered_map<Key, T, H2, P2, allocator_type>&& source );

  template<class H2, class P2>
  void merge( std::unordered_multimap<Key, T, H2, P2, allocator_type>& source );

  template<class H2, class P2>
  void merge( std::unordered_multimap<Key, T, H2, P2, allocator_type>&& source );

  std::size_t size() const
  {
    return super::size();
  }

  /// \returns An iterator over all keys.
  iterator begin();
  iterator end();

  /// \returns A const iterator over all keys.
  const_iterator begin() const;
  const_iterator end() const;

  /// \returns A const iterator over all keys.
  const_iterator cbegin() const;
  const_iterator cend() const;

  /// \returns True iff the set is empty.
  bool empty() const noexcept;

  /// \returns The amount of elements stored in this set.
  size_type max_size() const noexcept;
};

namespace utilities
{
/// \brief A unordered_map class in which aterms can be stored. 
template < class Key,
  class T,
  class Hash = std::hash<detail::reference_aterm<Key> >,
  class Pred = std::equal_to<detail::reference_aterm<Key> >,
  class Alloc = std::allocator< std::pair<const detail::reference_aterm<Key>, detail::reference_aterm<T> > >,
  bool ThreadSafe = false >

class unordered_map : public mcrl2::utilities::unordered_map< detail::reference_aterm<Key>, 
                                                              detail::reference_aterm<T>, Hash, Pred, Alloc, ThreadSafe, false >
{
  protected:
    using super = mcrl2::utilities::
        unordered_map<detail::reference_aterm<Key>, detail::reference_aterm<T>, Hash, Pred, Alloc, ThreadSafe, false>;

    detail::generic_aterm_container<mcrl2::utilities::unordered_map< detail::reference_aterm<Key>, detail::reference_aterm<T>, Hash, Pred, Alloc, ThreadSafe, false > > container_wrapper;

  public:

    /// Standard typedefs.
    using allocator_type = typename super::allocator_type;
    using value_type = typename super::value_type;
    using size_type = typename super::size_type;
    using reference = typename super::reference;
    using iterator = typename super::iterator;
    using const_iterator = typename super::const_iterator;

    /// \brief Default constructor.
    unordered_map()
      : super(),
      container_wrapper(*this)
    {}

    /// \brief Constructor.
    explicit unordered_map(const allocator_type& alloc)
      : super::unordered_map(alloc),
      container_wrapper(*this)
    {}

    /// \brief Constructor.
    explicit unordered_map(size_type n, const allocator_type& alloc = allocator_type())
      : super::unordered_map(n, alloc),
      container_wrapper(*this)
    {}

    unordered_map(size_type n, const value_type& val, const allocator_type& alloc = allocator_type())
      : super::unordered_map(n, detail::reference_aterm(val), alloc),
      container_wrapper(*this)
    {}

    /// \brief Constructor.
    template <class InputIterator>
    unordered_map(InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type())
      : super::unordered_map(first, last, alloc),
      container_wrapper(*this)
    {}

    /// \brief Constructor.
    unordered_map(const unordered_map& x)
      : super::unordered_map(x),
      container_wrapper(*this)
    {}

    /// \brief Constructor.
    unordered_map(const unordered_map& x, const allocator_type& alloc)
      : super::unordered_map(x, alloc),
      container_wrapper(*this)
    {}

    /// \brief Constructor.
    unordered_map(unordered_map&& x) noexcept
        : super::unordered_map(std::move(x)),
          container_wrapper(*this)
    {}

    /// \brief Constructor.
    unordered_map(unordered_map&& x, const allocator_type& alloc)
      : super::unordered_map(std::move(x), alloc),
      container_wrapper(*this)
    {}

    /// \brief Constructor. To be done later....
    unordered_map(std::initializer_list<value_type> il, const allocator_type& alloc = allocator_type())
      : super::unordered_map(il, alloc),
      container_wrapper(*this)
    {}

    /// \brief Standard assignment.
    unordered_map& operator=(const unordered_map& other) = default;

    /// \brief Standard move assignment.
    unordered_map& operator=(unordered_map&& other) = default;

    /// \brief Standard destructor.
    ~unordered_map() = default;

    void clear() noexcept;

    /// \brief Standard find function in a map.
    /// \returns Element with the specified key.
    template<typename ...Args>
    iterator find(const Args&... args);

    /// \brief Standard find function in a map.
    /// \returns Element with the specified key.
    template<typename ...Args>
    const_iterator find(const Args&... args) const;

    std::pair<iterator, bool> insert(const value_type& value);

    template< class P >
    std::pair<iterator, bool> insert(P&& value);

    iterator insert(const_iterator hint, value_type&& value);

    template< class P >
    iterator insert(const_iterator hint, P&& value);

    template< class InputIt >
    void insert(InputIt first, InputIt last);

    void insert(std::initializer_list<value_type> ilist);

    /*
    insert_return_type insert(node_type&& nh);

    iterator insert(const_iterator hint, node_type&& nh);
    */
    
    template <class M>
    std::pair<iterator, bool> insert_or_assign(const Key& k, M&& obj);

    template <class M>
    std::pair<iterator, bool> insert_or_assign(Key&& k, M&& obj);

    template <class M>
    iterator insert_or_assign(const_iterator hint, const Key& k, M&& obj);

    template <class M>
    iterator insert_or_assign(const_iterator hint, Key&& k, M&& obj);

    template< class... Args >
    std::pair<iterator, bool> emplace(Args&&... args);

    template <class... Args>
    iterator emplace_hint(const_iterator hint, Args&&... args);

    template< class... Args >
    std::pair<iterator, bool> try_emplace(const Key& k, Args&&... args);

    template< class... Args >
    std::pair<iterator, bool> try_emplace(Key&& k, Args&&... args);

    template< class... Args >
    iterator try_emplace(const_iterator hint, const Key& k, Args&&... args);

    template< class... Args >
    iterator try_emplace(const_iterator hint, Key&& k, Args&&... args);

    iterator erase(iterator pos);

    iterator erase(const_iterator pos);

    iterator erase(const_iterator first, const_iterator last);

    size_type erase(const Key& key);

    void swap(unordered_map& other) noexcept;

    /*
    node_type extract(const_iterator position);

    node_type extract(const Key& k);

    template<class H2, class P2>
    void merge(std::unordered_map<Key, T, H2, P2, allocator_type>& source);

    template<class H2, class P2>
    void merge(std::unordered_map<Key, T, H2, P2, allocator_type>&& source);

    template<class H2, class P2>
    void merge(std::unordered_multimap<Key, T, H2, P2, allocator_type>& source);

    template<class H2, class P2>
    void merge(std::unordered_multimap<Key, T, H2, P2, allocator_type>&& source);
    */

    std::size_t size() const
    {
      return super::size();
    }

    /// \returns An iterator over all keys.
    iterator begin();
    iterator end();

    /// \returns A const iterator over all keys.
    const_iterator begin() const;
    const_iterator end() const;

    /// \returns A const iterator over all keys.
    const_iterator cbegin() const;
    const_iterator cend() const;

    /// \returns True iff the set is empty.
    bool empty() const noexcept;

    /// \returns The amount of elements stored in this set.
    size_type max_size() const noexcept;

  protected:

    /// Function below is implemented in a .cpp file. 
    void rehash(std::size_t /* new_size */);
  
    void rehash_if_needed();
};
}

} // namespace atermpp
#endif // MCRL2_ATERMPP_STANDARD_CONTAINER_UNORDERED_MAP_H


