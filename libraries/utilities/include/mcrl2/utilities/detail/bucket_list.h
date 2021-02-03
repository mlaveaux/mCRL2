// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_DETAIL_BUCKETLIST_H_
#define MCRL2_UTILITIES_DETAIL_BUCKETLIST_H_

#include <cassert>
#include <atomic>
#include <cstddef>
#include <iterator>
#include <memory>

namespace mcrl2
{
namespace utilities
{
namespace detail
{

struct Sentinel{};

/// \brief A end of the iterator sentinel.
static constexpr Sentinel EndIterator = {};

/// \brief A compile time check for allocate_args in the given allocator, calls allocate(1) otherwise.
template<typename Allocator, typename ...Args>
inline auto allocate(Allocator& allocator, const Args&... args) -> decltype(allocator.allocate_args(args...))
{
  return allocator.allocate_args(args...);
}

template<typename Allocator, typename ...Args>
inline auto allocate(Allocator& allocator, const Args&...) -> decltype(allocator.allocate(1))
{
  return allocator.allocate(1);
}

/// \brief This essentially implements the std::forward_list, with the difference that
///        it does not own the nodes in the list. It just keeps track of the list
///        next pointers.
template<typename Key, typename Allocator>
class bucket_list
{
public:

  /// \brief The nodes of the bucket list without carrying any additional informations.
  ///        Used to make no different between the head and the tail of the list.
  class node_base
  {
  public:
    node_base() = default;
    node_base(node_base&& other)
    {
      m_next.store(other.m_next.load());
      other.m_next = nullptr;
    }

    /// \returns The next bucket in the linked list.
    node_base* next() const noexcept { return m_next.load(std::memory_order_relaxed); }

    /// \returns True if and only there is a next bucket in the linked list.
    bool has_next() const noexcept { return m_next.load(std::memory_order_relaxed) != nullptr; }

    /// \brief Set the next pointer to the given next pointer.
    void set_next(node_base* next) noexcept { m_next.store(next, std::memory_order_relaxed); }
    /// \returns True iff next has been replaced by value iff next is equal to expected.
    bool exchange(node_base*& expected, node_base* value) { return m_next.compare_exchange_weak(expected, value, std::memory_order_release, std::memory_order_relaxed); }
  protected:

    /// \brief Pointer to the next node.
    std::atomic<node_base*> m_next = nullptr;
  };

  /// \brief The nodes of the bucket list.
  class node final : public node_base
  {
  public:

    /// \brief Constructs a key by using the given arguments.
    template<typename ...Args>
    explicit node(Args&&... args)
      : m_key(std::forward<Args>(args)...)
    {}

    /// \returns A reference to the key stored in this node.
    Key& key() noexcept { return m_key; }
    const Key& key() const noexcept { return m_key; }

    /// \returns An implicit conversion to the underlying key.
    operator Key&() { return m_key; }
  private:
    /// \brief Store the actual key.
    Key m_key;
  };

  /// \brief Iterator over all keys in a bucket list.
  template<bool Constant = true>
  class key_iterator
  {
    friend class bucket_list;

  public:
    using value_type = typename std::conditional<Constant, const Key, Key>::type;
    using reference = typename std::conditional<Constant, const Key&, Key&>::type;
    using pointer = typename std::conditional<Constant, const Key*, Key*>::type;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    /// Default constructor.
    key_iterator() {}

    /// \brief Implicit conversion to const_iterator.
    operator key_iterator<true>() const
    {
      return key_iterator<true>(m_current_node);
    }

    key_iterator& operator++()
    {
      m_current_node = reinterpret_cast<node*>(m_current_node->next());
      return *this;
    }

    key_iterator operator++(int)
    {
      key_iterator copy(*this);
      ++(*this);
      return copy;
    }

    reference operator*() const
    {
      return m_current_node->key();
    }

    pointer operator->() const
    {
      return &m_current_node->key();
    }

    bool operator== (const key_iterator& it) const noexcept
    {
      return !(*this != it);
    }

    bool operator!= (const key_iterator& it) const noexcept
    {
      return m_current_node != it.m_current_node;
    }

    bool operator== (Sentinel) const noexcept
    {
      return !(*this != Sentinel());
    }

    bool operator!= (Sentinel) const noexcept
    {
      return m_current_node != nullptr;
    }

  private:
    explicit key_iterator(node_base* current)
      : m_current_node(reinterpret_cast<node*>(current))
    {}

    node* get_node() noexcept
    {
      return m_current_node;
    }

    const node* get_node() const noexcept
    {
      return m_current_node;
    }

    node* m_current_node = nullptr;
  };

public:
  using Bucket = bucket_list<Key, Allocator>;
  using iterator = key_iterator<false>;
  using const_iterator = key_iterator<true>;

  /// Rebind the passed to allocator to a bucket list node allocator.
  using NodeAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Bucket::node>;

  /// \returns The first element.
  Key& front() { return reinterpret_cast<node&>(*m_head.next()).key(); }
  const Key& front() const { return reinterpret_cast<const node&>(*m_head.next()).key(); }

  /// \returns An iterator over the keys of this bucket and successor buckets.
  iterator begin() { return iterator(m_head.next()); }
  iterator end() { return iterator(); }

  /// \return An iterator pointing to the before the head of the list.
  iterator before_begin() { return iterator(&m_head); }

  /// \returns A const iterator over the keys of this bucket and successor buckets.
  const_iterator begin() const { return const_iterator(m_head.next()); }
  const_iterator end() const { return const_iterator(); }

  const_iterator cbegin() const { return const_iterator(m_head.next()); }
  const_iterator cend() const { return const_iterator(); }

  /// \return A const iterator pointing to the before the head of the list.
  const_iterator before_begin() const { return const_iterator(const_cast<node_base*>(&m_head)); }

  /// \returns True iff this bucket has no elements.
  bool empty() const { return !m_head.has_next(); }

  /// \brief Empties the bucket list.
  void clear(NodeAllocator& allocator)
  {
    while(!empty())
    {
      erase_after(before_begin(), allocator);
    }
  }

  /// \brief Constructs an element using the allocator with the given arguments and insert it in the front.
  template<typename ...Args>
  void emplace_front(NodeAllocator& allocator, Args&& ...args)
  {
    // Allocate a new node.
    node* new_node = allocate(allocator, std::forward<Args>(args)...);
    std::allocator_traits<NodeAllocator>::construct(allocator, new_node, std::forward<Args>(args)...);

    // Ensure that the previous front is set behind this node.
    new_node->set_next(m_head.next());

    // Change the head to the newly allocated node.
    m_head.set_next(new_node);
  }

  /// \brief Constructs an element using the allocator with the given arguments and insert it in the front of the bucket iff it does not already exist.
  /// \returns True iff the insertion took place.
  /// \threadsafe
  template<typename ...Args,
           typename Equals>
  std::pair<iterator, bool> emplace_front_unique(NodeAllocator& allocator, const Equals& equals, Args&& ...args)
  {
    // Allocate a new node.
    node* new_node = allocate(allocator, std::forward<Args>(args)...);
    std::allocator_traits<NodeAllocator>::construct(allocator, new_node, std::forward<Args>(args)...);

    // This was the first and last when we started the operation.
    node_base* old_head = m_head.next();
    iterator old_end;

    // Change the head to the newly allocated node.
    do
    {
      // Ensure that the previous front is set behind this node.
      new_node->set_next(old_head);

      // Check whether the new node is not already contained in the bucket list.
      for (auto it = iterator(old_head); it != old_end; ++it)
      {
        if (equals(*it, std::forward<Args>(args)...))
        {
          // Clean up new node and leave bucket as is.          
          std::allocator_traits<NodeAllocator>::destroy(allocator, new_node);
          std::allocator_traits<NodeAllocator>::deallocate(allocator, new_node, 1);
          return std::make_pair(it, false);
        }
      }

      // Next iteration we only have to look at newly inserted nodes.
      old_end = iterator(old_head);
    }
    while(!m_head.exchange(old_head, new_node));

    return std::make_pair(iterator(new_node), true);
  }

  /// \brief Removes the element after the given iterator from the list. The returned iterator
  iterator erase_after(NodeAllocator& allocator, const_iterator it)
  {
    node* current_node = const_cast<node*>(it.get_node());
    assert(current_node->next() != nullptr); // Cannot erase after the last node.

    // Keep track of the node that we should remove.
    node* erased_node = reinterpret_cast<node*>(current_node->next());
    node_base* next_node = erased_node->next();

    // Update the next pointer of the current node.
    current_node->set_next(next_node);

    // Clean up the old node.
    std::allocator_traits<NodeAllocator>::destroy(allocator, erased_node);
    std::allocator_traits<NodeAllocator>::deallocate(allocator, erased_node, 1);

    return iterator(next_node);
  }

  /// \brief Moves the elements from other into this bucket after the given position.
  void splice_after(const_iterator pos, bucket_list& other)
  {
    if (!other.empty())
    {
      assert(begin() != other.begin());
      node* current_node = const_cast<node*>(pos.get_node());

      // Increment the position now as we are going to change the current_node.
      ++pos;

      // Sets the current position to be followed by the other bucket list.
      current_node->set_next(other.m_head.next());

      node* after_pos_node = const_cast<node*>(pos.get_node());
      if (after_pos_node != nullptr)
      {
        // Find the last node of the other list.
        iterator before_end;
        for (iterator it = other.begin(); it != other.end(); ++it)
        {
          before_end = it;
        }

        // Set the last element of other to the position after pos.
        node* last = before_end.get_node();
        last->set_next(after_pos_node);
      }

      // Make the other bucket empty, all elements belong to this bucket now.
      other.m_head.set_next(nullptr);
    }
  }

  /// \brief Moves the first node from the given bucket into this bucket after the given position.
  /// \details This is a non-standard addition to ensure efficient splicing (instead of splice_after(pos, other, other.begin(), ++other.begin()).
  void splice_front(const_iterator pos, bucket_list& other)
  {
    if (!other.empty())
    {
      // Sets the current position to be followed by the other bucket list.
      node* current_node = const_cast<node*>(pos.get_node());
      node_base* next_node = current_node->next();

      // Make the other head node the start of our bucket after pos.
      node_base* head_node = other.begin().get_node();
      current_node->set_next(head_node);

      // Make the next position the position after the head, but keep track of the current next node.
      node_base* next_head_node = head_node->next();
      head_node->set_next(next_node);

      // Make the other head point to the next element.
      other.m_head.set_next(next_head_node);
    }
  }

private:

  /// \returns True iff this bucket already contains the given node.
  bool contains(node* node)
  {
    node_base* element = &m_head;
    while(element->has_next())
    {
      if (node == element->next())
      {
        return true;
      }

      element = element->next();
    }

    return false;
  }

  /// \brief The first node in the bucket list.
  node_base m_head;
};

} // namespace detail
} // namespace utilities
} // namespace mcrl2

#endif // MCRL2_UTILITIES_DETAIL_BUCKETLIST_H_
