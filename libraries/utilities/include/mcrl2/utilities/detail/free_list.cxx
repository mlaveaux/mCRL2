// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_DETAIL_FREELIST_H_
#define MCRL2_UTILITIES_DETAIL_FREELIST_H_
MCRL2_MODULE;

#include <cassert>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <iterator>

#ifdef MCRL2_ENABLE_MODULES
  export module utilities:detail.free_list;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities::detail
{

/// \brief This essentially implements the std::forward_list, with the difference that
///        it does not own the nodes in the list. It just keeps track of the list
///        next pointers. Furthermore this assumes that the objects in the list will not
///        be accessed, so it stores the next pointer within the Element.
template<typename Element>
class free_list
{
public:

  /// \brief A sentinel value for the next of a slot to indicate that it belonged to the freelist.
  static constexpr std::uintptr_t FreeListSentinel = std::numeric_limits<std::uintptr_t>::max();

  /// \brief The nodes of the free list without.
  union slot
  {
  public:
    slot() {}
    ~slot() {}

    /// \returns The next bucket in the linked list.
    slot* next() const noexcept { return m_next; }

    /// \returns True if and only there is a next bucket in the linked list.
    bool has_next() const noexcept { return m_next != nullptr; }

    /// \brief Set the next pointer to the given next pointer.
    void next(slot* next) noexcept { m_next = next; }

    /// \returns A reference to the key stored in this node.
    Element& element() noexcept { return m_element; }

    /// \returns An implicit conversion to the underlying key.
    operator Element&() { return m_element; }

    /// \brief Mark this slot with a special value, destroys the slot.
    void mark() { m_next = reinterpret_cast<slot*>(FreeListSentinel); }

    /// \returns True if and only if this slot is marked.
    bool is_marked() const noexcept { return m_next == reinterpret_cast<slot*>(FreeListSentinel); }
  protected:

    /// \brief Pointer to the next node.
    slot* m_next = nullptr;

    /// \brief Store the actual element.
    Element m_element;
  };

  /// \brief Iterator over all keys in a bucket list.
  template<bool Constant = true>
  class slot_iterator : std::iterator_traits<Element>
  {
    struct Sentinel{};

    using tag = std::input_iterator_tag;

    friend class free_list<Element>;
    using Freelist = std::conditional_t<Constant, const free_list<Element>, free_list<Element>>;
    using reference = std::conditional_t<Constant, const Element&, Element&>;
    using slot_pointer = std::conditional_t<Constant, const slot*, slot*>;

  public:
    /// \brief A end of the iterator sentinel.
    static constexpr Sentinel EndIterator{};

    slot_iterator(slot_pointer slot)
      : m_slot(slot)
    {}

    slot_iterator(Freelist& list)
      : m_slot(list.head()->next())
    {}

    slot_iterator()
      : m_slot(nullptr)
    {}

    slot_iterator& operator++()
    {
      m_slot = m_slot->next();
      return *this;
    }

    template <bool Constant_ = Constant>
    reference operator*()
      requires(!Constant_)
    {
      return m_slot->element();
    }

    template <bool Constant_ = Constant>
    reference operator*() const
      requires(Constant_)
    {
      return m_slot->element();
    }

    bool operator != (const slot_iterator& it) const noexcept
    {
      return m_slot != it.m_slot;
    }

    bool operator != (Sentinel) const noexcept
    {
      return m_slot != nullptr;
    }

    slot_pointer get_slot()
    {
      return m_slot;
    }

  private:
    slot_pointer m_slot;
  };

public:
  using Freelist = free_list<Element>;
  using iterator = slot_iterator<false>;
  using const_iterator = slot_iterator<true>;

  free_list() = default;

  /// \returns An iterator over the keys of this bucket and successor buckets.
  iterator begin() { return iterator(head().next()); }
  iterator end() { return iterator(); }

  /// \return An iterator pointing to the before the head of the list.
  iterator before_begin() { return iterator(&head()); }

  /// \returns A const iterator over the keys of this bucket and successor buckets.
  const_iterator begin() const { return const_iterator(*this); }
  const_iterator end() const { return const_iterator(); }

  /// \return A const iterator pointing to the before the head of the list.
  const_iterator before_begin() const { return const_iterator(&m_head); }

  /// \brief Removes the element after the given iterator from the list. The returned iterator
  iterator erase_after(iterator it)
  {
    slot* current_node = it.get_node();
    assert(current_node->next() != nullptr); // Cannot erase after the last node.

    // Keep track of the node that we should remove.
    slot* erased_node = current_node->next();
    slot* next_node = erased_node->next();

    // Update the next pointer of the current node.
    current_node->next(next_node);
    return iterator(next_node);
  }

  /// \brief Empties the bucket list.
  void clear() { head().next(nullptr); }

  /// \returns True if the given pointer is already in the freelist.
  bool contains(Element* pointer)
  {
    union slot* slot = reinterpret_cast<union slot*>(pointer);
    for (auto it = begin(); it != end(); ++it)
    {
      if (slot == *it)
      {
        return true;
      }
    }

    return false;
  }

  /// \returns The length of the freelist.
  std::size_t count()
  {
    std::size_t length = 0;
    for (auto it = begin(); it != end(); ++it)
    {
      ++length;
    }
    return length;
  }

  /// \brief Mark all elements in this list with a special value, destroys the list.
  void destructive_mark()
  {
    for (auto it = begin(); it != end();)
    {
      slot* current_slot = it.get_slot();
      ++it;
      current_slot->mark();
    }

    head().next(nullptr);
  }

  /// \returns True if and only if this list is empty.
  bool empty() const noexcept
  {
    return !head().has_next();
  }

  /// \brief Inserts an element at the position after the iterator.
  iterator insert_after(iterator it, slot& element)
  {
    slot& current = *it.get_slot();
    element.next(current.next());
    current.next(&element);
    return iterator(&element);
  }

  /// \brief Removes the head of the list and returns a reference to this head slot.
  Element& pop_front()
  {
    slot* element = head().next();
    head().next(element->next());
    return element->element();
  }

  /// \brief Puts the given node before the head of the list.
  void push_front(slot& slot)
  {
    slot.next(head().next());
    head().next(&slot);
  }

  /// \brief Erases elements in the range (position, last], constant complexity since no destructors
  ///        will be called.
  void erase_after(iterator position, iterator last)
  {
    position.get_slot()->next(last.get_slot());
  }

private:
  /// \brief Provides access to the head as if it was an actual slot.
  slot& head() noexcept { return m_head; }
  const slot& head() const noexcept { return m_head; }

  /// \brief The first node in the bucket list.
  /// \details This wastes sizeof(Element) space, because slot can store both an Element or a next pointer.
  slot m_head;
};

} // namespace mcrl2::utilities::detail



#endif // MCRL2_UTILITIES_DETAIL_FREELIST_H_
