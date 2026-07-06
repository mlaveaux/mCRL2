// Author(s): David N. Jansen, Institute of Software, Chinese Academy of
// Sciences, Beijing, PR China
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/// \file lts/detail/simple_list.h
///
/// \brief Simple list implementation (with pool allocator)
///
/// \details This file supports partition refinement algorithms by implementing
/// a simple list data structure, possibly with a pool allocator.
///
/// The main difference to std::list<T> is that the simple list only stores
/// a pointer to the first list element; there is no sentinel.  The list
/// elements themselves form a doubly-linked list that is almost circular; only
/// the pointer from the last element to the next one is nullptr.  This allows
/// to find the last element in the list (namely as predecessor of the first).
/// A disadvantage is that the list does not allow to find the last element
/// through `std::prev(end())`; to alleviate this, we offer a separate function
/// `before_end()`.
///
/// Also, a simple list requires the stored elements to be trivially
/// destructible, so they can be stored in a pool allocator data structure.
///
/// The pool allocator uses larger blocks to allocate data items,
/// mostly list elements, but also other items can be stored, as long
/// as they are trivially destructible.  To access the pool, one uses
/// `simple_list<element type>::get_pool()`.  However, only elements of the
/// same size as `simple_list<element type>` entries can be deleted.  This is
/// a static pool that is available throughout the running time of the program.
/// (If you are really concerned about shaving off every small bit of time,
/// you might destroy the pool through a memory leak upon termination of the
/// program -- to that end, modify the destructor of `my_pool` to do nothing.)
///
/// \author David N. Jansen, Institute of Software, Chinese Academy of
/// Sciences, Beijing, PR China

#ifndef SIMPLE_LIST_H
#define SIMPLE_LIST_H

#include <cstddef>     // for std::size_t
#include <new>         // for placement new
#include <type_traits> // for std::is_trivially_destructible<class>
#include <functional>  // for std::less
#include <iterator>    // for std::forward_iterator_tag

// My provisional recommendation is to always use simple lists and pool
// allocators.  Using standard allocation and standard lists is 5-15% slower
// and uses perhaps 0.7% more memory.  Using standard allocation and simple
// lists is 10-20% slower and has no significant effect on memory use.  These
// numbers are based on a small set with not-so-large case studies for the JGKW
// bisimulation minimisation algorithm, none of which includes new bottom
// states.

#define USE_SIMPLE_LIST

#ifndef USE_SIMPLE_LIST
    #include <list>
#endif

#define USE_POOL_ALLOCATOR

namespace mcrl2::lts::detail
{





/* ************************************************************************* */
/*                                                                           */
/*                     M E M O R Y   M A N A G E M E N T                     */
/*                                                                           */
/* ************************************************************************* */





#ifdef USE_POOL_ALLOCATOR
    /// \cond INTERNAL_DOCS
    #define ONLY_IF_POOL_ALLOCATOR(...) __VA_ARGS__ // NOLINT(cppcoreguidelines-macro-usage)
    #ifndef USE_SIMPLE_LIST
        #error "Using the pool allocator also requires using the simple list"
    #endif

    /// \class pool
    /// \brief a pool allocator class
    /// \details This class allocates a large chunk of memory at once and hands
    /// out smaller parts.  It is supposed to be more efficient than calling
    /// new/delete, in particular because it assumes that T is trivially
    /// destructible, so it won't call destructors.  It allows to store
    /// elements of different sizes.
    ///
    /// Internally, it keeps a (single-linked) list of large chunks of size
    /// NR_ELEMENTS*sizeof(T)+sizeof(pointer). Each chunk contains a data area; for all chunks except the
    /// first one, this area is completely in use.
    ///
    /// There is a free list, a (single-linked) list of elements in the chunks
    /// that have been freed.  However, all elements in the free list have to
    /// have the same size as type T.
    template <class T, std::size_t NR_ELEMENTS = 4000>
    class my_pool
    {
      static_assert(std::is_trivially_destructible_v<T>);
    private:
      static_assert(sizeof(void*) <= sizeof(T));
      class pool_block_t
      {
      public:
        char data[NR_ELEMENTS * sizeof(T)]{}; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays) raw storage for the pool allocator
        pool_block_t* next_block;

        pool_block_t(pool_block_t* const new_next_block)
            : next_block(new_next_block)
        {}
        };                                                                      static_assert(sizeof(T) <= sizeof(pool_block_t::data));

        /// \brief first chunk in list of chunks
        /// \details All chunks except the first one are completely in use.
        pool_block_t* first_block;

        /// \brief start of part in the first chunk that is already in use
        void* begin_used_in_first_block;

        /// \brief first freed element
        void* first_free_T = nullptr;

        static void*& deref_void(void* addr)
        {
            return *static_cast<void**>(addr);
        }
      public:
        /// \brief constructor
        my_pool()
          : first_block(new pool_block_t(nullptr)),
            begin_used_in_first_block(
                                &first_block->data[sizeof(first_block->data)])
        {}


        /// \brief destructor
        ~my_pool()
        {
            pool_block_t* block(first_block);                                   assert(nullptr != block);
            do
            {
                pool_block_t* next_block = block->next_block;
                delete block;
                block = next_block;
            }
            while(nullptr != block);
        }


      private:
        /// \brief allocate and construct a new element of the same size as the
        /// free list
        template <class U, class... Args>
        U* construct_samesize(Args&&... args)
        {                                                                       static_assert(sizeof(T) == sizeof(U));
            void* new_el;                                                       assert(nullptr != first_block);
            if (first_block->data + sizeof(U) <= begin_used_in_first_block)
            {
                begin_used_in_first_block =
                    new_el = static_cast<char*>(begin_used_in_first_block) -
                                                                     sizeof(U);
            }
            else if (nullptr != first_free_T)
            {
                // free list is tested afterwards because I expect that there
                // won't be too many elements in the free list.
                new_el = first_free_T;
                first_free_T = deref_void(new_el);
            }
            else
            {
                first_block = new pool_block_t(first_block);
                begin_used_in_first_block =
                    new_el = &first_block->data[sizeof(first_block->data) -
                                                                    sizeof(U)];
            }
            return new(new_el) U(std::forward<Args>(args)...);
        }


        /// \brief allocate and construct a new element of a size that doesn't
        /// fit the free list
        template <class U, class... Args>
        U* construct_othersize(Args&&... args)
        {                                                                       static_assert(sizeof(U) != sizeof(T));
            void* new_el;                                                       assert(nullptr != first_block);
            if (first_block->data + sizeof(U) <= begin_used_in_first_block)
            {
                begin_used_in_first_block =
                    new_el = static_cast<char*>(begin_used_in_first_block) -
                                                                     sizeof(U);
            }
            else
            {
                if constexpr (sizeof(T) * 2 < sizeof(U))
                {
                    // There may be space for several T-elements,
                    // extend the free list accordingly
                    while (first_block->data + sizeof(T) <=
                                                     begin_used_in_first_block)
                    {
                        begin_used_in_first_block = static_cast<char*>
                                       (begin_used_in_first_block) - sizeof(T);
                        deref_void(begin_used_in_first_block) = first_free_T;
                        first_free_T = begin_used_in_first_block;
                    }
                }
                else if constexpr (sizeof(T) < sizeof(U))
                {
                    // There may be space for one T-element (but not more),
                    // extend the free list accordingly
                    if (first_block->data + sizeof(T) <=
                                                     begin_used_in_first_block)
                    {
                        begin_used_in_first_block = static_cast<char*>
                                       (begin_used_in_first_block) - sizeof(T);
                        deref_void(begin_used_in_first_block) = first_free_T;
                        first_free_T = begin_used_in_first_block;
                    }
                }                                                               assert(first_block->data + sizeof(T) > begin_used_in_first_block);
                first_block = new pool_block_t(first_block);
                begin_used_in_first_block =
                    new_el = &first_block->data[sizeof(first_block->data) -
                                                                    sizeof(U)];
            }
            return new(new_el) U(std::forward<Args>(args)...);
        }
      public:
        /// \brief allocate and construct a new element (of any type)
        template <class U, class... Args>
        U* construct(Args&&... args)
        {
          static_assert(std::is_trivially_destructible_v<U>);
          if constexpr (sizeof(U) == sizeof(T))
          {
            return construct_samesize<U>(std::forward<Args>(args)...);
            }
            else
            {                                                                   static_assert(sizeof(U) <= sizeof(first_block->data));
                return construct_othersize<U>(std::forward<Args>(args)...);
            }
        }


        /// \brief destroy and delete some element
        /// \details destroy() is only allowed if the destructor of U is
        /// trivial.  This ensures that in my_pool::~my_pool() we do not have
        /// to test whether some element has been freed before we destroy it.
        /// Also, the size of U has to be the same size as the size of T, so
        /// each entry in the free list has the same size.
        template <class U>
        void destroy(U* const old_el)
        {                                                                       static_assert(sizeof(T) == sizeof(U));
            old_el->~U();
            static_assert(std::is_trivially_destructible_v<U>);
#ifndef NDEBUG
                                                                                    // ensure that old_el points to an element in some block
                                                                                    static std::less<const void*> const total_order; // NOLINT(modernize-use-transparent-functors) explicit type required to compare differing pointer types
                                                                                    for (const pool_block_t* block(first_block);
                                                                                                assert(nullptr != block),
                                                                                                total_order(old_el, block->data) ||
                                                                                                total_order(&block->data[sizeof(block->data)], old_el + 1);
                                                                                                                                    block = block->next_block )
                                                                                    {  }
                                                                                #endif
            deref_void(old_el) = first_free_T;
            first_free_T = static_cast<void*>(old_el);
        }
    };
#else
    #define ONLY_IF_POOL_ALLOCATOR(...)
    /// \endcond
#endif // #define USE_POOL_ALLOCATOR

#ifdef USE_SIMPLE_LIST
    /// \class simple_list
    /// \brief a simple implementation of lists
    /// \details This class simplifies lists:  It assumes that list entries are
    /// trivially destructible, and it does not store the size of a list.
    /// Therefore, the destructor, erase() and splice() can be simplified.
    /// Also, the simple_list object itself is trivially destructible if the
    /// pool allocator is used; therefore, destroying a block_t object becomes
    /// trivial as well.
    template <class T>
    class simple_list
    {
      public:
        class const_iterator;

    #ifndef USE_POOL_ALLOCATOR
      private:
    #endif
        /// \brief list entry
        /// \details If the list is to use the pool allocator, its designated
        /// type must be `simple_list::entry` so elements can be erased.
        class entry
        {
          private:
            entry* next;
            entry* prev;
            T data;

            friend class simple_list;
            friend class const_iterator;
            friend class my_pool<entry>;

            template <class... Args>
            entry(entry* const new_next, entry* const new_prev, Args&&... args)
              : next(new_next),
                prev(new_prev),
                data(std::forward<Args>(args)...)
            {  }
        };

      private:
        /// \brief pointer to the beginning of the list
        entry* first;

      public:
        #ifdef USE_POOL_ALLOCATOR
            static my_pool<entry>& get_pool()
            {
                static my_pool<entry> pool;

                return pool;
            }
        #endif

        /// \brief constant iterator class for simple_list
        class const_iterator
        {
          public:
            using value_type = T;
            using pointer = T*;
            using reference = T&;
            using difference_type = std::ptrdiff_t;
            using iterator_category = std::forward_iterator_tag;
          protected:
            entry* ptr;

            const_iterator(const entry* const new_ptr)
              : ptr(const_cast<entry*>(new_ptr))
            {  }

            friend class simple_list;
          public:
            const_iterator() = default;
            const_iterator(const const_iterator& other) = default;
            const_iterator& operator=(const const_iterator& other) = default;

            const_iterator& operator++()
            {                                                                   assert(nullptr != ptr);
                ptr = ptr->next;
                return *this;
            }

            const_iterator& operator--()
            {                                                                   assert(nullptr != ptr);
                ptr = ptr->prev;                                                assert(nullptr != ptr->next);
                return *this;
            }

            const_iterator operator++(int)
            {                                                                   assert(nullptr != ptr);
                const_iterator temp(*this);
                ptr = ptr->next;
                return temp;
            }

            const_iterator operator--(int)
            {                                                                   assert(nullptr != ptr);
                const_iterator temp(*this);
                ptr = ptr->prev;                                                assert(nullptr != ptr->next);
                return temp;
            }

            const T& operator*() const
            {                                                                   assert(nullptr != ptr);
                return  ptr->data;
            }

            const T* operator->() const
            {                                                                   assert(nullptr != ptr);
                return &ptr->data;
            }

            bool operator==(const const_iterator& other) const
            {
                return ptr == other.ptr;
            }

            bool operator==(const T* const other) const
            {                                                                   assert(nullptr != other);
                // It is allowed to call this even if is_null().
                return ptr != nullptr && operator->() == other;
            }

        };

        /// \brief iterator class for simple_list
        class iterator : public const_iterator
        {
          public:
            using typename const_iterator::value_type;
            using typename const_iterator::pointer;
            using typename const_iterator::reference;
            using typename const_iterator::difference_type;
            using typename const_iterator::iterator_category;
          protected:
            iterator(entry* const new_ptr) : const_iterator(new_ptr)  {  }

            friend class simple_list;
          public:
            iterator() = default;

            iterator(const iterator& other) = default;

            iterator& operator=(const iterator& other) = default;

            iterator& operator++(){const_iterator::operator++(); return *this;}

            iterator& operator--(){const_iterator::operator--(); return *this;}

            iterator operator++(int)
            {
                iterator temp(*this);
                const_iterator::operator++();
                return temp;
            }

            iterator operator--(int)
            {
                iterator temp(*this);
                const_iterator::operator--();
                return temp;
            }

            T& operator*() const
            {
                return const_cast<T&>(const_iterator::operator*());
            }

            T* operator->() const
            {
                return const_cast<T*>(const_iterator::operator->());
            }
        };

        /// \brief class that stores either an iterator or a null value
        /// \details We cannot use C++14's "null forward iterators", as they
        /// are not guaranteed to compare unequal to valid iterators.  We also
        /// need to compare null iterators with non-null ones.
        class iterator_or_null : public iterator
        {
          public:
            using typename iterator::value_type;
            using typename iterator::pointer;
            using typename iterator::reference;
            using typename iterator::difference_type;
            using typename iterator::iterator_category;
            using iterator::operator*;
            using iterator::operator->;

            iterator_or_null() : iterator()  {  }

            iterator_or_null(std::nullptr_t) : iterator()
            {
                 const_iterator::ptr = nullptr;
            }

            iterator_or_null(const iterator& other) : iterator(other)  {  }

            bool is_null() const  {  return nullptr == const_iterator::ptr;  }

            iterator_or_null& operator=(std::nullptr_t)
            {
                const_iterator::ptr = nullptr;
                return *this;
            }
        };

        /// \brief constructor
        simple_list()
          : first(nullptr)
        {
          static_assert(std::is_trivially_destructible_v<entry>);
        }

        #ifndef USE_POOL_ALLOCATOR
            /// \brief destructor
            ~simple_list()
            {
                for (iterator iter = begin(); end() != iter; )
                {
                    iterator next = std::next(iter);
                    delete iter.ptr;
                    iter = next;
                }
            }
        #endif

        /// \brief return true iff the list is empty
        bool empty() const  {  return nullptr==first;  }

        /// \brief return an iterator to the first element of the list
        iterator begin()              {  return iterator(first);         }

        /// \brief return an iterator past the last element of the list
        static iterator end()         {  return iterator(nullptr);       }

        /// \brief return a constant iterator to the first element of the list
        const_iterator cbegin() const {  return const_iterator(first);   }

        /// \brief return a constant iterator past the last element of the list
        static const_iterator cend()  {  return end();                   }

        /// \brief return a constant iterator to the first element of the list
        const_iterator begin() const  {  return cbegin();                }

        /// \brief return an iterator to the last element of the list
        iterator before_end()
        {                                                                       assert(!empty());
            return iterator(first->prev);
        }

        const_iterator before_end() const
        {                                                                       assert(!empty());
            return const_iterator(first->prev);
        }

        /// \brief return a reference to the first element of the list
        T& front()
        {                                                                       assert(!empty());
            return first->data;
        }

        /// \brief return a reference to the last element of the list
        T& back()
        {                                                                       assert(!empty());
            return first->prev->data;
        }
                                                                                    [[nodiscard]]
                                                                                    bool check_linked_list() const
                                                                                    {
                                                                                        if (empty())
                                                                                        {
                                                                                            return true;
                                                                                        }
                                                                                        const_iterator i = first;
                                                                                        if (nullptr == i.ptr->prev)
                                                                                        {
                                                                                            return false;
                                                                                        }
                                                                                        while (nullptr != i.ptr->next)
                                                                                        {
                                                                                            if (i.ptr->next->prev != i.ptr)
                                                                                            {
                                                                                                return false;
                                                                                            }
                                                                                            ++i;
                                                                                            assert(i.ptr->prev->next == i.ptr);
                                                                                        }
                                                                                        return first->prev == i.ptr;
                                                                                    }
        /// \brief construct a new list entry before pos
        /// \details If pos==end(), construct a new list entry at the end
        template<class... Args>
        iterator emplace(iterator pos, Args&&... args)
        {                                                                       assert(end()==pos || !empty()); assert(end()==pos || nullptr!=pos.ptr->prev);
                                                                                #ifndef NDEBUG
                                                                                    assert(check_linked_list());
                                                                                    if (end() != pos)
                                                                                    { for(const_iterator i=begin(); i!=pos; ++i) { assert(end()!=i); } }
                                                                                #endif
            entry* const prev = end() == pos
                                ? (empty() ? nullptr : first->prev)
                                : pos.ptr->prev;
            entry* const new_entry(
                #ifdef USE_POOL_ALLOCATOR
                    get_pool().template construct<entry>
                #else
                    new entry
                #endif
                                 (pos.ptr, prev, std::forward<Args>(args)...));
            if (begin() == pos)
            {
                // we insert a new element before the current list begin, so
                // the begin should change.  This includes the case that the
                // list was empty before.
                first = new_entry;
            }
            else if (nullptr != prev)
            {
                // We insert an element not at the current list begin, so it
                // should be reachable from its predecessor.
                prev->next = new_entry;
            }
            if (end() != pos)
            {
                pos.ptr->prev = new_entry;
            }
            else
            {                                                                   assert(nullptr != first);
                first->prev = new_entry;
            }                                                                   assert(check_linked_list());
                                                                                #ifndef NDEBUG
                                                                                    assert((end() == pos ? before_end() : pos.ptr->prev) == new_entry);
                                                                                    for (const_iterator i=begin(); i != new_entry; ++i) { assert(end() != i); }
                                                                                #endif
            return iterator(new_entry);
        }


        /// Puts a new element at the end.
        template <class... Args>
        iterator emplace_back(Args&&... args)
        {
            return emplace(end(), std::forward<Args>(args)...);
        }


        /// \brief construct a new list entry after pos
        /// \details if pos==end(), the new list entry is created at the front.
        template<class... Args>
        iterator emplace_after(iterator pos, Args&&... args)
        {                                                                       assert(end()==pos || !empty()); assert(end()==pos || end()!=pos.ptr->prev);
                                                                                #ifndef NDEBUG
                                                                                    assert(check_linked_list());
                                                                                    if (end() != pos) {
                                                                                        for (const_iterator i = begin(); i != pos; ++i) { assert(end() != i); }
                                                                                    }
                                                                                #endif
            entry* const next = end() == pos ? begin().ptr : pos.ptr->next;
            entry* const new_entry(
                #ifdef USE_POOL_ALLOCATOR
                    get_pool().template construct<entry>
                #else
                    new entry
                #endif
                                 (next, pos.ptr, std::forward<Args>(args)...));
            if (end() == pos)
            {
                // we insert a new element before the current list begin, so
                // the begin should change.  This includes the case that the
                // list was empty before.
                new_entry->prev = empty() ? new_entry : before_end().ptr;
                first = new_entry;
            }
            else
            {
                pos.ptr->next = new_entry;
            }                                                                   assert(nullptr != first);
            if (nullptr == next)
            {
                first->prev = new_entry;
            }
            else
            {
                next->prev = new_entry;
            }                                                                   assert(check_linked_list());
                                                                                #ifndef NDEBUG
                                                                                    assert((end() == pos ? first : pos.ptr->next) == new_entry);
                                                                                    for (const_iterator i=begin(); i != new_entry; ++i) { assert(end() != i); }
                                                                                #endif
            return iterator(new_entry);
        }


        /// \brief construct a new list entry at the beginning
        template<class... Args>
        iterator emplace_front(Args&&... args)
        {
            return emplace_after(end(), std::forward<Args>(args)...);
        }


        /// The function moves the element pointed at by from_pos (that is in
        /// the list indicated by the 2nd parameter) just after position to_pos
        /// (that is in this list). If to_pos == end(), move the element to the
        /// beginning of this list.
        void splice_to_after(iterator const to_pos, simple_list<T>& from_list,
                                                       iterator const from_pos)
        {                                                                       assert(from_pos != to_pos);
                                                                                #ifndef NDEBUG
                                                                                    assert(check_linked_list());
                                                                                    if (end() != to_pos) {
                                                                                        assert(!empty());  assert(nullptr != to_pos.ptr->prev);
                                                                                        for(const_iterator i = begin(); i!=to_pos; ++i) { assert(end() != i); }
                                                                                    }
                                                                                    assert(end() != from_pos);  assert(nullptr != from_pos.ptr->prev);
                                                                                    assert(!from_list.empty());  assert(from_list.check_linked_list());
            /* remove element from_pos from its original list                */     for (const_iterator i = from_list.begin(); i != from_pos; ++i) {
                                                                                        assert(from_list.end() != i);
                                                                                    }
                                                                                #endif
            if (from_pos.ptr != from_list.first)
            {                                                                   assert(nullptr != from_list.first->next); // at least 2 elements in the list
                /* not the first element in from_list                        */ assert(from_pos == from_pos.ptr->prev->next);
                from_pos.ptr->prev->next = from_pos.ptr->next;
                if (nullptr != from_pos.ptr->next)
                {
                    /* not the last element in from_list                     */ assert(from_pos == from_pos.ptr->next->prev);
                    from_pos.ptr->next->prev = from_pos.ptr->prev;
                }
                else
                {
                    /* last element in from_list                             */ assert(from_pos == from_list.first->prev);
                    from_list.first->prev = from_pos.ptr->prev;
                }
            }
            else
            {
                /* first element in from_list                                */ assert(nullptr == from_pos.ptr->prev->next);
                from_list.first = from_pos.ptr->next;
                if (!from_list.empty())
                {
                    /* not the last element in from_list                     */ assert(from_pos == from_pos.ptr->next->prev);
                    from_pos.ptr->next->prev = from_pos.ptr->prev;
                }
            }
            // update the pointers of from_pos and insert from_pos into this
            // list
            entry* next;
            if (end() == to_pos)
            {
                // we insert the element before the current list begin, so
                // the begin should change.  This includes the case that the
                // list was empty before.
                if (!empty())
                {
                    from_pos.ptr->prev = before_end().ptr;
                }
                // else from_pos->prev = from_pos; -- will be set below.
                next = first;
                first = from_pos.ptr;
            }
            else
            {
                from_pos.ptr->prev = to_pos.ptr;
                next = to_pos.ptr->next;
                to_pos.ptr->next = from_pos.ptr;
            }                                                                   assert(nullptr != first);
            from_pos.ptr->next = next;
            if (nullptr == next)
            {
                first->prev = from_pos.ptr;
            }
            else
            {
                next->prev = from_pos.ptr;
            }                                                                   assert(check_linked_list());  assert(from_list.check_linked_list());
                                                                                #ifndef NDEBUG
                                                                                    assert(from_pos == (end() == to_pos ? first : to_pos.ptr->next));
                                                                                    for (const_iterator i=begin(); i!=from_pos; ++i) { assert(i!=end()); }
                                                                                    if (end() != to_pos) {
                                                                                        for (const_iterator i=begin(); i!=to_pos; ++i) { assert(end() != i); }
                                                                                    }
                                                                                #endif
        }


        /// \brief move a list entry from one position to another (possibly in
        /// a different list)
        /// The function moves the element pointed at by from_pos (that is in
        /// the list indicated by the 2nd parameter) just before position
        /// to_pos (that is in this list). If to_pos == end(), move the element
        /// to the end of this list.
        void splice(iterator const to_pos, simple_list<T>& from_list,
                                                       iterator const from_pos)
        {                                                                       assert(from_pos != to_pos);  assert(end() == to_pos || !empty());
                                                                                #ifndef NDEBUG
                                                                                    assert(end() == to_pos || nullptr != to_pos.ptr->prev);
                                                                                    assert(check_linked_list());
                                                                                    if (end() != to_pos) {
                                                                                        for (const_iterator i=begin(); i!=to_pos; ++i) { assert(end() != i); }
                                                                                    }
                                                                                    assert(from_list.end() != from_pos); assert(nullptr != from_pos.ptr->prev);
                                                                                    assert(!from_list.empty());  assert(from_list.check_linked_list());
            /* remove element from_pos from its original list                */     for (const_iterator i=from_list.begin(); i!=from_pos; ++i) {
                                                                                        assert(from_list.end() != i);
                                                                                    }
                                                                                #endif
            if (from_pos != from_list.first)
            {                                                                   assert(nullptr != from_list.first->next); // at least 2 elements in the list
                /* not the first element in from_list                        */ assert(from_pos.ptr == from_pos.ptr->prev->next);
                from_pos.ptr->prev->next = from_pos.ptr->next;
                if (nullptr != from_pos.ptr->next)
                {
                    /* not the last element in from_list                     */ assert(from_pos.ptr == from_pos.ptr->next->prev);
                    from_pos.ptr->next->prev = from_pos.ptr->prev;
                }
                else
                {
                    /* last element in from_list                             */ assert(from_pos.ptr == from_list.first->prev);
                    from_list.first->prev = from_pos.ptr->prev;
                }
            }
            else
            {
                /* first element in from_list                                */ assert(nullptr == from_pos.ptr->prev->next);
                from_list.first = from_pos.ptr->next;
                if (!from_list.empty())
                {
                    /* not the last element in from_list                     */ assert(from_pos.ptr == from_pos.ptr->next->prev);
                    from_pos.ptr->next->prev = from_pos.ptr->prev;
                }
            }
            // update the pointers of from_pos and insert from_pos into this
            // list
            from_pos.ptr->next = to_pos.ptr;
            if (end() == to_pos)
            {
                // from_pos becomes the last element in *this
                if (empty())
                {
                    from_pos.ptr->prev = from_pos.ptr;
                    first = from_pos.ptr;                                       assert(check_linked_list());  assert(from_list.check_linked_list());
                    return;
                }
                from_pos.ptr->prev = before_end().ptr;
                from_pos.ptr->prev->next = from_pos.ptr;
                first->prev = from_pos.ptr;
            }
            else
            {
                /* from_pos does not become the last element in *this        */ assert(!empty());
                from_pos.ptr->prev = to_pos.ptr->prev;
                to_pos.ptr->prev = from_pos.ptr;
                if (to_pos.ptr == first)
                {
                    // we insert a new element before the current list begin,
                    // so the begin should change.
                    first = from_pos.ptr;
                }
                else
                {
                    from_pos.ptr->prev->next = from_pos.ptr;
                }
            }                                                                   assert(check_linked_list());  assert(from_list.check_linked_list());
                                                                                #ifndef NDEBUG
                                                                                    assert((end() == to_pos ? before_end().ptr
                                                                                                            : iterator(to_pos.ptr->prev)) == from_pos.ptr);
                                                                                    for (const_iterator i=begin(); i != from_pos; ++i) { assert(end() != i); }
                                                                                    if (end() != to_pos) {
                                                                                        for (const_iterator i=begin(); i!=to_pos; ++i) { assert(end() != i); }
                                                                                    }
                                                                                #endif
        }


        /// \brief erase an element from a list
        void erase(iterator const pos)
        {                                                                       assert(end() != pos);  assert(nullptr != pos.ptr->prev);  assert(!empty());
                                                                                #ifndef NDEBUG
                                                                                    assert(check_linked_list());
                                                                                    for (const_iterator i = begin(); i != pos; ++i) { assert(end() != i); }
                                                                                #endif
            if (pos != first)
            {                                                                   assert(nullptr != first->next); // at least 2 elements in the list
                /* not the first element                                     */ assert(pos.ptr == pos.ptr->prev->next);
                pos.ptr->prev->next = pos.ptr->next;
                if (nullptr != pos.ptr->next)
                {
                    /* not the last element                                  */ assert(pos.ptr == pos.ptr->next->prev);
                    pos.ptr->next->prev = pos.ptr->prev;
                }
                else
                {
                    /* last element                                          */ assert(pos.ptr == first->prev);
                    first->prev = pos.ptr->prev;
                }
            }
            else
            {
                /* first element                                             */ assert(nullptr == pos.ptr->prev->next);
                first = pos.ptr->next;
                if (!empty())
                {
                    /* not the last element                                  */ assert(pos.ptr == pos.ptr->next->prev);
                    pos.ptr->next->prev = pos.ptr->prev;
                }
            }
            #ifdef USE_POOL_ALLOCATOR
                get_pool().destroy(pos.ptr);
            #else
                delete pos.ptr;
            #endif
        }


        /// The function computes the successor of pos in the list.  If pos is
        /// the last element of the list, it returns end().  It is an error if
        /// pos==end() or if pos is not in the list.
        #ifdef NDEBUG
            static // only in debug mode it accesses data of the list itself
        #endif
        iterator next(iterator pos)
                                                                                #ifndef NDEBUG
                                                                                    const // static functions cannot be const
                                                                                #endif
        {                                                                       assert(end() != pos);
                                                                                #ifndef NDEBUG
                                                                                    for (const_iterator i = begin(); i != pos; ++i) { assert(end() != i); }
                                                                                #endif
            return iterator(pos.ptr->next);
        }


        /// The function computes the successor of pos in the list.  If pos is
        /// the last element of the list, it returns end().  It is an error if
        /// pos==end() or if pos is not in the list.
        #ifdef NDEBUG
            static // only in debug mode it accesses data of the list itself
        #endif
        const_iterator next(const_iterator pos)
                                                                                #ifndef NDEBUG
                                                                                    const // static functions cannot be const
                                                                                #endif
        {                                                                       assert(end() != pos);
                                                                                #ifndef NDEBUG
                                                                                    for (const_iterator i = begin(); i != pos; ++i) { assert(end() != i); }
                                                                                #endif
            return const_iterator(pos.ptr->next);
        }


        /// The function computes the predecessor of pos in the list.  If pos is at
        /// the beginning of the list, it returns end().  It is an error if
        /// pos==end() or if pos is not in the list.
        iterator prev(iterator pos) const
        {                                                                       assert(pos!=end());
                                                                                #ifndef NDEBUG
                                                                                    for (const_iterator i = begin(); i != pos; ++i) { assert(end() != i); }
                                                                                #endif
            return begin() == pos ? end() : iterator(pos.ptr->prev);
        }


        /// The function computes the predecessor of pos in the list.  If pos is at
        /// the beginning of the list, it returns end().  It is an error if
        /// pos==end() or if pos is not in the list.
        const_iterator prev(const_iterator pos) const
        {                                                                       assert(end() != pos);
                                                                                #ifndef NDEBUG
                                                                                    for (const_iterator i = begin(); i != pos; ++i) { assert(end() != i); }
                                                                                #endif
            return begin() == pos ? end() : const_iterator(pos.ptr->prev);
        }


        bool operator==(const simple_list& other) const
        {
            const_iterator it = cbegin();
            const_iterator other_it = other.cbegin();
            while (cend() != it)
            {
                if (cend() == other_it || *it != *other_it)
                {
                    return false;
                }
                ++it;
                ++other_it;
            }
            return end()==other_it;
        }
    };

    template <class El>
    using iterator_or_null_t = typename simple_list<El>::iterator_or_null;
#else
    #define simple_list std::list

    template <class El>
    union iterator_or_null_t
    {
      public:
        typedef std::list<El>::iterator iterator;
        typedef std::list<El>::const_iterator const_iterator;
      private:
        const void* null;
        iterator iter;
      public:
        /// \brief Construct an uninitialized object
        iterator_or_null_t()
        {
            if constexpr (!std::is_trivially_destructible<iterator>::value)
            {
                // We still have to internally decide whether to construct
                // the iterator or not so the destructor knows what to do.
                null = nullptr;
            }
        }


        /// \brief Construct an object containing a null pointer
        explicit iterator_or_null_t(nullptr_t)
        {
            null = nullptr;
        }


        /// \brief Construct an object containing a valid iterator
        explicit iterator_or_null_t(const iterator other)
        {
            new (&iter) iterator(other);                                        assert(nullptr != null);
        }


        /// \brief Check whether the object contains a valid iterator
        bool is_null() const  {  return nullptr == null;  }


        /// \brief Destruct an object (whether it contains a valid iterator or
        // not)
        ~iterator_or_null_t()
        {
            if constexpr (!std::is_trivially_destructible<iterator>::value)
            {
                if (!is_null())  iter.~iterator();
            }
        }

        El* operator->()
        {                                                                       assert(nullptr != null);
            return iter.operator->();
        }
        El& operator*()
        {                                                                       assert(nullptr != null);
            return iter.operator*();
        }

        void operator=(nullptr_t)
        {
            if constexpr (!std::is_trivially_destructible<iterator>::value)
            {
                if (!is_null())  iter.~iterator();
            }
            null = nullptr;
        }


        explicit operator iterator() const
        {                                                                       assert(nullptr != null);
            return iter;
        }


        void operator=(const iterator& other)
        {
            if constexpr (!std::is_trivially_destructible<iterator>::value)
            {
                if (!is_null())  iter.~iterator();
            }
            new (&iter) iterator(other);                                        assert(nullptr != null);
        }

        /// \brief Compare the object with another iterator_or_null_t object
        /// \details The operator could be templated so that iterator_or_null_t
        /// objects of different types can be compared.
        bool operator==(const iterator_or_null_t other) const
        {
            if constexpr (sizeof(null) == sizeof(iter))
            {
                return &*iter == &*other.iter;
            }
            else
            {
                return (is_null() && other.is_null()) ||
                    (!is_null() && !other.is_null() && &*iter == &*other.iter);
            }
        }


        /// \brief Compare the object with another iterator_or_null_t object
        bool operator!=(const iterator_or_null_t other) const
        {
            return !operator==(other);
        }


        /// \brief Compare the object with an iterator
        /// \details If the object does not contain a valid iterator, it
        /// compares unequal with the iterator.
        bool operator==(const const_iterator other) const
        {                                                                       // assert(nullptr != &*other); -- generates a warning
            return (sizeof(null) == sizeof(iter) || !is_null()) &&
                                                             &*iter == &*other;
        }


        bool operator!=(const const_iterator other) const
        {
            return !operator==(other);
        }


        /// \brief Compare the object with a non-null pointer
        /// \details If the object does not contain a valid iterator, it
        /// compares unequal with the pointer.
        bool operator==(const El* const other) const
        {                                                                       assert(nullptr != other);
            return !is_null() && &*iter == other;
        }


        bool operator!=(const El* const other) const
        {
            return !operator==(other);
        }
    };
#endif

} // end namespace detail
// end namespace lts
// end namespace mcrl2

#endif // ifndef SIMPLE_LIST_H
