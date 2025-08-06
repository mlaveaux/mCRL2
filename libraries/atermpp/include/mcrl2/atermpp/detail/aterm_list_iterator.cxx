// Author(s): Jan Friso Groote, Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/detail/aterm_list_iterator.h
/// \brief Iterator for term_list.

#ifndef MCRL2_ATERMPP_ATERM_LIST_ITERATOR_H
#define MCRL2_ATERMPP_ATERM_LIST_ITERATOR_H
MCRL2_MODULE;

#include <cassert>
#include <cstddef>
#include <memory>
#include <iterator>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/global_aterm_pool.cxx"
#else
  export module atermpp:detail.aterm_list_iterator;

  import :aterm_core;
  import :detail.aterm_data;
  import :detail.aterm_core_data;
  import :detail.global_aterm_pool;
#endif

MCRL2_MODULE_EXPORT namespace atermpp
{

/// \cond INTERNAL_DOCS
namespace detail
{
  template <class Term>
  class _aterm_list;
}

/// \endcond

/// \brief Iterator for term_list.
template <typename Term>
class term_list_iterator 
{
    template<class T>
    friend class term_list;

  protected:
    const detail::_aterm_list<Term>* m_list;

    /// \brief Constructor from an aterm which must be a list.
    /// \param l A sequence of terms
    term_list_iterator(const detail::_aterm* l)
      : m_list(static_cast<const detail::_aterm_list<Term>*>(l))
    { 
      assert(l->function()==detail::g_term_pool().as_list()
             || l->function()==detail::g_term_pool().as_empty_list());
    } 

  public:
    using value_type = Term;
    using reference = Term&;
    using pointer = Term*;
    using difference_type = ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    /// \brief Default constructor.
    term_list_iterator()
      : m_list(nullptr)
    {}

    /// \brief Copy constructor.
    /// \param other A sequence of terms
    term_list_iterator(const term_list_iterator& other)
      : m_list(other.m_list)
    { 
    } 

    /// \brief Assignment
    /// \param other A sequence of terms
    term_list_iterator& operator=(const term_list_iterator& other) = default;

    /// \brief Dereference operator on an iterator
    const Term& operator*() const
    {
      assert(m_list->function()==detail::g_term_pool().as_list());
      return m_list->head();
    }

    /// Arrow operator on an iterator
    const Term* operator->() const
    {
      assert(m_list->function()==detail::g_term_pool().as_list());
      return &m_list->head();
    }
    
    /// \brief Prefix increment operator on iterator.
    term_list_iterator& operator++()
    {
      assert(m_list->function() == detail::g_term_pool().as_list());
      m_list = static_cast<detail::_aterm_list<Term>*>(detail::address(m_list->tail()));
      return *this;
    }

    /// \brief Postfix increment operator on iterator.
    term_list_iterator operator++(int)
    {
      assert(m_list->function() == detail::g_term_pool().as_list());
      const term_list_iterator temp = *this;
      m_list = static_cast<detail::_aterm_list<Term>*>(detail::address(m_list->tail()));
      return temp;
    }

    /// \brief Equality of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators point to the same term_list.
    bool operator ==(const term_list_iterator& other) const
    {
      return m_list == other.m_list;
    }

    /// \brief Inequality of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators do not point to the same term_list.
    bool operator !=(const term_list_iterator& other) const
    {
      return m_list != other.m_list;
    }

    /// \brief Comparison of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the pointer to this termlist is smaller than the other pointer.
    bool operator <(const term_list_iterator& other) const
    {
      return m_list < other.m_list;
    }

    /// \brief Comparison of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators point to the same term_list.
    bool operator <=(const term_list_iterator& other) const
    {
      return m_list <= other.m_list;
    }

    /// \brief Comparison of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators point to the same term_list.
    bool operator >(const term_list_iterator& other) const
    {
      return m_list > other.m_list;
    }

    /// \brief Comparison of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators point to the same term_list.
    bool operator >=(const term_list_iterator& other) const
    {
      return m_list >= other.m_list;
    }

};

/// \brief Reverse iterator for term_list.
template <typename Term>
class reverse_term_list_iterator 
{
    template<class T>
    friend class term_list;

  protected:
    std::size_t m_position;   // m_position refers one above the position to be deliverd. 
    std::unique_ptr<detail::_aterm_list<Term> const*[]> m_list_element_references;

    /// \brief Constructor from an aterm which must be a list.
    /// \param l A sequence of terms
    reverse_term_list_iterator(detail::_aterm const* l)
      : m_position(reinterpret_cast<const detail::_aterm_list<Term>*>(l)->size()),
        m_list_element_references((m_position==0?nullptr:new typename detail::_aterm_list<Term> const*[m_position]))
    { 
      assert(l->function()==detail::g_term_pool().as_list()
             || l->function()==detail::g_term_pool().as_empty_list());
      std::size_t j=0;
      for(detail::_aterm_list<Term> const* t=reinterpret_cast<detail::_aterm_list<Term> const*>(l); 
              t->function()==detail::g_term_pool().as_list(); 
              t=reinterpret_cast<detail::_aterm_list<Term> const*>(detail::address(t->tail())), j++)
      {
        m_list_element_references[j]=t;
      }
    } 

  public:
    using value_type = Term;
    using reference = Term&;
    using pointer = Term*;
    using difference_type = ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    /// \brief Default constructor.
    reverse_term_list_iterator()
      : m_position(0),
        m_list_element_references(nullptr)
    {
    }

    /// \brief The copy constructor is not available.
    /// \param other A sequence of terms
    reverse_term_list_iterator(const reverse_term_list_iterator& other) = delete;

    /// \brief Assignment is not available.
    /// \param other A sequence of terms
    reverse_term_list_iterator& operator=(const reverse_term_list_iterator& other) = delete;

    /// \brief Dereference operator on an iterator
    const Term& operator*() const
    {
      assert(m_list_element_references[m_position-1]->function()==detail::g_term_pool().as_list());
      return m_list_element_references[m_position-1]->head();
    }

    /// Arrow operator on an iterator
    const Term* operator->() const
    {
      assert(m_list_element_references[m_position-1]->function()==detail::g_term_pool().as_list());
      return &(m_list_element_references[m_position-1]->head());
    }
    
    /// \brief Prefix increment operator on iterator.
    reverse_term_list_iterator& operator++()
    {
      assert(m_list_element_references[m_position-1]->function() == detail::g_term_pool().as_list());
      m_position--;
      return *this;
    }

    /// \brief Postfix increment operator on iterator.
    void operator++(int)
    {
      assert(m_list_element_references[m_position-1]->function() == detail::g_term_pool().as_list());
      m_position--;
    }

    /// \brief Equality of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators point to the same term_list.
    bool operator ==(const reverse_term_list_iterator& other) const
    {
      return m_position == other.m_position;
    }

    /// \brief Inequality of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators do not point to the same term_list.
    bool operator !=(const reverse_term_list_iterator& other) const
    {
      return !(*this == other);
    }

    /// \brief Comparison of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the pointer to this termlist is smaller than the other pointer.
    bool operator <(const reverse_term_list_iterator& other) const
    {
      return m_position < other.m_position;
    }

    /// \brief Comparison of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators point to the same term_list.
    bool operator <=(const reverse_term_list_iterator& other) const
    {
      return m_position <= other.m_position;
    }

    /// \brief Comparison of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators point to the same term_list.
    bool operator >(const reverse_term_list_iterator& other) const
    {
      return m_position > other.m_position;
    }

    /// \brief Comparison of iterators.
    /// \param other The iterator with which this iterator is compared.
    /// \return true if the iterators point to the same term_list.
    bool operator >=(const reverse_term_list_iterator& other) const
    {
      return m_position >= other.m_position;
    }

};

} // namespace atermpp

#endif // MCRL2_ATERMPP_ATERM_LIST_ITERATOR_H
