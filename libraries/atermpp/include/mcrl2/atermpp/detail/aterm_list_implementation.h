// Author(s): Jan Friso Groote, Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_DETAIL_ATERM_LIST_IMPLEMENTATION_H
#define MCRL2_ATERMPP_DETAIL_ATERM_LIST_IMPLEMENTATION_H
#pragma once

#include "mcrl2/atermpp/aterm_appl.h"
#include "mcrl2/atermpp/aterm_list.h"

#include "mcrl2/utilities/exception.h"
#include "mcrl2/utilities/workarounds.h"

namespace atermpp
{

constexpr std::size_t LengthOfShortList = 10000;  /// \brief The length of a short list. If lists
                                                  ///        are short the stack can be used for temporary data.
                                                  ///        Otherwise the heap must be used to avoid stack overflow.
                                                  ///        The chosen value is rather arbitrary.

template <class Term>
void term_list<Term>::push_front(const Term& el)
{
   detail::g_thread_term_pool().create_appl(*this, detail::g_term_pool().as_list(), el, *this);
}

template <class Term>
template<typename ...Args>
void term_list<Term>::emplace_front(Args&&... arguments)
{
  detail::g_thread_term_pool().create_appl(*this, detail::g_term_pool().as_list(), Term(std::forward<Args>(arguments)...), *this);
}

template <typename Term>
inline
term_list<Term> push_back(const term_list<Term>& l, const Term& el)
{
  typedef typename term_list<Term>::const_iterator const_iterator;

  const std::size_t len = l.size();

  // The resulting list
  term_list<Term> result;
  result.push_front(el);

  if (len < LengthOfShortList)
  {
    // The list is short, use the stack for temporal storage.
    const_iterator* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(const_iterator, len);

    // Collect all elements of list in buffer.
    std::size_t j=0;
    for (const_iterator i = l.begin(); i != l.end(); ++i, ++j)
    {
      buffer[j]=i;
    }

    // Insert elements at the front of the list.
    while (j>0)
    {
      j=j-1;
      result.push_front(*buffer[j]);
    }
  }
  else
  {
    // The list is long. Use the heap to store intermediate data.
    std::vector<Term> buffer;
    buffer.reserve(len);

    for (const Term& t: l)
    {
      buffer.push_back(t);
    }

    // Insert elements at the front of the list
    for (typename std::vector<Term>::reverse_iterator i=buffer.rbegin(); i!=buffer.rend(); ++i)
    {
      result.push_front(*i);
    }
  }

  return result;
}


template <typename Term>
inline
term_list<Term> reverse(const term_list<Term>& l)
{
  term_list<Term> result;
  for(const Term& t: l)
  {
    result.push_front(t);
  }
  return result;
}

template <typename Term>
inline
term_list<Term> sort_list(const term_list<Term>& l, 
                             const std::function<bool(const Term&, const Term&)>& ordering 
                                  /* = [](const Term& t1, const Term& t2){ return t1<t2;}*/ )
{
  const std::size_t len = l.size();
  if (len<=1)
  {
    return l;
  }
 
  // The resulting list
  term_list<Term> result;

  if (len < LengthOfShortList)
  {
    // The list is short, use the stack for temporal storage.
    Term* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(Term, len);

    // Collect all elements of list in buffer.
    std::size_t j=0;
    for (const Term& t: l)
    {
      new (buffer+j) Term(t); // A mcrl2 stack allocator does not handle construction by default. 
      ++j;
    }
   
    std::sort(buffer, buffer+len, ordering);

    // Insert elements at the front of the list.
    while (j>0)
    {
      j=j-1;
      result.push_front(buffer[j]);
      buffer[j].~Term();    // Explicitly call the destructor, as an mCRL2 stack allocator does not do that itself. . 
    }
  }
  else
  {
    // The list is long. Use the heap to store intermediate data.
    std::vector<Term> buffer;
    buffer.reserve(len);

    for (const Term& t: l)
    {
      buffer.push_back(t);
    }

    // Sort using a standard algorithm.
    std::sort(buffer.begin(), buffer.end(), ordering);

    // Insert elements at the front of the list
    for (typename std::vector<Term>::reverse_iterator i=buffer.rbegin(); i!=buffer.rend(); ++i)
    {
      result.push_front(*i);
    }
  }
  return result;
}


template <typename Term>
inline
term_list<Term> remove_one_element(const term_list<Term>& list, const Term& t)
{
  typedef typename term_list<Term>::const_iterator const_iterator;

  std::size_t len=0;
  const_iterator i = list.begin();
  for( ; i!=list.end(); ++i, ++len)
  {
    if (*i==t)
    {
      break;
    }
  }

  if (i==list.end())
  {
    // Term t not found in list.
    return list;
  }

  const_iterator* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(const_iterator, len);

  term_list<Term> result = list;
  std::size_t k=0;
  for(const_iterator j = list.begin(); j != i; ++j, ++k)
  {
    buffer[k]=j;
    result.pop_front();
  }
  assert(len==k);
  assert(result.front()==t);
  result.pop_front(); // skip the element.

  while (k>0)
  {
    k=k-1;
    result.push_front(*buffer[k]);
  }

  return result;
}


template <typename Term1, typename Term2>
inline
typename std::conditional<std::is_convertible<Term2,Term1>::value,term_list<Term1>,term_list<Term2>>::type
operator+(const term_list<Term1>& l, const term_list<Term2>& m)
{
  static_assert(std::is_convertible< Term1, Term2 >::value ||
                std::is_convertible< Term2, Term1 >::value,
                       "Concatenated lists must be of convertible types. ");
  static_assert(sizeof(Term1) == sizeof(aterm),
                "aterm cast cannot be applied types derived from aterms where extra fields are added. ");
  static_assert(sizeof(Term2) == sizeof(aterm),
                "aterm cast cannot be applied types derived from aterms where extra fields are added. ");
  typedef typename std::conditional<std::is_convertible<Term2,Term1>::value,Term1,Term2>::type ResultType;
  typedef typename term_list<Term1>::const_iterator const_iterator;

  if (m.empty())
  {
    return reinterpret_cast<const term_list<ResultType>&>(l);
  }

  std::size_t len = l.size();

  if (len == 0)
  {
    return reinterpret_cast<const term_list<ResultType>&>(m);
  }

  term_list<ResultType> result = reinterpret_cast<const term_list<ResultType>&>(m);
  if (len < LengthOfShortList)
  {
    // The length is short. Use the stack for temporary storage.
    const_iterator* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(const_iterator, len);

    std::size_t j=0;
    for (const_iterator i = l.begin(); i != l.end(); ++i, ++j)
    {
      buffer[j]=i;
    }
    assert(j==len);

    // Insert elements at the front of the list
    while (j>0)
    {
      j=j-1;
      result.push_front(*buffer[j]);
    }
  }
  else
  {
    // The length of l is very long. Use the heap for temporary storage.
    std::vector<ResultType> buffer;
    buffer.reserve(len);

    for (const Term1& t: l)
    {
      buffer.push_back(t);
    }

    // Insert elements at the front of the list
    for(typename std::vector<ResultType>::const_reverse_iterator i=buffer.rbegin(); i!=buffer.rend(); ++i)
    {
      result.push_front(*i);
    }
  }
  return result;
}


namespace detail
{
  template <class Term, class Iter, class ATermConverter, class ATermFilter>
  inline aterm make_list_backward(Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter)
  {
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(aterm),"Term derived from an aterm must not have extra fields");

    aterm_list list;
    while (first != last)
    {
      --last;
      const Term t = convert_to_aterm(*last);
      if (aterm_filter(t))
      {
        list.push_front(t);
      }
    }

    return mcrl2::workaround::return_std_move(list);
  }

  template <class Term, class Iter, class ATermConverter>
  inline aterm make_list_backward(Iter first, Iter last, ATermConverter convert_to_aterm)
  {
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(aterm),"Term derived from an aterm must not have extra fields");

    aterm_list list;
    while (first != last)
    {
      --last;
      list.push_front(convert_to_aterm(*last));
    }

    return mcrl2::workaround::return_std_move(list);
  }

  // See the note at make_list_backwards for why there are two almost similar version of make_list_forward.
  template <class Term, class Iter, class ATermConverter, class ATermFilter>
  inline aterm make_list_forward(Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter)
  {
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(aterm),"Term derived from an aterm must not have extra fields");

    // The list that will be returned.
    aterm_list list;

    const std::size_t len = std::distance(first,last);
    if (len < LengthOfShortList)  // If the list is sufficiently short, use the stack.
    {
      Term* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(Term, len);
      Term *const buffer_begin=buffer;
      Term* i=buffer_begin;
      for(; first != last; ++first)
      {
        const Term t = convert_to_aterm(*first);
        if (aterm_filter(t))
        {
          // Placement new; The buffer is not properly initialised.
          new (i) Term(t);
          ++i;
        }
      }

      // Construct the list using the temporary array of elements.
      for( ; i != buffer_begin; )
      {
        --i;
        list.push_front(*i);
        (*i).~Term(); // Destroy the elements in the buffer explicitly.
      }
    }
    else
    {
      // The list is long. Therefore use the heap for temporary storage.
      std::vector<Term> buffer;
      buffer.reserve(len);
      for(; first != last; ++first)
      {
        const Term t = convert_to_aterm(*first);
        if (aterm_filter(t))
        {
          // Placement new; The buffer is not properly initialised.
          buffer.push_back(t);
        }
      }

      // Construct the list using the temporary array of elements.
      for(typename std::vector<Term>::const_reverse_iterator i = buffer.rbegin(); i != buffer.rend(); ++i)
      {
        list.push_front(*i);
      }
    }

    return mcrl2::workaround::return_std_move(list);
  }

  template <class Term, class Iter, class ATermConverter>
  inline aterm make_list_forward(Iter first, Iter last, ATermConverter convert_to_aterm)
  {
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(aterm),"Term derived from an aterm must not have extra fields");

    // The list that will be returned.
    aterm_list list;

    const std::size_t len = std::distance(first,last);
    if (len < LengthOfShortList) // If the list is sufficiently short, use the stack.
    {
      Term* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(Term, len);
      Term *const buffer_begin = buffer;
      Term* i = buffer_begin;
      for(; first != last; ++first,++i)
      {
        // Placement new; The buffer is not properly initialised.
        new (i) Term(convert_to_aterm(*first));
      }

      for( ; i != buffer_begin; )
      {
        --i;
        list.push_front(*i);
        (*i).~Term(); // Destroy the elements in the buffer explicitly.
      }
    }
    else
    {
      // The list is very long. Reserve memory on the heap.
      std::vector<Term> buffer;
      buffer.reserve(len);
      for(; first != last; ++first)
      {
        buffer.push_back(convert_to_aterm(*first));
      }

      for(typename std::vector<Term>::const_reverse_iterator i = buffer.rbegin(); i != buffer.rend(); ++i)
      {
        list.push_front(*i);
      }
    }

    return mcrl2::workaround::return_std_move(list);
  }
} // detail



} // namespace atermpp

#endif // MCRL2_ATERMPP_DETAIL_ATERM_LIST_IMPLEMENTATION_H
