// Author(s): Wieger Wesselink, Jan Friso Groote, Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/aterm_int.h
/// \brief Term containing an integer.

#ifndef MCRL2_ATERMPP_DETAIL_ATERM_LIST_H
#define MCRL2_ATERMPP_DETAIL_ATERM_LIST_H

#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/atermpp/concepts.h"
#include "mcrl2/utilities/stack_array.h"  

namespace atermpp
{

template <class Term>
class term_list;

namespace detail
{
  
constexpr std::size_t LengthOfShortList = 10000; /// \brief The length of a short list. If lists
                                                 ///        are short the stack can be used for temporary data.
                                                 ///        Otherwise the heap must be used to avoid stack overflow.
                                                 ///        The chosen value is rather arbitrary.

template <class Term>
struct do_not_convert_term
{
  void operator()(Term& result, const Term& t) const
  {
    result=t;
  }

  const Term& operator()(const Term& t) const
  {
    return t;
  }

  Term& operator()(Term& t) const
  {
    return t;
  } 
};

template <typename Term>
inline
void make_reverse(term_list<Term>& result, const term_list<Term>& l)
{
  make_term_list<Term>(result);
  for(const Term& t: l)
  {
    result.push_front(t);
  }
}

/// \brief Constructs a list starting from first to last. The iterators are traversed backwards and each element is
/// 		   converted using the TermConverter.
/// \details The functions make_list_backward and make_list_forward with three and four arguments are almost the same.
/// 			 The reason for this is that there is a 5% loss of speed of the toolset when merging these two functions.
///          This is caused by storing and protecting the intermediate value of the converted aterm. See Term t = convert_to_aterm(...).
template <IsATerm Term, typename Iter, typename ATermConverter>
inline aterm make_list_backward(Iter first, Iter last, ATermConverter convert_to_aterm)
{
  term_list<Term> result_list;
  make_list_backward<Term, Iter, ATermConverter>(result_list, first, last, convert_to_aterm);
  return result_list;
}

/// \brief Constructs a list starting from first to last where the result is put in result.
template <IsATerm Term, class Iter, class ATermConverter>
inline void make_list_backward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm)
{
    while (first != last)
    {
      --last;
      result.push_front(convert_to_aterm(*last));
    }
}

/// \brief Constructs a list starting from first to last. The iterators are traversed backwards and each element is
/// 		   converted using the TermConverter and inserted whenever TermFilter yields true for the converted element.
template <IsATerm Term, typename Iter, typename ATermConverter, typename ATermFilter>
inline aterm make_list_backward(Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter)
{
  term_list<Term> result_list;
  make_list_backward<Term, Iter, ATermConverter, ATermFilter>(result_list, first, last, convert_to_aterm, aterm_filter);
  return result_list;
}

/// \brief Construct a list iterating from the last to the first element. Result is put in the variable result.
template <IsATerm Term, class Iter, class ATermConverter, class ATermFilter>
inline void make_list_backward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter)
{
  Term t;
  while (first != last)
  {
    --last;
    t = convert_to_aterm(*last);
    if (aterm_filter(t))
    {
      result.push_front(t);
    }
  }
}

/// \brief Constructs a list starting from first to last. Each element is converted using the TermConverter.
template <IsATerm Term, class Iter, class ATermConverter>
aterm make_list_forward(Iter first, Iter last, ATermConverter convert_to_aterm)
{
  term_list<Term> result_list;
  make_list_forward<Term, Iter, ATermConverter>(result_list, first, last, convert_to_aterm);
  return result_list;
}

/// \brief Constructs a list starting from first to last. Each element is converted using the TermConverter.
template <IsATerm Term, class Iter, class ATermConverter>
inline void make_list_forward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm)
{
  while (first != last)
  {
    --last;
    result.push_front(convert_to_aterm(*last));
  }
}

/// \brief Constructs a list starting from first to last. Each element is converted using the TermConverter and inserted
/// 		   whenever TermFilter yields true for the converted element.
/// \details Will first store the converted elements in an array and then insert them into the list.
template <typename Term, class Iter, class ATermConverter, class ATermFilter>
aterm make_list_forward(Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter)
{
  term_list<Term> result_list;
  make_list_forward<Term, Iter, ATermConverter, ATermFilter>(result_list, first, last, convert_to_aterm, aterm_filter);
  return result_list;
}

template < class Term, typename ForwardTraversalIterator, class Transformer >
void make_list_forward_helper(term_list<Term>& result, ForwardTraversalIterator& p, const ForwardTraversalIterator last, Transformer transformer)
{
  assert(p!=last);
  make_term_appl(result, 
    detail::as_list(), 
    [&transformer, &p](Term& result)
      {
        if constexpr (mcrl2::utilities::is_applicable2<Transformer, Term&, const Term&>::value)   
        {
          transformer(reinterpret_cast<Term&>(result), *(p++));
        }
        else
        {
          reinterpret_cast<Term&>(result)=transformer(*(p++));
        }
      },
    [&transformer, &p, last](term_list<Term>& result)
      {
        if (p==last)
        {
          make_term_list(reinterpret_cast<term_list<Term>& >(result));
        }
        else 
        {
          make_list_forward_helper(reinterpret_cast<term_list<Term>& >(result), p, last, transformer);
        }
      });
}

/// \brief Constructs a list traversing the iterator from first to last, putting the result in place in the variable result. 
template <IsATerm Term, class Iter, class ATermConverter, class ATermFilter>
inline void make_list_forward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter)
{
  const std::size_t len = std::distance(first,last);
  if (first==last)
  {
    make_term_list(result); // Put the empty list in result.
    return;
  }
  else if (len < LengthOfShortList) // If the list is sufficiently short, use the stack.
  {
    make_list_forward_helper(result, first, last, convert_to_aterm);
  }
  else
  {
    // The list is very long. Reserve memory on the heap.
    std::vector<Term> buffer;
    buffer.reserve(len);
    for(; first != last; ++first)
    {
      if constexpr (mcrl2::utilities::is_applicable2<ATermConverter, Term&, const Term>::value)
      {
        buffer.emplace_back();
        convert_aterm(buffer.back(), *first);
      }
      else
      {
        buffer.emplace_back(convert_to_aterm(*first));
      }

    }

    make_term_list(result); // Put the empty list in result.
    for(typename std::vector<Term>::const_reverse_iterator i = buffer.rbegin(); i != buffer.rend(); ++i)
    {
      result.push_front(*i);
    }
  }
}

} // namespace detail

} // namespace atermpp

#endif // MCRL2_ATERMPP_DETAIL_ATERM_INT_H
