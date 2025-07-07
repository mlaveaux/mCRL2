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

namespace atermpp
{

template <class Term>
class term_list;

namespace detail
{

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

/// \brief Constructs a list starting from first to last. The iterators are traversed backwards and each element is
/// 		   converted using the TermConverter.
/// \details The functions make_list_backward and make_list_forward with three and four arguments are almost the same.
/// 			 The reason for this is that there is a 5% loss of speed of the toolset when merging these two functions.
///          This is caused by storing and protecting the intermediate value of the converted aterm. See Term t = convert_to_aterm(...).
template <typename Term, typename Iter, typename ATermConverter>
inline aterm make_list_backward(Iter first, Iter last, ATermConverter convert_to_aterm);

/// \brief Constructs a list starting from first to last where the result is put in result.
template <class Term, class Iter, class ATermConverter>
inline void make_list_backward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm);


/// \brief Constructs a list starting from first to last. The iterators are traversed backwards and each element is
/// 		   converted using the TermConverter and inserted whenever TermFilter yields true for the converted element.
template <typename Term, typename Iter, typename ATermConverter, typename ATermFilter>
inline aterm make_list_backward(Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter);

/// \brief Construct a list iterating from the last to the first element. Result is put in the variable result.
template <class Term, class Iter, class ATermConverter, class ATermFilter>
inline void make_list_backward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter);

/// \brief Constructs a list starting from first to last. Each element is converted using the TermConverter.
template <typename Term, class Iter, class ATermConverter>
aterm make_list_forward(Iter first, Iter last, ATermConverter convert_to_aterm);

/// \brief Constructs a list starting from first to last. Each element is converted using the TermConverter.
template <typename Term, class Iter, class ATermConverter>
inline void make_list_forward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm);

/// \brief Constructs a list starting from the first iterator element to the last. The result is put into the variable result. 
template <class Term, class Iter, class ATermConverter>
inline void make_list_backward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm);

/// \brief Constructs a list starting from first to last. Each element is converted using the TermConverter and inserted
/// 		   whenever TermFilter yields true for the converted element.
/// \details Will first store the converted elements in an array and then insert them into the list.
template <typename Term, class Iter, class ATermConverter, class ATermFilter>
aterm make_list_forward(Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter);

/// \brief Constructs a list traversing the iterator from first to last, putting the result in place in the variable result. 
template <class Term, class Iter, class ATermConverter, class ATermFilter>
inline void make_list_forward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm, ATermFilter aterm_filter);

} // namespace detail

} // namespace atermpp

#endif // MCRL2_ATERMPP_DETAIL_ATERM_INT_H
