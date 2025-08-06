// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/detail/algorithm_impl.h
/// \brief Implementations of algorithms.

#ifndef MCRL2_ATERMPP_DETAIL_ALGORITHM_IMPL_H
#define MCRL2_ATERMPP_DETAIL_ALGORITHM_IMPL_H

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm_list.cxx"
#else
  import atermpp;
#endif

namespace atermpp::detail
{
/// \brief Applies the function f to all children of a.
/// \param a A term
/// \param f A function on terms
/// \return The transformed term
template <typename Function>
aterm appl_apply(const aterm& a, const Function f)
{
  return aterm(a.function(), a.begin(), a.end(), f);
}

//--- find ----------------------------------------------------------------//

// we need to use our own traits classes to extract the value type from insert iterators
template <class Iterator>
struct iterator_value
{
  using type = typename std::iterator_traits<Iterator>::value_type;
};

template <class Container>
struct iterator_value<std::insert_iterator<Container> >
{
  using type = typename Container::value_type;
};

template <class Container>
struct iterator_value<std::back_insert_iterator<Container> >
{
  using type = typename Container::value_type;
};

template <class Container>
struct iterator_value<std::front_insert_iterator<Container> >
{
  using type = typename Container::value_type;
};

/// \internal
// used to abort the recursive find
struct found_term_exception
{
  aterm t;

  found_term_exception(const aterm& t_)
    : t(t_)
  {}
};

/// \brief Implements the for_each algorithm
/// \param t A term
/// \param op A unary function on terms
/// \return The result of the algorithm
template <typename UnaryFunction>
UnaryFunction for_each_impl(aterm t, UnaryFunction op)
{
  if (t.type_is_list())
  {
    const aterm_list& l = down_cast<aterm_list>(t);
    for (const aterm& x: l)
    {
      for_each_impl(x, op);
    }
  }
  else if (t.type_is_appl())
  {
    if (op(t))
    {
      for (const aterm& x: t)
      {
        for_each_impl(x, op);
      }
    }
  }
  return op;
}

/// \brief Implements the find_if algorithm
/// If the term t is found, it is stored in output and true is returned
/// \param t A term
/// \param match A predicate function on terms
/// \param output The variable to store the match in
/// \return true if a match was found, false otherwise
template <typename MatchPredicate>
bool find_if_impl(const aterm& t, MatchPredicate match, aterm& output)
{
  if (t.type_is_appl())
  {
    if (match(t))
    {
      output = t;
      return true;
    }
    for (const aterm& x: t)
    {
      if (find_if_impl(x, match, output))
      {
        return true;
      }
    }
  }
  else if (t.type_is_list())
  {
    const aterm_list& l = down_cast<aterm_list>(t);
    for (const aterm& x: l)
    {
      if (find_if_impl(x, match, output))
      {
        return true;
      }
    }
  }
  return false;
}

/// \brief Implements the find_all_if algorithm
/// \param t A term
/// \param op A predicate function on terms
/// \param destBegin The beginning of a range to where the results are written
template <typename MatchPredicate, typename OutputIterator>
void find_all_if_impl(const aterm& t, MatchPredicate op, OutputIterator& destBegin)
{
  using value_type = typename iterator_value<OutputIterator>::type;

  if (t.type_is_list())
  {
    const aterm_list& l = down_cast<aterm_list>(t);
    for (const aterm& x: l)
    {
      find_all_if_impl<MatchPredicate>(x, op, destBegin);
    }
  }
  else if (t.type_is_appl())
  {
    if (op(t))
    {
      *destBegin++ = vertical_cast<value_type>(t);
    }
    for (const aterm& x: t)
    {
      find_all_if_impl<MatchPredicate>(x, op, destBegin);
    }
  }
  else
  {
    return;
  }
}

//--- partial find --------------------------------------------------------//

/// \brief Implements the partial_find_if_impl algorithm
/// \param t A term
/// \param match A predicate function on terms
/// \param stop A predicate function on terms
template <typename MatchPredicate, typename StopPredicate>
aterm partial_find_if_impl(const aterm& t, MatchPredicate match, StopPredicate stop)
{
  if (t.type_is_appl())
  {
    if (match(t))
    {
      return t; // report the match
    }
    if (stop(t))
    {
      return aterm(); // nothing was found
    }
    for (const aterm& x: t)
    {
      aterm result = partial_find_if_impl<MatchPredicate, StopPredicate>(x, match, stop);
      if (result != aterm())
      {
        return result;
      }
    }
  }

  if (t.type_is_list())
  {
    const aterm_list& l = down_cast<aterm_list>(t);
    for (const aterm& x: l)
    {
      aterm result = partial_find_if_impl<MatchPredicate, StopPredicate>(x, match, stop);
      if (result != aterm())
      {
        return result;
      }
    }
  }
  return aterm();
}

/// \brief Implements the partial_find_all_if algorithm
/// \param t A term
/// \param match A predicate function on terms
/// \param stop A predicate function on terms
/// \param destBegin The beginning of a range to where the results are written
template <typename MatchPredicate, typename StopPredicate, typename OutputIterator>
void partial_find_all_if_impl(const aterm& t, MatchPredicate match, StopPredicate stop, OutputIterator& destBegin)
{
  if (t.type_is_appl())
  {
    if (match(t))
    {
      *destBegin++ = t;
    }
    if (stop(t))
    {
      return;
    }
    for (const aterm& x: t)
    {
      partial_find_all_if_impl<MatchPredicate, StopPredicate>(x, match, stop, destBegin);
    }
  }

  if (t.type_is_list())
  {
    const aterm_list& l = down_cast<aterm_list>(t);
    for (const aterm& x: l)
    {
      partial_find_all_if_impl<MatchPredicate, StopPredicate>(x, match, stop, destBegin);
    }
  }
}

} // namespace atermpp::detail



#endif // MCRL2_ATERMPP_DETAIL_ALGORITHM_IMPL_H
