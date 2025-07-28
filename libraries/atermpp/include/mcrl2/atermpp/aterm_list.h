// Author(s): Wieger Wesselink, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_LIST_H
#define MCRL2_ATERMPP_ATERM_LIST_H

#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/atermpp/concepts.h"
#include "mcrl2/atermpp/detail/aterm_list.h"
#include "mcrl2/atermpp/detail/aterm_list_iterator.h"
#include "mcrl2/atermpp/function_symbol.h"
#include "mcrl2/atermpp/detail/type_traits_impl.h"

#include <mcrl3_ffi.h>

#include <ranges>
#include <type_traits>

namespace atermpp
{

/// \brief A list of aterm objects.
template <typename Term>
class term_list : public aterm
{
protected:
  /// \brief Constructor for term lists from internally constructed terms delivered as reference.
  explicit term_list(mcrl3::ffi::unprotected_aterm_t t) noexcept
      : aterm(t)
  {
    assert(!defined() || type_is_list());
  }

public:
  /// The type of object, T stored in the term_list.
  using value_type = Term;

  /// Pointer to T.
  using pointer = Term*;

  /// Reference to T.
  using reference = Term&;

  /// Const reference to T.
  using const_reference = const Term&;

  /// An unsigned integral type.
  using size_type = std::size_t;

  /// A signed integral type.
  using difference_type = ptrdiff_t;

  /// Iterator used to iterate through an term_list.
  using iterator = term_list_iterator<Term>;

  /// Const iterator used to iterate through an term_list.
  using const_iterator = term_list_iterator<Term>;

  /// Const iterator used to iterate through an term_list.
  using const_reverse_iterator = reverse_term_list_iterator<Term>;

  /// \brief Default constructor. Creates an empty list.
  term_list() noexcept
      : aterm(mcrl3::ffi::term_empty_list())
  {}

  /// \brief Constructor from an aterm.
  /// \param t A list.
  explicit term_list(const aterm& t) noexcept
      : aterm(t)
  {
    // assert(!defined() || type_is_list());
    assert(type_is_list()); // A list should not be a default aterm.
  }

  /// \brief Copy constructor.
  /// \param t A list.
  term_list(const term_list<Term>& t) noexcept
      : aterm(t)
  {
    assert(!defined() || type_is_list());
  }

  /// \brief Move constructor.
  /// \param t A list.
  term_list(term_list<Term>&& t) noexcept
      : aterm(std::move(t))
  {
    assert(!defined() || type_is_list());
  }

  /// This class has user-declared copy constructor so declare copy and move assignment.
  term_list& operator=(const term_list& other) noexcept = default;
  term_list& operator=(term_list&& other) noexcept = default;

  /// \brief Creates a term_list with the elements from first to last.
  /// \details It is assumed that the range can be traversed from last to first.
  /// \param first The start of a range of elements.
  /// \param last The end of a range of elements.
  template <class Iter>
  explicit term_list(Iter first,
      Iter last,
      typename std::enable_if<std::is_base_of<std::bidirectional_iterator_tag,
          typename std::iterator_traits<Iter>::iterator_category>::value>::type* = nullptr)
      : aterm(detail::make_list_backward<Term, Iter, detail::do_not_convert_term<Term>>(first,
            last,
            detail::do_not_convert_term<Term>()))
  {
    assert(!defined() || type_is_list());
  }

  /// \brief Creates a term_list with the elements from first to last converting the elements before inserting.
  /// \details It is assumed that the range can be traversed from last to first. The operator () in the class
  ///          ATermConverter is applied to each element before inserting it in the list.
  /// \param first The start of a range of elements.
  /// \param last The end of a range of elements.
  /// \param convert_to_aterm A class with a () operation, which is applied to each element
  ///                   before it is put into the list.
  template <class Iter, class ATermConverter>
  explicit term_list(Iter first,
      Iter last,
      const ATermConverter& convert_to_aterm,
      typename std::enable_if<std::is_base_of<std::bidirectional_iterator_tag,
          typename std::iterator_traits<Iter>::iterator_category>::value>::type* = 0)
      : aterm(detail::make_list_backward<Term, Iter, ATermConverter>(first, last, convert_to_aterm))
  {
    assert(!defined() || type_is_list());
  }

  /// \brief Creates a term_list with the elements from first to last, converting and filtering the list.
  /// \details It is assumed that the range can be traversed from last to first. The operator () in the class
  ///          ATermConverter is applied to each element before inserting it in the list. Elements are only
  ///          inserted if the operator () of the class ATermFilter yields true when applied to such an element.
  /// \param first The start of a range of elements.
  /// \param last The end of a range of elements.
  /// \param convert_to_aterm A class with a () operation, which is applied to each element
  ///                   before it is put into the list.
  /// \param aterm_filter A class with an operator () that is used to determine whether elements can be inserted in the
  /// list.
  template <class Iter, class ATermConverter, class ATermFilter>
  explicit term_list(Iter first,
      Iter last,
      const ATermConverter& convert_to_aterm,
      const ATermFilter& aterm_filter,
      typename std::enable_if<std::is_base_of<std::bidirectional_iterator_tag,
          typename std::iterator_traits<Iter>::iterator_category>::value>::type* = 0)
      : aterm(detail::make_list_backward<Term, Iter, ATermConverter, ATermFilter>(first,
            last,
            convert_to_aterm,
            aterm_filter))
  {
    assert(!defined() || type_is_list());
  }

  /// \brief Creates a term_list from the elements from first to last.
  /// \details The range is traversed from first to last. This requires
  ///           to copy the elements internally, which is less efficient
  ///           than this function with random access iterators as arguments.
  /// \param first The start of a range of elements.
  /// \param last The end of a range of elements.
  template <class Iter>
  explicit term_list(Iter first,
      Iter last,
      typename std::enable_if<!std::is_base_of<std::bidirectional_iterator_tag,
          typename std::iterator_traits<Iter>::iterator_category>::value>::type* = nullptr)
      : aterm(detail::make_list_forward<Term, Iter, detail::do_not_convert_term<Term>>(first,
            last,
            detail::do_not_convert_term<Term>()))
  {
    assert(!defined() || type_is_list());
  }

  /// \brief Creates a term_list from the elements from first to last converting the elements before inserting.
  /// \details The range is traversed from first to last. This requires
  ///           to copy the elements internally, which is less efficient
  ///           than this function with random access iterators as arguments.
  ///           The operator () in the class
  ///           ATermConverter is applied to each element before inserting it in the list.
  /// \param first The start of a range of elements.
  /// \param last The end of a range of elements.
  /// \param convert_to_aterm A class with a () operation, whic is applied to each element
  ///                      before it is put into the list.
  template <class Iter, class ATermConverter>
  explicit term_list(Iter first,
      Iter last,
      const ATermConverter& convert_to_aterm,
      typename std::enable_if<!std::is_base_of<std::bidirectional_iterator_tag,
          typename std::iterator_traits<Iter>::iterator_category>::value>::type* = nullptr)
      : aterm(detail::make_list_forward<Term, Iter, ATermConverter>(first, last, convert_to_aterm))
  {
    assert(!defined() || type_is_list());
  }

  /// \brief Creates a term_list from the elements from first to last converting and filtering the elements before
  /// inserting. \details The range is traversed from first to last. This requires
  ///           to copy the elements internally, which is less efficient
  ///           than this function with random access iterators as arguments.
  ///           The operator () in the class ATermConverter is applied to
  ///           each element before inserting it in the list. Elements are only
  ///           inserted if the operator () of the class ATermFilter yields true when applied to such an element.
  /// \param first The start of a range of elements.
  /// \param last The end of a range of elements.
  /// \param convert_to_aterm A class with a () operation, whic is applied to each element
  ///                      before it is put into the list.
  /// \param aterm_filter A class with an operator () that is used to determine whether elements can be inserted in the
  /// list.
  template <class Iter, class ATermConverter, class ATermFilter>
  explicit term_list(Iter first,
      Iter last,
      const ATermConverter& convert_to_aterm,
      const ATermFilter& aterm_filter,
      typename std::enable_if<!std::is_base_of<std::random_access_iterator_tag,
          typename std::iterator_traits<Iter>::iterator_category>::value>::type* = nullptr)
      : aterm(detail::make_list_forward<Term, Iter, ATermConverter>(first, last, convert_to_aterm, aterm_filter))
  {
    assert(!defined() || type_is_list());
  }

  /// \brief Creates a term_list from the elements in the range.
  template <std::ranges::range R>
    requires std::is_convertible_v<std::ranges::range_value_t<R>, Term>
  explicit term_list(R&& r)
      : aterm(detail::make_list_forward<Term, std::ranges::iterator_t<R>, detail::do_not_convert_term<Term>>(
            std::ranges::begin(r),
            std::ranges::end(r),
            detail::do_not_convert_term<Term>()))
  {
    assert(!defined() || type_is_list());
  }

  /// \brief A constructor based on an initializer list.
  /// \details This constructor is not made explicit to conform to initializer lists in standard containers.
  /// \param init The initialiser list.
  term_list(std::initializer_list<Term> init)
      : aterm(detail::make_list_backward<Term,
            typename std::initializer_list<Term>::const_iterator,
            detail::do_not_convert_term<Term>>(init.begin(), init.end(), detail::do_not_convert_term<Term>()))
  {
    assert(!defined() || type_is_list());
  }

  /// \brief Returns the tail of the list.
  /// \return The tail of the list.
  const term_list<Term>& tail() const
  {
    assert(!empty());
    return m_term[1];
  }

  /// \brief Removes the first element of the list.
  void pop_front() { *this = m_term[1]; }

  /// \brief Returns the first element of the list.
  /// \return The term at the head of the list.
  const Term& front() const { return static_cast<const detail::_aterm_list<Term>&>(*m_term).head(); }

  /// \brief Inserts a new element at the beginning of the current list.
  /// \param el The term that is added.
  void push_front(const Term& el) {}

  /// \brief Construct and insert a new element at the beginning of the current list.
  /// \param el The term that is added.
  template <typename... Args>
  void emplace_front(Args&&... arguments)
  {}

  /// \brief Returns the size of the term_list.
  /// \details The complexity of this function is linear in the size of the list.
  /// \return The size of the list.
  [[nodiscard]] size_type size() const
  {
    std::size_t size = 0;
    for (const_iterator i = begin(); i != end(); ++i)
    {
      ++size;
    }
    return size;
  }

  /// \brief Returns true if the list's size is 0.
  /// \return True iff the list is empty.
  [[nodiscard]] bool empty() const { return type_is_empty_list(); }

  /// \brief Returns a const_iterator pointing to the beginning of the term_list.
  /// \return The beginning of the list.
  const_iterator begin() const { return const_iterator(*this); }

  /// \brief Returns a const_iterator pointing to the end of the term_list.
  /// \return The end of the list.
  const_iterator end() const { return const_iterator(unprotected_aterm(mcrl3::ffi::term_empty_list())); }

  /// \brief Returns a const_reverse_iterator pointing to the end of the term_list.
  /// \details This operator requires linear time and memory in the size of the list to yield the iterator.
  /// \return The end of the list.
  const_reverse_iterator rbegin() const { return const_reverse_iterator(m_term); }

  /// \brief Returns a const_iterator pointing to the end of the term_list.
  /// \return The end of the list.
  const_reverse_iterator rend() const { return const_reverse_iterator(); }

  /// \brief Returns the largest possible size of the term_list.
  /// \return The largest possible size of the list.
  [[nodiscard]] size_type max_size() const { return std::numeric_limits<std::size_t>::max(); }
};

/// \brief Make an empty list and put it in target;
/// \param target The variable to which the empty list is assigned.
template <class Term>
void make_term_list(term_list<Term>& target)
{
  target = atermpp::down_cast<term_list<Term>>(mcrl3::ffi::term_empty_list());
}

/// \brief Creates a term_list with the elements from first to last.
/// \details It is assumed that the range can be traversed from last to first.
/// \param target The variable to which the list is assigned.
/// \param first The start of a range of elements.
/// \param last The end of a range of elements.
template <class Term, class Iter>
void make_term_list(term_list<Term>& target,
    Iter first,
    Iter last,
    typename std::enable_if<std::is_base_of<std::bidirectional_iterator_tag,
        typename std::iterator_traits<Iter>::iterator_category>::value>::type* = nullptr)
{
  detail::make_list_backward<Term, Iter, detail::do_not_convert_term<Term>>(target,
      first,
      last,
      detail::do_not_convert_term<Term>());
  assert(!target.defined() || target.type_is_list());
}

/// \brief Creates a term_list with the elements from first to last converting the elements before inserting.
/// \details It is assumed that the range can be traversed from last to first. The operator () in the class
///          ATermConverter is applied to each element before inserting it in the list.
/// \param target The variable to which the list is assigned.
/// \param first The start of a range of elements.
/// \param last The end of a range of elements.
/// \param convert_to_aterm A class with a () operation, which is applied to each element
///                   before it is put into the list.
template <class Term, class Iter, class ATermConverter>
void make_term_list(term_list<Term>& target,
    Iter first,
    Iter last,
    const ATermConverter& convert_to_aterm,
    typename std::enable_if<std::is_base_of<std::bidirectional_iterator_tag,
        typename std::iterator_traits<Iter>::iterator_category>::value>::type* = 0)
{
  detail::make_list_backward<Term, Iter, ATermConverter>(target, first, last, convert_to_aterm);
  assert(!target.defined() || target.type_is_list());
}

/// \brief Creates a term_list with the elements from first to last, converting and filtering the list.
/// \details It is assumed that the range can be traversed from last to first. The operator () in the class
///          ATermConverter is applied to each element before inserting it in the list. Elements are only
///          inserted if the operator () of the class ATermFilter yields true when applied to such an element.
/// \param target The variable to which the list is assigned.
/// \param first The start of a range of elements.
/// \param last The end of a range of elements.
/// \param convert_to_aterm A class with a () operation, which is applied to each element
///                   before it is put into the list.
/// \param aterm_filter A class with an operator () that is used to determine whether elements can be inserted in the
/// list.
template <class Term, class Iter, class ATermConverter, class ATermFilter>
void make_term_list(term_list<Term>& target,
    Iter first,
    Iter last,
    const ATermConverter& convert_to_aterm,
    const ATermFilter& aterm_filter,
    typename std::enable_if<std::is_base_of<std::bidirectional_iterator_tag,
        typename std::iterator_traits<Iter>::iterator_category>::value>::type* = 0)
{
  detail::make_list_backward<Term, Iter, ATermConverter, ATermFilter>(target,
      first,
      last,
      convert_to_aterm,
      aterm_filter);
  assert(!target.defined() || target.type_is_list());
}

/// \brief Creates a term_list from the elements from first to last.
/// \details The range is traversed from first to last. This requires
///           to copy the elements internally, which is less efficient
///           than this function with random access iterators as arguments.
/// \param target The variable to which the list is assigned.
/// \param first The start of a range of elements.
/// \param last The end of a range of elements.
template <class Term, class Iter>
void make_term_list(term_list<Term>& target,
    Iter first,
    Iter last,
    typename std::enable_if<!std::is_base_of<std::bidirectional_iterator_tag,
        typename std::iterator_traits<Iter>::iterator_category>::value>::type* = nullptr)
{
  detail::make_list_forward<Term, Iter, detail::do_not_convert_term<Term>>(target,
      first,
      last,
      detail::do_not_convert_term<Term>());
  assert(!target.defined() || target.type_is_list());
}

/// \brief Creates a term_list from the elements from first to last converting the elements before inserting.
/// \details The range is traversed from first to last. This requires
///           to copy the elements internally, which is less efficient
///           than this function with random access iterators as arguments.
///           The operator () in the class
///           ATermConverter is applied to each element before inserting it in the list.
/// \param target The variable to which the list is assigned.
/// \param first The start of a range of elements.
/// \param last The end of a range of elements.
/// \param convert_to_aterm A class with a () operation, which is applied to each element
///                      before it is put into the list.
template <class Term, class Iter, class ATermConverter>
void make_term_list(term_list<Term>& target,
    Iter first,
    Iter last,
    const ATermConverter& convert_to_aterm,
    typename std::enable_if<!std::is_base_of<std::bidirectional_iterator_tag,
        typename std::iterator_traits<Iter>::iterator_category>::value>::type* = nullptr)
{
  detail::make_list_forward<Term, Iter, ATermConverter>(target, first, last, convert_to_aterm);
  assert(!target.defined() || target.type_is_list());
}

/// \brief Creates a term_list from the elements from first to last converting and filtering the elements before
/// inserting. \details The range is traversed from first to last. This requires
///           to copy the elements internally, which is less efficient
///           than this function with random access iterators as arguments.
///           The operator () in the class ATermConverter is applied to
///           each element before inserting it in the list. Elements are only
///           inserted if the operator () of the class ATermFilter yields true when applied to such an element.
/// \param target The variable to which the list is assigned.
/// \param first The start of a range of elements.
/// \param last The end of a range of elements.
/// \param convert_to_aterm A class with a () operation, which is applied to each element
///                      before it is put into the list.
/// \param aterm_filter A class with an operator () that is used to determine whether elements can be inserted in the
/// list.
template <class Term, class Iter, class ATermConverter, class ATermFilter>
void make_term_list(term_list<Term>& target,
    Iter first,
    Iter last,
    const ATermConverter& convert_to_aterm,
    const ATermFilter& aterm_filter,
    typename std::enable_if<!std::is_base_of<std::random_access_iterator_tag,
        typename std::iterator_traits<Iter>::iterator_category>::value>::type* = nullptr)
{
  detail::make_list_forward<Term, Iter, ATermConverter>(target, first, last, convert_to_aterm, aterm_filter);
  assert(!target.defined() || target.type_is_list());
}

/// \brief A constructor based on an initializer list.
/// \details This constructor is not made explicit to conform to initializer lists in standard containers.
/// \param target The variable to which the list is assigned.
/// \param init The initialiser list.
template <class Term>
void make_term_list(term_list<Term>& target, std::initializer_list<Term> init)
{
  target = detail::make_list_backward<Term,
      typename std::initializer_list<Term>::const_iterator,
      detail::do_not_convert_term<Term>>(init.begin(), init.end(), detail::do_not_convert_term<Term>());
  assert(!target.defined() || target.type_is_list());
}

/// \cond INTERNAL_DOCS
namespace detail
{

/// \brief Template specialization to make a term_list recognizable as a container type (see
///        type_traits.h and detail/type_traits_impl.h).
template <typename T>
struct is_container_impl<atermpp::term_list<T>> : public std::true_type
{};

} // namespace detail

/// \brief A term_list with elements of type aterm.
using aterm_list = term_list<aterm>;

/// \brief Returns the list with the elements in reversed order.
/// \param l A list.
/// \details This operator is linear in the size of the list.
/// \return The reversed list.
template <typename Term>
inline term_list<Term> reverse(const term_list<Term>& l)
{
  if (l.size() < 2)
  {
    return l;
  }
  term_list<Term> result;
  for (const Term& t : l)
  {
    result.push_front(t);
  }
  return result;
}

/// \brief Returns the list with the elements sorted according to given ordering which is by default the ordering of
/// addresses of terms. \param l A list. \param ordering An total orderings relation on Term, by default the ordering
/// relation on Terms. \details This operator has complexity nlog n where n is the size of the list. \return The sorted
/// list.
template <IsATerm Term>
inline term_list<Term> sort_list(
    const term_list<Term>& l,
    const std::function<bool(const Term&, const Term&)>& ordering
    = [](const Term& t1, const Term& t2) { return t1 < t2; })
{
  const std::size_t len = l.size();
  if (len <= 1)
  {
    return l;
  }

  // The resulting list
  term_list<Term> result;

  if (len < detail::LengthOfShortList)
  {
    // The list is short, use the stack for temporal storage.
    Term* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(Term, len);

    // Collect all elements of list in buffer.
    std::size_t j = 0;
    for (const Term& t : l)
    {
      new (buffer + j) Term(t); // A mcrl2 stack allocator does not handle construction by default.
      ++j;
    }

    std::sort(buffer, buffer + len, ordering);

    // Insert elements at the front of the list.
    while (j > 0)
    {
      j = j - 1;
      result.push_front(buffer[j]);
      buffer[j].~Term(); // Explicitly call the destructor, as an mCRL2 stack allocator does not do that itself. .
    }
  }
  else
  {
    // The list is long. Use the heap to store intermediate data.
    std::vector<Term> buffer;
    buffer.reserve(len);

    for (const Term& t : l)
    {
      buffer.push_back(t);
    }

    // Sort using a standard algorithm.
    std::sort(buffer.begin(), buffer.end(), ordering);

    // Insert elements at the front of the list
    for (typename std::vector<Term>::reverse_iterator i = buffer.rbegin(); i != buffer.rend(); ++i)
    {
      result.push_front(*i);
    }
  }
  return result;
}

/// \brief Returns the merged list sorted according to the <-operator, which is by default the ordering of addresses of
/// terms. \param l1 An ordered list. \param l2 Another ordered list. \param ordering An total orderings relation on
/// Term, by default the ordering relation on Terms. \details This operator is linear in the cumulative length of l1 and
/// l2. In debug mode it checks whether l1 and l2 are ordered. \return The sorted list.
template <IsATerm Term>
inline term_list<Term> merge_lists(
    const term_list<Term>& l1,
    const term_list<Term>& l2,
    const std::function<bool(const Term&, const Term&)>& ordering
    = [](const Term& t1, const Term& t2) { return t1 < t2; })
{
  const std::size_t len1 = l1.size();
  const std::size_t len2 = l2.size();
  if (len1 == 0)
  {
    assert(l2 == sort_list(l2, ordering));
    return l2;
  }
  if (len2 == 0)
  {
    assert(l1 == sort_list(l1, ordering));
    return l1;
  }

  // The resulting list
  term_list<Term> result;

  term_list<Term> i1 = l1;
  term_list<Term> i2 = l2;
  if (len1 + len2 < detail::LengthOfShortList)
  {
    // The list is short, use the stack for temporal storage.
    Term* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(Term, len1 + len2);

    // Collect all elements of list in buffer.
    std::size_t j = 0;
    while (!i1.empty() && !i2.empty())
    {
      if (ordering(i1.front(), i2.front()))
      {
        new (buffer + j) Term(i1.front());
        i1.pop_front();
      }
      else
      {
        new (buffer + j) Term(i2.front());
        i2.pop_front();
      }
      ++j;
    }
    if (i1.empty())
    {
      result = i2;
    }
    else
    {
      assert(i2.empty());
      result = i1;
    }

    // Insert elements at the front of the list.
    while (j > 0)
    {
      j = j - 1;
      result.push_front(buffer[j]);
      buffer[j].~Term(); // Explicitly call the destructor, as an mCRL2 stack allocator does not do that itself. .
    }
  }
  else
  {
    // The list is long. Use the heap to store intermediate data.
    std::vector<Term> buffer;
    buffer.reserve(len1 + len2);

    while (!i1.empty() && !i2.empty())
    {
      if (ordering(i1.front(), i2.front()))
      {
        buffer.push_back(i1.front());
        i1.pop_front();
      }
      else
      {
        buffer.push_back(i2.front());
        i2.pop_front();
      }
    }
    if (i1.empty())
    {
      result = i2;
    }
    else
    {
      assert(i2.empty());
      result = i1;
    }

    // Insert elements at the front of the list
    for (typename std::vector<Term>::reverse_iterator i = buffer.rbegin(); i != buffer.rend(); ++i)
    {
      result.push_front(*i);
    }
  }
  assert(result.size() == len1 + len2);
  assert(result == sort_list(result, ordering));
  return result;
}

/// \brief Returns the concatenation of two lists with convertible element types.
///  \details The type of the result is either the type of l, if the elements of m
///           can be converted implicitly to the type of the elements of l. Otherwise if the
///           elements of l can be converted implicitly to the type of the elements
///           of m, the result type is that or m.
/// \param l A list.
/// \param m A list.
/// \details The complexity of this operator is linear in the length of l.
/// \return The concatenation of the lists l followed by m.

template <IsATerm Term1, IsATerm Term2>
inline typename std::conditional<std::is_convertible<Term2, Term1>::value, term_list<Term1>, term_list<Term2>>::type
operator+(const term_list<Term1>& l, const term_list<Term2>& m)
{
  static_assert(std::is_convertible<Term1, Term2>::value || std::is_convertible<Term2, Term1>::value,
      "Concatenated lists must be of convertible types. ");
  typedef typename std::conditional<std::is_convertible<Term2, Term1>::value, Term1, Term2>::type ResultType;
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
  if (len < detail::LengthOfShortList)
  {
    // The length is short. Use the stack for temporary storage.
    const_iterator* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(const_iterator, len);

    std::size_t j = 0;
    for (const_iterator i = l.begin(); i != l.end(); ++i, ++j)
    {
      buffer[j] = i;
    }
    assert(j == len);

    // Insert elements at the front of the list
    while (j > 0)
    {
      j = j - 1;
      result.push_front(*buffer[j]);
    }
  }
  else
  {
    // The length of l is very long. Use the heap for temporary storage.
    std::vector<ResultType> buffer;
    buffer.reserve(len);

    for (const Term1& t : l)
    {
      buffer.push_back(t);
    }

    // Insert elements at the front of the list
    for (typename std::vector<ResultType>::const_reverse_iterator i = buffer.rbegin(); i != buffer.rend(); ++i)
    {
      result.push_front(*i);
    }
  }
  return result;
}

/// \brief Appends a new element at the end of the list. Note
///        that the complexity of this function is O(n), with n the number of
///        elements in the list!!!
/// \param l The list to which the term is appended.
/// \param el A term.
/// \return The list l with elem appended at the end.
template <typename Term>
inline term_list<Term> push_back(const term_list<Term>& l, const Term& el)
{
  typedef typename term_list<Term>::const_iterator const_iterator;

  const std::size_t len = l.size();

  // The resulting list
  term_list<Term> result;
  result.push_front(el);

  if (len < detail::LengthOfShortList)
  {
    // The list is short, use the stack for temporal storage.
    const_iterator* buffer = MCRL2_SPECIFIC_STACK_ALLOCATOR(const_iterator, len);

    // Collect all elements of list in buffer.
    std::size_t j = 0;
    for (const_iterator i = l.begin(); i != l.end(); ++i, ++j)
    {
      buffer[j] = i;
    }

    // Insert elements at the front of the list.
    while (j > 0)
    {
      j = j - 1;
      result.push_front(*buffer[j]);
    }
  }
  else
  {
    // The list is long. Use the heap to store intermediate data.
    std::vector<Term> buffer;
    buffer.reserve(len);

    for (const Term& t : l)
    {
      buffer.push_back(t);
    }

    // Insert elements at the front of the list
    for (typename std::vector<Term>::reverse_iterator i = buffer.rbegin(); i != buffer.rend(); ++i)
    {
      result.push_front(*i);
    }
  }

  return result;
}

/// \brief Converts the given term list to a vector.
template <typename T>
std::vector<T> as_vector(const atermpp::term_list<T>& x)
{
  return std::vector<T>(x.begin(), x.end());
}

/// \brief Converts the given term list to a set.
template <typename T>
std::set<T> as_set(const atermpp::term_list<T>& x)
{
  return std::set<T>(x.begin(), x.end());
}

} // namespace atermpp

namespace std
{
//
/// \brief Swaps two term_lists.
/// \details This operation is more efficient than exchanging terms by an assignment,
///          as swapping does not require to change the protection of terms.
/// \param t1 The first term
/// \param t2 The second term
template <class T>
inline void swap(atermpp::term_list<T>& t1, atermpp::term_list<T>& t2) noexcept
{
  t1.swap(t2);
}

/// \brief The standard hash class.
template <class Term>
struct hash<atermpp::term_list<Term>>
{
  /// \brief A specialization of the standard std::hash function.
  /// \param l The list for which a hash value is calculated.
  /// \return A hash value for l.
  std::size_t operator()(const atermpp::term_list<Term>& l) const
  {
    std::hash<atermpp::aterm> hasher;
    return hasher(l);
  }
};

} // namespace std

#endif // MCRL2_ATERMPP_ATERM_LIST_H
