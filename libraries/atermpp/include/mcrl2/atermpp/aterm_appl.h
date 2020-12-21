// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/aterm_appl.h
/// \brief The term_appl class represents function application.

#ifndef MCRL2_ATERMPP_ATERM_APPL_H
#define MCRL2_ATERMPP_ATERM_APPL_H

#include "mcrl2/atermpp/detail/aterm_appl_iterator.h"
#include "mcrl2/atermpp/detail/aterm_list.h"
#include "mcrl2/atermpp/detail/global_aterm_pool.h"
#include "mcrl2/utilities/type_traits.h"

namespace atermpp
{

template <class Term>
class term_appl: public aterm
{
protected:
  /// \brief Constructor.
  /// \param t A pointer internal data structure from which the term is constructed.
  /// \details This function is explicitly protected such that is not used in common code. 
  explicit term_appl(detail::_term_appl *t)
   : aterm(reinterpret_cast<detail::_aterm*>(t))
  {
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(std::size_t),"Term derived from an aterm must not have extra fields");
  }

  /// \brief Explicit constructor from an aterm.
  /// \param t The aterm from which the term is constructed.
  /// \details This function is explicitly protected, to avoid its use in general code. 
  explicit term_appl(const aterm& t) 
   : aterm(t)
  {
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(std::size_t),"Term derived from an aterm must not have extra fields");
  } 

public:
  /// The type of object, T stored in the term_appl.
  typedef Term value_type;

  /// Pointer to T.
  typedef Term* pointer;

  /// Reference to T.
  typedef Term& reference;

  /// Const reference to T.
  typedef const Term const_reference;

  /// An unsigned integral type.
  typedef std::size_t size_type;

  /// A signed integral type.
  typedef ptrdiff_t difference_type;

  /// Iterator used to iterate through an term_appl.
  typedef term_appl_iterator<Term> iterator;

  /// Const iterator used to iterate through an term_appl.
  typedef term_appl_iterator<Term> const_iterator;

  /// \brief Default constructor.
  term_appl():aterm()
  {}

  /// This class has user-declared copy constructor so declare default copy and move operators.
  term_appl(const term_appl& other) noexcept = default;
  term_appl& operator=(const term_appl& other) noexcept = default;
  term_appl(term_appl&& other) noexcept = default;
  term_appl& operator=(term_appl&& other) noexcept = default;

  /// \details The iterator range is traversed more than once. If only one traversal is required
  ///          use term_appl with a TermConverter argument. But this function
  ///          is substantially less efficient.
  ///          The length of the iterator range should must match the arity of the function symbol.
  /// \param sym A function symbol.
  /// \param begin The start of a range of elements.
  /// \param end The end of a range of elements.
  template <class ForwardIterator,
            typename std::enable_if<mcrl2::utilities::is_iterator<ForwardIterator>::value>::type* = nullptr,
            typename std::enable_if<!std::is_same<typename ForwardIterator::iterator_category, std::input_iterator_tag>::value>::type* = nullptr,
            typename std::enable_if<!std::is_same<typename ForwardIterator::iterator_category, std::output_iterator_tag>::value>::type* = nullptr>
  term_appl(const function_symbol& sym,
            ForwardIterator begin,
            ForwardIterator end)
  {
    detail::g_term_pool().create_appl_dynamic(*this, sym, begin, end);
    static_assert((std::is_base_of<aterm, Term>::value),"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(std::size_t),"Term derived from an aterm must not have extra fields");
    static_assert(!std::is_same<typename ForwardIterator::iterator_category, std::input_iterator_tag>::value,
                  "A forward iterator has more requirements than an input iterator.");
    static_assert(!std::is_same<typename ForwardIterator::iterator_category, std::output_iterator_tag>::value,
                  "A forward iterator has more requirements than an output iterator.");
  }

  /// \details The given iterator is traversed only once. So it can be used with an input iterator.
  ///          This means that the TermConverter is applied exactly once to each element.
  ///          The length of the iterator range must be equal to the arity of the function symbol.
  /// \param sym A function symbol.
  /// \param begin The start of a range of elements.
  /// \param end The end of a range of elements.
  template <class InputIterator,
            typename std::enable_if<mcrl2::utilities::is_iterator<InputIterator>::value>::type* = nullptr,
            typename std::enable_if<std::is_same<typename InputIterator::iterator_category, std::input_iterator_tag>::value>::type* = nullptr>
  term_appl(const function_symbol& sym,
            InputIterator begin,
            InputIterator end)
    : term_appl(sym, begin, end, [](const Term& term) -> const Term& { return term; } )
  {
    static_assert((std::is_base_of<aterm, Term>::value),"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(std::size_t),"Term derived from an aterm must not have extra fields");
    static_assert(std::is_same<typename InputIterator::iterator_category, std::input_iterator_tag>::value,
                  "The InputIterator is missing the input iterator tag.");
  }

  /// \details The given iterator is traversed only once. So it can be used with an input iterator.
  ///          This means that the TermConverter is applied exactly once to each element.
  ///          The length of the iterator range must be equal to the arity of the function symbol.
  /// \param sym A function symbol.
  /// \param begin The start of a range of elements.
  /// \param end The end of a range of elements.
  /// \param converter An class or lambda term containing an operator Term operator()(const Term& t) which is
  ///        applied to each each element in the iterator range before it becomes an argument of this term.
  template <class InputIterator,
            class TermConverter,
            typename std::enable_if<mcrl2::utilities::is_iterator<InputIterator>::value>::type* = nullptr>
  term_appl(const function_symbol& sym,
            InputIterator begin,
            InputIterator end,
            TermConverter converter)
  {
    detail::g_term_pool().create_appl_dynamic(*this, sym, converter, begin, end);
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(std::size_t),"Term derived from an aterm must not have extra fields");
    static_assert(!std::is_same<typename InputIterator::iterator_category, std::output_iterator_tag>::value,
                  "The InputIterator has the output iterator tag.");
  }

  /// \brief Constructor.
  /// \param sym A function symbol.
  term_appl(const function_symbol& sym)
  {
    detail::g_term_pool().create_term(*this, sym);
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(std::size_t),"Term derived from an aterm must not have extra fields");
  }

  /// \brief Constructor for n-arity function application.
  /// \param symbol A function symbol.
  /// \param arguments The arguments of the function application.
  template<typename ...Terms>
  term_appl(const function_symbol& symbol, const Terms& ...arguments)
  {
    detail::g_term_pool().create_appl(*this, symbol, arguments...);
    static_assert(detail::are_terms<Terms...>::value, "Arguments of function application should be terms.");
    static_assert(std::is_base_of<aterm, Term>::value,"Term must be derived from an aterm");
    static_assert(sizeof(Term)==sizeof(std::size_t),"Term derived from an aterm must not have extra fields");
  }

  /// \brief Returns the function symbol belonging to an aterm_appl.
  /// \return The function symbol of this term.
  const function_symbol& function() const
  {
    return m_term->function();
  }

  /// \brief Returns the number of arguments of this term.
  /// \return The number of arguments of this term.
  size_type size() const
  {
    return m_term->function().arity();
  }

  /// \brief Returns true if the term has no arguments.
  /// \return True if this term has no arguments.
  bool empty() const
  {
    return size()==0;
  }

  /// \brief Returns an iterator pointing to the first argument.
  /// \return An iterator pointing to the first argument.
  const_iterator begin() const
  {
    return const_iterator(reinterpret_cast<const Term*>(&(reinterpret_cast<const detail::_term_appl*>(m_term)->arg(0))));
  }

  /// \brief Returns a const_iterator pointing past the last argument.
  /// \return A const_iterator pointing past the last argument.
  const_iterator end() const
  {
    return const_iterator(reinterpret_cast<const Term*>(&reinterpret_cast<const detail::_term_appl*>(m_term)->arg(size())));
  }

  /// \brief Returns the largest possible number of arguments.
  /// \return The largest possible number of arguments.
  constexpr size_type max_size() const
  {
    return std::numeric_limits<size_type>::max();
  }

  /// \brief Returns the i-th argument.
  /// \param i A positive integer.
  /// \return The argument with the given index.
  const Term& operator[](const size_type i) const
  {
    assert(i < size()); // Check the bounds.
    return reinterpret_cast<const detail::_term_appl*>(m_term)->arg(i);
  }
};

typedef term_appl<aterm> aterm_appl;
} // namespace atermpp

namespace std
{

/// \brief Swaps two term_applss.
/// \details This operation is more efficient than exchanging terms by an assignment,
///          as swapping does not require to change the protection of terms.
/// \param t1 The first term.
/// \param t2 The second term.
template <class T>
inline void swap(atermpp::term_appl<T>& t1, atermpp::term_appl<T>& t2) noexcept
{
  t1.swap(t2);
}

/// \brief Standard hash function.
template<class T>
struct hash<atermpp::term_appl<T> >
{
  std::size_t operator()(const atermpp::term_appl<T>& t) const
  {
    return std::hash<atermpp::aterm>()(t);
  }
};
} // namespace std

#endif // MCRL2_ATERMPP_ATERM_APPL_H
