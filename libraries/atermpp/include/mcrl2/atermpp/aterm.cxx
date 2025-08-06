// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/aterm.h
/// \brief The term_appl class represents function application.

#ifndef MCRL2_ATERMPP_ATERM_APPL_H
#define MCRL2_ATERMPP_ATERM_APPL_H
MCRL2_MODULE;

#include <cassert>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm_core.cxx"
  #include "mcrl2/atermpp/concepts.cxx"
  #include "mcrl2/atermpp/detail/aterm_appl_iterator.cxx"
  #include "mcrl2/atermpp/detail/aterm_core.cxx"
  #include "mcrl2/atermpp/detail/aterm_list_data.cxx"
  #include "mcrl2/atermpp/detail/global_aterm_pool.cxx"
  #include "mcrl2/utilities/type_traits.cxx"
#else
  export module atermpp:aterm;

  import :aterm_core;
  import :concepts;
  import :function_symbol;
  import :detail.thread_aterm_pool;
  import :detail.aterm_appl_iterator;
  import :detail.aterm_data;
  import :detail.aterm_core_data;
  import :detail.aterm_list_data; 

  import utilities;
#endif

MCRL2_MODULE_EXPORT namespace atermpp
{

class aterm : public aterm_core
{
protected:
  /// \brief Constructor.
  /// \param t A pointer internal data structure from which the term is constructed.
  /// \details This function is explicitly protected such that is not used in common code.
  explicit aterm(detail::_term_appl* t)
      : aterm_core(reinterpret_cast<detail::_aterm*>(t))
  {}

public:
  /// An unsigned integral type.
  using size_type = std::size_t;

  /// A signed integral type.
  using difference_type = ptrdiff_t;

  /// Iterator used to iterate through an term_appl.
  using iterator = term_appl_iterator<aterm>;

  /// Const iterator used to iterate through an term_appl.
  using const_iterator = term_appl_iterator<aterm>;

  /// \brief Default constructor.
  aterm()
      : aterm_core()
  {}

  /// This class has user-declared copy constructor so declare default copy and move operators.
  aterm(const aterm& other) noexcept = default;
  aterm& operator=(const aterm& other) noexcept = default;
  aterm(aterm&& other) noexcept = default;
  aterm& operator=(aterm&& other) noexcept = default;

  /// \brief Constructor that provides an aterm based on a function symbol and forward iterator providing the arguments.
  /// \details The iterator range is traversed more than once. If only one traversal is required
  ///          use term_appl with a TermConverter argument. But this function
  ///          is substantially less efficient.
  ///          The length of the iterator range must match the arity of the function symbol.
  /// \param sym A function symbol.
  /// \param begin The start of a range of elements.
  /// \param end The end of a range of elements.
  template <class ForwardIterator>
  aterm(const function_symbol& sym, ForwardIterator begin, ForwardIterator end)
    requires (mcrl2::utilities::is_iterator<ForwardIterator>::value
          && !std::is_same_v<typename ForwardIterator::iterator_category, std::input_iterator_tag>
          && !std::is_same_v<typename ForwardIterator::iterator_category, std::output_iterator_tag>)
  {
    detail::g_thread_term_pool().create_appl_dynamic(*this, sym, begin, end);
  }

  /// \brief Constructor that provides an aterm based on a function symbol and an input iterator providing the
  /// arguments. \details The given iterator is traversed only once. So it can be used with an input iterator.
  ///          This means that the TermConverter is applied exactly once to each element.
  ///          The length of the iterator range must be equal to the arity of the function symbol.
  /// \param sym A function symbol.
  /// \param begin The start of a range of elements.
  /// \param end The end of a range of elements.
  template <class InputIterator>
  aterm(const function_symbol& sym, InputIterator begin, InputIterator end)
    requires mcrl2::utilities::is_iterator<InputIterator>::value
          && std::is_same_v<typename InputIterator::iterator_category, std::input_iterator_tag>
      : aterm(sym, begin, end, [](const unprotected_aterm_core& term) -> const unprotected_aterm_core& { return term; })
  {
    static_assert(std::is_same_v<typename InputIterator::iterator_category, std::input_iterator_tag>,
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
  template <class InputIterator, class TermConverter>
  aterm(const function_symbol& sym, InputIterator begin, InputIterator end, TermConverter converter)
    requires mcrl2::utilities::is_iterator<InputIterator>::value
  {
    detail::g_thread_term_pool().create_appl_dynamic(*this, sym, converter, begin, end);
    static_assert(!std::is_same_v<typename InputIterator::iterator_category, std::output_iterator_tag>,
        "The InputIterator has the output iterator tag.");
  }

  /// \brief Constructor.
  /// \param sym A function symbol.
  aterm(const function_symbol& sym) { detail::g_thread_term_pool().create_term(*this, sym); }

  /// \brief Constructor for n-arity function application.
  /// \param symbol A function symbol.
  /// \param arguments The arguments of the function application.
  template <typename... Terms>
  aterm(const function_symbol& symbol, const Terms&... arguments)
  {
    detail::g_thread_term_pool().create_appl(*this, symbol, arguments...);
  }

  /// \brief Returns the function symbol belonging to an aterm.
  /// \return The function symbol of this term.
  const function_symbol& function() const { return m_term->function(); }

  /// \brief Returns the number of arguments of this term.
  /// \return The number of arguments of this term.
  size_type size() const { return m_term->function().arity(); }

  /// \brief Returns true if the term has no arguments.
  /// \return True if this term has no arguments.
  bool empty() const { return size() == 0; }

  /// \brief Returns an iterator pointing to the first argument.
  /// \return An iterator pointing to the first argument.
  const_iterator begin() const
  {
    return const_iterator(&static_cast<const aterm&>(reinterpret_cast<const detail::_term_appl*>(m_term)->arg(0)));
  }

  /// \brief Returns a const_iterator pointing past the last argument.
  /// \return A const_iterator pointing past the last argument.
  const_iterator end() const
  {
    return const_iterator(&static_cast<const aterm&>(reinterpret_cast<const detail::_term_appl*>(m_term)->arg(size())));
  }

  /// \brief Returns the largest possible number of arguments.
  /// \return The largest possible number of arguments.
  constexpr size_type max_size() const { return std::numeric_limits<size_type>::max(); }

  /// \brief Returns the i-th argument.
  /// \param i A positive integer.
  /// \return The argument with the given index.
  const aterm& operator[](const size_type i) const
  {
    assert(i < size()); // Check the bounds.
    return static_cast<const aterm&>(reinterpret_cast<const detail::_term_appl*>(m_term)->arg(i));
  }
};

using term_callback = void (*)(const aterm&);

extern void add_deletion_hook(const function_symbol&, term_callback);

/// \brief Constructor an aterm in a variable based on a function symbol and an forward iterator providing the
/// arguments. \details The iterator range is traversed more than once. If only one traversal is required
///          use term_appl with a TermConverter argument. But this function
///          is substantially less efficient.
///          The length of the iterator range must match the arity of the function symbol.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param sym A function symbol.
/// \param begin The start of a range of elements.
/// \param end The end of a range of elements.
template <class Term,
    class ForwardIterator,
    std::enable_if_t<mcrl2::utilities::is_iterator<ForwardIterator>::value>* = nullptr,
    std::enable_if_t<!std::is_same_v<typename ForwardIterator::iterator_category, std::input_iterator_tag>>* = nullptr>
void make_term_appl(Term& target, const function_symbol& sym, ForwardIterator begin, ForwardIterator end)
  requires (mcrl2::utilities::is_iterator<ForwardIterator>::value
        && !std::is_same_v<typename ForwardIterator::iterator_category, std::input_iterator_tag>
        && !std::is_same_v<typename ForwardIterator::iterator_category, std::output_iterator_tag>)
{
  detail::g_thread_term_pool().create_appl_dynamic(target, sym, begin, end);

  static_assert((std::is_base_of_v<aterm, Term>), "Term must be derived from an aterm");
  static_assert(sizeof(Term) == sizeof(std::size_t), "Term derived from an aterm must not have extra fields");
}

/// \brief Constructor an aterm in a variable based on a function symbol and an input iterator providing the arguments.
/// \details The given iterator is traversed only once. So it can be used with an input iterator.
///          This means that the TermConverter is applied exactly once to each element.
///          The length of the iterator range must be equal to the arity of the function symbol.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param sym A function symbol.
/// \param begin The start of a range of elements.
/// \param end The end of a range of elements.
template <class Term, class InputIterator>
void make_term_appl(Term& target, const function_symbol& sym, InputIterator begin, InputIterator end)
  requires mcrl2::utilities::is_iterator<InputIterator>::value
        && std::is_same_v<typename InputIterator::iterator_category, std::input_iterator_tag>
{
  make_term_appl(target, sym, begin, end, [](const Term& term) -> const Term& { return term; });

  static_assert((std::is_base_of_v<aterm, Term>), "Term must be derived from an aterm");
  static_assert(sizeof(Term) == sizeof(std::size_t), "Term derived from an aterm must not have extra fields");
  static_assert(std::is_same_v<typename InputIterator::iterator_category, std::input_iterator_tag>,
      "The InputIterator is missing the input iterator tag.");
}

/// \brief Constructor an aterm in a variable based on a function symbol and an forward iterator providing the
/// arguments. \details The given iterator is traversed only once. So it can be used with an input iterator.
///          This means that the TermConverter is applied exactly once to each element.
///          The length of the iterator range must be equal to the arity of the function symbol.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param sym A function symbol.
/// \param begin The start of a range of elements.
/// \param end The end of a range of elements.
/// \param converter An class or lambda term containing an operator Term operator()(const Term& t) which is
///        applied to each each element in the iterator range before it becomes an argument of this term.
template <class Term, class InputIterator, class TermConverter>
void make_term_appl(Term& target,
    const function_symbol& sym,
    InputIterator begin,
    InputIterator end,
    TermConverter converter)
  requires mcrl2::utilities::is_iterator<InputIterator>::value
{
  detail::g_thread_term_pool().create_appl_dynamic(target, sym, converter, begin, end);

  static_assert(std::is_base_of_v<aterm, Term>, "Term must be derived from an aterm");
  static_assert(sizeof(Term) == sizeof(std::size_t), "Term derived from an aterm must not have extra fields");
}

/// \brief Make an term_appl consisting of a single function symbol.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param sym A function symbol.
template <class Term>
void make_term_appl(Term& target, const function_symbol& sym)
{
  detail::g_thread_term_pool().create_term(target, sym);

  static_assert(std::is_base_of_v<aterm, Term>, "Term must be derived from an aterm");
  static_assert(sizeof(Term) == sizeof(std::size_t), "Term derived from an aterm must not have extra fields");
}

/// \brief Make an aterm application for n-arity function application.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param symbol A function symbol.
/// \param arguments The arguments of the function application.
template <class Term, typename... Terms>
void make_term_appl(Term& target, const function_symbol& symbol, const Terms&... arguments)
{
  detail::g_thread_term_pool().create_appl(target, symbol, arguments...);
}

/// \brief Constructor for n-arity function application with an index.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param symbol A function symbol.
/// \param arguments The arguments of the function application.
template <class Term, class INDEX_TYPE, typename... Terms>
void make_term_appl_with_index(aterm& target, const function_symbol& symbol, const Terms&... arguments)
{
  detail::g_thread_term_pool().create_appl_index<Term, INDEX_TYPE>(target, symbol, arguments...);
}

/// \brief A universal cast from an aterm to another aterm that is convertible in either direction. Less strict than vertical_cast.
/// \param  t A term of a type inheriting from an aterm.
/// \return  A term of type const Derived&.
template <IsATerm Derived, IsATerm Base>
  requires std::is_convertible_v<std::remove_reference_t<Base>, std::remove_reference_t<Derived>>
        || std::is_convertible_v<std::remove_reference_t<Derived>, std::remove_reference_t<Base>>
const Derived& down_cast(const Base& t)
{
  // Runtime check that the cast is valid.
  assert(Derived(static_cast<const aterm&>(t)) != aterm());

  // UB: Only allowed when we constructed an actual Derived type
  return reinterpret_cast<const Derived&>(reinterpret_cast<const detail::_aterm&>(t));
}

/// \brief A cast form an aterm derived class to a class that inherits in (possibly multiple steps) from this class.
/// \details The derived class is not allowed to contain extra fields. This conversion does not require runtime
/// computation
///          effort. Also see down_cast.
/// \param t The term that is converted.
/// \return A term of type Derived.
template <IsATerm Derived, IsATerm Base>
  requires std::is_base_of_v<std::remove_reference_t<Base>, std::remove_reference_t<Derived>>
const Derived& vertical_cast(const Base& t)
{
  // Runtime check that the cast is valid.
  assert(Derived(static_cast<const aterm&>(t)) != aterm());
  
  return reinterpret_cast<const Derived&>(reinterpret_cast<const detail::_aterm&>(t));
}

/// \brief Send the term in textual form to the ostream.
/// \param out The stream to which the term is sent.
/// \param t   The term that is printed to the stream.
/// \return The stream to which the term is written.
std::ostream& operator<<(std::ostream& out, const atermpp::aterm& t);

/// \brief Transform an aterm to an ascii string.
/// \param t The input aterm.
/// \return A string representation of the given term derived from an aterm.
inline std::string pp(const atermpp::aterm& t)
{
  std::ostringstream oss;
  oss << t;
  return oss.str();
}


namespace detail
{
  
/// This function can be used to translate an term_appl_iterator of one sort into another.
template <IsATerm Derived, IsATerm Base>
term_appl_iterator<Derived> aterm_appl_iterator_cast(term_appl_iterator<Base> a)
  requires std::is_base_of_v<aterm, Base> && std::is_base_of_v<aterm, Derived>
{
  return term_appl_iterator<Derived>(reinterpret_cast<const Derived*>(a.m_term));
}

}

} // namespace atermpp

namespace std
{

/// \brief Swaps two aterms.
/// \details This operation is more efficient than exchanging terms by an assignment,
///          as swapping does not require to change the protection of terms.
///          In order to be used in the standard containers, the declaration must
///          be preceded by an empty template declaration. This swap function is
///          not used for classes that derive from the aterm class. A specific
///          swap function must be provided for derived classes.
/// \param t1 The first term
/// \param t2 The second term
template <>
inline void swap(atermpp::unprotected_aterm_core& t1, atermpp::unprotected_aterm_core& t2) noexcept
{
  t1.swap(t2);
}
} // namespace std
namespace std
{

/// \brief Swaps two term_applss.
/// \details This operation is more efficient than exchanging terms by an assignment,
///          as swapping does not require to change the protection of terms.
/// \param t1 The first term.
/// \param t2 The second term.
inline void swap(atermpp::aterm& t1, atermpp::aterm& t2) noexcept
{
  t1.swap(t2);
}

/// \brief Standard hash function.
template <>
struct hash<atermpp::aterm>
{
  std::size_t operator()(const atermpp::aterm& t) const { return std::hash<atermpp::aterm_core>()(t); }
};

} // namespace std

#endif // MCRL2_ATERMPP_ATERM_APPL_H
