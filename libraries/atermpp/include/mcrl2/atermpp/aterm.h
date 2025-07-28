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

#include "mcrl2/atermpp/concepts.h"
#include "mcrl2/atermpp/function_symbol.h"
#include "mcrl2/atermpp/unprotected_aterm.h"
#include "mcrl2/utilities/type_traits.h"

#include <cassert>
#include <sstream>
#include <type_traits>

#include <mcrl3_ffi.h>

namespace atermpp
{
class aterm : public unprotected_aterm
{
public:
  /// An unsigned integral type.
  using size_type = std::size_t;

  /// A signed integral type.
  using difference_type = ptrdiff_t;

  /// \brief Default constructor.
  aterm()
      : unprotected_aterm()
  {}

  /// \brief Copy constructor.
  /// \param other The term that is copied.
  /// \details  This class has a non-trivial destructor so explicitly define the copy and move operators.
  aterm(const aterm& other) noexcept { m_term = other.m_term; m_root = mcrl3::ffi::term_protect(other.m_term); }

  /// Construct a term from the internal representation of an integer term.
  explicit aterm(mcrl3::ffi::aterm_t term)
    : unprotected_aterm(term.term),
      m_root(term.root)
  {}

  /// Construct an aterm from an unprotected term, does NOT change the protection of the term.
  explicit aterm(mcrl3::ffi::unprotected_aterm_t term)
    : unprotected_aterm(term)
  {
    assert(defined());
  }

  /// Desstructor of the term.
  ~aterm() { if (defined()) { mcrl3::ffi::term_unprotect(m_root); } }

  /// \brief Assignment operator.
  /// \param other The term that is assigned.
  /// \details This class has a non-trivial destructor so explicitly define the copy and move assignments.
  aterm& operator=(const aterm& other) noexcept
  {
    if (this != &other)
    {
      if (defined())
      {
        mcrl3::ffi::term_unprotect(m_root);
      }
    }

    m_term = other.m_term;
    m_root = mcrl3::ffi::term_protect(other.m_term);
    return *this;
  }

  aterm(aterm&& other) noexcept = default;
  aterm& operator=(aterm&& other) noexcept = default;

  /// \brief Assignment operator, to be used when the busy flags do not need to be set.
  /// \details This is only safe in the parallel context when the busy flag is already
  ///          known to be set. This is also checked by an assert. This can be used for
  ///          instance in a lambda function that is passed in a make_.... function, as
  ///          this unprotected assign will only be called when a term is constructed.
  /// \param other The aterm_core that will be assigned.
  template <bool CHECK_BUSY_FLAG = true>
  aterm& unprotected_assign(const aterm& other) noexcept
  {
    if constexpr (CHECK_BUSY_FLAG)
    {
      assert(mcrl3::ffi::term_pool_is_busy_set() && "The busy flag must be set before unprotected_assign is called.");
    }

    m_term = other.m_term;
    m_root = mcrl3::ffi::term_protect(other.m_term);
    return *this;
  }

  /// \brief Constructor that provides an aterm based on a function symbol and forward iterator providing the arguments.
  /// \details The iterator range is traversed more than once. If only one traversal is required
  ///          use term_appl with a TermConverter argument. But this function
  ///          is substantially less efficient.
  ///          The length of the iterator range must match the arity of the function symbol.
  /// \param sym A function symbol.
  /// \param begin The start of a range of elements.
  /// \param end The end of a range of elements.
  template <class ForwardIterator,
      typename std::enable_if<mcrl2::utilities::is_iterator<ForwardIterator>::value>::type* = nullptr,
      typename std::enable_if<
          !std::is_same<typename ForwardIterator::iterator_category, std::input_iterator_tag>::value>::type* = nullptr,
      typename std::enable_if<
          !std::is_same<typename ForwardIterator::iterator_category, std::output_iterator_tag>::value>::type* = nullptr>
  aterm(const function_symbol& sym, ForwardIterator begin, ForwardIterator end)
  {
    static_assert(!std::is_same<typename ForwardIterator::iterator_category, std::input_iterator_tag>::value,
        "A forward iterator has more requirements than an input iterator.");
    static_assert(!std::is_same<typename ForwardIterator::iterator_category, std::output_iterator_tag>::value,
        "A forward iterator has more requirements than an output iterator.");
  }

  /// \brief Constructor that provides an aterm based on a function symbol and an input iterator providing the
  /// arguments. \details The given iterator is traversed only once. So it can be used with an input iterator.
  ///          This means that the TermConverter is applied exactly once to each element.
  ///          The length of the iterator range must be equal to the arity of the function symbol.
  /// \param sym A function symbol.
  /// \param begin The start of a range of elements.
  /// \param end The end of a range of elements.
  template <class InputIterator,
      typename std::enable_if<mcrl2::utilities::is_iterator<InputIterator>::value>::type* = nullptr,
      typename std::enable_if<
          std::is_same<typename InputIterator::iterator_category, std::input_iterator_tag>::value>::type* = nullptr>
  aterm(const function_symbol& sym, InputIterator begin, InputIterator end)
      : aterm(sym, begin, end, [](const unprotected_aterm& term) -> const unprotected_aterm& { return term; })
  {
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
  aterm(const function_symbol& sym, InputIterator begin, InputIterator end, TermConverter converter)
  {
    static_assert(!std::is_same<typename InputIterator::iterator_category, std::output_iterator_tag>::value,
        "The InputIterator has the output iterator tag.");
  }

  /// \brief Constructor.
  /// \param sym A function symbol.
  aterm(const function_symbol& sym) {}

  /// \brief Constructor for n-arity function application.
  /// \param symbol A function symbol.
  /// \param arguments The arguments of the function application.
  template <IsATerm... Terms>
  aterm(const function_symbol& symbol, const Terms&... arguments)
    : aterm(mcrl3::ffi::term_create_appl(symbol.get(), convert_to_array(arguments...).data(), sizeof...(Terms)))
  {}

  /// \brief Returns the function symbol belonging to an aterm.
  /// \return The function symbol of this term.
  [[nodiscard]] const function_symbol function() const
  {
    return reinterpret_cast<const function_symbol&>(mcrl3::ffi::term_get_function_symbol(m_term).ptr);
  }

  /// \brief Returns the i-th argument.
  /// \param i A positive integer.
  /// \return The argument with the given index.
  const aterm& operator[](const std::size_t i) const
  {
    assert(i < size()); // Check the bounds.
    // return reinterpret_cast<const aterm&>(mcrl3::ffi::term_get_argument(m_term, i).ptr);
  }

private:
  template<typename... Args>
  static std::array<mcrl3::ffi::unprotected_aterm_t, sizeof...(Args)> convert_to_array(const Args&... args)
  {
    return {static_cast<mcrl3::ffi::unprotected_aterm_t>(args)...};
  }

  /// \brief Returns true if the term has no arguments.
  /// \return True if this term has no arguments.
  [[nodiscard]] bool empty() const { return size() == 0; }
private:
  mcrl3::ffi::root_index_t m_root;
};

using term_callback = void (*)(const aterm&);

inline
void add_deletion_hook(const function_symbol& symbol, term_callback hook)
{
  // mcrl3::ffi::function_symbol_t func = symbol.get();
  // mcrl3::ffi::register_deletion_hook(&func, remap(hook);
}

/// \brief Constructor an aterm in a variable based on a function symbol and an forward iterator providing the
/// arguments. \details The iterator range is traversed more than once. If only one traversal is required
///          use term_appl with a TermConverter argument. But this function
///          is substantially less efficient.
///          The length of the iterator range must match the arity of the function symbol.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param sym A function symbol.
/// \param begin The start of a range of elements.
/// \param end The end of a range of elements.
template <IsATerm Term,
    class ForwardIterator,
    typename std::enable_if<mcrl2::utilities::is_iterator<ForwardIterator>::value>::type* = nullptr,
    typename std::enable_if<
        !std::is_same<typename ForwardIterator::iterator_category, std::input_iterator_tag>::value>::type* = nullptr,
    typename std::enable_if<
        !std::is_same<typename ForwardIterator::iterator_category, std::output_iterator_tag>::value>::type* = nullptr>
void make_term_appl(Term& target, const function_symbol& sym, ForwardIterator begin, ForwardIterator end)
{
  static_assert(!std::is_same<typename ForwardIterator::iterator_category, std::input_iterator_tag>::value,
      "A forward iterator has more requirements than an input iterator.");
  static_assert(!std::is_same<typename ForwardIterator::iterator_category, std::output_iterator_tag>::value,
      "A forward iterator has more requirements than an output iterator.");
}

/// \brief Constructor an aterm in a variable based on a function symbol and an input iterator providing the arguments.
/// \details The given iterator is traversed only once. So it can be used with an input iterator.
///          This means that the TermConverter is applied exactly once to each element.
///          The length of the iterator range must be equal to the arity of the function symbol.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param sym A function symbol.
/// \param begin The start of a range of elements.
/// \param end The end of a range of elements.
template <IsATerm Term,
    class InputIterator,
    typename std::enable_if<mcrl2::utilities::is_iterator<InputIterator>::value>::type* = nullptr,
    typename std::enable_if<
        std::is_same<typename InputIterator::iterator_category, std::input_iterator_tag>::value>::type* = nullptr>
void make_term_appl(Term& target, const function_symbol& sym, InputIterator begin, InputIterator end)
{
  make_term_appl(target, sym, begin, end, [](const Term& term) -> const Term& { return term; });

  static_assert(std::is_same<typename InputIterator::iterator_category, std::input_iterator_tag>::value,
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
template <IsATerm Term,
    class InputIterator,
    class TermConverter,
    typename std::enable_if<mcrl2::utilities::is_iterator<InputIterator>::value>::type* = nullptr>
void make_term_appl(Term& target,
    const function_symbol& sym,
    InputIterator begin,
    InputIterator end,
    TermConverter converter)
{
  static_assert(!std::is_same<typename InputIterator::iterator_category, std::output_iterator_tag>::value,
      "The InputIterator has the output iterator tag.");
}

/// \brief Make an term_appl consisting of a single function symbol.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param sym A function symbol.
template <IsATerm Term>
void make_term_appl(Term& target, const function_symbol& sym)
{

}

/// \brief Make an aterm application for n-arity function application.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param symbol A function symbol.
/// \param arguments The arguments of the function application.
template <class Term, typename... Terms>
void make_term_appl(Term& target, const function_symbol& symbol, const Terms&... arguments)
{}

/// \brief Constructor for n-arity function application with an index.
/// \param target The variable in which the result will be put. This variable may be used for scratch purposes.
/// \param symbol A function symbol.
/// \param arguments The arguments of the function application.
template <class Term, class INDEX_TYPE, typename... Terms>
void make_term_appl_with_index(aterm& target, const function_symbol& symbol, const Terms&... arguments)
{}

/// \brief A universal cast from an aterm to another aterm that is convertible in either direction. Less strict than
/// vertical_cast. \param  t A term of a type inheriting from an aterm. \return  A term of type const Derived&.
template <IsATerm Derived, IsATerm Base>
  requires std::is_convertible_v<std::remove_reference_t<Base>, std::remove_reference_t<Derived>>
           || std::is_convertible_v<std::remove_reference_t<Derived>, std::remove_reference_t<Base>>
const Derived& down_cast(const Base& t)
{
  // Runtime check that the cast is valid.
  assert(Derived(static_cast<const aterm&>(t)) != aterm());

  // UB: Only allowed when we constructed an actual Derived type
  return reinterpret_cast<const Derived&>(t);
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

  return reinterpret_cast<const Derived&>(t);
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

} // namespace atermpp

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
  std::size_t operator()(const atermpp::aterm& t) const { return std::hash<atermpp::unprotected_aterm>()(t); }
};

} // namespace std

#endif // MCRL2_ATERMPP_ATERM_APPL_H
