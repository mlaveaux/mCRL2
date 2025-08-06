// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_DETAIL_ATERM_APPL_H
#define MCRL2_ATERMPP_DETAIL_ATERM_APPL_H
MCRL2_MODULE;

#include <cassert>
#include <memory>
#include <array>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm_core.cxx"
  #include "mcrl2/atermpp/concepts.cxx"
#else
  export module atermpp:detail.aterm_data;

  import :aterm_core;
  import :concepts;
  import :function_symbol;
  import :detail.aterm_core_data;

  import utilities;
#endif

MCRL2_MODULE_EXPORT namespace atermpp
{

// Forward declarations
template <typename Term>
class term_appl_iterator;

class aterm;

namespace detail
{

/// \brief This class stores a term followed by N arguments. Where N should be equal to
///        the arity of the function symbol. These arguments do have room reserved for them
///        during the creation of the _aterm_appl.
template<std::size_t N = 1>
class _aterm_appl : public _aterm
{
public:

  /// \brief Constructs a term application with the given symbol and arguments.
  template<typename... Terms>
  _aterm_appl(const function_symbol& sym, const Terms&... arguments)
    requires (sizeof...(Terms) == N && are_terms<Terms...>::value)
      : _aterm(sym),
        m_arguments{{arguments...}}
  {
    assert(N == sym.arity()); // The arity of the function symbol matches.
  }

  /// \brief constructs a term application with the given symbol and an iterator where the number
  ///        of elements is equal to the template parameter N.
  template<typename Iterator>
  _aterm_appl(const function_symbol& sym, Iterator it, [[maybe_unused]] Iterator end)
    requires (mcrl2::utilities::is_iterator<Iterator>::value)
      : _aterm(sym)
  {
    for (std::size_t i = 0; i < N; ++i)
    {
      assert(it != end);
      m_arguments[i] = *it;
      ++it;
    }

    assert(it == end);
  }

  /// \brief Constructs a term application with the given symbol and arguments.
  _aterm_appl(const function_symbol& sym, std::array<unprotected_aterm_core, N> arguments)
    : _aterm(sym),
      m_arguments(arguments)
  {
    assert(sym.arity() == N);
  }

  /// \brief constructs a term application with the given symbol and its arguments from the iterator.
  template<typename Iterator>
  _aterm_appl(const function_symbol& symbol, Iterator it, [[maybe_unused]] Iterator end, bool)
    requires (mcrl2::utilities::is_iterator<Iterator>::value)
      : _aterm(symbol)
  {
    for (std::size_t i = 0; i < symbol.arity(); ++i)
    {
      // Prevent bound checking, the allocator must make sure that symbol.arity() arguments fit.
      assert(it != end);
      m_arguments.data()[i] = *it;
      ++it;
    }

    assert(it == end);
  }

  /// \returns A const reference to the arguments at the ith position.
  const aterm_core& arg(std::size_t index) const
  {
    return static_cast<const aterm_core&>(m_arguments.data()[index]);
  }

  /// \brief Convert any known number of arguments aterm<N> to the default _aterm_appl.
  explicit operator _aterm_appl<1>& ()
  {
    return reinterpret_cast<_aterm_appl<1>&>(*this);
  }

private:
  std::array<unprotected_aterm_core, N> m_arguments; /// \brief Array of arguments.
};

/// A default instantiation for the underlying term application.
using _term_appl = _aterm_appl<>;

/// \brief This class allocates _aterm_appl objects where the size is based on the arity of
///        the function symbol.
/// \details The template T is required to be an object that implicitly converts to an _aterm_appl.
template<typename T = _term_appl>
class _aterm_appl_allocator
{
private:
  /// \returns The size (in bytes) of a class T with arity number of arguments placed at the end.
  constexpr static std::size_t term_appl_size(std::size_t arity)
  {
    return sizeof(T) + (arity - 1) * sizeof(aterm_core);
  }

public:
  using size_type = std::size_t;
  using pointer = T*;
  using value_type = T;

  template <class U>
  struct rebind
  {
    using other = _aterm_appl_allocator<U>;
  };

  /// \brief Allocates space for an _aterm_appl where the arity is given by the function symbol.
  template<typename ForwardIterator>
  T* allocate_args(const function_symbol& symbol, ForwardIterator, ForwardIterator)
  {
    // We assume that object T contains the _aterm_appl<aterm, 1> at the end and reserve extra space for parameters.
    char* newTerm = m_packed_allocator.allocate(term_appl_size(symbol.arity()));
    return reinterpret_cast<T*>(newTerm);
  }

  /// \brief Allocates space for an _aterm_appl where the arity is given by the function symbol.
  /// \details Assumes that arguments contains symbol.arity() number of terms.
  T* allocate_args(const function_symbol& symbol, unprotected_aterm_core*)
  {
    // We assume that object T contains the _aterm_appl<aterm, 1> at the end and reserve extra space for parameters.
    char* newTerm = m_packed_allocator.allocate(term_appl_size(symbol.arity()));
    return reinterpret_cast<T*>(newTerm);
  }

  /// \brief Constructs an _aterm_appl with arguments taken from begin, the arity is given by the function symbol.
  template<typename ForwardIterator>
  void construct(T* element, const function_symbol& symbol, ForwardIterator begin, ForwardIterator end)
  {
    new (element) T(symbol, begin, end, true);
  }

  /// \brief Specialize destroy for _aterm_appl to only destroy the function symbol. The reference count for the aterm does not have to be decreased.
  void destroy(T* element)
  {
    assert(element != nullptr);

    // Only destroy the function symbol.
    _term_appl& term = static_cast<_term_appl&>(*element);
    term.function().~function_symbol();
  }

  void deallocate(T* element, std::size_t)
  {
    assert(element != nullptr);

    // Deallocate the memory of this aterm.
    _term_appl& term = static_cast<_term_appl&>(*element);
    m_packed_allocator.deallocate(reinterpret_cast<char*>(element), term_appl_size(term.function().arity()));
  }

  // These member functions are to ensure parity with the memory_pool.
  constexpr std::size_t capacity() const { return 0; }
  constexpr std::size_t consolidate() const noexcept { return 0; }
  constexpr bool has_free_slots() const noexcept { return false; }

private:
  std::allocator<char> m_packed_allocator;
};

static_assert(sizeof(_term_appl) == sizeof(_aterm) + sizeof(aterm_core), "Sanity check: aterm size");

template <IsATerm Derived, IsATerm Base>
term_appl_iterator<Derived> aterm_appl_iterator_cast(term_appl_iterator<Base> a)
  requires std::is_base_of_v<aterm, Base> && std::is_base_of_v<aterm, Derived>;

} // namespace detail
} // namespace atermpp

#endif // MCRL2_ATERMPP_DETAIL_ATERM_APPL_H
