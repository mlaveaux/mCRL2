// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_DETAIL_ATERM_HASH_H_
#define MCRL2_ATERMPP_DETAIL_ATERM_HASH_H_
MCRL2_MODULE;

#include <cstdint>
#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/aterm_data.cxx"
  #include "mcrl2/atermpp/detail/aterm_int_data.cxx"
  #include "mcrl2/atermpp/detail/function_symbol_hash.cxx"
  #include "mcrl2/utilities/hash_utility.cxx"
  #include "mcrl2/utilities/detail/memory_utility.cxx"
#else
  export module atermpp:detail.aterm_hash;

  import :aterm_core;
  import :function_symbol;
  import :detail.aterm_core_data;
  import :detail.function_symbol_hash;
  import :detail.aterm_data;
  import :detail.aterm_int_data;

  import utilities;
#endif

MCRL2_MODULE_EXPORT namespace std
{

/// \brief specialization of the standard std::hash function.
template<>
struct hash<atermpp::detail::_aterm*>
{
  std::size_t operator()(const atermpp::detail::_aterm* term) const
  {
    // The hash function given here shifts a term four positions to the right.
    // This is very effective in the toolset, as the aterms are often stored 
    // in consecutive positions and this means they are stored in consecutive positions
    // in hash tables based on chaining. However, this has function is very disadvantageous
    // for other hash tables, especially those with open addressing. In that case
    // it is much better to use the hash function below, which unfortunately can lead
    // to a performance drop of 1/3, compared to the ">> 4" hash below. 
    // The reason for the performance drop is most likely that with better hashing
    // the items are better dispersed in hash tables, leading to worse cache behaviour. 
    //   std::hash<const void*> hasher;
    //   return hasher(reinterpret_cast<const void*>(term));

    return reinterpret_cast<std::uintptr_t>(term) >> 4;
    
  }
};

/// \brief specialization of the standard std::hash function.
template<>
struct hash<atermpp::unprotected_aterm_core>
{
  std::size_t operator()(const atermpp::unprotected_aterm_core& term) const
  {
    const std::hash<atermpp::detail::_aterm*> hasher;
    return hasher(atermpp::detail::address(term));
  }
};

/// \brief specialization of the standard std::hash function.
template<>
struct hash<atermpp::aterm_core>
{
  std::size_t operator()(const atermpp::aterm_core& t) const
  {
    const std::hash<atermpp::detail::_aterm*> aterm_hasher;
    return aterm_hasher(atermpp::detail::address(t));
  }
};

} // namespace std


namespace atermpp::detail
{

/// \brief Indicates that the number of arguments is not known at compile time.
constexpr std::size_t DynamicNumberOfArguments = std::numeric_limits<std::size_t>::max();

/// \brief Computes the hash of the given term.
/// \details Can be optimized with loop unrolling when template parameter N is provided.
///          However, this assumes that every passed to term has arity equal to N.
template<std::size_t N = DynamicNumberOfArguments>
struct aterm_hasher
{
  using is_transparent = void;

  std::size_t operator()(const _aterm& term) const noexcept;
  std::size_t operator()(const function_symbol& symbol) const noexcept;
  std::size_t operator()(const function_symbol& symbol, unprotected_aterm_core arguments[]) const noexcept;

  template <typename ForwardIterator>
  std::size_t operator()(const function_symbol& symbol, ForwardIterator begin, ForwardIterator end) const noexcept 
    requires mcrl2::utilities::is_iterator<ForwardIterator>::value;
};

/// \brief Computes the hash of the given term.
/// \details This version only works whenever N is a compile time constant.
template<std::size_t N = 0>
struct aterm_hasher_finite : public aterm_hasher<N>
{
  using is_transparent = void;

  using aterm_hasher<N>::operator();
  std::size_t operator()(const function_symbol& symbol, std::array<unprotected_aterm_core, N> key) const noexcept;

  template<typename ...Args>
  std::size_t operator()(const function_symbol& symbol, const Args&... args) const noexcept;
};

/// \brief Computes the hash of the integral term arguments.
struct aterm_int_hasher
{
  using is_transparent = void;

  inline std::size_t operator()(const _aterm_int& term) const noexcept;
  inline std::size_t operator()(std::size_t value) const noexcept;
};

/// \brief Returns true iff first and second are value-equivalent.
/// \details Can be optimized with loop unrolling when template parameter N is provided.
///          However, this assumes that every passed to term has arity equal to N.
template<std::size_t N = DynamicNumberOfArguments>
struct aterm_equals
{
  using is_transparent = void;

  bool operator()(const _aterm& first, const _aterm& second) const noexcept;
  bool operator()(const _aterm& term, const function_symbol& symbol) const noexcept;
  bool operator()(const _aterm& term, const function_symbol& symbol, unprotected_aterm_core arguments[]) const noexcept;

  template <typename ForwardIterator>
  bool operator()(const _aterm& term,
      const function_symbol& symbol,
      ForwardIterator begin,
      ForwardIterator end) const noexcept
    requires mcrl2::utilities::is_iterator<ForwardIterator>::value;
};

template<std::size_t N>
struct aterm_equals_finite : public aterm_equals<N>
{
  using aterm_equals<N>::operator();
  bool operator()(const _aterm& term, const function_symbol& symbol, std::array<unprotected_aterm_core, N> key) const noexcept;

  template<typename ...Args>
  bool operator()(const _aterm& term, const function_symbol& symbol, const Args&... args) const noexcept;
};

/// \brief Returns true iff the given term(s) or value are equivalent.
struct aterm_int_equals
{
  using is_transparent = void;

  inline bool operator()(const _aterm_int& first, const _aterm_int& second) const noexcept;
  inline bool operator()(const _aterm_int& term, std::size_t value) const noexcept;
};

/// \brief Auxiliary function to combine hnr with aterms.
inline
std::size_t combine(const std::size_t hnr, const unprotected_aterm_core& term)
{
  const std::hash<unprotected_aterm_core> hasher;
  return mcrl2::utilities::detail::hash_combine_cheap(hnr, hasher(term));
}

/// Implementation

template<std::size_t N>
std::size_t aterm_hasher<N>::operator()(const _aterm& term) const noexcept
{
  const function_symbol& f = term.function();
  std::size_t hnr = operator()(f);

  // The arity is defined by the function symbol iff N is unchanged and the arity is N otherwise.
  const std::size_t arity = (N == DynamicNumberOfArguments) ? f.arity() : N;

  // This is a function application with arguments, hash each argument and combine the result.
  const _aterm_appl<>& term_appl = static_cast<const _aterm_appl<>&>(term);
  for (std::size_t i = 0; i < arity; ++i)
  {
    hnr = combine(hnr, term_appl.arg(i));
  }

  return hnr;
}
template<std::size_t N>
std::size_t aterm_hasher<N>::operator()(const function_symbol& symbol) const noexcept
{
  const std::hash<function_symbol> function_hasher;
  return function_hasher(symbol);
}

template<std::size_t N>
std::size_t aterm_hasher<N>::operator()(const function_symbol& symbol, unprotected_aterm_core arguments[]) const noexcept
{
  // The arity is defined by the function symbol iff N is unchanged and the arity is N otherwise.
  const std::size_t arity = (N == DynamicNumberOfArguments) ? symbol.arity() : N;

  // This is a function application with arguments, hash each argument and combine the result.
  std::size_t hnr = operator()(symbol);
  for (std::size_t i = 0; i < arity; ++i)
  {
    hnr = combine(hnr, arguments[i]);
  }

  return hnr;
}

template <std::size_t N>
template <typename ForwardIterator>
inline std::size_t aterm_hasher<N>::operator()(const function_symbol& symbol,
    ForwardIterator it,
    [[maybe_unused]] ForwardIterator end) const noexcept
  requires mcrl2::utilities::is_iterator<ForwardIterator>::value
{
  // The arity is defined by the function symbol iff N is unchanged and the arity is N otherwise.
  const std::size_t arity = (N == DynamicNumberOfArguments) ? symbol.arity() : N;

  // This is a function application with arguments, hash each argument and combine the result.
  std::size_t hnr = operator()(symbol);
  for (std::size_t i = 0; i < arity; ++i)
  {
    assert(it != end);
    hnr = combine(hnr, *it);
    ++it;
  }

  assert(it == end);
  return hnr;
}

template<std::size_t N>
std::size_t aterm_hasher_finite<N>::operator()(const function_symbol& symbol, std::array<unprotected_aterm_core, N> arguments) const noexcept
{
  std::size_t hnr = operator()(symbol);

  // This is a function application with arguments, hash each argument and combine the result.
  for (std::size_t i = 0; i < N; ++i)
  {
    hnr = combine(hnr, arguments[i]);
  }

  return hnr;
}

template <std::size_t I = 0, typename... Tp>
std::size_t combine_args(std::size_t seed, const Tp&...)
  requires(I == sizeof...(Tp))
{
  return seed;
}

template <std::size_t I = 0, typename... Tp>
std::size_t combine_args(std::size_t seed, const Tp&... t)
  requires(I < sizeof...(Tp))
{
  return combine_args<I+1>(combine(seed, std::get<I>(std::forward_as_tuple(t...))), t...);
}

template<std::size_t N>
template<typename ...Args>
std::size_t aterm_hasher_finite<N>::operator()(const function_symbol& symbol, const Args&... args) const noexcept
{
  std::size_t hnr = operator()(symbol);

  // This is a function application with arguments, hash each argument and combine the result.
  return combine_args(hnr, args...);
}

std::size_t aterm_int_hasher::operator()(const _aterm_int& term) const noexcept
{
  return aterm_int_hasher()(term.value());
}

// The size_t hashfunction below has been taken from a note on stackoverflow
// by Wolfgang Brehm. 
static inline std::size_t xorshift(const std::size_t n, const std::size_t i)
{
  return n^(n>>i);
}

std::size_t aterm_int_hasher::operator()(std::size_t value) const noexcept
{
  const std::size_t p = 0x5555555555555555ull; // pattern of alternating 0 and 1
  const std::size_t c = 17316035218449499591ull;// random odd integer constant; 
  return c*xorshift(p*xorshift(value,32),32);
}

template<std::size_t N>
bool aterm_equals<N>::operator()(const _aterm& first, const _aterm& second) const noexcept
{
  if (&first == &second)
  {
    // If the pointers are equal they match by definition
    return true;
  }

  // The arity is defined by the function symbol iff N is unchanged and the arity is N otherwise.
  const std::size_t arity = (N == DynamicNumberOfArguments) ? first.function().arity() : N;

  // Check whether the remaining arguments match
  for (std::size_t i = 0; i < arity; ++i)
  {
    if (static_cast<const _term_appl&>(first).arg(i)
          != static_cast<const _term_appl&>(second).arg(i))
    {
      return false;
    }
  }

  return first.function() == second.function();
}

template<std::size_t N>
bool aterm_equals<N>::operator()(const _aterm& term, const function_symbol& symbol) const noexcept
{
  return term.function() == symbol;
}

template<std::size_t N>
bool aterm_equals<N>::operator()(const _aterm& term, const function_symbol& symbol, unprotected_aterm_core arguments[]) const noexcept
{
  // Each argument should be equal.
  for (std::size_t i = 0; i < symbol.arity(); ++i)
  {
    if (static_cast<const _term_appl&>(term).arg(i) != arguments[i])
    {
      return false;
    }
  }

  return term.function() == symbol;
}

template <std::size_t N>
template <typename ForwardIterator>
inline bool aterm_equals<N>::operator()(const _aterm& term,
    const function_symbol& symbol,
    ForwardIterator it,
    [[maybe_unused]] ForwardIterator end) const noexcept
  requires mcrl2::utilities::is_iterator<ForwardIterator>::value
{
  const std::size_t arity = (N == DynamicNumberOfArguments) ? symbol.arity() : N;

  // Each argument should be equal.
  for (std::size_t i = 0; i < arity; ++i)
  {
    assert(it != end);
    if (static_cast<const _term_appl&>(term).arg(i) != (*it))
    {
      return false;
    }
    ++it;
  }

  assert(it == end);
  return term.function() == symbol;
}

template<std::size_t N>
bool aterm_equals_finite<N>::operator()(const _aterm& term, const function_symbol& symbol, std::array<unprotected_aterm_core, N> arguments) const noexcept
{
  // Each argument should be equal.
  for (std::size_t i = 0; i < N; ++i)
  {
    if (static_cast<const _aterm_appl<N>&>(term).arg(i) != arguments[i])
    {
      return false;
    }
  }

  return term.function() == symbol;
}

template <std::size_t I = 0, typename... Tp>
bool equal_args(const _aterm_appl<8>&, const Tp&...)
  requires(I == sizeof...(Tp))
{
  return true;
}

template <std::size_t I = 0, typename... Tp>
bool equal_args(const _aterm_appl<8>& term, const Tp&... t)
  requires(I < sizeof...(Tp))
{
  return term.arg(I) == std::get<I>(std::forward_as_tuple(t...)) && equal_args<I+1>(term, t...);
}

template<std::size_t N>
template<typename ...Args>
bool aterm_equals_finite<N>::operator()(const _aterm& term, const function_symbol& symbol, const Args&... args) const noexcept
{
  return term.function() == symbol && equal_args(static_cast<const _aterm_appl<8>&>(term), args...);
}

bool aterm_int_equals::operator()(const _aterm_int& first, const _aterm_int& second) const noexcept
{
  return (first.value() == second.value());
}

bool aterm_int_equals::operator()(const _aterm_int& term, std::size_t value) const noexcept
{
  return (term.value() == value);
}

} // namespace atermpp::detail


#endif // MCRL2_ATERMPP_DETAIL_ATERM_HASH_H_
