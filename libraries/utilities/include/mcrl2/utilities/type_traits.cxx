// Author(s): Wieger Wesselink, Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_TYPE_TRAITS_H
#define MCRL2_UTILITIES_TYPE_TRAITS_H
MCRL2_MODULE;

#include <functional>
#include <type_traits>

#ifdef MCRL2_ENABLE_MODULES
  export module utilities:type_traits;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities
{

/// \brief A typetrait that is std::true_type iff std::begin() and std::end() can be called on type T
template<typename T, typename = void>
struct is_iterable : std::false_type
{};

template<typename T>
struct is_iterable<T, std::void_t<decltype(std::begin(std::declval<T>())), decltype(std::end(std::declval<T>()))>>
    : std::true_type
{};

/// \brief A concept for whether type T has iterator traits defined.
template <typename T>
concept HasIteratorTraits = requires {
    typename std::iterator_traits<T>::value_type;
};

/// \brief A typetrait that is std::true_type iff the given type has the iterator traits.
template<typename T>
struct is_iterator : std::bool_constant<HasIteratorTraits<T>>
{};

/// \brief Checks whether condition holds for all types passed as variadic template.
template<typename... Conds>
struct forall : std::true_type
{};

template<typename Cond, typename... Conds>
struct forall<Cond, Conds...> : std::conditional_t<Cond::value, forall<Conds...>, std::false_type>
{};

// Helpers
template<typename T>
constexpr bool is_iterable_v = is_iterable<T>::value;

template<typename T>
constexpr bool is_iterator_v = is_iterator<T>::value;

template<typename FunctionType, typename ResultType, typename = void>
struct is_constant_function_yielding : public std::false_type
{};

template<typename FunctionType, typename ResultType>
struct is_constant_function_yielding<FunctionType,
    ResultType,
    typename std::is_convertible<typename std::invoke_result_t<FunctionType>::type, ResultType>> : public std::true_type
{};

template<typename FunctionType, typename ArgumentType, typename = void>
struct is_applicable : public std::false_type
{};

template<typename FunctionType, typename ArgumentType>
struct is_applicable<FunctionType, ArgumentType, typename std::invoke_result_t<FunctionType, ArgumentType>>
    : public std::true_type
{};

template<typename FunctionType, typename ArgumentType1, typename ArgumentType2, typename = void>
struct is_applicable2 : public std::false_type
{};

template<typename FunctionType, typename ArgumentType1, typename ArgumentType2>
struct is_applicable2<FunctionType,
    ArgumentType1,
    ArgumentType2,
    typename std::invoke_result_t<FunctionType, ArgumentType1, ArgumentType2>> : public std::true_type
{};

template<typename T>
struct function_traits : public function_traits<decltype(&T::operator())>
{};

// For generic types, directly use the result of the signature of its 'operator()'

template<typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType (ClassType::*)(Args...) const>
// we specialize for pointers to member function
{
  // arity is the number of arguments.
  enum
  {
    arity = sizeof...(Args)
  };

  using result_type = ReturnType;

  template<size_t i>
  struct arg
  {
    using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
    // the i-th argument is equivalent to the i-th tuple element of a tuple
    // composed of those arguments.
  };
};

} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_TYPE_TRAITS_H
