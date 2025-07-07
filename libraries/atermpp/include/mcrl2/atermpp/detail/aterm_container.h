// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_DETAIL_ATERM_CONTAINER_H
#define MCRL2_ATERMPP_DETAIL_ATERM_CONTAINER_H

#include "mcrl2/atermpp/concepts.h"
#include "mcrl2/atermpp/unprotected_aterm.h"
#include "mcrl2/utilities/noncopyable.h"

#include <mcrl3_ffi.h>

#include <functional>
#include <type_traits>



namespace atermpp::detail
{

/// \brief Provides safe storage of unprotected_aterm instances in a container by marking
///        them during garbage collection.
///
/// \details Can not be inherited since it is being registered during construction and the
///          vptr is being updated by inherited classes otherwise.
class aterm_container final : private mcrl2::utilities::noncopyable
{
public:
  aterm_container(std::function<void()> mark_func, std::function<std::size_t()> size_func)
      : m_root(mcrl3::ffi::container_protect()),
        mark_func(std::move(mark_func)),
        size_func(std::move(size_func))
  {}
  
  ~aterm_container()
  {
    // Unprotect the root term of the container.
    mcrl3::ffi::container_unprotect(m_root);
  }

  /// \brief Ensure that all the terms in the containers.
  void mark() const { mark_func(); }

  [[nodiscard]] inline std::size_t size() const { return size_func(); }

private:
  mcrl3::ffi::root_index_t m_root;
  
  std::function<void()> mark_func;
  std::function<std::size_t()> size_func;
};

template <class T, typename Type = void>
class reference_aterm;

template <class T>
struct is_pair_helper : public std::false_type
{};

template <class T, class U>
struct is_pair_helper<std::pair<T, U>> : public std::true_type
{};

template <class T>
struct is_pair : public is_pair_helper<typename std::decay<T>::type>
{};

template <class T>
struct is_reference_aterm_helper : public std::false_type
{};

template <class T>
struct is_reference_aterm_helper<reference_aterm<T>> : public std::true_type
{};

template <class T>
struct is_reference_aterm : public is_reference_aterm_helper<typename std::decay<T>::type>
{};

/// \brief Base class that should not be used.
template <class T, typename Type>
class reference_aterm
{
protected:
  using T_type = typename std::decay<T>::type;

  T_type m_t;
public:
  reference_aterm() = default;

  reference_aterm(const T_type& other) noexcept
      : m_t(other)
  {}

  template <class... Args>
  reference_aterm(Args&&... args) noexcept
      : m_t(std::forward<Args>(args)...)
  {}

  template <class... Args>
  reference_aterm(const Args&... args) noexcept
      : m_t(args...)
  {}

  reference_aterm(T_type&& other) noexcept
      : m_t(std::forward(other))
  {}

  const reference_aterm& operator=(const T_type& other) noexcept
  {
    static_assert(std::is_base_of<unprotected_aterm, T_type>::value);
    m_t = other;
    return m_t;
  }

  const reference_aterm& operator=(T_type&& other) noexcept
  {
    static_assert(std::is_base_of<unprotected_aterm, T_type>::value);
    m_t = std::forward(other);
    return m_t;
  }

  /// Converts implicitly to a protected term of type T.
  operator T_type&() { return m_t; }

  operator const T_type&() const { return m_t; }

  /// For types that are not a std::pair, or a type convertible to an unprotected_aterm
  /// it is necessary that a dedicated mark function is provided that calls mark_term
  /// on all unprotected_aterm types in the class T, when this class is stored in an atermpp container.
  /// See below for an example, where the code is given for pairs, aterms and built in types.
  /// The container is traversed during garbage collection, such that all terms in these
  /// containers are protected individually, without putting them all explicitly in
  /// protection sets.
  void mark() const { m_t.mark(); }
};

/// \brief A reference unprotected_aterm applied to fundamental types, such as int, bool. Nothing needs to happen with such
///       terms. But a special class is needed, because such types are not classes, and we cannot derive from it.
template <class T>
class reference_aterm<T, typename std::enable_if<std::is_fundamental<typename std::decay<T>::type>::value>::type>
{
protected:
  using T_type = typename std::decay<T>::type;
  T_type m_t;

public:

  /// \brief Default constructor.
  reference_aterm() noexcept = default;

  reference_aterm(const T& other) noexcept
      : m_t(other)
  {}

  reference_aterm(T_type&& other) noexcept
      : m_t(std::move(other))
  {}

  const T& operator=(const T& other) noexcept
  {
    m_t = other;
    return m_t;
  }

  const T& operator=(T&& other) noexcept
  {
    m_t = std::move(other);
    return m_t;
  }

  /// Converts implicitly to a protected term of type T.
  operator T&() { return m_t; }

  operator const T&() const { return m_t; }

  void mark(/* todo */) const { /* Do nothing */ }

  bool operator==(const reference_aterm& other) const { return m_t == other.m_t; }
};

/// \brief An unprotected term that is stored inside an aterm_container.
template <IsATerm T>
class reference_aterm<T, typename std::enable_if<std::is_base_of<unprotected_aterm, T>::value>::type>
    : public unprotected_aterm
{
public:
  /// \brief Default constructor.
  reference_aterm() noexcept = default;

  explicit reference_aterm(const unprotected_aterm& other) noexcept { m_term = detail::address(other); }

  reference_aterm(const T& other) noexcept
      : unprotected_aterm(detail::address(other))
  {}

  reference_aterm(unprotected_aterm&& other) noexcept
      : unprotected_aterm(detail::address(other))
  {}

  const reference_aterm& operator=(const unprotected_aterm& other) noexcept;
  const reference_aterm& operator=(unprotected_aterm&& other) noexcept;

  /// Converts implicitly to a protected term of type T.
  operator T&()
  {
    return reinterpret_cast<T&>(*this);
  }

  operator const T&() const
  {
    return reinterpret_cast<const T&>(*this);
  }

  void mark() const
{
    unprotected_aterm::mark();
  }
};

template <typename T>
typename std::pair<typename std::conditional<is_reference_aterm<typename T::first_type>::value,
                       typename T::first_type,
                       reference_aterm<typename T::first_type>>::type,
    typename std::conditional<is_reference_aterm<typename T::second_type>::value,
        typename T::second_type,
        reference_aterm<typename T::second_type>>::type>
reference_aterm_pair_constructor_helper(const T& other)
{
  if constexpr (is_reference_aterm<typename T::first_type>::value && is_reference_aterm<typename T::second_type>::value)
  {
    return other;
  }
  else if constexpr (is_reference_aterm<typename T::first_type>::value
                     && !is_reference_aterm<typename T::second_type>::value)
  {
    return std::pair(other.first, reference_aterm<typename T::second_type>(other.second));
  }
  else if constexpr (!is_reference_aterm<typename T::first_type>::value
                     && is_reference_aterm<typename T::second_type>::value)
  {
    return std::pair(reference_aterm<typename T::first_type>(other.first), other.second);
  }
  else
  {
    static_assert(!is_reference_aterm<typename T::first_type>::value
                      && !is_reference_aterm<typename T::second_type>::value,
        "Logic error");

    return std::pair(reference_aterm<typename T::first_type>(other.first),
        reference_aterm<typename T::second_type>(other.second));
  }
}

/// \brief A pair that is stored into an atermpp container. This class takes care that all aterms that occur
/// (recursively) inside
///        such a pair are marked, whereas non-aterm types are not marked.
template <typename T>
class reference_aterm<T, typename std::enable_if<is_pair<T>::value>::type>
    : public std::pair<typename std::conditional<is_reference_aterm<typename T::first_type>::value,
                           typename T::first_type,
                           reference_aterm<typename T::first_type>>::type,
          typename std::conditional<is_reference_aterm<typename T::second_type>::value,
              typename T::second_type,
              reference_aterm<typename T::second_type>>::type>
{
protected:
  using super = std::pair<typename std::conditional<is_reference_aterm<typename T::first_type>::value,
                              typename T::first_type,
                              reference_aterm<typename T::first_type>>::type,
      typename std::conditional<is_reference_aterm<typename T::second_type>::value,
          typename T::second_type,
          reference_aterm<typename T::second_type>>::type>;
  using std_pair = T;

public:
  /// \brief Default constructor.
  reference_aterm() = default;

  reference_aterm(const reference_aterm& other)
      : super()
  {
    *this = other;
  }

  reference_aterm(const std_pair& other)
      : super(reference_aterm_pair_constructor_helper(other))
  {}

  reference_aterm(std_pair&& other)
      : super(reference_aterm<typename T::first_type>(std::move(other.first)),
            reference_aterm<typename T::second_type>(std::move(other.second)))
  {}

  reference_aterm& operator=(const reference_aterm& other)
  {
    super::first = other.first;
    super::second = other.second;
    return *this;
  }

  const reference_aterm& operator=(const std_pair& other)
  {
    super::first = other.first;
    super::second = other.second;
    return *this;
  }

  reference_aterm& operator=(reference_aterm&& other)
  {
    super::first = std::move(other.first);
    super::second = std::move(other.second);
    return *this;
  }

  const reference_aterm& operator=(std_pair&& other)
  {
    super::first = other.first;
    super::second = other.second;
    return *this;
  }

  /// Converts implicitly to a protected term of type std::pair<T,U>..
  operator std_pair&() { return reinterpret_cast<std_pair>(*this); }

  operator const std_pair&() const { return *reinterpret_cast<std_pair const*>(this); }

  void mark() const
  {
    if constexpr (is_reference_aterm<typename T::first_type>::value)
    {
      super::first.mark();
    }
    else
    {
      reference_aterm<typename T::first_type>(super::first).mark();
    }

    if constexpr (is_reference_aterm<typename T::second_type>::value)
    {
      super::second.mark();
    }
    else
    {
      reference_aterm<typename T::second_type>(super::second).mark();
    }
  }
};

template <typename T, typename Allocator>
class aterm_allocator
{
public:
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  // static_assert(std::is_same_v<value_type, typename Allocator::value_type>, "Types should be equal");

  template <class U>
  struct rebind
  {
    using other = aterm_allocator<U, typename Allocator::template rebind<U>::other>;
  };

  aterm_allocator() = default;

  /// \details The unused parameter is to make the interface equivalent
  ///          to the allocator.
  T* allocate(size_type n, const void* hint = nullptr) { return m_allocator.allocate(n, hint); }

  /// \details The unused parameter is to make the interface equivalent
  ///          to the allocator.
  void deallocate(T* p, size_type n);

  // Move assignment and construction is possible.
  aterm_allocator(aterm_allocator&&) = default;
  aterm_allocator& operator=(aterm_allocator&&) = default;

private:
  Allocator m_allocator;
};

template <typename Container>
class generic_aterm_container
{
public:
  /// \brief Constructor
  generic_aterm_container(const Container& container)
      : m_container(
            [&container]()
            {
              // Marking contained terms.
              for (const typename Container::value_type& element : container)
              {
                static_assert(is_reference_aterm<reference_aterm<typename Container::value_type>>::value);
                if constexpr (is_reference_aterm<typename Container::value_type>::value)
                {
                  static_assert(is_reference_aterm<typename Container::value_type>::value);
                  element.mark();
                }
                else
                {
                  static_assert(!is_reference_aterm<typename Container::value_type>::value);
                  reference_aterm<typename Container::value_type>(element).mark();
                }
              }
            },
            [&container]() -> std::size_t
            {
              // Return the number of elements in the container.
              return container.size();
            })
  {}

  // Container is a reference so unclear what to do in these cases.
  generic_aterm_container(const generic_aterm_container&) = delete;
  generic_aterm_container(generic_aterm_container&&) = delete;

  // It is fine here if the container gets updated, but the functions stay the same.
  generic_aterm_container& operator=(const generic_aterm_container&) { return *this; };

  generic_aterm_container& operator=(generic_aterm_container&) { return *this; }

protected:
  aterm_container m_container;
};

} // namespace atermpp::detail


namespace std
{

/// \brief specialization of the standard std::hash function.
template <class T>
struct hash<atermpp::detail::reference_aterm<T>>
{
  std::size_t operator()(const atermpp::detail::reference_aterm<T>& t) const { return std::hash<T>()(t); }
};

} // namespace std

#endif // MCRL2_ATERMPP_DETAIL_ATERM_CONTAINER_H
