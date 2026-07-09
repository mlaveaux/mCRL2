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

#include <compare>
#include <concepts>
#include <stack>
#include <type_traits>
#include <utility>
#include <vector>
#include "mcrl2/atermpp/aterm_core.h"
#include "mcrl2/atermpp/detail/aterm_core.h"
#include "mcrl2/atermpp/detail/aterm_pool_storage_implementation.h"
#include "mcrl2/utilities/noncopyable.h"

namespace atermpp
{

using term_mark_stack = std::stack<std::reference_wrapper<detail::_aterm>>;

inline void mark_term(const aterm_core& t, term_mark_stack& todo)
{
  if (t.defined())
  {
    detail::mark_term(*atermpp::detail::address(t), todo);
  }
}

namespace detail
{

class thread_aterm_pool;
thread_aterm_pool& g_thread_term_pool();

/// \brief Acquires a shared lock on the global thread aterm pool.
/// \details This free function is declared here but defined in aterm_implementation.cpp, where
///          thread_aterm_pool is a complete type. It allows reference_aterm::operator= below to be
///          defined inline, which is required to work around an MSVC limitation regarding out-of-line
///          definitions of members of constrained partial specializations, without this header needing
///          the complete definition of thread_aterm_pool, which Clang requires when the body is inline.
mcrl2::utilities::shared_guard lock_shared_aterm_pool();

/// \brief Provides safe storage of unprotected_aterm_core instances in a container by marking
///        them during garbage collection.
///
/// \details Can not be inherited since it is being registered during construction and the
///          vptr is being updated by inherited classes otherwise.
class aterm_container final : private mcrl2::utilities::noncopyable
{
public:
  aterm_container(std::function<void(term_mark_stack&)> mark_func, std::function<std::size_t()> size_func);
  ~aterm_container();

  /// \brief Mark all terms held by this container.
  void mark(term_mark_stack& todo) const
  {
    mark_func(todo);
  }

  std::size_t size() const
  {
    return size_func();
  }

private:
  std::function<void(term_mark_stack&)> mark_func;
  std::function<std::size_t()> size_func;
};

// Forward declaration.
template<typename T>
struct markable_aterm;

/// \brief Concept for types that support GC marking.
/// \details A Markable type provides void mark(term_mark_stack&) const that recursively marks
///          all aterm_core instances it contains.
template<typename T>
concept Markable = requires(const T& t, term_mark_stack& todo)
{
  { t.mark(todo) } -> std::same_as<void>;
};

/// \brief Type trait to detect markable_aterm instantiations.
template<typename T>
struct is_markable_aterm : std::false_type
{};

template<typename T>
struct is_markable_aterm<markable_aterm<T>> : std::true_type
{};

/// \brief Returns T if T is already a markable_aterm<U>, otherwise markable_aterm<T>.
/// \details Prevents double-wrapping when T is already a markable_aterm.
template<typename T>
using markable_t = std::conditional_t<is_markable_aterm<std::decay_t<T>>::value, T, markable_aterm<T>>;

// ---- Specialization: fundamental types ----

/// \brief Wrapper for fundamental types (int, bool, etc.) — mark is a no-op.
template<typename T>
  requires std::is_fundamental_v<std::decay_t<T>>
struct markable_aterm<T>
{
private:
  using stored_type = std::decay_t<T>;
  stored_type m_t{};

public:
  markable_aterm() noexcept = default;

  markable_aterm(const T& other) noexcept
    : m_t(other)
  {}

  markable_aterm(stored_type&& other) noexcept
    : m_t(std::move(other))
  {}

  markable_aterm& operator=(const T& other) noexcept
  {
    m_t = other;
    return *this;
  }

  markable_aterm& operator=(T&& other) noexcept
  {
    m_t = std::move(other);
    return *this;
  }

  operator T&() noexcept { return m_t; }
  operator const T&() const noexcept { return m_t; }

  void mark(term_mark_stack& /* todo */) const noexcept {}

  bool operator==(const markable_aterm& other) const { return m_t == other.m_t; }

  auto operator<=>(const markable_aterm& other) const { return m_t <=> other.m_t; }
};

// ---- Specialization: aterm_core-derived types ----

/// \brief Unprotected wrapper for aterm_core-derived types stored inside an aterm container.
/// \details Stores a raw _aterm* without registering it in the GC variable root set.
///          The owning container registers once via aterm_container instead.
template<typename T>
  requires std::is_base_of_v<aterm_core, T>
struct markable_aterm<T> : unprotected_aterm_core
{
  /// \brief Default constructor.
  markable_aterm() noexcept = default;

  explicit markable_aterm(const unprotected_aterm_core& other) noexcept
  {
    m_term = detail::address(other);
  }

  markable_aterm(const T& other) noexcept
    : unprotected_aterm_core(detail::address(other))
  {}

  markable_aterm(unprotected_aterm_core&& other) noexcept
    : unprotected_aterm_core(detail::address(other))
  {}

  markable_aterm& operator=(const unprotected_aterm_core& other) noexcept;
  markable_aterm& operator=(unprotected_aterm_core&& other) noexcept;

  operator T&() noexcept
  {
    static_assert(std::is_base_of_v<aterm_core, T>, "Term must be derived from an aterm_core");
    static_assert(sizeof(T) == sizeof(std::size_t), "Term derived from an aterm_core must not have extra fields");
    return reinterpret_cast<T&>(*this);
  }

  operator const T&() const noexcept
  {
    static_assert(std::is_base_of_v<aterm_core, T>, "Term must be derived from an aterm_core");
    static_assert(sizeof(T) == sizeof(std::size_t), "Term derived from an aterm_core must not have extra fields");
    return reinterpret_cast<const T&>(*this);
  }

  void mark(term_mark_stack& todo) const
  {
    if (defined())
    {
      mark_term(*m_term, todo);
    }
  }
};

// ---- Specialization: std::pair<F, S> ----

/// \brief Wrapper for std::pair — recursively marks both elements.
/// \details Each element is stored as markable_t<F> / markable_t<S> to avoid double-wrapping.
///          Because markable_t<T> always provides mark(), the mark() implementations below
///          call it directly without any if-constexpr dispatch.
template<typename F, typename S>
struct markable_aterm<std::pair<F, S>> : std::pair<markable_t<F>, markable_t<S>>
{
private:
  using std_pair = std::pair<F, S>;
  using super = std::pair<markable_t<F>, markable_t<S>>;

public:
  markable_aterm() = default;
  markable_aterm(const markable_aterm&) = default;
  markable_aterm(markable_aterm&&) = default;
  markable_aterm& operator=(const markable_aterm&) = default;
  markable_aterm& operator=(markable_aterm&&) = default;

  markable_aterm(const std_pair& other)
    : super(markable_t<F>(other.first), markable_t<S>(other.second))
  {}

  markable_aterm(std_pair&& other)
    : super(markable_t<F>(std::move(other.first)), markable_t<S>(std::move(other.second)))
  {}

  markable_aterm& operator=(const std_pair& other)
  {
    super::first = other.first;
    super::second = other.second;
    return *this;
  }

  markable_aterm& operator=(std_pair&& other)
  {
    super::first = std::move(other.first);
    super::second = std::move(other.second);
    return *this;
  }

  operator std_pair&()
  {
    return *reinterpret_cast<std_pair*>(this);
  }

  operator const std_pair&() const
  {
    return *reinterpret_cast<const std_pair*>(this);
  }

  void mark(term_mark_stack& todo) const
  {
    super::first.mark(todo);
    super::second.mark(todo);
  }
};

// ---- Specialization: std::vector<T, Alloc> ----

/// \brief Wrapper for std::vector — recursively marks all elements.
/// \details Elements are stored as markable_t<T> to avoid double-wrapping.
template<typename T, typename Alloc>
struct markable_aterm<std::vector<T, Alloc>>
  : std::vector<markable_t<T>, typename std::allocator_traits<Alloc>::template rebind_alloc<markable_t<T>>>
{
private:
  using stored_type = markable_t<T>;
  using stored_alloc = typename std::allocator_traits<Alloc>::template rebind_alloc<stored_type>;
  using std_vector = std::vector<T, Alloc>;
  using super = std::vector<stored_type, stored_alloc>;

public:
  markable_aterm() = default;

  markable_aterm(const std_vector& other)
  {
    this->reserve(other.size());
    for (const auto& v : other)
      this->emplace_back(v);
  }

  markable_aterm(std_vector&& other)
  {
    this->reserve(other.size());
    for (auto& v : other)
      this->emplace_back(std::move(v));
  }

  markable_aterm& operator=(const std_vector& other)
  {
    this->clear();
    this->reserve(other.size());
    for (const auto& v : other)
      this->emplace_back(v);
    return *this;
  }

  markable_aterm& operator=(std_vector&& other)
  {
    this->clear();
    this->reserve(other.size());
    for (auto& v : other)
      this->emplace_back(std::move(v));
    return *this;
  }

  operator std_vector&()
  {
    return *reinterpret_cast<std_vector*>(this);
  }

  operator const std_vector&() const
  {
    return *reinterpret_cast<const std_vector*>(this);
  }

  void mark(term_mark_stack& todo) const
  {
    for (const auto& element : *this)
      element.mark(todo);
  }
};

// ---- Primary template: user-defined Markable types ----

/// \brief Wrapper for user-defined types that provide their own mark() implementation.
/// \details T must satisfy the Markable concept by providing void mark(term_mark_stack&) const.
///          This is the extension point for types that contain aterm_core instances but are
///          not fundamental, aterm-derived, std::pair, or std::vector.
template<typename T>
struct markable_aterm
{
private:
  using stored_type = std::decay_t<T>;
  stored_type m_t;

public:
  markable_aterm() = default;

  markable_aterm(const markable_aterm&) = default;
  markable_aterm(markable_aterm&&) = default;
  markable_aterm& operator=(const markable_aterm&) = default;
  markable_aterm& operator=(markable_aterm&&) = default;
  ~markable_aterm() = default;

  markable_aterm(const stored_type& other) noexcept
    : m_t(other)
  {}

  /// \details Constrained such that it does not hijack the copy/move constructors above.
  template<typename... Args>
    requires (sizeof...(Args) != 1 || (!std::is_same_v<std::decay_t<Args>, markable_aterm> && ...))
  markable_aterm(Args&&... args) noexcept
    : m_t(std::forward<Args>(args)...)
  {}

  markable_aterm& operator=(const stored_type& other) noexcept
  {
    m_t = other;
    return *this;
  }

  markable_aterm& operator=(stored_type&& other) noexcept
  {
    m_t = std::move(other);
    return *this;
  }

  operator stored_type&() { return m_t; }
  operator const stored_type&() const { return m_t; }

  bool operator==(const markable_aterm& other) const
    requires std::equality_comparable<stored_type>
  {
    return m_t == other.m_t;
  }

  auto operator<=>(const markable_aterm& other) const
    requires std::three_way_comparable<stored_type>
  {
    return m_t <=> other.m_t;
  }

  void mark(term_mark_stack& todo) const
  {
    static_assert(Markable<stored_type>,
      "T must be a fundamental type, an aterm_core-derived type, a std::pair, "
      "a std::vector, or provide void mark(term_mark_stack&) const");
    m_t.mark(todo);
  }
};

/// \brief Backward-compatible alias for markable_aterm.
template<typename T>
using reference_aterm = markable_aterm<T>;

template<typename T, typename Allocator>
class aterm_allocator
{
public:
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  template <class U>
  struct rebind
  {
    using other = aterm_allocator<U, typename Allocator::template rebind<U>::other>;
  };

  aterm_allocator() = default;

  /// \details The unused parameter is to make the interface equivalent to the allocator.
  T* allocate(size_type n, const void* hint = nullptr)
  {
    return m_allocator.allocate(n, hint);
  }

  /// \details The unused parameter is to make the interface equivalent to the allocator.
  void deallocate(T* p, size_type n);

  // Move assignment and construction is possible.
  aterm_allocator(aterm_allocator&&) = default;
  aterm_allocator& operator=(aterm_allocator&&) = default;

private:
  Allocator m_allocator;
};

template<typename Container>
class generic_aterm_container
{
public:
  /// \brief Constructor. Registers the container with the GC.
  generic_aterm_container(const Container& container)
    : m_container(
        [&container](term_mark_stack& todo) {
          for (const auto& element : container)
          {
            if constexpr (Markable<typename Container::value_type>)
            {
              element.mark(todo);
            }
            else
            {
              markable_aterm<typename Container::value_type>(element).mark(todo);
            }
          }
        },
        [&container]() -> std::size_t { return container.size(); })
  {}

  // Container is a reference so copy/move construction is not meaningful.
  generic_aterm_container(const generic_aterm_container&) = delete;
  generic_aterm_container(generic_aterm_container&&) = delete;

  // The captured reference stays valid; assignment is a no-op.
  generic_aterm_container& operator=(const generic_aterm_container&) { return *this; }
  generic_aterm_container& operator=(generic_aterm_container&&) { return *this; }

protected:
  aterm_container m_container;
};

} // namespace detail
} // namespace atermpp

namespace std
{

/// \brief specialization of the standard std::hash function.
template<typename T>
struct hash<atermpp::detail::markable_aterm<T>>
{
  std::size_t operator()(const atermpp::detail::markable_aterm<T>& t) const
  {
    return std::hash<T>()(t);
  }
};

} // namespace std

#endif // MCRL2_ATERMPP_DETAIL_ATERM_CONTAINER_H
