// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_SHARED_REFERENCE_H_
#define MCRL2_UTILITIES_SHARED_REFERENCE_H_
MCRL2_MODULE;

#include <cassert>
#include <atomic>
#include <compare>
#include <type_traits>



#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/tagged_pointer.cxx"
#else
  export module utilities:shared_reference;

  import :tagged_pointer;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities
{

/// \brief Enable to count the number of reference count changes.
constexpr bool EnableReferenceCountMetrics = false;

/// \brief Stores a reference count that can be incremented and decremented.
/// \details The template variable is such that reference_count_changes corresponds
///          to the amount of times the reference count of type T changed.
template<typename T, bool ThreadSafe = false>
class shared_reference_counted
{
public:
  shared_reference_counted()
    : m_reference_count(0)
  {}

  /// \brief Obtain the reference count.
  std::size_t reference_count() const
  {
    return m_reference_count;
  }

  /// \brief Increment the reference count by one.
  void increment_reference_count() const
  {
    if constexpr (ThreadSafe)
    {
      m_reference_count.fetch_add(1, std::memory_order_relaxed);
    }
    else
    {
      ++m_reference_count;
    }
    count_reference_count_changes();
  }

  /// \brief Decrement the reference count by one.
  void decrement_reference_count() const
  {
    if constexpr (ThreadSafe)
    {
      m_reference_count.fetch_sub(1, std::memory_order_release);
    }
    else
    {
      --m_reference_count;
    }
    count_reference_count_changes();
  }

  /// \brief Obtain the number of times that this reference count has changed.
  static std::atomic<std::size_t>& reference_count_changes()
  {
    static std::atomic<std::size_t> g_reference_count_changes;
    return g_reference_count_changes;
  }

  /// \brief Increment the number of reference count changes.
  static void count_reference_count_changes()
  {
    if constexpr (EnableReferenceCountMetrics)
    {
      ++reference_count_changes();
    }
  }

protected:
  using SizeType = std::conditional_t<ThreadSafe, std::atomic<std::size_t>, std::size_t>;

  // The underlying reference counter can always be changed.
  mutable SizeType m_reference_count;
};

/// \brief A reference counted reference to a shared_reference_counted object.
/// \details Similar to a shared_ptr except that reference counts are only atomic when
///          thread safety is desired and that it stores the reference count in the
///          inherited object.
template<typename T>
class shared_reference
{
public:

  /// \brief The default constructor.
  shared_reference() noexcept
    : m_reference(nullptr)
  {}

  /// \brief Takes ownership of the passed reference, which means
  ///        that its reference count is incremented.
  shared_reference(T* reference) noexcept
    : m_reference(reference)
  {
    assert(defined());
    m_reference->increment_reference_count();
  }

  /// \brief Copy constructor.
  shared_reference(const shared_reference<T>& other) noexcept
    : m_reference(other.m_reference)
  {
    if (defined())
    {
      m_reference->increment_reference_count();
    }
  }

  /// \brief Move constructor.
  shared_reference(shared_reference<T>&& other) noexcept
    : m_reference(other.m_reference)
  {
    other.m_reference = nullptr;
  }

  ~shared_reference()
  {
    if (defined())
    {
      m_reference->decrement_reference_count();
    }
  }

  /// \brief Copy assignment constructor.
  shared_reference<T>& operator=(const shared_reference<T>& other) noexcept
  {
    // Increment first to prevent the same reference from getting a reference count of zero temporarily.
    if (other.defined())
    {
      other.m_reference->increment_reference_count();
    }

    // Decrement the reference from the reference that is currently referred to.
    if (defined())
    {
      m_reference->decrement_reference_count();
    }

    m_reference=other.m_reference;
    return *this;
  }

  /// \brief Move assignment constructor.
  shared_reference<T>& operator=(shared_reference<T>&& other) noexcept
  {
    if (defined())
    {
      m_reference->decrement_reference_count();
    }

    m_reference = other.m_reference;
    other.m_reference = nullptr;
    return *this;
  }

  /// \brief Check whether the shared_reference has a valid reference.
  bool defined() const
  {
    return m_reference.get() != nullptr;
  }

  /// Address equality operator.
  bool operator ==(const shared_reference<T>& other) const noexcept
  {
    return m_reference == other.m_reference;
  }

  /// Address comparison operator.
  std::weak_ordering operator <=>(const shared_reference<T>& other) const noexcept
  {
    return m_reference <=> other.m_reference;
  }

  T* operator->() const noexcept
  {
    assert(defined());
    return m_reference.get();
  }

  T* get() const noexcept
  {
    assert(defined());
    return m_reference.get();
  }

  T& operator*() const noexcept
  {
    assert(defined());
    return *m_reference;
  }

  /// \brief Swaps *this with the other shared reference.
  void swap(shared_reference<T>& other) noexcept { m_reference.swap(other.m_reference); }

  bool tagged() const noexcept
  {
    return m_reference.tagged();
  }

  void tag() const
  {
    m_reference.tag();
  }

  void untag() const
  {
    m_reference.untag();
  }

private:
  mutable utilities::tagged_pointer<T> m_reference;
};

} // namespace mcrl2::utilities


namespace std
{

template<typename T>
void swap(mcrl2::utilities::shared_reference<T>& a, mcrl2::utilities::shared_reference<T>& b) noexcept
{
  a.swap(b);
}

} // namespace std

#endif // MCRL2_UTILITIES_SHARED_REFERENCE_H_
