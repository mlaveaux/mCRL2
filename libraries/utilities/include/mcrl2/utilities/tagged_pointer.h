// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_TAGGED_POINTER_H_
#define MCRL2_UTILITIES_TAGGED_POINTER_H_

#include <compare>
#include <functional>
#include <type_traits>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/configuration.h"
  #include "mcrl2/utilities/detail/atomic_wrapper.h"
#else
  import utilities;
#endif

namespace mcrl2::utilities
{

/// \brief Applies a tag to a pointer.
/// \returns A pointer where one bit that is unused is set to one.
template<typename T>
T* tag(const T* p)
{
  return reinterpret_cast<T*>(reinterpret_cast<std::size_t>(p) | 1);
}

template<typename T>
T* tag(const detail::atomic_wrapper<T*>& p)
{
  return reinterpret_cast<T*>(reinterpret_cast<std::size_t>(p.load(std::memory_order_relaxed)) | 1);
}

/// \returns The original pointer that can be deferenced.
template<typename T>
T* pointer(const T* p)
{
  return reinterpret_cast<T*>(reinterpret_cast<std::size_t>(p) & ~static_cast<std::size_t>(1));
}

template<typename T>
T* pointer(const detail::atomic_wrapper<T*>& p)
{
  return reinterpret_cast<T*>(reinterpret_cast<std::size_t>(p.load(std::memory_order_relaxed)) & ~static_cast<std::size_t>(1));
}

/// \returns True iff this pointer has been tagged.
template<typename T>
bool tagged(const T* p)
{
  return reinterpret_cast<std::size_t>(p) & 1;
}

template<typename T>
bool tagged(const detail::atomic_wrapper<T*>& p)
{
  return reinterpret_cast<std::size_t>(p.load(std::memory_order_relaxed)) & 1;
}

/// \brief A pointer storage object that uses a least significant bit as a mark.
///        Can be used by objects that are 8 bytes aligned in memory.
template<typename T>
class tagged_pointer
{
public:
  tagged_pointer() = default;
  explicit tagged_pointer(T* p)
    : m_pointer(p)
  {}

  /// \returns True iff this pointer has been tagged,
  bool tagged() const
  {
    return mcrl2::utilities::tagged(m_pointer);
  }

  /// \brief Apply a tag to the pointer that can be checked with tagged().
  void tag() const
  {
    m_pointer = mcrl2::utilities::tag(m_pointer);
  }

  /// \brief Remove the tag.
  void untag() const
  {
    m_pointer = const_cast<T*>(get());
  }

  bool defined() const
  {
    return get() != nullptr;
  }

  tagged_pointer& operator=(std::nullptr_t)
  {
    m_pointer = mcrl2::utilities::tag(static_cast<T*>(nullptr));
    return *this;
  }

  bool operator==(std::nullptr_t) const
  {
    return get() == nullptr;
  }

  bool operator==(const tagged_pointer& other) const
  {
    return get() == other.get();
  }

  std::weak_ordering operator<=>(const tagged_pointer& other) const
  {
    return this->get() <=> other.get();
  }

  const T& operator*() const
  {
    return *get();
  }
  
  T& operator*()
  {
    return *get();
  }

  const T* operator->() const
  {
    return get();
  }

  T* operator->()
  {
    return get();
  }

  const T* get() const
  {
    return pointer(m_pointer);
  }

  T* get()
  {
    return pointer(m_pointer);
  }

  void swap(tagged_pointer<T>& other) noexcept
  {
    // This is not atomic, so swaps are only allowed when no tags are being applied concurrently.
    auto tmp = other.m_pointer;
    other.m_pointer = m_pointer;
    m_pointer = tmp;
  }

private:
  mutable std::conditional_t<detail::GlobalThreadSafe, detail::atomic_wrapper<T*>, T*> m_pointer = nullptr;
};

} // namespace mcrl2::utilities

namespace std
{

/// \brief specialization of the standard std::hash function.
template<typename T>
struct hash<mcrl2::utilities::tagged_pointer<T>>
{
  std::size_t operator()(const mcrl2::utilities::tagged_pointer<T>& p) const
  {
    const std::hash<const T*> hasher;
    return hasher(p.get());
  }
};

} // namespace std

#endif // MCRL2_UTILITIES_TAGGED_POINTER_H
