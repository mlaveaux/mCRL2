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

#include <functional>
#include <type_traits>

namespace mcrl2::utilities
{

/// \brief Given a pointer p to an 8-byte aligned memory applies a tag.
/// \returns A pointer where one bit is set to one that is unused.
template<typename T>
T* tag(const T* p)
{
  return reinterpret_cast<T*>(reinterpret_cast<std::size_t>(p) | 1);
}

/// \returns The original pointer that can be deferenced.
template<typename T>
T* pointer(const T* p)
{
  return reinterpret_cast<T*>(reinterpret_cast<std::size_t>(p) & ~static_cast<std::size_t>(1));
}

/// \returns True iff this pointer has been tagged.
template<typename T>
bool tagged(const T* p)
{
  return reinterpret_cast<std::size_t>(p) & 1;
}

/// \brief A pointer storage object that uses a least significant bit as a mark.
///        Can be used by objects that are 8 bytes aligned in memory.
template<typename T>
class tagged_pointer
{
public:
  tagged_pointer() {}
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

  void operator=(std::nullptr_t)
  {
    m_pointer = mcrl2::utilities::tag(static_cast<T*>(nullptr));
  }

  bool operator==(std::nullptr_t) const
  {
    return get() == nullptr;
  }

  bool operator!=(std::nullptr_t) const
  {
    return !(*this == nullptr);
  }

  bool operator==(const tagged_pointer& other) const
  {
    return get() == other.get();
  }

  bool operator!=(const tagged_pointer& other) const
  {
    return !(*this == other);
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

private:
  mutable T* m_pointer = nullptr;
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
