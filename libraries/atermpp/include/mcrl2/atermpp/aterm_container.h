// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_CONTAINER_H
#define MCRL2_ATERMPP_ATERM_CONTAINER_H

#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/atermpp/detail/aterm.h"
#include "mcrl2/atermpp/detail/aterm_pool_storage_implementation.h"

#include <stack>

namespace atermpp
{

/// \brief Provides safe storage of unprotected_aterm instances in a container by marking
///        them during garbage collection.
class aterm_container
{
public:
  aterm_container();
  virtual ~aterm_container();

  /// \brief Ensure that all the terms in the containers.
  virtual void mark(std::stack<std::reference_wrapper<detail::_aterm>>& todo) const = 0;
};

/// \brief An unprotected term that is stored inside an aterm_container.
template<typename T>
class reference_aterm : public unprotected_aterm
{
public:
  reference_aterm() noexcept = default;
  reference_aterm(const unprotected_aterm& other) noexcept
  {
    m_term = detail::address(other);
  }

  reference_aterm(unprotected_aterm&& other) noexcept
  {
    m_term = detail::address(other);
  }

  /// Converts implicitly to a protected term of type T.
  operator T&()
  {
    return reinterpret_cast<T&>(*this);
  }

  operator const T&() const
  {
    return reinterpret_cast<const T&>(*this);
  }
};

template<typename Container>
class generic_aterm_container : public aterm_container
{
public:
  /// \brief Provides access to the underlying container.
  Container& container() { return m_container; }
  const Container& container() const { return m_container; }

  void mark(std::stack<std::reference_wrapper<detail::_aterm>>& todo) const override
  {
    for (auto it = m_container.begin(); it != m_container.end(); ++it)
    {
      const aterm& element = *it;

      // Mark all terms (and their subterms) that are reachable, i.e the root set.
      detail::_aterm* term = detail::address(element);
      if (element.defined() && !term->is_marked())
      {
        // Mark the term itself as reachable.
        term->mark();

        // This variable is not a default term and that term has not been marked.
        mark_term(*term, todo);
      }
    }
  }

private:
  Container m_container;
};

} // namespace atermpp

namespace std
{

/// \brief specialization of the standard std::hash function.
template<class T>
struct hash<atermpp::reference_aterm<T>>
{
  std::size_t operator()(const atermpp::reference_aterm<T>& t) const
  {
    return std::hash<atermpp::aterm>()(t);
  }
};

} // namespace std

#endif // MCRL2_ATERMPP_ATERM_CONTAINER_H
