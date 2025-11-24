// Author(s): Maurice Laveaux and Menno Bartels
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_PBES_DETAIL_ANY_VIEW_H
#define MCRL2_PBES_DETAIL_ANY_VIEW_H

#include <concepts>
#include <iterator>
#include <memory>
#include <ranges>

namespace mcrl2::pbes_system::detail
{

/// A type-erasing iterator wrapper, this would be so much easier with dyn traits.
template<typename T>
class any_iterator
{
private:
  struct iterator_base
  {
    virtual ~iterator_base() = default;
    virtual std::unique_ptr<iterator_base> clone() const = 0;
    virtual T& operator*() = 0;
    virtual any_iterator& operator++() = 0;
    virtual bool equal(const iterator_base& other) const = 0;
  };

  template<typename Iterator>
  struct iterator_derived : iterator_base
  {
    Iterator it;

    explicit iterator_derived(Iterator&& iter)
      : it(std::forward<Iterator>(iter))
    {}

    std::unique_ptr<iterator_base> clone() const override { return std::make_unique<iterator_derived>(it); }

    T& operator*() override { return *it; }

    any_iterator& operator++() override
    {
      ++it;
      return static_cast<any_iterator&>(*this);
    }

    bool equal(const iterator_base& other) const override
    {
      if (auto* derived_other = dynamic_cast<const iterator_derived*>(&other))
      {
        return it == derived_other->it;
      }
      return false;
    }
  };

  std::unique_ptr<iterator_base> impl;

public:
  using iterator_category = std::input_iterator_tag;
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using pointer = T*;
  using reference = T&;

  any_iterator() = default;

  template<typename Iterator>
  any_iterator(Iterator&& it)
    : impl(std::make_unique<iterator_derived<std::decay_t<Iterator>>>(std::forward<Iterator>(it)))
  {}

  any_iterator(const any_iterator& other)
    : impl(other.impl ? other.impl->clone() : nullptr)
  {}

  any_iterator& operator=(const any_iterator& other)
  {
    if (this != &other)
    {
      impl = other.impl ? other.impl->clone() : nullptr;
    }
    return *this;
  }

  any_iterator(any_iterator&&) = default;
  T& operator*() { return impl->operator*(); }
  any_iterator& operator++() { return impl->operator++(); }
  any_iterator operator++(int)
  {
    any_iterator tmp(*this);
    ++(*this);
    return tmp;
  }

  bool operator==(const any_iterator& other) const { return impl && other.impl && impl->equal(*other.impl); }
};

/// A type-erasing wrapper for ranges that can hold any range type
template<typename T>
class any_view
{
private:
  /// Interface of a range for type erasure
  struct base
  {
    virtual ~base() = default;
    virtual std::unique_ptr<base> clone() const = 0;
    virtual any_iterator<T> begin() = 0;
    virtual any_iterator<T> end() = 0;
  };

  /// Implementation of the range interface for a specific range type
  template<typename Range>
  struct derived : base
  {
    Range range;

    explicit derived(Range&& r)
      : range(std::forward<Range>(r))
    {}

    std::unique_ptr<base> clone() const override { return std::make_unique<derived>(range); }

    any_iterator<T> begin() override { return any_iterator<T>(std::ranges::begin(range)); }

    any_iterator<T> end() override { return any_iterator<T>(std::ranges::end(range)); }
  };

  std::unique_ptr<base> impl;

public:
  template<typename Range>
    requires std::ranges::range<Range> && std::same_as<std::ranges::range_value_t<Range>, T>
  any_view(Range&& range)
    : impl(std::make_unique<derived<std::decay_t<Range>>>(std::forward<Range>(range)))
  {}

  any_view(const any_view& other)
    : impl(other.impl ? other.impl->clone() : nullptr)
  {}
  any_view& operator=(const any_view& other)
  {
    if (this != &other)
    {
      impl = other.impl ? other.impl->clone() : nullptr;
    }
    return *this;
  }

  any_view(any_view&&) = default;
  any_view& operator=(any_view&&) = default;

  auto begin() { return impl->begin(); }
  auto end() { return impl->end(); }
};

// static_assert(std::input_iterator<any_iterator<int>>, "any_iterator should satisfy the input iterator concept");
static_assert(std::ranges::range<any_view<int>>, "any_view should be a range");

} // namespace mcrl2::pbes_system::detail

#endif // MCRL2_PBES_DETAIL_ANY_VIEW_H