// Author(s): Maurice Laveaux and Menno Bartels
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <ranges>

namespace mcrl2::pbes_system::detail
{

/// Cartesian product of two ranges, should be available in C++23.
template<typename Range1, typename Range2>
class cartesian_product_view : public std::ranges::view_interface<cartesian_product_view<Range1, Range2>>
{
public:
  using value_type
    = std::pair<typename std::ranges::range_value_t<Range1>, typename std::ranges::range_value_t<Range2>>;

  cartesian_product_view() = default;

  cartesian_product_view(Range1 r1, Range2 r2)
    : m_range1(r1),
      m_range2(r2)
  {}

  class iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using iterator_concept = std::forward_iterator_tag;
    using value_type = cartesian_product_view::value_type;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type;

    iterator() = default;

    iterator(auto it1, auto it1_end, auto it2, auto it2_begin, auto it2_end)
      : m_it1(it1),
        m_it1_end(it1_end),
        m_it2(it2),
        m_it2_begin(it2_begin),
        m_it2_end(it2_end)
    {
      if (m_it1 != m_it1_end && m_it2 == m_it2_end)
      {
        m_it1 = m_it1_end; // End iterator
      }
    }

    value_type operator*() const { return std::make_pair(*m_it1, *m_it2); }

    iterator& operator++()
    {
      ++m_it2;
      if (m_it2 == m_it2_end)
      {
        ++m_it1;
        if (m_it1 != m_it1_end)
        {
          m_it2 = m_it2_begin;
        }
      }
      return *this;
    }

    iterator operator++(int)
    {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const iterator& other) const
    {
      return m_it1 == other.m_it1 && (m_it1 == m_it1_end || m_it2 == other.m_it2);
    }

    bool operator!=(const iterator& other) const
    {
      return !(*this == other);
    }

  private:
    decltype(std::ranges::begin(std::declval<const Range1&>())) m_it1, m_it1_end;
    decltype(std::ranges::begin(std::declval<const Range2&>())) m_it2, m_it2_begin, m_it2_end;
  };

  iterator begin() const
  {
    auto it1 = std::ranges::begin(m_range1);
    auto it1_end = std::ranges::end(m_range1);
    auto it2_begin = std::ranges::begin(m_range2);
    auto it2_end = std::ranges::end(m_range2);
    return iterator(it1, it1_end, it2_begin, it2_begin, it2_end);
  }

  iterator end() const
  {
    auto it1_end = std::ranges::end(m_range1);
    auto it2_begin = std::ranges::begin(m_range2);
    auto it2_end = std::ranges::end(m_range2);
    return iterator(it1_end, it1_end, it2_end, it2_begin, it2_end);
  }

private:
  Range1 m_range1;
  Range2 m_range2;
};

template<typename Range1, typename Range2>
auto cartesian_product(Range1 r1, Range2 r2)
{
  return cartesian_product_view<Range1, Range2>(r1, r2);
}

} // namespace mcrl2::pbes_system::detail
