// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file standard_containers_test.cpp
/// \brief Regression tests for the GC-aware standard container wrappers, which
///        expose the standard API by composition instead of inheritance.

#define BOOST_TEST_MODULE standard_containers_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/atermpp/aterm_int.h"
#include "mcrl2/atermpp/detail/global_aterm_pool.h"
#include "mcrl2/atermpp/standard_containers/deque.h"
#include "mcrl2/atermpp/standard_containers/indexed_set.h"
#include "mcrl2/atermpp/standard_containers/stack.h"
#include "mcrl2/atermpp/standard_containers/unordered_map.h"
#include "mcrl2/atermpp/standard_containers/unordered_set.h"
#include "mcrl2/atermpp/standard_containers/vector.h"

#include <algorithm>

using namespace atermpp;

namespace
{

aterm term(std::size_t i)
{
  return aterm_int(i);
}

/// \brief Forces garbage collection; the contents of all registered containers must survive.
void collect()
{
  detail::g_thread_term_pool().collect();
}

} // namespace

BOOST_AUTO_TEST_CASE(deque_basic_operations)
{
  atermpp::deque<aterm> d;
  BOOST_CHECK(d.empty());

  d.push_back(term(1));
  d.push_front(term(0));
  d.emplace_back(term(2));
  d.emplace_front(term(3));

  // Contents: 3 0 1 2
  BOOST_CHECK_EQUAL(d.size(), 4);
  BOOST_CHECK(d.front() == term(3));
  BOOST_CHECK(d.back() == term(2));
  BOOST_CHECK(d[1] == term(0));
  BOOST_CHECK(d.at(2) == term(1));

  collect();

  // pop_front was missing from the wrapper before the composition rewrite.
  d.pop_front();
  BOOST_CHECK(d.front() == term(0));
  d.pop_back();
  BOOST_CHECK(d.back() == term(1));
  BOOST_CHECK_EQUAL(d.size(), 2);

  // Iteration in both directions.
  BOOST_CHECK(std::distance(d.begin(), d.end()) == 2);
  BOOST_CHECK(std::distance(d.rbegin(), d.rend()) == 2);
  BOOST_CHECK(*d.cbegin() == term(0));
  BOOST_CHECK(*d.crbegin() == term(1));
}

BOOST_AUTO_TEST_CASE(deque_assign_resize_clear)
{
  atermpp::deque<aterm> d;

  // assign was missing from the wrapper before the composition rewrite.
  d.assign(10, term(42));
  BOOST_CHECK_EQUAL(d.size(), 10);
  BOOST_CHECK(std::all_of(d.begin(), d.end(), [](const aterm& t) { return t == term(42); }));

  std::vector<aterm> values = {term(1), term(2), term(3)};
  d.assign(values.begin(), values.end());
  BOOST_CHECK_EQUAL(d.size(), 3);
  BOOST_CHECK(d.front() == term(1));

  collect();
  BOOST_CHECK(d.front() == term(1));

  d.resize(5);
  BOOST_CHECK_EQUAL(d.size(), 5);
  d.resize(2, term(9));
  BOOST_CHECK_EQUAL(d.size(), 2);

  d.clear();
  BOOST_CHECK(d.empty());
}

BOOST_AUTO_TEST_CASE(deque_insert_erase_swap_compare)
{
  atermpp::deque<aterm> d = {term(1), term(2)};
  atermpp::deque<aterm> e(d);
  BOOST_CHECK(d == e);

  e.push_back(term(3));
  BOOST_CHECK(d != e);

  // Three-way comparison is consistent with equality.
  BOOST_CHECK((d <=> d) == 0);
  BOOST_CHECK((d <=> e) != 0);

  auto it = d.insert(d.begin(), term(0));
  BOOST_CHECK(*it == term(0));
  it = d.erase(d.begin());
  BOOST_CHECK(*it == term(1));

  d.swap(e);
  BOOST_CHECK_EQUAL(d.size(), 3);
  BOOST_CHECK_EQUAL(e.size(), 2);

  // Copy and move assignment run under the GC guard.
  e = d;
  BOOST_CHECK(e == d);
  atermpp::deque<aterm> f;
  f = std::move(e);
  BOOST_CHECK(f == d);
}

BOOST_AUTO_TEST_CASE(vector_basic_operations)
{
  atermpp::vector<aterm> v;
  BOOST_CHECK(v.empty());

  for (std::size_t i = 0; i < 100; ++i)
  {
    v.push_back(term(i));
  }
  BOOST_CHECK_EQUAL(v.size(), 100);
  BOOST_CHECK(v.front() == term(0));
  BOOST_CHECK(v.back() == term(99));
  BOOST_CHECK(v[10] == term(10));
  BOOST_CHECK(v.at(20) == term(20));

  collect();
  BOOST_CHECK(v[10] == term(10));

  v.pop_back();
  BOOST_CHECK(v.back() == term(98));

  v.reserve(1000);
  BOOST_CHECK(v.capacity() >= 1000);
  v.shrink_to_fit();

  // assign was missing from the wrapper before the composition rewrite.
  v.assign(4, term(7));
  BOOST_CHECK_EQUAL(v.size(), 4);
  BOOST_CHECK(v.front() == term(7));

  v.assign({term(1), term(2)});
  BOOST_CHECK_EQUAL(v.size(), 2);
  BOOST_CHECK(v.back() == term(2));

  atermpp::vector<aterm> w(v);
  BOOST_CHECK(v == w);
  w.emplace_back(term(3));
  BOOST_CHECK(v != w);
  BOOST_CHECK((v <=> v) == 0);

  v = w;
  BOOST_CHECK(v == w);

  v.clear();
  BOOST_CHECK(v.empty());
}

BOOST_AUTO_TEST_CASE(vector_thread_safe_variant)
{
  // The ThreadSafe variant acquires the exclusive lock for modifications.
  atermpp::vector<aterm, std::allocator<detail::markable_aterm<aterm>>, true> v;
  v.push_back(term(1));
  v.emplace_back(term(2));
  BOOST_CHECK_EQUAL(v.size(), 2);
  v.pop_back();
  v.clear();
  BOOST_CHECK(v.empty());
}

BOOST_AUTO_TEST_CASE(unordered_set_basic_operations)
{
  atermpp::unordered_set<aterm> s;
  BOOST_CHECK(s.empty());

  for (std::size_t i = 0; i < 50; ++i)
  {
    auto [it, inserted] = s.insert(term(i));
    BOOST_CHECK(inserted);
  }
  BOOST_CHECK_EQUAL(s.size(), 50);

  collect();

  // find/count/contains were inherited without GC protection before the rewrite.
  BOOST_CHECK(s.contains(term(25)));
  BOOST_CHECK(!s.contains(term(1000)));
  BOOST_CHECK(s.find(term(10)) != s.end());
  BOOST_CHECK(s.find(term(1000)) == s.end());
  BOOST_CHECK_EQUAL(s.count(term(10)), 1);

  BOOST_CHECK_EQUAL(s.erase(term(10)), 1);
  BOOST_CHECK(!s.contains(term(10)));
  BOOST_CHECK_EQUAL(s.size(), 49);

  auto [it2, inserted2] = s.emplace(term(10));
  BOOST_CHECK(inserted2);

  atermpp::unordered_set<aterm> t(s);
  BOOST_CHECK(t == s);
  t = s;
  BOOST_CHECK(t == s);

  t.clear();
  BOOST_CHECK(t.empty());

  s.swap(t);
  BOOST_CHECK(s.empty());
  BOOST_CHECK_EQUAL(t.size(), 50);
}

BOOST_AUTO_TEST_CASE(utilities_unordered_set_basic_operations)
{
  atermpp::utilities::unordered_set<aterm> s;
  BOOST_CHECK(s.empty());

  for (std::size_t i = 0; i < 50; ++i)
  {
    auto [it, inserted] = s.emplace(term(i));
    BOOST_CHECK(inserted);
  }
  BOOST_CHECK_EQUAL(s.size(), 50);

  collect();

  BOOST_CHECK(s.contains(term(25)));
  BOOST_CHECK(!s.contains(term(1000)));
  BOOST_CHECK(s.find(term(10)) != s.end());
  BOOST_CHECK_EQUAL(s.count(term(10)), 1);

  s.erase(term(10));
  BOOST_CHECK(!s.contains(term(10)));

  std::size_t count = 0;
  for (const aterm& t: s)
  {
    static_cast<void>(t);
    count++;
  }
  BOOST_CHECK_EQUAL(count, s.size());

  s.clear();
  BOOST_CHECK(s.empty());
}

BOOST_AUTO_TEST_CASE(unordered_map_basic_operations)
{
  atermpp::unordered_map<aterm, aterm> m;
  BOOST_CHECK(m.empty());

  // operator[] was inherited without GC protection before the rewrite.
  m[term(1)] = term(10);
  m[term(2)] = term(20);
  BOOST_CHECK_EQUAL(m.size(), 2);

  collect();

  BOOST_CHECK(m.at(term(1)) == term(10));
  BOOST_CHECK(m.contains(term(2)));
  BOOST_CHECK(!m.contains(term(3)));
  BOOST_CHECK_EQUAL(m.count(term(1)), 1);
  BOOST_CHECK(m.find(term(2)) != m.end());
  BOOST_CHECK(m.find(term(3)) == m.end());

  auto [it, inserted] = m.try_emplace(term(3), term(30));
  BOOST_CHECK(inserted);
  auto [it2, inserted2] = m.try_emplace(term(3), term(31));
  BOOST_CHECK(!inserted2);
  BOOST_CHECK(m.at(term(3)) == term(30));

  auto [it3, inserted3] = m.insert_or_assign(term(3), term(33));
  BOOST_CHECK(!inserted3);
  BOOST_CHECK(m.at(term(3)) == term(33));

  BOOST_CHECK_EQUAL(m.erase(term(3)), 1);
  BOOST_CHECK(!m.contains(term(3)));

  atermpp::unordered_map<aterm, aterm> n(m);
  BOOST_CHECK(n == m);
  n = m;
  BOOST_CHECK(n == m);

  std::size_t count = 0;
  for (const auto& [key, value]: m)
  {
    static_cast<void>(key);
    static_cast<void>(value);
    count++;
  }
  BOOST_CHECK_EQUAL(count, m.size());

  m.clear();
  BOOST_CHECK(m.empty());
}

BOOST_AUTO_TEST_CASE(unordered_map_fundamental_mapped_type)
{
  // Fundamental mapped types are stored in the fundamental markable_aterm specialization.
  atermpp::unordered_map<aterm, std::size_t> m;
  m[term(1)] = 42;

  std::size_t value = m.at(term(1));
  BOOST_CHECK_EQUAL(value, 42);

  collect();

  m[term(1)] = 43;
  value = m[term(1)];
  BOOST_CHECK_EQUAL(value, 43);
}

BOOST_AUTO_TEST_CASE(utilities_unordered_map_basic_operations)
{
  atermpp::utilities::unordered_map<aterm, aterm> m;
  BOOST_CHECK(m.empty());

  for (std::size_t i = 0; i < 100; ++i)
  {
    m.emplace(term(i), term(i + 1000));
  }
  BOOST_CHECK_EQUAL(m.size(), 100);

  collect();

  BOOST_CHECK(m.find(term(10)) != m.end());
  BOOST_CHECK((*m.find(term(10))).second == term(1010));
  BOOST_CHECK_EQUAL(m.count(term(10)), 1);

  // operator[] inserts a default value for absent keys.
  aterm& value = m[term(1000)];
  BOOST_CHECK(!value.defined());
  BOOST_CHECK_EQUAL(m.size(), 101);
  value = term(0);
  BOOST_CHECK(m.at(term(1000)) == term(0));

  m.erase(term(1000));
  BOOST_CHECK(m.find(term(1000)) == m.end());

  std::size_t count = 0;
  for (const auto& element: m)
  {
    static_cast<void>(element);
    count++;
  }
  BOOST_CHECK_EQUAL(count, m.size());

  m.clear();
  BOOST_CHECK(m.empty());
}

BOOST_AUTO_TEST_CASE(stack_basic_operations)
{
  atermpp::stack<aterm> s;
  BOOST_CHECK(s.empty());

  s.push(term(1));
  s.push(term(2));

  // Since C++17, emplace returns a reference to the created element.
  aterm& top = s.emplace(term(3));
  BOOST_CHECK(top == term(3));

  BOOST_CHECK_EQUAL(s.size(), 3);
  BOOST_CHECK(s.top() == term(3));

  collect();
  BOOST_CHECK(s.top() == term(3));

  s.pop();
  BOOST_CHECK(s.top() == term(2));

  atermpp::stack<aterm> t;
  t.swap(s);
  BOOST_CHECK(s.empty());
  BOOST_CHECK_EQUAL(t.size(), 2);
}

BOOST_AUTO_TEST_CASE(indexed_set_basic_operations)
{
  atermpp::indexed_set<aterm> s;

  auto [index1, inserted1] = s.insert(term(1));
  auto [index2, inserted2] = s.insert(term(2));
  BOOST_CHECK(inserted1);
  BOOST_CHECK(inserted2);
  BOOST_CHECK(index1 != index2);

  auto [index3, inserted3] = s.insert(term(1));
  BOOST_CHECK(!inserted3);
  BOOST_CHECK_EQUAL(index1, index3);

  collect();

  BOOST_CHECK_EQUAL(s.size(), 2);
  BOOST_CHECK_EQUAL(s.index(term(1)), index1);
  BOOST_CHECK_EQUAL(s.index(term(42)), atermpp::indexed_set<aterm>::npos);

  // contains performs the lookup under a single guard.
  BOOST_CHECK(s.contains(term(1)));
  BOOST_CHECK(!s.contains(term(42)));

  BOOST_CHECK(s.at(index1) == term(1));
  BOOST_CHECK(s[index2] == term(2));
  BOOST_CHECK(s.find(term(2)) != s.end());
  BOOST_CHECK(s.find(term(42)) == s.end());

  // The free contains function used by generic algorithms.
  BOOST_CHECK(mcrl2::utilities::detail::contains(s, term(1)));
  BOOST_CHECK(!mcrl2::utilities::detail::contains(s, term(42)));

  s.clear();
  BOOST_CHECK_EQUAL(s.size(), 0);
}

BOOST_AUTO_TEST_CASE(containers_survive_gc)
{
  // Fill all containers, force garbage collection repeatedly, and check contents.
  atermpp::deque<aterm> d;
  atermpp::vector<aterm> v;
  atermpp::unordered_set<aterm> s;
  atermpp::unordered_map<aterm, aterm> m;

  for (std::size_t i = 0; i < 100; ++i)
  {
    d.push_back(term(i));
    v.push_back(term(i));
    auto result = s.insert(term(i));
    BOOST_CHECK(result.second);
    m[term(i)] = term(i + 1000);
  }

  for (std::size_t i = 0; i < 10; ++i)
  {
    collect();
  }

  for (std::size_t i = 0; i < 100; ++i)
  {
    BOOST_CHECK(d[i] == term(i));
    BOOST_CHECK(v[i] == term(i));
    BOOST_CHECK(s.contains(term(i)));
    BOOST_CHECK(m.at(term(i)) == term(i + 1000));
  }
}
