// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#define BOOST_TEST_MODULE parallel_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/utilities/configuration.h"
#include "mcrl2/atermpp/aterm_int.h"
#include "mcrl2/atermpp/standard_containers/deque.h"
#include "mcrl2/atermpp/standard_containers/unordered_map.h"
#include "mcrl2/atermpp/standard_containers/unordered_set.h"
#include "mcrl2/atermpp/standard_containers/vector.h"
#include "mcrl2/atermpp/detail/global_aterm_pool.h"

#include <thread>

using namespace atermpp;

namespace
{

/// \brief Runs the given container mutation loop in a thread while the main
///        thread performs garbage collection extensively.
template<typename F>
void run_with_gc_stress(F&& f)
{
  std::thread local(std::forward<F>(f));

  while (true)
  {
    for (std::size_t i = 0; i < 1000; ++i)
    {
      atermpp::detail::g_thread_term_pool().collect();
    }

    if (local.joinable())
    {
      local.join();
      break;
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(parallel_vector)
{
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
  {
    // One thread continuously modifies a local atermpp::vector of aterms while the main thread performs garbage collection extensively.
    run_with_gc_stress([]()
    {
      // Push a constant to avoid term creations.
      atermpp::aterm_int value(0);

      for (std::size_t i = 0; i < 1000; ++i)
      {
        atermpp::vector<atermpp::aterm> vector;
        for (std::size_t i = 0; i < 1000; ++i)
        {
          vector.push_back(value);
        }
      }
    });
  }
}

BOOST_AUTO_TEST_CASE(parallel_deque_pop_front)
{
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
  {
    // pop_front and assign are structural modifications that must hold the GC guard;
    // they were unguarded before the composition rewrite.
    run_with_gc_stress([]()
    {
      atermpp::aterm_int value(0);

      for (std::size_t i = 0; i < 500; ++i)
      {
        atermpp::deque<atermpp::aterm> deque;
        deque.assign(500, value);

        while (!deque.empty())
        {
          deque.pop_front();
        }
      }
    });
  }
}

BOOST_AUTO_TEST_CASE(parallel_unordered_set)
{
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
  {
    run_with_gc_stress([]()
    {
      atermpp::aterm_int value(0);

      for (std::size_t i = 0; i < 500; ++i)
      {
        atermpp::unordered_set<atermpp::aterm> set;
        for (std::size_t j = 0; j < 100; ++j)
        {
          auto result = set.insert(atermpp::aterm_int(j));
          BOOST_CHECK(result.second);
        }
        for (std::size_t j = 0; j < 100; ++j)
        {
          set.erase(atermpp::aterm_int(j));
        }
        static_cast<void>(value);
      }
    });
  }
}

BOOST_AUTO_TEST_CASE(parallel_unordered_map_subscript)
{
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
  {
    // operator[] inserts elements and can trigger a rehash; it was unguarded
    // before the composition rewrite.
    run_with_gc_stress([]()
    {
      for (std::size_t i = 0; i < 500; ++i)
      {
        atermpp::unordered_map<atermpp::aterm, atermpp::aterm> map;
        for (std::size_t j = 0; j < 100; ++j)
        {
          map[atermpp::aterm_int(j)] = atermpp::aterm_int(j + 1);
        }
        map.clear();
      }
    });
  }
}

BOOST_AUTO_TEST_CASE(parallel_container_assignment)
{
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
  {
    // Copy and move assignment are structural modifications that must hold the
    // GC guard; they were defaulted and unguarded before the composition rewrite.
    run_with_gc_stress([]()
    {
      atermpp::aterm_int value(0);

      atermpp::vector<atermpp::aterm> source(100, value);
      for (std::size_t i = 0; i < 1000; ++i)
      {
        atermpp::vector<atermpp::aterm> copy;
        copy = source;
        atermpp::vector<atermpp::aterm> moved;
        moved = std::move(copy);
      }
    });
  }
}

