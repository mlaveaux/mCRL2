// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "benchmark_shared.h"

#include "mcrl2/utilities/stopwatch.h"

int main(int, char*[])
{
  std::size_t amount = 1500000;
  std::size_t iterations = 50;

  detail::g_term_pool().enable_garbage_collection(false);

  // Keep track of the total time spend on these parts.
  long long creation_time = 0;
  long long garbage_collect_time = 0;

  // Generate various function applications.
  aterm_appl f1 = create_nested_function("f", "c", 1, 1 * amount);
  aterm_appl f2 = create_nested_function("g", "d", 2, 4 * amount);
  aterm_appl f4 = create_nested_function("h", "e", 4, amount/4);
  aterm_appl f7 = create_nested_function("i", "f", 7, amount/8);

  for (std::size_t i = 0; i < iterations; ++i)
  {
    mcrl2::utilities::stopwatch stopwatch;
    {
      // These will be cleaned up every iteration.
      aterm_appl f2 = create_nested_function("j", "g", 2, 4 * amount);
    }
    creation_time += stopwatch.time();

    // Trigger garbage collection.
    stopwatch.reset();
    detail::g_term_pool().collect();
    garbage_collect_time += stopwatch.time();
  }

  std::cerr << "Creating terms took " << creation_time << " milliseconds.\n";
  std::cerr << "Garbage collection took " << garbage_collect_time << " milliseconds.\n";
}
