// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "benchmark_shared.h"

#include "mcrl2/atermpp/aterm_int.h"
#include "mcrl2/utilities/stopwatch.h"

using namespace atermpp;

int main(int, char*[])
{
  std::size_t amount = 3000000;
  std::size_t iterations = 1000;
  std::size_t number_of_threads = 1;

  // Define a function that repeatedly creates integers.
  auto create_integers = [amount, iterations, number_of_threads](void) -> void
    {
      // Store the integers in a vector to prevent them from being deleted.
      std::vector<aterm_int> integers(amount);

      // Track the time that the first iteration (when the terms are created) takes.
      mcrl2::utilities::stopwatch stopwatch;
      bool first_run = true;

      for (std::size_t i = 0; i < iterations / number_of_threads; ++i)
      {
        for (std::size_t j = 0; j < amount; ++j)
        {
          integers[i] = aterm_int(i);
        }

        if (first_run)
        {
          first_run = false;
          std::cerr << "Creating " << amount << " integers took " << stopwatch.time() << " milliseconds.\n";
        }
      }
    };

  detail::g_term_pool().enable_garbage_collection(false);
  benchmark_threads(number_of_threads, create_integers);

  return 0;
}
