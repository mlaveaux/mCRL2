// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "benchmark_shared.h"

#include "mcrl2/atermpp/aterm_list.h"
#include "mcrl2/utilities/stopwatch.h"

using namespace atermpp;

int main(int, char*[])
{
  std::size_t length = 1000000;
  std::size_t iterations = 1000;
  std::size_t number_of_threads = 1;

  // Defines a function that repeatedly creates list terms.
  auto create_list = [number_of_threads, iterations, length](void) -> void
    {
      function_symbol c("c", 0);
      aterm_appl c_term(c);

      // Track the time that the first iteration (when the term is created) takes.
      mcrl2::utilities::stopwatch stopwatch;
      bool first_run = true;

      aterm_list list;
      for (std::size_t i = 0; i < iterations / number_of_threads; ++i)
      {
        aterm_list tmp_list;
        for (std::size_t j = 0; j < length; ++j)
        {
          tmp_list.push_front(c_term);
        }
        list = tmp_list;

        if (first_run)
        {
          first_run = false;
          std::cerr << "Creating list of length " << length << " took " << stopwatch.time() << " milliseconds.\n";
        }
      }

    };

  detail::g_term_pool().enable_garbage_collection(false);
  benchmark_threads(number_of_threads, create_list);

  return 0;
}
