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

using namespace atermpp;

int main(int argc, char* argv[])
{
  std::size_t number_of_arguments = 0;
  std::size_t size = 2000000;
  std::size_t iterations = 1000;
  std::size_t number_of_threads = 1;

  // Accept one argument for the number of arguments.
  if (argc > 1)
  {
    // The first argument is the path of the executable.
    number_of_arguments = static_cast<std::size_t>(std::stoi(argv[1]));
  }

  // Define a function that repeatedly creates nested function applications with a converter.
  auto nested_function = [iterations, number_of_arguments, size, number_of_threads](void) -> void
    {
      // Track the time that the first iteration (when the term is created) takes.
      mcrl2::utilities::stopwatch stopwatch;
      bool first_run = true;

      aterm_appl f;
      for (std::size_t i = 0; i < iterations / number_of_threads; ++i)
      {
        f = create_nested_function<true>("f", "c", number_of_arguments, size / (number_of_arguments + 1));

        if (first_run)
        {
          first_run = false;
          std::cerr << "Creating nested function application with " << number_of_arguments << " arguments and " << (size / (number_of_arguments + 1)) << " depth with converter took " << stopwatch.time() << " milliseconds.\n";
        }
      }
    };

  detail::g_term_pool().enable_garbage_collection(false);
  benchmark_threads(number_of_threads, nested_function);

  return 0;
}
