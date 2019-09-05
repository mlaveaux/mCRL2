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

/// \brief Benchmark the creation of function symbols
int main(int, char*[])
{
  std::size_t amount = 50000;
  std::size_t iterations = 1000;
  std::size_t number_of_threads = 1;

  // Define a function that repeatedly creates function symbols.
  auto create_function_symbols = [amount, iterations, number_of_threads](void) -> void
    {
      // Store them in a vector to prevent them from being deleted.
      std::vector<function_symbol> symbols(amount);

      // Track the time that the first iteration (when the terms are created) takes.
      mcrl2::utilities::stopwatch stopwatch;
      bool first_run = true;

      // Generate function symbols f + suffix, where suffix from 0 to amount.
      std::string name("f");
      for (std::size_t k = 0; k < iterations / number_of_threads; ++k)
      {
        for (std::size_t i = 0; i < amount; ++i)
        {
          symbols[i] = function_symbol(name + std::to_string(i), 0);
        }

        if (first_run)
        {
          first_run = false;
          std::cerr << "Creating " << amount << " function symbols took " << stopwatch.time() << " milliseconds.\n";
        }
      }
    };

  benchmark_threads(number_of_threads, create_function_symbols);
}
