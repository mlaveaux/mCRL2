// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "benchmark_shared.h"

using namespace atermpp;

int main(int argc, char* argv[])
{
  detail::g_term_pool().enable_garbage_collection(false);
  std::size_t number_of_arguments = 1;
  std::size_t number_of_threads = 1;

  // Accept one argument for the number of arguments.
  if (argc > 1)
  {
    // The first argument is the path of the executable.
    number_of_arguments = static_cast<std::size_t>(std::stoi(argv[1]));
  }

  if (argc > 2)
  {
    // The first argument is the path of the executable.
    number_of_threads = static_cast<std::size_t>(std::stoi(argv[2]));
  }

  std::size_t size = 400000;
  std::size_t iterations = 1000;

  std::vector<aterm_appl> terms;
  for (int id = 0; id < number_of_threads; ++id)
  {
    terms.push_back(create_nested_function("f", std::to_string(id), number_of_arguments, size));
  }

  // Define a function that repeatedly creates nested function applications.
  auto nested_function = [&](int id) -> void
    {
      aterm_appl f;
      for (std::size_t i = 0; i < iterations / number_of_threads; ++i)
      {
        f = create_nested_function("f", std::to_string(id), number_of_arguments, size);
      }
    };

  benchmark_threads(number_of_threads, nested_function);
  return 0;
}
