// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_STOPWATCH_H_
#define MCRL2_UTILITIES_STOPWATCH_H_

#include <chrono>

namespace mcrl2
{
namespace utilities
{

/// \brief Implements a simple stopwatch that starts on construction.
class stopwatch
{
public:
  stopwatch()
  {
    reset();
  }

  /// \brief Reset the stopwatch to count from this moment onwards.
  void reset()
  {
    m_timestamp = std::chrono::system_clock::now();
  }

  /// \returns The time in milliseconds since the last reset.
  long long time()
  {
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - m_timestamp).count();
  }

private:
  std::chrono::time_point<std::chrono::system_clock> m_timestamp;
};

} // namespace utilities
} // namespace mcrl2

#endif // MCRL2_UTILITIES_STOPWATCH_H_
