// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/utilities/cache_metric.h"
#include <sstream>

using namespace mcrl2::utilities;

std::string cache_metric::message() const
{
  std::stringstream str;
  std::size_t total_count = m_hit_count + m_miss_count;
  str << m_hit_count << " times found out of " << total_count << " calls (" << (static_cast<double>(m_hit_count) / static_cast<double>(total_count)) * 100.0 << " %)";
  return str.str();
}
