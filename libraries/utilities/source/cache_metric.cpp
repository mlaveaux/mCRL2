// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
MCRL2_MODULE;

#include <sstream>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/cache_metric.cxx"
#else
  module utilities;

  import :cache_metric;
#endif

using namespace mcrl2::utilities;

std::string cache_metric::message() const
{
  std::stringstream str;
  std::size_t total_count = m_hit_count + m_miss_count;
  str << m_hit_count << " times found out of " << total_count << " calls (" << static_cast<double>(m_hit_count) / static_cast<double>(m_miss_count) * 100 << " %)";
  return str.str();
}
