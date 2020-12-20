// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_CONFIGURATION_H
#define MCRL2_ATERMPP_ATERM_CONFIGURATION_H

namespace atermpp
{
namespace detail
{

/// \brief Enables thread safety for the global term and function symbol pools.
constexpr static bool GlobalThreadSafe = false;

/// \brief Enable to print garbage collection statistics.
constexpr static bool EnableGarbageCollectionMetrics = false;

/// \brief Enable to print hashtable collision, size and number of buckets.
constexpr static bool EnableHashtableMetrics = false;

/// \brief Enable to obtain the percentage of terms found compared to allocated.
constexpr static bool EnableCreationMetrics = false;

/// \brief Enable garbage collection.
constexpr static bool EnableGarbageCollection = true && !GlobalThreadSafe;

} // namespace detail
} // namespace atermpp

#endif // MCRL2_ATERMPP_ATERM_CONFIGURATION_H
