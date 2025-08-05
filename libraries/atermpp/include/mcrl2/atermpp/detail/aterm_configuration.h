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

namespace atermpp::detail
{

/// \brief Enable garbage collection.
constexpr static bool EnableGarbageCollection = true;

/// \brief Enable the block allocator for terms.
constexpr static bool EnableBlockAllocator = true;

/// \brief Enable to print garbage collection statistics.
constexpr static bool EnableGarbageCollectionMetrics = false;

/// Performs garbage collection intensively for testing purposes.
constexpr static bool EnableAggressiveGarbageCollection = false;

/// \brief Enable to print hashtable collision, size and number of buckets.
constexpr static bool EnableHashtableMetrics = false;

/// \brief Enable to obtain the percentage of terms found compared to allocated.
constexpr static bool EnableCreationMetrics = false;

/// \brief Keep track of the number of variables registered.
constexpr static bool EnableVariableRegistrationMetrics = false;

} // namespace atermpp::detail


#endif // MCRL2_ATERMPP_ATERM_CONFIGURATION_H
