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
MCRL2_MODULE;

#ifdef MCRL2_ENABLE_MODULES
  export module atermpp:detail.aterm_configuration;
#endif

MCRL2_MODULE_EXPORT namespace atermpp::detail
{

/// \brief Enable garbage collection.
constexpr bool EnableGarbageCollection = true;

/// \brief Enable the block allocator for terms.
constexpr bool EnableBlockAllocator = true;

/// \brief Enable to print garbage collection statistics.
constexpr bool EnableGarbageCollectionMetrics = false;

/// Performs garbage collection intensively for testing purposes.
constexpr bool EnableAggressiveGarbageCollection = false;

/// \brief Enable to print hashtable collision, size and number of buckets.
constexpr bool EnableHashtableMetrics = false;

/// \brief Enable to obtain the percentage of terms found compared to allocated.
constexpr bool EnableCreationMetrics = false;

/// \brief Keep track of the number of variables registered.
constexpr bool EnableVariableRegistrationMetrics = false;

} // namespace atermpp::detail


#endif // MCRL2_ATERMPP_ATERM_CONFIGURATION_H
