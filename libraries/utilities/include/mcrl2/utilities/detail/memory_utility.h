// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/detail/memory_utility.h
/// \brief This file contains a workaround that allows to
/// assign a small array of elements on the stack. For most
/// systems this can be done by the function alloca, or by:
///
///    type var[n];
///
/// for variable n, but for the MSVC compilers, the
/// special function _alloca must be used.

#ifndef MCRL2_UTILITIES_DETAIL_MEMORY_UTILITY_H
#define MCRL2_UTILITIES_DETAIL_MEMORY_UTILITY_H

#include "mcrl2/utilities/platform.h"

// Reserve a local array of type TYPE and SIZE elements (where SIZE
// is not necessarily a constant value). These will be allocated on the stack,
// but the memory will not be initialised, nor will a destructor be called
// on these memory positions when the reserved data is freed.

#ifdef MCRL2_PLATFORM_WINDOWS
  #include <malloc.h>
  #define MCRL2_SPECIFIC_STACK_ALLOCATOR(TYPE, SIZE) reinterpret_cast<TYPE*>(_malloca((SIZE)*sizeof(TYPE)))
#elif (MCRL2_PLATFORM_LINUX || MCRL2_PLATFORM_MAC)
  #include <alloca.h>
  #define MCRL2_SPECIFIC_STACK_ALLOCATOR(TYPE, SIZE) reinterpret_cast<TYPE*>(alloca((SIZE)*sizeof(TYPE)))
#elif MCRL2_PLATFORM_FREEBSD
  #include <cstdlib>
 #define MCRL2_SPECIFIC_STACK_ALLOCATOR(TYPE, SIZE) reinterpret_cast<TYPE*>(alloca((SIZE)*sizeof(TYPE)))
#else
  static_assert(false, "MCRL2_SPECIFICATION_STACK_ALLOCATOR has not yet been defined for your platform.");
#endif

#endif // MCRL2_UTILITIES_DETAIL_MEMORY_UTILITY_H

