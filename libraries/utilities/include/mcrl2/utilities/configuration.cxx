// Copyright: see the accompanying file COPYING or copy at
// Author(s): Jan Friso Groote
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_CONFIGURATION_H
#define MCRL2_UTILITIES_CONFIGURATION_H
MCRL2_MODULE;

#ifdef MCRL2_ENABLE_MODULES
  export module utilities:configuration;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities::detail
{

/// Enables thread safety for the whole toolset.
  constexpr bool GlobalThreadSafe =     
#ifdef MCRL2_ENABLE_MULTITHREADING 
    true;
#else
    false;
#endif 
}

#endif // MCRL2_UTILITIES_CONFIGURATION_H