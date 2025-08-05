// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file logger.h

#ifndef MCRL2_UTILITIES_LOGGER_H
#define MCRL2_UTILITIES_LOGGER_H

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/logger.cxx"
#endif

/// \brief mCRL2log(LEVEL) provides the stream used to log.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define mCRL2log(LEVEL) if (mcrl2::log::mCRL2logEnabled(LEVEL)) mcrl2::log::logger(LEVEL).get()

#endif // MCRL2_UTILITIES_LOGGER_H
