// Author(s): Frank Stappers
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/basename.h
/// \brief functions to get the executable's basename.

#ifndef MCRL2_BASENAME_H
#define MCRL2_BASENAME_H

#include "mcrl2/utilities/platform.h"

#include <cstdio>
#include <string>

#ifdef MCRL2_PLATFORM_LINUX
#include <unistd.h>
#endif

#ifdef MCRL2_PLATFORM_WINDOWS
#include <windows.h>
#endif

#ifdef MCRL2_PLATFORM_MAC
#include <mach-o/dyld.h>
#endif

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::utilities
{

    /// \brief Returns the basename of a tool.
    /// \return A string
    inline std::string get_executable_basename()
    {
      std::string path;
#ifdef MCRL2_PLATFORM_LINUX
      path = "";
      pid_t pid = getpid();
      char buf[10];
      sprintf(buf,"%d",pid);
      std::string _link = "/proc/";
      _link.append(buf);
      _link.append("/exe");
      char proc[512];
      int ch = readlink(_link.c_str(),proc,512);
      if (ch != -1)
      {
        proc[ch] = 0;
        path = proc;
        std::string::size_type t = path.find_last_of('/');
        path = path.substr(0,t);
      }
#endif // MCRL2_PLATFORM_LINUX

#ifdef MCRL2_PLATFORM_MAC
      char* pathbuf = NULL;
      uint32_t bufsize = 0;
      _NSGetExecutablePath(pathbuf, &bufsize);
      pathbuf = new char[bufsize];
      if (_NSGetExecutablePath(pathbuf, &bufsize) != 0)
      {
        throw mcrl2::runtime_error("Could not retrieve path to main executable (_NSGetExecutablePath returned nonzero).");
      }
      path = pathbuf;
      delete[] pathbuf;
      std::string::size_type t = path.find_last_of("/");
      path = path.substr(0,t);
#endif //MCRL2_PLATFORM_MAC

#ifdef MCRL2_PLATFORM_WINDOWS
      char buffer[MAX_PATH];  //always use MAX_PATH for filepaths
      GetModuleFileName(NULL, buffer, sizeof(buffer));
      path = buffer;
      std::string::size_type t = path.find_last_of("\\");
      path = path.substr(0,t);
#endif // MCRL2_PLATFORM_WINDOWS

      return path;
    }
}

#endif //MCRL2_BASENAME_H
