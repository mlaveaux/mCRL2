// Author(s): Jeroen Keiren, Sjoerd Cranen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/dynamiclibrary.h


#ifndef MCRL2_UTILITIES_DYNAMIC_LIBRARY_H
#define MCRL2_UTILITIES_DYNAMIC_LIBRARY_H

#include "mcrl2/utilities/logger.h"
#include "mcrl2/utilities/platform.h"

#ifdef MCRL2_PLATFORM_WINDOWS
  #include <windows.h>
  typedef HMODULE library_handle;
  typedef FARPROC library_proc;

  inline
  library_handle get_module_handle(const std::string& fname)
  {
    return LoadLibrary(fname.c_str());
  }

  inline
  library_proc get_proc_address(library_handle handle, const std::string& procname)
  {
    return GetProcAddress(handle, procname.c_str());
  }

  inline
  bool close_module_handle(library_handle handle)
  {
    return FreeLibrary(handle);
  }

  inline
  std::string get_last_error()
  {
    std::string result;
    char *buffer = 0;
    DWORD last = GetLastError();
    DWORD x = FormatMessage(
      FORMAT_MESSAGE_ALLOCATE_BUFFER | // Let Windows allocate string
      FORMAT_MESSAGE_FROM_SYSTEM,      // Retrieve a system error message
        nullptr,               // No source needed
        last,               // Requested error message
        0,             // Default language
        reinterpret_cast<LPSTR>(&buffer), // Output buffer
        0,             // Minimum allocation size
        nullptr             // No arguments here
      );
    if (x)
    {
      result = buffer;
    }
    else
    {
      result = "Unknown error.";
    }
    LocalFree(buffer);
    return result;
  }
#else
  #include <dlfcn.h>
using library_handle = void*;
using library_proc = void*;

inline library_handle get_module_handle(const std::string& fname)
{
  return dlopen(fname.c_str(), RTLD_NOW);
  }

  inline
  library_proc get_proc_address(library_handle handle, const std::string& procname)
  {
    return dlsym(handle, procname.c_str());
  }

  inline
  bool close_module_handle(library_handle handle)
  {
    return (dlclose(handle) == 0);
  }

  inline
  std::string get_last_error()
  {
    return std::string(dlerror());
  }
#endif

class dynamic_library 
{
  protected:
    library_handle m_library = nullptr;
    std::string m_filename;

    void load() 
    {
      if (m_library == nullptr)
      {
        m_library = get_module_handle(m_filename.c_str());
        if (m_library == nullptr)
        {
          throw std::runtime_error("Could not load library (" + m_filename + "): " + get_last_error());
        }
      }
    }
  
    void unload() 
    {
      if (m_library!=nullptr)
      {
        if (!close_module_handle(m_library))
        {
          throw std::runtime_error("Could not close library (" + m_filename + "): " + get_last_error());
        }
        m_library = nullptr;
      }
    }
  
  public:
    dynamic_library(const std::string& filename = std::string())
        : m_filename(filename)
    {
    }

    virtual ~dynamic_library() 
    {
      try
      {
        unload();
      }
      catch(std::runtime_error& error)
      {
        mCRL2log(mcrl2::log::error) << "Error while unloading dynamic library: " << error.what() << std::endl;
      }
    }
  
    library_proc proc_address(const std::string& name) 
    {
      if (m_library == nullptr)
      {
        load();
      }
      library_proc result = get_proc_address(m_library, name.c_str());
      if (result == nullptr)
      {
        throw std::runtime_error("Could not find proc address (" + m_filename + ":" + name + "): " + get_last_error());
      }
      return result;
    }
};

#endif // MCRL2_UTILITIES_DYNAMIC_LIBRARY_H
