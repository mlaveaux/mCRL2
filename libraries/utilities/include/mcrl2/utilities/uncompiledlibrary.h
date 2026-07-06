// Author(s): Jeroen Keiren, Sjoerd Cranen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/uncompiledlibrary.h
/// \brief Extends the dynamic_library from dynamiclibrary.h to be able to compile a
///        source file and load the resulting library.

#ifndef MCRL2_UTILITY_UNCOMPILED_LIBRARY_H
#define MCRL2_UTILITY_UNCOMPILED_LIBRARY_H

/*
 * Extends the dynamic_library from dynamiclibrary.h to be able to compile a
 * source file and load the resulting library.
 *
 * Usage:
 *
 *   uncompiled_library mylib;
 *   mylib.compile(source_filename);
 *   myfunc = mylib.proc_address("myfunc");
 *   myfunc(10);
 *
 * Remarks:
 *
 * The source is compiled using a script that must take two string arguments.
 * The first argument is the source file, the second is the destination file.
 * After (successful) termination, only the source and destination files must
 * remain on disk -- it is the responsibility of the script to remove any
 * temporary files.
 *
 */

#include "mcrl2/utilities/platform.h"

#ifndef MCRL2_PLATFORM_WINDOWS

#include <array>
#include <cerrno>
#include <list>
#include "mcrl2/utilities/dynamiclibrary.h"
#include "mcrl2/utilities/file_utility.h"

class uncompiled_library : public dynamic_library
{
  protected:
    std::list<std::string> m_tempfiles;
    std::string m_compile_script;

  public:
    uncompiled_library() = delete;

    uncompiled_library(const std::string& script)
     : m_compile_script(script) 
    {}

    void compile(const std::string& filename) 
    {
      std::stringstream commandline;
      commandline << '"' << m_compile_script << "\" " << filename << " " << " 2>&1";
      
      // Execute script.
      FILE* stream = popen(commandline.str().c_str(), "r");
      if (stream == nullptr)
      {
        throw std::runtime_error("Could not execute compile script.");
      }

      // Script produces one file per line. Last file is the shared library,
      // preceding files are temporary files that should be removed when the
      // library is unloaded.
      std::array<char, 1024> buf;
      while(!feof(stream))
      {
        if (fgets(buf.data(), buf.size(), stream) != nullptr)
        {
          std::string line(buf.data());
          assert(*line.rbegin() == '\n');
          line.erase(line.size() - 1);
          mCRL2log(mcrl2::log::debug) << "  Read line: " << line << std::endl;
          
          // Check that reported file exists. If not, produce error message and
          // flush script output to the log.
          if (!mcrl2::utilities::file_exists(line))
          {
            mCRL2log(mcrl2::log::error) << "Compile script " << m_compile_script << " produced unexpected output:\n";
            mCRL2log(mcrl2::log::error) << line << std::endl;
            while (fgets(buf.data(), buf.size(), stream) != nullptr)
            {
              mCRL2log(mcrl2::log::error) << std::string(buf.data());
            }
            pclose(stream);
            throw std::runtime_error("Compile script failed.");
          }
          m_tempfiles.push_back(line);
        }
        else if(ferror(stream) && errno == EINTR)
        {
            // On OSX, interrupts sometimes arrive during the call to read(), which
            // is called by fgets. If an interrupt arrives, we just ignore it
            // an clear the error status of the stream, and try again.
            mCRL2log(mcrl2::log::debug) << "Reading was interrupted. Clearing error status and retrying" << std::endl;
            perror("Error according to errno");
            clearerr(stream);
        }
      }
      
      if (ferror(stream))
      {
        pclose(stream);
        throw std::runtime_error("There was a problem reading the output of the compile script.");
      }
      
      pclose(stream);

      m_filename = m_tempfiles.back();
    }

    void leave_files()
    {
      m_tempfiles.clear();
    }

    void cleanup() 
    {
      for (auto& tempfile: m_tempfiles)
      {
        if (remove(tempfile.c_str()))
        {
          std::stringstream s;
          s << "Could not remove file: " << tempfile;
          throw std::runtime_error(s.str());
        }
        else
        {
          mCRL2log(mcrl2::log::debug) << "Temporary file '" << tempfile << "' deleted." << std::endl;
        }
      }
    }

    ~uncompiled_library() override
    {
#ifndef NDEBUG // In debug mode, the compiled rewriter has not been removed directly after loading, 
               // and we still have to remove it.
      try
      {
        cleanup();
      }
      catch (std::runtime_error &error)
      {
        mCRL2log(mcrl2::log::error) << "Could not cleanup temporary files: " << error.what() << std::endl;
      }
#endif
    }

};

#endif // MCRL2_PLATFORM_WINDOWS

#endif // MCRL2_UTILITY_UNCOMPILED_LIBRARY_H
