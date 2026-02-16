// Author(s): anonymous
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/io.h
/// \brief IO routines for boolean equation systems.

#ifndef MCRL2_PBES_IO_H
#define MCRL2_PBES_IO_H

#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/extended_pbes.h"
#include "mcrl2/utilities/file_utility.h"

namespace mcrl2::pbes_system
{

void save_bes_pgsolver(const pbes& bes, std::ostream& stream, bool maxpg = true);

const std::vector<utilities::file_format>& pbes_file_formats();

inline
bool is_pbes_file_format(const utilities::file_format& format)
{
  for (const utilities::file_format& i: pbes_file_formats())
  {
    if (i == format)
    {
      return true;
    }
  }
  return false;
}

inline
const utilities::file_format& pbes_format_internal() { return pbes_file_formats()[0]; }
inline
const utilities::file_format& pbes_format_text() { return pbes_file_formats()[1]; }
inline
const utilities::file_format& pbes_format_internal_bes() { return pbes_file_formats()[2]; }
inline
const utilities::file_format& pbes_format_pgsolver() { return pbes_file_formats()[3]; }

inline
utilities::file_format guess_format(const std::string& filename)
{
  for (const utilities::file_format& it : pbes_file_formats())
  {
    if (it.matches(filename))
    {
      return it;
    }
  }
  return utilities::file_format();
}

/// \brief Save a PBES in the format specified.
/// \param pbes The PBES to be stored
/// \param stream The stream to which the output is saved.
/// \param format Determines the format in which the result is written. If unspecified, or
///        pbes_file_unknown is specified, then a default format is chosen.
void save_pbes(const pbes& pbes,
               std::ostream& stream,
               utilities::file_format format = utilities::file_format());

/// \brief Load a PBES from file.
/// \param pbes The PBES to which the result is loaded.
/// \param stream The stream from which to load the PBES.
/// \param format The format that should be assumed for the file in infilename. If unspecified, or
///        pbes_file_unknown is specified, then a default format is chosen.
/// \param source The source from which the stream originates. Used for error messages.
void load_pbes(pbes& pbes, std::istream& stream, utilities::file_format format, const std::string& source = "");

/// \brief save_pbes Saves a PBES to a file.
/// \param pbes The PBES to save.
/// \param filename The file to save the PBES in.
/// \param format The format in which to save the PBES.
/// \param welltypedness_check If set to false, skips checking whether pbes is well typed before
///                            saving it to file.
///
/// The format of the file in infilename is guessed if format is not given or if it is equal to
/// utilities::file_format().
void save_pbes(const pbes& pbes, const std::string& filename,
               utilities::file_format format = utilities::file_format(),
               bool welltypedness_check = true);

/// \brief Load pbes from file.
/// \param pbes The pbes to which the result is loaded.
/// \param filename The file from which to load the PBES.
/// \param format The format in which the PBES is stored in the file.
///
/// The format of the file in infilename is guessed if format is not given or if it is equal to
/// utilities::file_format().
void load_pbes(pbes& pbes,
               const std::string& filename,
               utilities::file_format format = utilities::file_format());

void save_extended_pbes(const pbes& pbes,
               std::ostream& stream,
               utilities::file_format format = utilities::file_format());

void load_extended_pbes(extended_pbes& pbes, std::istream& stream, utilities::file_format format, const std::string& source = "");

void save_extended_pbes(const extended_pbes& pbes, const std::string& filename,
               utilities::file_format format = utilities::file_format(),
               bool welltypedness_check = true);

void load_extended_pbes(extended_pbes& pbes,
               const std::string& filename,
               utilities::file_format format = utilities::file_format());

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_IO_H
