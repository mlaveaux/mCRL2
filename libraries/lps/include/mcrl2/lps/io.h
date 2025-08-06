// Author(s): Sjoerd Cranen; simplified by Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/io.h
/// \brief IO routines for linear process specifications.

#ifndef MCRL2_LPS_IO_H
#define MCRL2_LPS_IO_H

#include <iostream>
#include <fstream>

#include "mcrl2/atermpp/aterm_io_binary.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/stochastic_specification.h"


namespace mcrl2::lps
{

/// \brief Writes LPS to the stream.
atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const specification& spec);
atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const stochastic_specification& spec);

/// \brief Reads LPS from the stream.
atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, specification& spec);
atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, stochastic_specification& spec);

/// \brief Save an LPS in the format specified.
/// \param spec The LPS to be stored
/// \param stream The stream to which the output is written.
/// \param target A string indicating the target stream, used to print information.
///        utilities::file_format::unknown() is specified, then a default format is chosen.
///
template <typename Specification>
void save_lps(const Specification& spec, std::ostream& stream, const std::string& target = "")
{
  mCRL2log(log::debug) << "Saving LPS" << (target.empty() ? "" : " to " + target) << ".\n";
  atermpp::binary_aterm_ostream(stream) << spec;
}

/// \brief Load LPS from file.
/// \param spec The LPS to which the result is loaded.
/// \param stream The stream from which to load the LPS (it is assumed to have been opened in the
///               right mode).
/// \param source The source from which the stream originates. Used for error messages.
template <typename Specification>
void load_lps(Specification& spec, std::istream& stream, const std::string& source = "")
{
  mCRL2log(log::debug) << "Loading LPS" << (source.empty() ? "" : " from " + source) << ".\n";

  if constexpr (std::is_same_v<Specification, specification>)
  {
    stochastic_specification stoch_spec;
    atermpp::binary_aterm_istream(stream) >> stoch_spec;
    spec=remove_stochastic_operators(stoch_spec);
  }
  else
  {
    atermpp::binary_aterm_istream(stream) >> spec;
  }
}

/// \brief Saves an LPS to a file.
/// \param spec The LPS to save.
/// \param filename The file to save the LPS in. If empty, the file is written to stdout.
template <typename Specification>
void save_lps(const Specification& spec, const std::string& filename)
{
  if (filename.empty() || filename == "-")
  {
    save_lps(spec, std::cout, "standard output");
    return;
  }

  std::ofstream ofs(filename, std::ios_base::binary);
  if (!ofs.good())
  {
    throw mcrl2::runtime_error("Could not open file " + filename + ".");
  }
  save_lps(spec, ofs, filename);
}

/// \brief Load LPS from file.
/// \param spec The LPS to which the result is loaded.
/// \param filename The file from which to load the LPS. If empty, the file is read from stdin.

template <typename Specification>
void load_lps(Specification& spec, const std::string& filename)
{
  if (filename.empty() || filename == "-")
  {
    load_lps(spec, std::cin, "standard input");
    return;
  }

  std::ifstream ifs(filename, std::ios_base::binary);
  if (!ifs.good())
  {
    throw mcrl2::runtime_error("Could not open file " + filename + ".");
  }
  load_lps(spec, ifs, filename);
}

} // namespace mcrl2::lps

#endif // MCRL2_LPS_IO_H
