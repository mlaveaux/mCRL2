// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_PBES_DETAIL_PBES_IO_H
#define MCRL2_PBES_DETAIL_PBES_IO_H

#include "mcrl2/pbes/io.h"
#include "mcrl2/pbes/extended_pbes.h"
#include "mcrl2/pbes/pbes.h"

#include <fstream>

namespace mcrl2::pbes_system
{

/// \brief Writes the pbes to a stream.
atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const pbes& pbes);

/// \brief Reads a pbes from a stream.
atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, pbes& pbes);

/// \brief Writes the pbes to a stream.
atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const extended_pbes& pbes);

/// \brief Reads a pbes from a stream.
atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, extended_pbes& pbes);

namespace detail
{

/// \brief Loads a PBES from filename, or from stdin if filename equals "".
pbes load_pbes(const std::string& filename);

/// \brief Saves an PBES to filename, or to stdout if filename equals "".
void save_pbes(const pbes& pbesspec, const std::string& filename);

} // namespace detail

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_DETAIL_PBES_IO_H
