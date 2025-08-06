// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_IO_BINARY_H
#define MCRL2_ATERMPP_ATERM_IO_BINARY_H

#include <memory>

#include "mcrl2/atermpp/aterm_io.h"

#include "mcrl2/utilities/bitstream.h"
#include "mcrl2/atermpp/standard_containers/deque.h"
#include "mcrl2/atermpp/standard_containers/indexed_set.h"

namespace atermpp
{

/// \brief Writes terms in a streamable binary aterm format to an output stream.
/// \details The streamable aterm format:
///
///          Aterms (and function symbols) are written as packets (with an identifier in the header) and their
///          indices are derived from the number of aterms, resp. symbols, that occur before them in this stream. For each term
///          we first ensure that its arguments and symbol are written to the stream (avoiding duplicates). Then its
///          symbol index followed by a number of indices (depending on the arity) for its argments are written as integers.
///          Packet headers also contain a special value to indicate that the read term should be visible as output as opposed to
///          being only a subterm.
///          The start of the stream is a zero followed by a header and a version and a term with function symbol index zero
///          indicates the end of the stream.
///
class binary_aterm_ostream final : public aterm_ostream
{
public:
  /// \brief Provide the output stream to which the terms are written.
  binary_aterm_ostream(std::ostream& os);
  binary_aterm_ostream(std::shared_ptr<mcrl2::utilities::obitstream> stream);

  ~binary_aterm_ostream() override;

  /// \brief Writes an aterm in a compact binary format where subterms are shared. The term that is
  ///        written itself is not shared whenever it occurs as the argument of another term.
  void put(const aterm &term) override;

private:
  /// \brief Write a function symbol to the output stream.
  std::size_t write_function_symbol(const function_symbol& symbol);

  /// \returns The number of bits needed to index terms.
  unsigned int term_index_width();

  /// \returns The number of bits needed to index function symbols.
  unsigned int function_symbol_index_width();

  std::shared_ptr<mcrl2::utilities::obitstream> m_stream;

  unsigned int m_term_index_width = 0U;       ///< caches the result of term_index_width().
  unsigned int m_function_symbol_index_width; ///< caches the result of function_symbol_index_width().

  atermpp::indexed_set<aterm> m_terms; ///< An index of already written terms.
  mcrl2::utilities::indexed_set<function_symbol> m_function_symbols; ///< An index of already written function symbols.
};

/// \brief Reads terms from a stream in the steamable binary aterm format.
class binary_aterm_istream final : public aterm_istream
{
public:
  /// \brief Provide the input stream from which terms are read.
  binary_aterm_istream(std::istream& is);
  binary_aterm_istream(std::shared_ptr<mcrl2::utilities::ibitstream> stream);

  void get(aterm& t) override;

private:
  /// \returns The number of bits needed to index terms.
  unsigned int term_index_width();

  /// \returns The number of bits needed to index function symbols.
  unsigned int function_symbol_index_width();

  std::shared_ptr<mcrl2::utilities::ibitstream> m_stream;

  unsigned int m_term_index_width = 0U;       ///< caches the result of term_index_width().
  unsigned int m_function_symbol_index_width; ///< caches the result of function_symbol_index_width().

  atermpp::deque<aterm> m_terms; ///< An index of read terms.
  std::deque<function_symbol> m_function_symbols; ///< An index of read function symbols.
};

} // namespace atermpp

bool is_a_binary_aterm(std::istream& is);

#endif // MCRL2_ATERMPP_ATERM_IO_BINARY_H
