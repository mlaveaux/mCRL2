// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_SYMBOLIC_LDD_IOSTREAM_H
#define MCRL2_SYMBOLIC_LDD_IOSTREAM_H

#ifdef MCRL2_ENABLE_SYLVAN

#include "mcrl2/utilities/bitstream.h"

#include <sylvan_ldd.hpp>

#include <deque>
#include <iostream>
#include <memory>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/indexed_set.cxx"
#else
  import utilities;
#endif

namespace std
{

/// \brief specialization of the standard std::hash function.
template<>
struct hash<sylvan::ldds::ldd>
{
  std::size_t operator()(const sylvan::ldds::ldd& ldd) const
  {
    std::hash<std::uint64_t> hasher;
    return hasher(ldd.get());
  }
};

} // namespace std

namespace mcrl2::symbolic
{

/// \brief Writes ldds in a streamable binary format to an output stream.
/// \details The streamable ldd format:
///
///          Every LDD is traversed in order and assigned a unique number.
///          Whenever traversal encounters an LDD of which all children have
///          been visited it is written to the stream as 0:[value, down_index,
///          right_index]. An output LDD (as returned by
///          binary_ldd_istream::get()) is written as 1:index.
class binary_ldd_ostream
{
public:
  binary_ldd_ostream(std::ostream& os);
  binary_ldd_ostream(std::shared_ptr<mcrl2::utilities::obitstream> stream);

  ~binary_ldd_ostream();

  void put(const sylvan::ldds::ldd& U);

private:
  /// \returns The number of bits needed to index ldds.
  unsigned int ldd_index_width();
  
  std::shared_ptr<mcrl2::utilities::obitstream> m_stream;
  mcrl2::utilities::indexed_set<sylvan::ldds::ldd> m_nodes; ///< An index of already written ldds.

};

class binary_ldd_istream
{
public:
  binary_ldd_istream(std::istream& os);
  binary_ldd_istream(std::shared_ptr<mcrl2::utilities::ibitstream> stream);

  sylvan::ldds::ldd get();

private:
  /// \returns The number of bits needed to index ldds.
  unsigned int ldd_index_width(bool input = false);
  
  std::shared_ptr<mcrl2::utilities::ibitstream> m_stream;
  std::deque<sylvan::ldds::ldd> m_nodes; ///< An index of read ldds.
};

/// \brief Write the given term to the stream.
inline binary_ldd_ostream& operator<<(binary_ldd_ostream& stream, const sylvan::ldds::ldd& term) { stream.put(term); return stream; }

/// \brief Read the given term from the stream, but for aterm_list we want to use a specific one that performs validation (defined below).
inline binary_ldd_istream& operator>>(binary_ldd_istream& stream, sylvan::ldds::ldd& term) { term = stream.get(); return stream; }

} // namespace mcrl2::symbolic

#endif // MCRL2_ENABLE_SYLVAN

#endif // MCRL2_SYMBOLIC_LDD_IOSTREAM_H