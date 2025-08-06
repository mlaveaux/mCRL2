// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_IO_H
#define MCRL2_ATERMPP_ATERM_IO_H

#include <cassert>
#include <string>

#include "mcrl2/atermpp/aterm_int.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/concepts.cxx"
  #include "mcrl2/atermpp/aterm_list.cxx"
#else
  import atermpp;
#endif

namespace atermpp
{

/// \brief A function that is applied to all terms. The resulting term should only use
///        a subset of the original arguments (i.e. not introduce new terms).
/// \details Typical usage is removing the index traits from function symbols that represent operators.
using aterm_transformer = aterm(const aterm&);

/// \brief The default transformer that maps each term to itself.
inline aterm identity(const aterm& x)
{
  return x;
}

/// \brief The general aterm stream interface, which enables the use of a transformer to
///        change the written/read terms.
class aterm_stream
{
public:
  virtual ~aterm_stream();

  /// \brief Sets the given transformer to be applied to following writes.
  void set_transformer(aterm_transformer transformer) { m_transformer = transformer; }

  /// \returns The currently assigned transformer function.
  aterm_transformer* get_transformer() const { return m_transformer; }

protected:
  aterm_transformer* m_transformer = identity;
};

/// \brief The interface for a class that writes aterm to a stream.
///        Every written term is retrieved by the corresponding aterm_istream::get() call.
class aterm_ostream : public aterm_stream
{
public:
  ~aterm_ostream() override;

  /// \brief Write the given term to the stream.
  virtual void put(const aterm& term) = 0;
};

/// \brief The interface for a class that reads aterm from a stream.
///        The default constructed term aterm() indicates the end of the stream.
class aterm_istream : public aterm_stream
{
public:
  ~aterm_istream() override;

  /// \brief Reads an aterm from this stream.
  virtual void get(aterm& t) = 0;
};

// These free functions provide input/output operators for these streams.

/// \brief Sets the given transformer to be applied to following reads.
inline aterm_istream& operator>>(aterm_istream& stream, aterm_transformer transformer)
{
  stream.set_transformer(transformer);
  return stream;
}
inline aterm_ostream& operator<<(aterm_ostream& stream, aterm_transformer transformer)
{
  stream.set_transformer(transformer);
  return stream;
}

/// \brief Write the given term to the stream.
inline aterm_ostream& operator<<(aterm_ostream& stream, const aterm& term)
{
  stream.put(term);
  return stream;
}

/// \brief Read the given term from the stream, but for aterm_list we want to use a specific one that performs
/// validation (defined below).
inline aterm_istream& operator>>(aterm_istream& stream, aterm& term)
{
  stream.get(term);
  return stream;
}

// Utility functions

/// \brief A helper class to restore the state of the aterm_{i,o}stream objects upon destruction. Currently, onlt
///        preserves the transformer object.
class aterm_stream_state
{
public:
  aterm_stream_state(aterm_stream& stream)
      : m_stream(stream)
  {
    m_transformer = stream.get_transformer();
  }

  ~aterm_stream_state() { m_stream.set_transformer(m_transformer); }

private:
  aterm_stream& m_stream;
  aterm_transformer* m_transformer;
};

/// \brief Write any container (that is not an aterm itself) to the stream.
template <typename T>
  requires std::ranges::range<T> && (!IsATerm<T>)
inline aterm_ostream& operator<<(aterm_ostream& stream, const T& container)
{
  // Write the number of elements, followed by each element in the container.
  stream << aterm_int(std::distance(container.begin(), container.end()));

  for (const auto& element : container)
  {
    stream << element;
  }

  return stream;
}

/// \brief Read any container (that is not an aterm itself) from the stream.
template <typename T>
  requires std::ranges::range<T> && (!IsATerm<T>)
inline aterm_istream& operator>>(aterm_istream& stream, T& container)
{
  // Insert the next nof_elements into the container.
  aterm_int nof_elements;
  stream >> nof_elements;

  auto it = std::inserter(container, container.end());
  for (std::size_t i = 0; i < nof_elements.value(); ++i)
  {
    typename T::value_type element;
    stream >> element;
    it = element;
  }

  return stream;
}

// r-value forwarding to the streaming operators, used for convenience.
template<typename T>
inline aterm_ostream&& operator<<(aterm_ostream&& stream, const T& t)
{  
  static_cast<aterm_ostream&>(stream) << t;
  return std::move(stream);
}

template<typename T>
inline aterm_istream&& operator>>(aterm_istream&& stream, T& t)
{
  static_cast<aterm_istream&>(stream) >> t;
  return std::move(stream);
}

/// \brief Sends the name of a function symbol to an ostream.
/// \param out The out stream.
/// \param f The function symbol to be output.
/// \return The stream.
inline std::ostream& operator<<(std::ostream& out, const function_symbol& f)
{
  return out << f.name();
}

/// \brief Prints the name of a function symbol as a string.
/// \param f The function symbol.
/// \return The string representation of r.
inline const std::string& pp(const function_symbol& f)
{
  return f.name();
}

/// \brief Writes term t to a stream in binary aterm format.
void write_term_to_binary_stream(const aterm& t, std::ostream& os);

/// \brief Reads a term from a stream in binary aterm format.
void read_term_from_binary_stream(std::istream& is, aterm& t);

/// \brief Writes term t to a stream in textual format.
void write_term_to_text_stream(const aterm& t, std::ostream& os);

/// \brief Reads a term from a stream which contains the term in textual format.
void read_term_from_text_stream(std::istream& is, aterm& t);

/// \brief Reads an aterm from a string. The string can be in either binary or text format.
aterm read_term_from_string(const std::string& s);

/// \brief Reads an aterm_list from a string. The string can be in either binary or text format.
/// \details If the input is not a string, an aterm is returned of the wrong type.
/// \return The term corresponding to the string.
inline aterm_list read_list_from_string(const std::string& s)
{
  const aterm_list l = down_cast<aterm_list>(read_term_from_string(s));
  assert(l.type_is_list());
  return l;
}

/// \brief Reads an aterm_int from a string. The string can be in either binary or text format.
/// \details If the input is not an int, an aterm is returned of the wrong type.
/// \return The aterm_int corresponding to the string.
inline aterm_int read_int_from_string(const std::string& s)
{
  const aterm_int n = down_cast<aterm_int>(read_term_from_string(s));
  assert(n.type_is_int());
  return n;
}

/// \brief Reads an aterm from a string. The string can be in either binary or text format.
/// \details If the input is not an aterm, an aterm is returned of the wrong type.
/// \return The term corresponding to the string.
inline aterm read_appl_from_string(const std::string& s)
{
  const aterm a = read_term_from_string(s);
  assert(a.type_is_appl());
  return a;
}

} // namespace atermpp

#endif // MCRL2_ATERMPP_ATERM_IO_H
