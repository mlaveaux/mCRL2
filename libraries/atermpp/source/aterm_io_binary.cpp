// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/atermpp/aterm_io_binary.h"
#include "mcrl2/atermpp/standard_containers/stack.h"

namespace atermpp
{
using namespace mcrl2::utilities;

/// \brief The magic value for a binary aterm format stream.
/// \details As of version 0x8305 the magic and version are written as 2 bytes not encoded as variable-width integers. To ensure
///          compatibility with older formats the previously variable-width encoding is mimicked by prefixing them with 1000 (0x8).
static constexpr std::uint16_t BAF_MAGIC = 0x8baf;

/// \brief The BAF_VERSION constant is the version number of the ATerms written in BAF
///        format. As of 29 August 2013 this version number is used by the mCRL2
///        toolset. Whenever the file format of mCRL2 files is changed, the BAF_VERSION
///        has to be increased.
///
/// \details History:
///
/// before 2013       : version 0x0300
/// 29 August 2013    : version changed to 0x0301
/// 23 November 2013  : version changed to 0x0302 (introduction of index for variable types)
/// 24 September 2014 : version changed to 0x0303 (introduction of stochastic distribution)
///  2 April 2017     : version changed to 0x0304 (removed a few superfluous fields in the format)
/// 19 July 2019      : version changed to 0x8305 (introduction of the streamable aterm format)
/// 28 February 2020  : version changed to 0x8306 (added ability to stream aterm_int, implemented structured streaming for all objects)
/// 24 January 2023   : version changed to 0x8307 (removed NoIndex from Variables, Boolean variables. Made the .lts format more
///                                                compact by not storing states with a default probability 1)
/// 6  August 2024    : version changed to 0x8308 (introduced machine numbers)
static constexpr std::uint16_t BAF_VERSION = 0x8308;

/// \brief Each packet has a header consisting of a type.
/// \details Either indicates a function symbol, a term (either shared or output) or an arbitrary integer.
enum class packet_type
{
  function_symbol = 0,
  aterm,
  aterm_output,
  aterm_int_output
};

/// \brief The number of bits needed to store an element of packet_type.
static constexpr unsigned int packet_bits = 2;

binary_aterm_ostream::binary_aterm_ostream(std::shared_ptr<mcrl2::utilities::obitstream> stream)
  : m_stream(stream)
{
  // The term with function symbol index 0 indicates the end of the stream, its actual value does not matter.
  m_function_symbols.insert(function_symbol("end_of_stream", 0));
  m_function_symbol_index_width = 1;

  // Write the header of the binary aterm format.
  m_stream->write_bits(0, 8);
  m_stream->write_bits(BAF_MAGIC, 16);
  m_stream->write_bits(BAF_VERSION, 16);  
}

binary_aterm_ostream::binary_aterm_ostream(std::ostream& stream)
  : binary_aterm_ostream(std::make_shared<mcrl2::utilities::obitstream>(stream))
{}

binary_aterm_ostream::~binary_aterm_ostream()
{
  // Write the end of the stream.
  m_stream->write_bits(static_cast<std::size_t>(packet_type::aterm), packet_bits);
  m_stream->write_bits(0, function_symbol_index_width());
}

/// \brief Keep track of whether the term can be written to the stream.
struct write_todo
{
  unprotected_aterm term;
  bool write = false;

  write_todo(unprotected_aterm term)
   : term(term)
  {}

  void mark() const
  {
    term.mark();
  }
};

void binary_aterm_ostream::put(const aterm& term)
{
  // Traverse the term bottom up and store the subterms (and function symbol) before the actual term.
  atermpp::stack<write_todo> stack;
  stack.emplace(static_cast<const aterm&>(term));

  do
  {
    write_todo& current = stack.top();
    aterm transformed = m_transformer(static_cast<const aterm&>(static_cast<const aterm&>(current.term)));

    // Indicates that this term is output and not a subterm, these should always be written.
    bool is_output = stack.size() == 1;
    if (m_terms.index(down_cast<aterm>(current.term)) >= m_terms.size() || is_output)
    {
      if (current.write)
      {
        if (transformed.type_is_int())
        {
          if (is_output)
          {
            // If the integer is output, write the header and just an integer
            m_stream->write_bits(static_cast<std::size_t>(packet_type::aterm_int_output), packet_bits);
            m_stream->write_integer(static_cast<const aterm_int&>(static_cast<const aterm&>(current.term)).value());
          }
          else
          {
            std::size_t symbol_index = write_function_symbol(transformed.function());

            // Write the packet identifier of an aterm_int followed by its value.
            m_stream->write_bits(static_cast<std::size_t>(packet_type::aterm), packet_bits);
            m_stream->write_bits(symbol_index, function_symbol_index_width());
            m_stream->write_integer(static_cast<const aterm_int&>(static_cast<const aterm&>(current.term)).value());
          }
        }
        else
        {
          std::size_t symbol_index = write_function_symbol(transformed.function());

          // Write the packet identifier, followed by the indices of its function symbol and arguments.
          m_stream->write_bits(static_cast<std::size_t>(is_output ? packet_type::aterm_output : packet_type::aterm), packet_bits);
          m_stream->write_bits(symbol_index, function_symbol_index_width());

          for (const aterm& argument : transformed)
          {
            std::size_t index = m_terms.index(argument);
            assert(index < m_terms.size()); // Every argument must already be written.
            m_stream->write_bits(index, term_index_width());
          }
        }

        if (!is_output)
        {
          // Only regular terms (not output) are shared and as such need a unique index.
          [[maybe_unused]]
          bool assigned = m_terms.insert(down_cast<aterm>(current.term)).second;
          assert(assigned); // This term should have a new index assigned.
          m_term_index_width = static_cast<std::uint8_t>(std::log2(m_terms.size()) + 1);
        }

        stack.pop();
      }
      else
      {
        // Add all the arguments to the stack; to be processed first.
        for (const unprotected_aterm& argument : transformed)
        {
          const aterm& term = static_cast<const aterm&>(argument);
          if (m_terms.index(term) >= m_terms.size())
          {
            // Only add arguments that have not been written before.
            stack.emplace(term);
          }
        }

        current.write = true;
      }
    }
    else
    {
     stack.pop(); // This term was already written and as such should be skipped. This can happen if
                  // one term has two equal subterms.
    }
  }
  while (!stack.empty());
}

unsigned int binary_aterm_ostream::term_index_width()
{
  assert(m_term_index_width == static_cast<unsigned int>(std::log2(m_terms.size()) + 1));
  return m_term_index_width;
}

unsigned int binary_aterm_ostream::function_symbol_index_width()
{
  assert(m_function_symbol_index_width == static_cast<unsigned int>(std::log2(m_function_symbols.size()) + 1));
  return m_function_symbol_index_width;
}

binary_aterm_istream::binary_aterm_istream(std::shared_ptr<mcrl2::utilities::ibitstream> stream)
  : m_stream(stream)
{
  // The term with function symbol index 0 indicates the end of the stream.
  m_function_symbols.emplace_back();
  m_function_symbol_index_width = 1;

  // Read the binary aterm format header.
  if (m_stream->read_bits(8) != 0 || m_stream->read_bits(16) != BAF_MAGIC)
  {
    throw mcrl2::runtime_error("Error while reading: missing the BAF_MAGIC control sequence.");
  }

  std::size_t version = m_stream->read_bits(16);
  if (version != BAF_VERSION)
  {
    throw mcrl2::runtime_error("The BAF version (" + std::to_string(version) + ") of the input file is incompatible with the version (" + std::to_string(BAF_VERSION) +
                               ") of this tool. The input file must be regenerated. ");
  }
}

binary_aterm_istream::binary_aterm_istream(std::istream& is)
  : binary_aterm_istream(std::make_shared<mcrl2::utilities::ibitstream>(is))
{}

std::size_t binary_aterm_ostream::write_function_symbol(const function_symbol& symbol)
{
  std::size_t result = m_function_symbols.index(symbol);

  if (result < m_function_symbols.size())
  {
    return result;
  }
  else
  {
    // The function symbol has not been written yet, write it now and return its index.
    m_stream->write_bits(static_cast<std::size_t>(packet_type::function_symbol), packet_bits);
    m_stream->write_string(symbol.name());
    m_stream->write_integer(symbol.arity());

    auto result = m_function_symbols.insert(symbol);
    m_function_symbol_index_width = static_cast<unsigned int>(std::log2(m_function_symbols.size()) + 1);

    return result.first;
  }
}

void binary_aterm_istream::get(aterm& t)
{
  while(true)
  {
    // Determine the type of the next packet.
    std::size_t header = m_stream->read_bits(packet_bits);
    packet_type packet = static_cast<packet_type>(header);

    if (packet == packet_type::function_symbol)
    {
      // Read a single function symbol and insert it into the already read function symbols.
      std::string name = m_stream->read_string();
      std::size_t arity = m_stream->read_integer();
      m_function_symbols.emplace_back(name, arity);
      m_function_symbol_index_width = static_cast<unsigned int>(std::log2(m_function_symbols.size()) + 1);
    }
    else if (packet == packet_type::aterm_int_output)
    {
      // Read the integer from the stream and construct an aterm_int.
      std::size_t value = m_stream->read_integer();
      make_aterm_int(static_cast<aterm_int&>(t), value);
      return;
    }
    else if (packet == packet_type::aterm || packet == packet_type::aterm_output)
    {
      // First read the function symbol of the following term.
      function_symbol symbol = m_function_symbols[m_stream->read_bits(function_symbol_index_width())];

      if (!symbol.defined())
      {
        // The term with function symbol zero marks the end of the stream.
        t=aterm();
        return;
      }
      else if (symbol.type_is_int())
      {
        // Read the integer from the stream and construct a shared aterm_int in the index.
        std::size_t value = m_stream->read_integer();
        m_terms.emplace_back(aterm_int(value));
        m_term_index_width = static_cast<unsigned int>(std::log2(m_terms.size()) + 1);
      }
      else
      {
        // Read arity number of arguments from the stream and search them in the already defined set of terms.
        std::vector<aterm> arguments(symbol.arity());
        for (std::size_t argument = 0; argument < symbol.arity(); ++argument)
        {
          arguments[argument] = m_terms[m_stream->read_bits(term_index_width())];
        }

        // Transform the resulting term.
        aterm transformed = m_transformer(aterm(symbol, arguments.begin(), arguments.end()));

        if (packet == packet_type::aterm_output)
        {
          // This aterm was marked as output in the file so return it.
          t=transformed;
          return;
        }
        else
        {
          // Construct the term appl from the function symbol and the already read arguments and insert it.
          m_terms.emplace_back(transformed);
          m_term_index_width = static_cast<unsigned int>(std::log2(m_terms.size()) + 1);
        }
      }
    }
  }
}

unsigned int binary_aterm_istream::term_index_width()
{
  assert(m_term_index_width == static_cast<unsigned int>(std::log2(m_terms.size()) + 1));
  return m_term_index_width;
}

unsigned int binary_aterm_istream::function_symbol_index_width()
{
  assert(m_function_symbol_index_width == static_cast<unsigned int>(std::log2(m_function_symbols.size()) + 1));
  return m_function_symbol_index_width;
}

void write_term_to_binary_stream(const aterm& t, std::ostream& os)
{
  binary_aterm_ostream(os) << t;
}

void read_term_from_binary_stream(std::istream& is, aterm& t)
{
  binary_aterm_istream(is).get(t);
}

} // namespace atermpp
