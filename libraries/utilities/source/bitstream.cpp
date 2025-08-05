// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/utilities/bitstream.h"

#include <cassert>
#include <iostream>

#ifdef MCRL2_PLATFORM_WINDOWS
#include <io.h>
#include <fcntl.h>
#endif

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
  #include "mcrl2/utilities/power_of_two.cxx"
#else
  import utilities;
#endif

using namespace mcrl2::utilities;

/// \brief Encodes an unsigned variable-length integer using the most significant bit (MSB) algorithm.
///        This function assumes that the value is stored as little endian.
/// \param value The input value. Any standard integer type is allowed.
/// \param output A pointer to a piece of reserved memory. Must have a minimum size dependent on the input size (32 bit = 5 bytes, 64 bit = 10 bytes).
/// \returns The number of bytes used in the output.
/// \details Implementation taken from https://techoverflow.net/2013/01/25/efficiently-encoding-variable-length-integers-in-cc/
template<typename int_t = std::size_t>
static size_t encode_variablesize_int(int_t value, uint8_t* output)
{
  size_t outputSize = 0;

  // While more than 7 bits of data are left, occupy the last output byte
  // and set the next byte flag.
  while (value > 127)
  {
    // | 128: Sets the next byte flag.
    output[outputSize] = (static_cast<uint8_t>(value & 127)) | 128;

    // Remove the seven bits we just wrote from value.
    value >>= 7;
    outputSize++;
  }

  output[outputSize] = (static_cast<uint8_t>(value));
  return outputSize + 1;
}

/// \brief Decodes an unsigned variable-length integer using the MSB algorithm.
/// \param stream The stream from which the bytes for this value are read.
/// \details Implementation taken from https://techoverflow.net/2013/01/25/efficiently-encoding-variable-length-integers-in-cc/, 
///          but converted to a streaming method.
template<typename int_t = std::size_t>
int_t decode_variablesize_int(ibitstream& stream)
{
  int_t value = 0;
  for (size_t i = 0; i < integer_encoding_size<int_t>(); i++)
  {
    // Read the next byte from the stream.
    std::size_t byte = stream.read_bits(8);

    // Take 7 bits (mask 0x01111111) from byte and shift it before the bits already written to value.
    value |= (static_cast<int_t>(byte) & 127) << (7 * i);

    if (!(byte & 128))
    {
      // If the next-byte flag is not set then we are finished.
      break;
    }
    else if (i >= integer_encoding_size<int_t>() - 1)
    {
      // The next-byte flag was set, but we cannot represent it using int_t.
      throw std::runtime_error("Fail to read an int from the input");
    }
  }

  return value;
}

/// \brief Change the current stream to binary mode (no handle of newline characters),
static void set_stream_binary([[maybe_unused]] const std::string& name, [[maybe_unused]] FILE* handle)
{
#ifdef MCRL2_PLATFORM_WINDOWS
  if (_setmode(_fileno(handle), _O_BINARY) == -1)
  {
    mCRL2log(mcrl2::log::warning) << "Cannot set " << name << " to binary mode.\n";
  }
  else
  {
    mCRL2log(mcrl2::log::debug) << "Converted " << name << " to binary mode.\n";
  }

  fflush(stderr);
#endif // MCRL2_PLATFORM_WINDOWS
}

obitstream::obitstream(std::ostream& stream)
  : stream(stream)
{
  // Ensures that the given stream is changed to binary mode.
  if (stream.rdbuf() == std::cout.rdbuf())
  {
    set_stream_binary("cout", stdout);
  }
  else if (stream.rdbuf() == std::cerr.rdbuf())
  {
    set_stream_binary("cerr", stderr);
  }
}

void obitstream::write_bits(std::size_t value, unsigned int number_of_bits)
{
  // Add val to the buffer by masking out additional bits and put them at left-most position free in the buffer.
  assert(number_of_bits <= std::numeric_limits<std::size_t>::digits);
  if (number_of_bits < std::numeric_limits<std::size_t>::digits) // the check is needed to avoid undefined behavior
  {
    value &= (static_cast<std::size_t>(1) << number_of_bits) - 1;
  }

  write_buffer |= std::bitset<128>(value) << ((128 - bits_in_buffer) - number_of_bits);
  bits_in_buffer += number_of_bits;

  // Write 8 bytes if available
  if (bits_in_buffer >= 64)
  {
    unsigned long long write_value = (write_buffer >> 64).to_ullong();
    write_buffer <<= 64;
    bits_in_buffer -= 64;

    for (int32_t i = 7; i >= 0; --i)
    {
      // Write the 8 * i most significant bits and mask out the other values.
      stream.put(static_cast<char>((write_value >> (8 * i)) & 255));

      if (stream.fail())
      {
        throw mcrl2::runtime_error("Failed to write bytes to the output file/stream.");
      }
    }
  }
}

void obitstream::write_string(const std::string& string)
{
  // Write length.
  write_integer(string.size());

  // Write actual string.
  write(reinterpret_cast<const std::uint8_t*>(string.c_str()), string.size());
}

void obitstream::write_integer(std::size_t val)
{
  std::size_t nr_bytes = encode_variablesize_int(val, integer_buffer);

  write(integer_buffer, nr_bytes);
}

ibitstream::ibitstream(std::istream& stream)
  : stream(stream)
{
  // Ensures that the given stream is changed to binary mode.
  if (stream.rdbuf() == std::cin.rdbuf())
  {
    set_stream_binary("cin", stdin);
  }
}

const char* ibitstream::read_string()
{
  std::size_t length;

  // Get length of string.
  length = read_integer();

  // Assure buffer can hold the string.
  if (m_text_buffer.size() < (length + 1))
  {
    m_text_buffer.resize(round_up_to_power_of_two(length + 1));
  }

  // Read the actual string.
  read(length, reinterpret_cast<std::uint8_t*>(m_text_buffer.data()));
  m_text_buffer[length] = '\0';

  return m_text_buffer.data();
}

std::size_t ibitstream::read_bits(unsigned int number_of_bits)
{
  // Read at most the number of bits of a std::size_t.
  assert(number_of_bits <= std::numeric_limits<std::size_t>::digits);

  while (bits_in_buffer < number_of_bits)
  {
    // Read bytes until the buffer is sufficiently full.
    int byte = stream.get();

    if (stream.eof())
    {
      throw mcrl2::runtime_error("Unexpected end-of-file reached in the input file/stream.");
    }
    else if (stream.fail())
    {
      throw mcrl2::runtime_error("Failed to read bytes from the input file/stream.");
    }

    // Shift the 8 bits to the first free (120 - bits_in_buffer) position in the buffer.
    read_buffer |= std::bitset<128>(static_cast<std::size_t>(byte)) << (56 + 64 - bits_in_buffer);
    bits_in_buffer += 8;
  }

  // Read nr_bits from the buffer by shifting them to the least significant bits and masking out the remaining bits.
  std::size_t value = (read_buffer >> (128 - number_of_bits)).to_ullong();

  // Shift the first bit to the first position in the buffer.
  read_buffer <<= number_of_bits;
  bits_in_buffer -= number_of_bits;

  return value;
}

std::size_t ibitstream::read_integer()
{
  return decode_variablesize_int(*this);
}

// Private functions

void obitstream::flush()
{
  // Writing the buffer full to 64 bits should flush it internally, this also guarantees that the unnecessary bits are zeroed out.
  write_bits(0, 64 - bits_in_buffer);
  assert(bits_in_buffer == 0);

  stream.flush();
  if (stream.fail())
  {
    throw mcrl2::runtime_error("Failed to write the last byte to the output file/stream.");
  }
}

void obitstream::write(const uint8_t* buffer, std::size_t size)
{
  for (std::size_t index = 0; index < size; ++index)
  {
    // Write a single byte for every entry in the buffer that was filled (size).
    write_bits(buffer[index], 8);
  }
}

void ibitstream::read(std::size_t size, std::uint8_t* buffer)
{
  for (std::size_t index = 0; index < size; ++index)
  {
    // Read a single byte for every entry into the buffer that was filled (size).
    std::size_t value = read_bits(8);
    buffer[index] = static_cast<std::uint8_t>(value);
  }
}
