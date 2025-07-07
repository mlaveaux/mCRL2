// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_BITSTREAM_H
#define MCRL2_UTILITIES_BITSTREAM_H

#include <bitset>
#include <cstdint>
#include <vector>


namespace mcrl2::utilities
{

/// \returns The number of bits needed to represent a value of type T in most significant bit encoding.
template<typename T>
constexpr std::size_t integer_encoding_size()
{
  return ((sizeof(T) + 1) * 8) / 7;
}

/// \brief A bitstream provides per bit writing of data to any stream (including stdout).
/// \details Internally uses bitpacking and buffering for compact and efficient IO.
class obitstream
{
public:
  /// \brief Provides the stream on which the write function operate.
  obitstream(std::ostream& stream);
  ~obitstream() { flush(); }

  /// \brief Write the num_of_bits least significant bits in descending order from value.
  /// @param value Variable that contains the bits.
  /// @param num_of_bits Number of bits to write to the output stream.
  void write_bits(std::size_t value, unsigned int num_of_bits);

  /// \brief Write the given string to the output stream.
  /// \details Encoded in bits using <length, string>
  void write_string(std::string_view string);

  /// \brief Write the given value to the output stream.
  /// \details Uses most significant bit encoding.
  void write_integer(std::size_t value);

private:
  /// \brief Flush the remaining bits in the buffer to the output stream.
  /// \details Note that this aligns it to the next byte, e.g. when bits_in_buffer is 6 then two zero bits are added redundantly.
  void flush();

  /// \brief Writes size bytes from the given buffer.
  void write(const std::uint8_t* buffer, std::size_t size);

  std::ostream& stream;

  /// \brief Buffer that is filled starting from bit 127 when writing
  std::bitset<128> write_buffer = 0;

  unsigned int bits_in_buffer = 0; ///< how many bits in are used in the buffer.

  std::uint8_t integer_buffer[integer_encoding_size<std::size_t>()]; ///< Reserved space to store an n byte integer.
};

/// \brief The counterpart of obitstream, guarantees that the same data is read as has been written when calling the read operators
///        in the same sequence as the corresponding write operators.
class ibitstream
{
public:
  /// \brief Provides the stream on which the read function operate.
  ibitstream(std::istream& stream);

  /// \brief Reads an num_of_bits bits from the input stream and stores them in the least significant part (in descending order) of the return value.
  /// \param num_of_bits Number of bits to read from the input stream.
  std::size_t read_bits(unsigned int num_of_bits);

  /// \returns A pointer to the read string.
  /// \details Remains valid until the next call to read_string.
  const char* read_string();

  /// \returns A natural number that was read from the binary stream encoded in most significant bit encoding.
  std::size_t read_integer();

private:
  /// \brief Read size bytes into the provided buffer.
  void read(std::size_t size, std::uint8_t* buffer);

  std::istream& stream;

  /// \brief Buffer that is filled starting from bit 127 when reading.
  std::bitset<128> read_buffer = 0;

  unsigned int bits_in_buffer = 0; ///< how many bits in the buffer are used.

  std::vector<char> m_text_buffer; ///< A temporary buffer to store char array strings.
};

} // namespace mcrl2::utilities


#endif // MCRL2_UTILITIES_BITSTREAM_H
