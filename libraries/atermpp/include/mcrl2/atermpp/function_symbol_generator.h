// Author(s): Wieger Wesselink, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/function_symbol_generator.h
/// \brief add your file description here.

#ifndef MCRL2_ATERMPP_FUNCTION_SYMBOL_GENERATOR_H
#define MCRL2_ATERMPP_FUNCTION_SYMBOL_GENERATOR_H

#include "mcrl2/atermpp/function_symbol.h"
#include "mcrl2/utilities/configuration.h"
#include "mcrl2/utilities/text_utility.h"

#include <mcrl3/shared_counter.h>

#include <atomic>
#include <cassert>
#include <cctype>
#include <string>

namespace atermpp
{

static inline std::atomic<std::size_t>& generator_sequence_number()
{
  static std::atomic<size_t> n = 0;
  return n;
}

/// \brief Generates unique function symbols with a given prefix.
class function_symbol_generator // : private mcrl2::utilities::noncopyable
{
protected:
  std::string m_prefix;

  /// \brief Cache the value that is set in the constructor
  std::size_t m_initial_index;

  /// \brief A local string cache to prevent allocating new strings for every function symbol generated.
  std::string m_string_buffer;

  /// \brief A reference to the index as present in the function symbol generator.
  std::size_t m_index;

  /// \brief The address of the central index for this prefix.
  mcrl3::shared_counter m_central_index;

public:
  /// \brief Constructor
  /// \param[in] prefix The prefix of the generated generated strings.
  /// \pre The prefix may not be empty, and it may not have trailing digits
  function_symbol_generator(const std::string& prefix)
  {
    std::size_t index = generator_sequence_number().fetch_add(1);

    m_prefix = prefix + (index > 0 ? std::to_string(index) + "_" : "");
    m_string_buffer = m_prefix;
    assert(!prefix.empty() && !(std::isdigit(*prefix.rbegin())));

    // Obtain a reference to the first index possible.
    m_central_index
        = mcrl3::shared_counter(mcrl3::ffi::function_symbol_register_prefix(m_prefix.c_str(), m_prefix.length()));
    m_index = *m_central_index;

    m_initial_index = m_index;
  }

  /// \brief Restores the index back to the value that was initially assigned in the constructor.

  void clear() { m_index = m_initial_index; }

  ~function_symbol_generator() { mcrl3::ffi::function_symbol_deregister_prefix(m_prefix.c_str(), m_prefix.length()); }

  /// \brief Generates a unique function symbol with the given prefix followed by a number.
  function_symbol operator()(std::size_t arity = 0)
  {
    // Check whether in the meantime other variables have been used with the same prefix.
    if (m_index <= *m_central_index)
    {
      m_index = *m_central_index;
      m_initial_index = m_index;
    }
    // Put the number m_index after the prefix in the string buffer.
    mcrl2::utilities::number2string(m_index, m_string_buffer, m_prefix.size());

    // Increase the index.
    ++m_index;

    // Generate a new function symbol with prefix + index.
    function_symbol f(m_string_buffer, arity, false);
    return f;
  }
};

} // namespace atermpp

#endif // MCRL2_ATERMPP_FUNCTION_SYMBOL_GENERATOR_H
