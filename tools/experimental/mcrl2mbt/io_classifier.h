// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <istream>
#include <set>
#include <string>

namespace mcrl2
{

/// \brief Classifies action names as input, output, or internal.
///
/// The classification file uses the following format:
///   input
///     action_name_1
///     action_name_2
///   output
///     action_name_3
///
/// Any action not listed is treated as internal.
class io_classifier
{
public:
  enum class action_type
  {
    input,
    output,
    internal
  };

  /// \brief Load classifications from a file.
  /// \throws std::runtime_error on I/O or parse error.
  void load(const std::string& filename);

  /// \brief Read classifications from an input stream.
  /// \throws std::runtime_error on I/O or parse error.
  void read(std::istream& input);

  /// \brief Classify an action by its label name.
  action_type classify(const std::string& action_name) const;

  bool is_input(const std::string& action_name) const { return classify(action_name) == action_type::input; }

  bool is_output(const std::string& action_name) const { return classify(action_name) == action_type::output; }

private:
  std::set<std::string> m_inputs;
  std::set<std::string> m_outputs;
};

} // namespace mcrl2
