// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <boost/json.hpp>

namespace json = boost::json;

namespace mcrl2::mbt_protocol
{

using message_id = std::uint64_t;

/// \brief Configuration negotiated during hello handshake.
struct session_config
{
  std::size_t tau_closure_depth = 1000;
  std::uint64_t early_output_timeout_ms = 1000;
  std::uint64_t heartbeat_interval_ms = 5000;
  std::uint64_t heartbeat_timeout_ms = 15000;
};

enum class error_code
{
  malformed_message,
  unknown_type,
  protocol_mismatch,
  unsupported_config,
  not_ready,
  input_not_enabled,
  output_unprocessed,
  quiescence_unexpected,
  lps_unavailable,
  internal_error
};

std::string_view to_string(error_code code);

/// \brief First message sent by the tool after connecting.
json::object make_hello(const std::string& lps_identifier, const std::string& lps_hash);

/// \brief Keepalive ping.
json::object make_heartbeat();

/// \brief Positive acknowledgement for input, output, or quiescence.
/// \param kind  One of "input", "output", "quiescence".
json::object make_ack(message_id in_reply_to, std::string_view kind);

/// \brief Non-fatal warning (e.g. early output received).
json::object make_warning(message_id in_reply_to, std::string_view code, std::string_view detail = "");

/// \brief Fatal error with a message id context.
json::object make_error(message_id in_reply_to, error_code code, std::string_view detail = "");

/// \brief Fatal error without a message id (used during handshake).
json::object make_error_no_id(error_code code, std::string_view detail = "");

/// \brief Reply to get_enabled.
json::object make_enabled(message_id in_reply_to,
  const std::vector<std::string>& inputs,
  const std::vector<std::string>& outputs,
  bool quiescence);

/// \brief Reply to reset.
json::object make_reset_ack(message_id in_reply_to);

/// \brief Initiate graceful close.
json::object make_close(std::string_view reason = "");

// ─── Incoming (adapter → tool) ────────────────────────────────────────────────

struct incoming_message
{
  std::string type;
  std::optional<message_id> id;
  std::optional<message_id> in_reply_to;
  json::object payload; ///< Full parsed object for field access.
};

/// \brief Parse a text WebSocket frame into an incoming_message.
/// \throws std::runtime_error on malformed JSON or missing "type" field.
incoming_message parse_incoming(std::string_view text);

/// \brief Parse session_config from the "config" sub-object of a hello response.
/// Missing fields fall back to defaults.
session_config parse_session_config(const json::object& config_obj);

} // namespace mcrl2::mbt_protocol
