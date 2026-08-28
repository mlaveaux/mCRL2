// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mbt_protocol.h"

#include <atomic>
#include <stdexcept>

namespace mcrl2::mbt_protocol
{

// \brief Message counter used for unique message identifiers.
static std::atomic<message_id> s_next_id{1};

static message_id next_id()
{
  return s_next_id.fetch_add(1, std::memory_order_relaxed);
}

std::string_view to_string(error_code code)
{
  switch (code)
  {
  case error_code::malformed_message:
    return "malformed_message";
  case error_code::unknown_type:
    return "unknown_type";
  case error_code::protocol_mismatch:
    return "protocol_mismatch";
  case error_code::unsupported_config:
    return "unsupported_config";
  case error_code::not_ready:
    return "not_ready";
  case error_code::input_not_enabled:
    return "input_not_enabled";
  case error_code::output_unprocessed:
    return "output_unprocessed";
  case error_code::quiescence_unexpected:
    return "quiescence_unexpected";
  case error_code::lps_unavailable:
    return "lps_unavailable";
  case error_code::internal_error:
    return "internal_error";
  }
  return "unknown";
}

json::object make_hello(const std::string& lps_identifier, const std::string& lps_hash)
{
  return {{"type", "hello"},
    {"role", "mbt_tool"},
    {"id", next_id()},
    {"protocol_version", "0.2"},
    // TODO: mCRL2 version
    {"tool", {{"name", "mbt_tool"}, {"version", "0.2"}}},
    {"lps", {{"identifier", lps_identifier}, {"hash", lps_hash}}}
  };
}

json::object make_heartbeat()
{
  return {{"type", "heartbeat"}, {"id", next_id()}};
}

json::object make_ack(message_id in_reply_to, std::string_view kind)
{
  return {{"type", "ack"}, {"id", next_id()}, {"in_reply_to", in_reply_to}, {"kind", std::string(kind)}};
}

json::object make_warning(message_id in_reply_to, std::string_view code, std::string_view detail)
{
  json::object obj{{"type", "warning"}, {"id", next_id()}, {"in_reply_to", in_reply_to}, {"code", std::string(code)}};
  if (!detail.empty())
  {
    obj["detail"] = std::string(detail);
  }
  return obj;
}

json::object make_error(message_id in_reply_to, error_code code, std::string_view detail)
{
  json::object obj{{"type", "error"},
    {"id", next_id()},
    {"in_reply_to", in_reply_to},
    {"code", std::string(to_string(code))}};
  if (!detail.empty())
  {
    obj["detail"] = std::string(detail);
  }
  return obj;
}

json::object make_error_no_id(error_code code, std::string_view detail)
{
  json::object obj{{"type", "error"}, {"id", next_id()}, {"code", std::string(to_string(code))}};
  if (!detail.empty())
  {
    obj["detail"] = std::string(detail);
  }
  return obj;
}

json::object make_enabled(message_id in_reply_to,
  const std::vector<std::string>& inputs,
  const std::vector<std::string>& outputs,
  bool quiescence)
{
  json::array input_arr;
  for (const auto& s: inputs)
  {
    input_arr.push_back(json::value(s));
  }
  json::array output_arr;
  for (const auto& s: outputs)
  {
    output_arr.push_back(json::value(s));
  }

  return {{"type", "enabled"},
    {"id", next_id()},
    {"in_reply_to", in_reply_to},
    {"inputs", std::move(input_arr)},
    {"outputs", std::move(output_arr)},
    {"quiescence", quiescence}};
}

json::object make_reset_ack(message_id in_reply_to)
{
  return {{"type", "reset_ack"}, {"id", next_id()}, {"in_reply_to", in_reply_to}};
}

json::object make_close(std::string_view reason)
{
  json::object obj{{"type", "close"}, {"id", next_id()}};
  if (!reason.empty())
  {
    obj["reason"] = std::string(reason);
  }
  return obj;
}

incoming_message parse_incoming(std::string_view text)
{
  json::value val;
  try
  {
    val = json::parse(text);
  }
  catch (const std::exception& e)
  {
    throw std::runtime_error(std::string("malformed JSON: ") + e.what());
  }

  if (!val.is_object())
  {
    throw std::runtime_error("expected a JSON object");
  }

  json::object obj = val.as_object();

  const json::value* type_val = obj.if_contains("type");
  if (type_val == nullptr || !type_val->is_string())
  {
    throw std::runtime_error("missing or non-string 'type' field");
  }

  incoming_message msg;
  msg.type = json::value_to<std::string>(*type_val);
  msg.payload = obj;

  if (const json::value* id_val = obj.if_contains("id"))
  {
    if (id_val->is_number())
    {
      msg.id = json::value_to<message_id>(*id_val);
    }
  }

  if (const json::value* irt_val = obj.if_contains("in_reply_to"))
  {
    if (irt_val->is_number())
    {
      msg.in_reply_to = json::value_to<message_id>(*irt_val);
    }
  }

  return msg;
}

namespace
{

template<typename T>
T get_or(const json::object& obj, std::string_view key, T default_value)
{
  const json::value* val = obj.if_contains(key);
  if (val != nullptr && val->is_number())
  {
    return json::value_to<T>(*val);
  }
  return default_value;
}

} // namespace

session_config parse_session_config(const json::object& config_obj)
{
  session_config cfg;
  cfg.tau_closure_depth = get_or<std::size_t>(config_obj, "tau_closure_depth", cfg.tau_closure_depth);
  cfg.early_output_timeout_ms
    = get_or<std::uint64_t>(config_obj, "early_output_timeout_ms", cfg.early_output_timeout_ms);
  cfg.heartbeat_interval_ms = get_or<std::uint64_t>(config_obj, "heartbeat_interval_ms", cfg.heartbeat_interval_ms);
  cfg.heartbeat_timeout_ms = get_or<std::uint64_t>(config_obj, "heartbeat_timeout_ms", cfg.heartbeat_timeout_ms);
  return cfg;
}

} // namespace mcrl2::mbt_protocol
