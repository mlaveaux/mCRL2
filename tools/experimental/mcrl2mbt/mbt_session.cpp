// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mbt_session.h"

#include <algorithm>
#include <chrono>
#include <stdexcept>

#include "mcrl2/utilities/logger.h"

namespace mcrl2
{

mbt_session::mbt_session(const lps::specification& spec,
  const lps::explorer_options& options,
  const data::rewriter& rewr,
  const io_classifier& classifier,
  const std::string& lps_identifier,
  const std::string& lps_hash)
  : m_client(m_ioc),
    m_states(spec, options, rewr),
    m_classifier(classifier),
    m_lps_identifier(lps_identifier),
    m_lps_hash(lps_hash),
    m_send_timer(m_ioc),
    m_recv_watchdog(m_ioc)
{
  m_states.reset_to_initial();

  m_client.set_message_handler([this](const mbt_protocol::incoming_message& msg) { on_message(msg); });

  m_client.set_error_handler(
    [this](const std::string& desc)
    {
      mCRL2log(log::error) << "connection error: " << desc << std::endl;
      m_had_error = true;
      m_status = status::closed;
    });
}

void mbt_session::connect(const std::string& url)
{
  m_client.connect(url);
  m_status = status::handshaking;

  // Send hello immediately after connecting.
  m_client.send(mbt_protocol::make_hello(m_lps_identifier, m_lps_hash));
  arm_heartbeat_recv_watchdog();
}

void mbt_session::run()
{
  m_client.run();
}

// ─── message dispatch ─────────────────────────────────────────────────────────

void mbt_session::on_message(const mbt_protocol::incoming_message& msg)
{
  arm_heartbeat_recv_watchdog(); // reset watchdog on every received message

  if (msg.type == "heartbeat")
  {
    on_heartbeat(msg);
    return;
  }

  if (msg.type == "close")
  {
    on_close(msg);
    return;
  }

  if (m_status == status::handshaking)
  {
    if (msg.type == "hello")
    {
      on_hello_response(msg);
    }
    else
    {
      const auto err
        = mbt_protocol::make_error_no_id(mbt_protocol::error_code::protocol_mismatch, "expected hello response");
      m_client.send(err);
      m_client.close("protocol mismatch");
    }
    return;
  }

  if (m_status != status::ready)
  {
    return;
  }

  if (msg.type == "input")
  {
    on_input(msg);
  }
  else if (msg.type == "output")
  {
    on_output(msg);
  }
  else if (msg.type == "quiescence")
  {
    on_quiescence(msg);
  }
  else if (msg.type == "get_enabled")
  {
    on_get_enabled(msg);
  }
  else if (msg.type == "reset")
  {
    on_reset(msg);
  }
  else
  {
    if (msg.id)
    {
      m_client.send(mbt_protocol::make_error(*msg.id, mbt_protocol::error_code::unknown_type, msg.type));
    }
    else
    {
      m_client.send(mbt_protocol::make_error_no_id(mbt_protocol::error_code::unknown_type, msg.type));
    }
  }
}

void mbt_session::on_hello_response(const mbt_protocol::incoming_message& msg)
{
  try
  {
    // Validate protocol version before accepting the session configuration.
    const json::value* ver_val = msg.payload.if_contains("protocol_version");
    if (ver_val == nullptr || !ver_val->is_string())
    {
      m_had_error = true;
      m_client.send(mbt_protocol::make_error_no_id(mbt_protocol::error_code::protocol_mismatch,
        "missing or non-string protocol_version"));
      m_client.close("protocol mismatch");
      return;
    }

    const std::string role = json::value_to<std::string>(*msg.payload.if_contains("role"));
    if (role != "adapter")
    {
      m_had_error = true;
      m_client.send(mbt_protocol::make_error_no_id(mbt_protocol::error_code::protocol_mismatch,
        "unexpected role: " + role));
      m_client.close("protocol mismatch");
      return;
    }

    const std::string ver = json::value_to<std::string>(*ver_val);
    if (ver != "0.2")
    {
      m_had_error = true;
      m_client.send(mbt_protocol::make_error_no_id(mbt_protocol::error_code::protocol_mismatch,
        "unsupported protocol version: " + ver));
      m_client.close("protocol mismatch");
      return;
    }

    const json::value* config_val = msg.payload.if_contains("config");
    if (config_val != nullptr && config_val->is_object())
    {
      m_config = mbt_protocol::parse_session_config(config_val->as_object());
    }

    m_status = status::ready;
    arm_heartbeat_send_timer();

    mCRL2log(log::verbose) << "hello handshake complete; session ready" << std::endl;
  }
  catch (const std::exception& e)
  {
    m_had_error = true;
    m_client.send(mbt_protocol::make_error_no_id(mbt_protocol::error_code::internal_error, e.what()));
    m_client.close("internal error during hello");
  }
}

void mbt_session::on_input(const mbt_protocol::incoming_message& msg)
{
  if (!msg.id)
  {
    m_client.send(
      mbt_protocol::make_error_no_id(mbt_protocol::error_code::malformed_message, "input message missing id"));
    return;
  }

  std::string label;
  try
  {
    label = extract_action(msg);
  }
  catch (const std::exception& e)
  {
    m_client.send(mbt_protocol::make_error(*msg.id, mbt_protocol::error_code::malformed_message, e.what()));
    return;
  }

  if (!m_states.apply_action(label, m_config.tau_closure_depth))
  {
    m_client.send(mbt_protocol::make_error(*msg.id, mbt_protocol::error_code::input_not_enabled, label));
    return;
  }

  reevaluate_early_set();
  m_client.send(mbt_protocol::make_ack(*msg.id, "input"));
}

void mbt_session::on_output(const mbt_protocol::incoming_message& msg)
{
  if (!msg.id)
  {
    m_client.send(
      mbt_protocol::make_error_no_id(mbt_protocol::error_code::malformed_message, "output message missing id"));
    return;
  }

  std::string label;
  try
  {
    label = extract_action(msg);
  }
  catch (const std::exception& e)
  {
    m_client.send(mbt_protocol::make_error(*msg.id, mbt_protocol::error_code::malformed_message, e.what()));
    return;
  }
  const mbt_protocol::message_id reply_to = *msg.id;

  if (m_states.apply_action(label, m_config.tau_closure_depth))
  {
    reevaluate_early_set();
    m_client.send(mbt_protocol::make_ack(reply_to, "output"));
    return;
  }

  // Output not currently enabled — put it in the early set.
  m_client.send(mbt_protocol::make_warning(reply_to, "output_early"));

  if (m_config.early_output_timeout_ms == 0)
  {
    // Timeout zero: immediately expire (fire error inline, no timer needed).
    m_client.send(mbt_protocol::make_error(reply_to, mbt_protocol::error_code::output_unprocessed, label));
    return;
  }

  m_early_outputs.push_back({label, reply_to, net::steady_timer{m_ioc}});
  auto it = std::prev(m_early_outputs.end());
  it->timer.expires_after(std::chrono::milliseconds(m_config.early_output_timeout_ms));
  // Capture by value so this handler is safe even if the list is reorganised.
  const mbt_protocol::message_id id_copy = reply_to;
  const std::string label_copy = label;
  it->timer.async_wait(
    [this, id_copy, label_copy](const boost::system::error_code& ec)
    {
      if (ec == boost::asio::error::operation_aborted)
      {
        return;
      }

      // Find the entry by stable reply-to id (on_reset may have cleared the list already).
      auto found = std::find_if(m_early_outputs.begin(),
        m_early_outputs.end(),
        [id_copy](const early_entry& e) { return e.reply_to == id_copy; });
      if (found == m_early_outputs.end())
      {
        return; // Removed by on_reset before this handler ran.
      }
      m_client.send(mbt_protocol::make_error(id_copy, mbt_protocol::error_code::output_unprocessed, label_copy));
      m_early_outputs.erase(found);
    });
}

void mbt_session::on_quiescence(const mbt_protocol::incoming_message& msg)
{
  if (!msg.id)
  {
    m_client.send(
      mbt_protocol::make_error_no_id(mbt_protocol::error_code::malformed_message, "quiescence message missing id"));
    return;
  }

  try
  {
    if (!m_states.refine_to_quiescent(m_config.tau_closure_depth, m_classifier))
    {
      m_client.send(mbt_protocol::make_error(*msg.id,
        mbt_protocol::error_code::quiescence_unexpected,
        "no quiescent states reachable"));
      return;
    }

    reevaluate_early_set();
    m_client.send(mbt_protocol::make_ack(*msg.id, "quiescence"));
  }
  catch (const std::exception& e)
  {
    m_client.send(mbt_protocol::make_error(*msg.id, mbt_protocol::error_code::internal_error, e.what()));
  }
}

void mbt_session::on_get_enabled(const mbt_protocol::incoming_message& msg)
{
  if (!msg.id)
  {
    m_client.send(
      mbt_protocol::make_error_no_id(mbt_protocol::error_code::malformed_message, "get_enabled message missing id"));
    return;
  }

  try
  {
    const auto enabled = m_states.get_enabled(m_config.tau_closure_depth, m_classifier);
    m_client.send(mbt_protocol::make_enabled(*msg.id, enabled.inputs, enabled.outputs, enabled.quiescence));
  }
  catch (const std::exception& e)
  {
    m_client.send(mbt_protocol::make_error(*msg.id, mbt_protocol::error_code::internal_error, e.what()));
  }
}

void mbt_session::on_reset(const mbt_protocol::incoming_message& msg)
{
  if (!msg.id)
  {
    m_client.send(
      mbt_protocol::make_error_no_id(mbt_protocol::error_code::malformed_message, "reset message missing id"));
    return;
  }

  try
  {
    m_states.reset_to_initial();

    // Cancel all early output timers without sending error/ack.
    for (auto& entry: m_early_outputs)
    {
      entry.timer.cancel();
    }
    m_early_outputs.clear();

    m_client.send(mbt_protocol::make_reset_ack(*msg.id));
  }
  catch (const std::exception& e)
  {
    m_client.send(mbt_protocol::make_error(*msg.id, mbt_protocol::error_code::internal_error, e.what()));
  }
}

void mbt_session::on_close(const mbt_protocol::incoming_message& msg)
{
  (void)msg;
  m_status = status::closed;
  m_send_timer.cancel();
  m_recv_watchdog.cancel();
  m_client.close();
}

void mbt_session::on_heartbeat(const mbt_protocol::incoming_message& /*msg*/)
{
  // Nothing to reply; watchdog was reset by on_message.
}

void mbt_session::arm_heartbeat_send_timer()
{
  m_send_timer.expires_after(std::chrono::milliseconds(m_config.heartbeat_interval_ms));
  m_send_timer.async_wait(
    [this](const boost::system::error_code& ec)
    {
      if (ec == boost::asio::error::operation_aborted)
      {
        return;
      }
      m_client.send(mbt_protocol::make_heartbeat());
      arm_heartbeat_send_timer();
    });
}

void mbt_session::arm_heartbeat_recv_watchdog()
{
  m_recv_watchdog.expires_after(std::chrono::milliseconds(m_config.heartbeat_timeout_ms));
  m_recv_watchdog.async_wait(
    [this](const boost::system::error_code& ec)
    {
      if (ec == boost::asio::error::operation_aborted)
      {
        return;
      }
      mCRL2log(log::error) << "heartbeat timeout: no message received from adapter" << std::endl;
      m_had_error = true;
      m_status = status::closed;
      m_client.close("heartbeat timeout");
    });
}

void mbt_session::reevaluate_early_set()
{
  bool found = true;
  while (found)
  {
    found = false;
    for (auto it = m_early_outputs.begin(); it != m_early_outputs.end(); ++it)
    {
      if (m_states.apply_action(it->label, m_config.tau_closure_depth))
      {
        it->timer.cancel();
        m_client.send(mbt_protocol::make_ack(it->reply_to, "output"));
        m_early_outputs.erase(it);
        found = true;
        break; // restart scan from the beginning
      }
    }
  }
}

std::string mbt_session::extract_action(const mbt_protocol::incoming_message& msg)
{
  const json::value* val = msg.payload.if_contains("multi_action");
  if (val == nullptr || !val->is_string())
  {
    throw std::runtime_error("message missing 'multi_action' field");
  }
  return json::value_to<std::string>(*val);
}

} // namespace mcrl2
