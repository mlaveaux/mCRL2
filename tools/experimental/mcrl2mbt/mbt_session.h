// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <list>
#include <string>

#include <boost/asio/io_context.hpp>
#include <boost/asio/steady_timer.hpp>

#include "mcrl2/data/rewriter.h"
#include "mcrl2/lps/explorer_options.h"
#include "mcrl2/lps/specification.h"

#include "io_classifier.h"
#include "mbt_client.h"
#include "mbt_protocol.h"
#include "state_set.h"

namespace net = boost::asio;

namespace mcrl2
{

/// \brief Implements the MBT protocol session over a WebSocket connection.
///
/// The session connects to an adapter, exchanges hello handshakes, and then
/// handles input/output/quiescence messages according to IOCO semantics.
class mbt_session
{
public:
  enum class status
  {
    disconnected,
    handshaking,
    ready,
    closed
  };

  mbt_session(const lps::specification& spec,
    const lps::explorer_options& options,
    const data::rewriter& rewr,
    const io_classifier& classifier,
    const std::string& lps_identifier,
    const std::string& lps_hash);

  /// \brief Connect to the adapter WebSocket server.
  void connect(const std::string& url);

  /// \brief Run until the session ends (blocking — drives io_context).
  void run();

  status current_status() const { return m_status; }

  /// \brief True if the session ended due to a connection or protocol error.
  bool had_error() const { return m_had_error; }

private:
  void on_message(const mbt_protocol::incoming_message& msg);
  void on_hello_response(const mbt_protocol::incoming_message& msg);
  void on_input(const mbt_protocol::incoming_message& msg);
  void on_output(const mbt_protocol::incoming_message& msg);
  void on_quiescence(const mbt_protocol::incoming_message& msg);
  void on_get_enabled(const mbt_protocol::incoming_message& msg);
  void on_reset(const mbt_protocol::incoming_message& msg);
  void on_close(const mbt_protocol::incoming_message& msg);
  void on_heartbeat(const mbt_protocol::incoming_message& msg);

  void arm_heartbeat_send_timer();
  void arm_heartbeat_recv_watchdog();

  /// \brief Process early outputs that may now be enabled after S changed.
  void reevaluate_early_set();

  /// \brief Parse "multi_action" from a message payload as a wire representation.
  static mbt_protocol::wire_multi_action extract_action(const mbt_protocol::incoming_message& msg);

  net::io_context m_ioc{1}; // must be declared before m_client and the timers
  mbt_client m_client;
  state_set m_states;
  io_classifier m_classifier;
  mbt_protocol::session_config m_config;
  status m_status = status::disconnected;

  bool m_had_error = false;
  std::string m_lps_identifier;
  std::string m_lps_hash;

  net::steady_timer m_send_timer;
  net::steady_timer m_recv_watchdog;

  struct early_entry
  {
    mbt_protocol::wire_multi_action label;
    mbt_protocol::message_id reply_to;
    net::steady_timer timer;
  };
  std::list<early_entry> m_early_outputs;
};

} // namespace mcrl2
