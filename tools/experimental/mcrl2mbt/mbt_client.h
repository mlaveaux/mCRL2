// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <deque>
#include <functional>
#include <string>
#include <string_view>

#include <boost/asio/ip/tcp.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>
#include <boost/json.hpp>

#include "mbt_protocol.h"

namespace beast = boost::beast;
namespace websocket = beast::websocket;
namespace net = boost::asio;
namespace json = boost::json;
using tcp = net::ip::tcp;

namespace mcrl2
{

/// \brief Async WebSocket client used by the MBT session.
///
/// All callbacks are invoked on the io_context strand (single-threaded).
/// Use send() to queue outgoing messages; they are flushed one at a time.
class mbt_client
{
public:
  using message_handler = std::function<void(const mbt_protocol::incoming_message&)>;
  using error_handler = std::function<void(const std::string& description)>;

  explicit mbt_client(net::io_context& ioc);

  /// \brief Connect and perform WebSocket handshake (synchronous convenience wrapper).
  /// \throws std::runtime_error on failure.
  void connect(const std::string& url);

  /// \brief Queue a message for async send. Thread-safe from the io_context thread.
  void send(const json::object& message);

  /// \brief Set callback for every complete incoming message.
  void set_message_handler(message_handler handler) { m_message_handler = std::move(handler); }

  /// \brief Set callback for unrecoverable errors (called before close).
  void set_error_handler(error_handler handler) { m_error_handler = std::move(handler); }

  /// \brief Run the io_context until the session ends (blocking).
  void run();

  /// \brief Initiate graceful close.
  void close(std::string_view reason = "");

  bool is_connected() const { return m_connected; }

private:
  void schedule_read();
  void on_read(beast::error_code ec, std::size_t bytes);
  void flush_write_queue();
  void on_write(beast::error_code ec, std::size_t bytes);
  void handle_error(const std::string& description);

  net::io_context& m_ioc;
  websocket::stream<tcp::socket> m_ws;
  beast::flat_buffer m_read_buffer;
  std::deque<std::string> m_write_queue;
  bool m_writing = false;
  bool m_connected = false;
  bool m_closing = false;
  message_handler m_message_handler;
  error_handler m_error_handler;
};

} // namespace mcrl2
