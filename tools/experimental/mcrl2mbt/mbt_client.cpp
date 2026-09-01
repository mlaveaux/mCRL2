// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mbt_client.h"

#include <stdexcept>

#include <boost/asio/connect.hpp>

#include "mcrl2/utilities/exception.h"
#include "mcrl2/utilities/logger.h"

namespace mcrl2
{

namespace
{

struct websocket_url
{
  std::string host;
  std::string port;
  std::string target;
};

websocket_url parse_websocket_url(const std::string& url)
{
  const std::string prefix = "ws://";
  if (url.compare(0, prefix.size(), prefix) != 0)
  {
    throw mcrl2::runtime_error("only ws:// URLs are supported");
  }

  const std::string remainder = url.substr(prefix.size());
  if (remainder.empty())
  {
    throw mcrl2::runtime_error("invalid URL: missing host");
  }

  std::string authority = remainder;
  std::string target = "/";
  const std::size_t slash = remainder.find('/');
  if (slash != std::string::npos)
  {
    authority = remainder.substr(0, slash);
    target = remainder.substr(slash);
  }

  if (authority.empty())
  {
    throw mcrl2::runtime_error("invalid URL: missing host");
  }

  std::string host = authority;
  std::string port = "80";
  const std::size_t colon = authority.rfind(':');
  if (colon != std::string::npos)
  {
    host = authority.substr(0, colon);
    port = authority.substr(colon + 1);
    if (host.empty() || port.empty())
    {
      throw mcrl2::runtime_error("invalid URL: malformed host:port");
    }
  }

  return {host, port, target};
}

} // namespace

mbt_client::mbt_client(net::io_context& ioc)
  : m_ioc(ioc),
    m_ws(ioc)
{}

void mbt_client::connect(const std::string& url)
{
  const websocket_url endpoint = parse_websocket_url(url);

  tcp::resolver resolver{m_ioc};
  const auto results = resolver.resolve(endpoint.host, endpoint.port);
  net::connect(m_ws.next_layer(), results.begin(), results.end());
  m_ws.handshake(endpoint.host + ":" + endpoint.port, endpoint.target);
  m_connected = true;

  schedule_read();
}

void mbt_client::send(const json::object& message)
{
  const std::string serialised = json::serialize(message);
  m_write_queue.push_back(serialised);
  if (!m_writing)
  {
    flush_write_queue();
  }
}

void mbt_client::run()
{
  m_ioc.run();
}

void mbt_client::close(std::string_view reason)
{
  if (m_closing || !m_connected)
  {
    return;
  }
  m_closing = true;
  websocket::close_reason cr{websocket::close_code::normal};
  cr.reason.assign(reason.data(), reason.size());
  m_ws.async_close(cr,
    [](beast::error_code ec)
    {
      if (ec && ec != websocket::error::closed)
      {
        mCRL2log(log::warning) << "WebSocket close error: " << ec.message() << std::endl;
      }
    });
}

void mbt_client::schedule_read()
{
  m_ws.async_read(m_read_buffer, [this](beast::error_code ec, std::size_t bytes) { on_read(ec, bytes); });
}

void mbt_client::on_read(beast::error_code ec, std::size_t bytes)
{
  if (ec == websocket::error::closed || ec == net::error::eof)
  {
    m_connected = false;
    m_ioc.stop();
    return;
  }

  if (ec)
  {
    handle_error("read error: " + ec.message());
    return;
  }

  const std::string text = beast::buffers_to_string(m_read_buffer.data());
  m_read_buffer.consume(bytes);

  try
  {
    const auto msg = mbt_protocol::parse_incoming(text);
    if (m_message_handler)
    {
      m_message_handler(msg);
    }
  }
  catch (const std::exception& e)
  {
    mCRL2log(log::warning) << "failed to parse incoming message: " << e.what() << std::endl;
  }

  schedule_read();
}

void mbt_client::flush_write_queue()
{
  if (m_write_queue.empty() || m_closing)
  {
    m_writing = false;
    return;
  }

  m_writing = true;
  const std::string& front = m_write_queue.front();
  m_ws.text(true);
  m_ws.async_write(net::buffer(front), [this](beast::error_code ec, std::size_t bytes) { on_write(ec, bytes); });
}

void mbt_client::on_write(beast::error_code ec, std::size_t bytes)
{
  (void)bytes;
  if (ec)
  {
    handle_error("write error: " + ec.message());
    return;
  }

  m_write_queue.pop_front();
  flush_write_queue();
}

void mbt_client::handle_error(const std::string& description)
{
  mCRL2log(log::error) << "mbt_client: " << description << std::endl;
  if (m_error_handler)
  {
    m_error_handler(description);
  }
  m_connected = false;
  m_ioc.stop();
}

} // namespace mcrl2
