// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE mock_adapter_test
#include <boost/test/included/unit_test.hpp>

#include <cstdio>
#include <string>
#include <thread>

#include "boost/json/src.hpp" // NOLINT — single TU definition of Boost.JSON
#include <boost/asio/ip/tcp.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>

#include "mcrl2/data/rewriter.h"
#include "mcrl2/lps/explorer_options.h"
#include "mcrl2/lps/parse.h"

#include "../io_classifier.cpp" // NOLINT
#include "../mbt_client.cpp" // NOLINT
#include "../mbt_protocol.cpp" // NOLINT
#include "../mbt_session.cpp" // NOLINT
#include "../state_set.cpp" // NOLINT

namespace beast = boost::beast;
namespace websocket = beast::websocket;
namespace net = boost::asio;
namespace json = boost::json;
using tcp = net::ip::tcp;

using namespace mcrl2;

// Single output action looping; "out" is classified as output.
static const std::string SPEC_TEXT = "act out;\n"
                                     "proc P = out.P;\n"
                                     "init P;\n";

/// Describes the sequence of assertions and replies the mock adapter performs.
struct mock_adapter_result
{
  bool hello_received = false;
  bool get_enabled_acked = false;
  bool reset_acked = false;
  bool close_sent = false;
};

/// \brief Send a serialised JSON object over the WebSocket.
static void ws_send(websocket::stream<tcp::socket>& ws, const json::object& msg)
{
  ws.write(net::buffer(json::serialize(msg)));
}

/// \brief Read one message and parse it as a JSON object.
static json::object ws_recv(websocket::stream<tcp::socket>& ws)
{
  beast::flat_buffer buf;
  ws.read(buf);
  json::value v = json::parse(beast::buffers_to_string(buf.data()));
  BOOST_REQUIRE(v.is_object());
  return v.as_object();
}

static mock_adapter_result run_mock_adapter(tcp::acceptor& acceptor)
{
  mock_adapter_result result;

  net::io_context ioc;
  tcp::socket sock{ioc};
  acceptor.accept(sock);

  websocket::stream<tcp::socket> ws{std::move(sock)};
  ws.accept();

  // 1. Receive the hello message from the tool.
  {
    const auto msg = ws_recv(ws);
    BOOST_REQUIRE(msg.contains("type"));
    BOOST_REQUIRE_EQUAL(json::value_to<std::string>(msg.at("type")), "hello");
    BOOST_REQUIRE(msg.contains("protocol_version"));
    BOOST_REQUIRE_EQUAL(json::value_to<std::string>(msg.at("protocol_version")), "0.2");
    BOOST_REQUIRE(msg.contains("lps"));
    BOOST_REQUIRE(msg.at("lps").is_object());
    BOOST_REQUIRE(msg.at("lps").as_object().contains("hash"));
    result.hello_received = true;
  }

  // 2. Reply with hello + minimal config (disabling heartbeat timeouts).
  {
    json::object reply;
    reply["type"] = "hello";
    reply["role"] = "adapter";
    reply["protocol_version"] = "0.2";
    reply["config"] = json::object{{"tau_closure_depth", json::value(static_cast<std::int64_t>(0))},
      {"early_output_timeout_ms", json::value(static_cast<std::int64_t>(0))},
      {"heartbeat_interval_ms", json::value(static_cast<std::int64_t>(3600000))},
      {"heartbeat_timeout_ms", json::value(static_cast<std::int64_t>(3600000))}};
    ws_send(ws, reply);
  }

  // 3. Send get_enabled — expect an enabled response back.
  {
    json::object req;
    req["type"] = "get_enabled";
    req["id"] = json::value(static_cast<std::int64_t>(1));
    ws_send(ws, req);

    const auto resp = ws_recv(ws);
    BOOST_REQUIRE(resp.contains("type"));
    BOOST_REQUIRE_EQUAL(json::value_to<std::string>(resp.at("type")), "enabled");
    BOOST_REQUIRE(resp.contains("outputs"));
    BOOST_REQUIRE(resp.at("outputs").is_array());
    BOOST_REQUIRE(resp.at("outputs").as_array().size() >= 1u);
    // Each output is a wire multi-action array.
    const auto& ma = resp.at("outputs").as_array()[0].as_array();
    BOOST_REQUIRE(ma.size() >= 1u);
    BOOST_CHECK_EQUAL(json::value_to<std::string>(ma[0].as_object().at("name")), "out");
    result.get_enabled_acked = true;
  }

  // 4. Send reset — expect reset_ack.
  {
    json::object req;
    req["type"] = "reset";
    req["id"] = json::value(static_cast<std::int64_t>(2));
    ws_send(ws, req);

    const auto resp = ws_recv(ws);
    BOOST_REQUIRE(resp.contains("type"));
    BOOST_REQUIRE_EQUAL(json::value_to<std::string>(resp.at("type")), "reset_ack");
    result.reset_acked = true;
  }

  // 5. Send close.
  {
    json::object req;
    req["type"] = "close";
    req["reason"] = "test complete";
    ws_send(ws, req);
    result.close_sent = true;
  }

  return result;
}

BOOST_AUTO_TEST_CASE(mock_adapter_round_trip)
{
  // Bind a listener on an OS-assigned port.
  net::io_context listener_ioc;
  tcp::acceptor acceptor{listener_ioc, tcp::endpoint(tcp::v4(), 0)};
  const unsigned short port = acceptor.local_endpoint().port();
  const std::string url = "ws://127.0.0.1:" + std::to_string(port);

  // Run the mock adapter (blocking) in a background thread.
  mock_adapter_result adapter_result;
  std::exception_ptr adapter_exc;
  std::thread adapter_thread(
    [&]()
    {
      try
      {
        adapter_result = run_mock_adapter(acceptor);
      }
      catch (...)
      {
        adapter_exc = std::current_exception();
      }
    });

  // Build the session.
  lps::specification spec = lps::parse_linear_process_specification(SPEC_TEXT);
  data::rewriter rewr(spec.data());
  lps::explorer_options opts;

  // Create a classifier: "out" is an output action.
  std::string tmpl = "/tmp/mbt_ma_test_XXXXXX";
  const int fd = mkstemp(tmpl.data());
  BOOST_REQUIRE(fd >= 0);
  const std::string io_content = "input\noutput\n  out\n";
  ::write(fd, io_content.data(), io_content.size());
  ::close(fd);
  io_classifier cls;
  cls.load(tmpl);
  std::remove(tmpl.c_str());

  mbt_session session(spec, opts, rewr, cls, "test_spec", "deadbeef");
  session.connect(url);
  session.run(); // blocks until adapter closes the connection

  adapter_thread.join();

  // Re-throw any exception from the adapter thread.
  if (adapter_exc)
  {
    std::rethrow_exception(adapter_exc);
  }

  BOOST_CHECK(adapter_result.hello_received);
  BOOST_CHECK(adapter_result.get_enabled_acked);
  BOOST_CHECK(adapter_result.reset_acked);
  BOOST_CHECK(adapter_result.close_sent);
  BOOST_CHECK(!session.had_error());
}
