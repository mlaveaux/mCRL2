// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE mbt_protocol_test
#include <boost/test/included/unit_test.hpp>

// Provide Boost.JSON implementation in this translation unit.
#include <boost/json/src.hpp>

#include "../mbt_protocol.cpp" // NOLINT

using namespace mcrl2::mbt_protocol;
namespace json = boost::json;

static incoming_message round_trip(const json::object& msg)
{
  return parse_incoming(json::serialize(msg));
}

BOOST_AUTO_TEST_CASE(parse_incoming_basic)
{
  const auto msg = round_trip(make_ack(5, "output"));
  BOOST_CHECK_EQUAL(msg.type, "ack");
  BOOST_REQUIRE(msg.id.has_value());
  BOOST_REQUIRE(msg.in_reply_to.has_value());
  BOOST_CHECK_EQUAL(*msg.in_reply_to, 5u);
}

BOOST_AUTO_TEST_CASE(parse_incoming_malformed_json)
{
  BOOST_CHECK_THROW(parse_incoming("not json"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(parse_incoming_missing_type)
{
  BOOST_CHECK_THROW(parse_incoming("{\"id\":1}"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(parse_incoming_non_object)
{
  BOOST_CHECK_THROW(parse_incoming("[1,2,3]"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(session_config_overrides)
{
  const json::object obj{{"tau_closure_depth", json::value(500)},
    {"early_output_timeout_ms", json::value(200)},
    {"heartbeat_interval_ms", json::value(1000)},
    {"heartbeat_timeout_ms", json::value(3000)}};
  const session_config cfg = parse_session_config(obj);
  BOOST_CHECK_EQUAL(cfg.tau_closure_depth, 500u);
  BOOST_CHECK_EQUAL(cfg.early_output_timeout_ms, 200u);
  BOOST_CHECK_EQUAL(cfg.heartbeat_interval_ms, 1000u);
  BOOST_CHECK_EQUAL(cfg.heartbeat_timeout_ms, 3000u);
}

BOOST_AUTO_TEST_CASE(ids_are_unique)
{
  const auto a = make_heartbeat();
  const auto b = make_heartbeat();
  BOOST_CHECK(json::value_to<message_id>(a.at("id")) != json::value_to<message_id>(b.at("id")));
}
