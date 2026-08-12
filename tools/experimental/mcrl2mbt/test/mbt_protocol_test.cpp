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

// ─── helpers ─────────────────────────────────────────────────────────────────

static incoming_message round_trip(const json::object& msg)
{
  return parse_incoming(json::serialize(msg));
}

// ─── make_hello ──────────────────────────────────────────────────────────────

BOOST_AUTO_TEST_CASE(hello_fields)
{
  const auto msg = make_hello("my.lps", "deadbeef");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("type")), "hello");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("lps_identifier")), "my.lps");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("lps_hash")), "deadbeef");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("protocol_version")), "0.2");
  BOOST_CHECK(msg.contains("id"));
}

// ─── make_heartbeat ──────────────────────────────────────────────────────────

BOOST_AUTO_TEST_CASE(heartbeat_fields)
{
  const auto msg = make_heartbeat();
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("type")), "heartbeat");
  BOOST_CHECK(msg.contains("id"));
}

// ─── make_ack ────────────────────────────────────────────────────────────────

BOOST_AUTO_TEST_CASE(ack_fields)
{
  const auto msg = make_ack(42, "input");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("type")), "ack");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("kind")), "input");
  BOOST_CHECK_EQUAL(json::value_to<message_id>(msg.at("in_reply_to")), 42u);
}

// ─── make_error / make_error_no_id ───────────────────────────────────────────

BOOST_AUTO_TEST_CASE(error_fields)
{
  const auto msg = make_error(7, error_code::input_not_enabled, "act_a");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("type")), "error");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("code")), "input_not_enabled");
  BOOST_CHECK_EQUAL(json::value_to<message_id>(msg.at("in_reply_to")), 7u);
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("detail")), "act_a");
}

BOOST_AUTO_TEST_CASE(error_no_id_no_detail)
{
  const auto msg = make_error_no_id(error_code::internal_error);
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("type")), "error");
  BOOST_CHECK(!msg.contains("in_reply_to"));
  BOOST_CHECK(!msg.contains("detail"));
}

// ─── make_enabled ────────────────────────────────────────────────────────────

BOOST_AUTO_TEST_CASE(enabled_fields)
{
  const auto msg = make_enabled(3, {"a", "b"}, {"c"}, true);
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("type")), "enabled");
  BOOST_CHECK_EQUAL(msg.at("inputs").as_array().size(), 2u);
  BOOST_CHECK_EQUAL(msg.at("outputs").as_array().size(), 1u);
  BOOST_CHECK_EQUAL(msg.at("quiescence").as_bool(), true);
}

// ─── make_reset_ack ──────────────────────────────────────────────────────────

BOOST_AUTO_TEST_CASE(reset_ack_fields)
{
  const auto msg = make_reset_ack(99);
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("type")), "reset_ack");
  BOOST_CHECK_EQUAL(json::value_to<message_id>(msg.at("in_reply_to")), 99u);
}

// ─── make_close ──────────────────────────────────────────────────────────────

BOOST_AUTO_TEST_CASE(close_with_reason)
{
  const auto msg = make_close("done");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("type")), "close");
  BOOST_CHECK_EQUAL(json::value_to<std::string>(msg.at("reason")), "done");
}

BOOST_AUTO_TEST_CASE(close_without_reason)
{
  const auto msg = make_close();
  BOOST_CHECK(!msg.contains("reason"));
}

// ─── parse_incoming ──────────────────────────────────────────────────────────

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

// ─── parse_session_config ────────────────────────────────────────────────────

BOOST_AUTO_TEST_CASE(session_config_defaults)
{
  const session_config cfg = parse_session_config(json::object{});
  BOOST_CHECK_EQUAL(cfg.tau_closure_depth, 1000u);
  BOOST_CHECK_EQUAL(cfg.early_output_timeout_ms, 1000u);
  BOOST_CHECK_EQUAL(cfg.heartbeat_interval_ms, 5000u);
  BOOST_CHECK_EQUAL(cfg.heartbeat_timeout_ms, 15000u);
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

// ─── message id monotonically increases ──────────────────────────────────────

BOOST_AUTO_TEST_CASE(ids_are_unique)
{
  const auto a = make_heartbeat();
  const auto b = make_heartbeat();
  BOOST_CHECK(json::value_to<message_id>(a.at("id")) != json::value_to<message_id>(b.at("id")));
}
