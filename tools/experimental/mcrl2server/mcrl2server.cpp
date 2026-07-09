// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2server.cpp
/// \brief A websocket server that exposes the state space of an LPS.
///
/// The server loads a linear process specification (LPS) and uses the lps
/// explorer to compute the outgoing transitions of states on demand. Clients
/// connect over a websocket and request the initial state and the successors
/// of previously discovered states. Discovered states are assigned stable
/// integer identifiers so that they can be referred to across messages.

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/asio/ip/tcp.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>
#include <boost/json.hpp>
// Compile Boost.JSON in this translation unit so that no separate library has
// to be linked. This header must be included in exactly one source file.
#include <boost/json/src.hpp>

#include "mcrl2/atermpp/standard_containers/indexed_set.h"
#include "mcrl2/data/rewriter.h"
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/lps/explorer.h"
#include "mcrl2/lps/explorer_options.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/state.h"
#include "mcrl2/utilities/input_tool.h"
#include "mcrl2/utilities/logger.h"

namespace beast = boost::beast;
namespace websocket = beast::websocket;
namespace net = boost::asio;
namespace json = boost::json;
using tcp = net::ip::tcp;

using namespace mcrl2;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;
using mcrl2::data::tools::rewriter_tool;

namespace
{

/// \brief Wraps the lps explorer and assigns stable integer identifiers to the
///        states that have been discovered so far.
///
/// \details Only the non-stochastic, non-timed case is handled for now. The
///          explorer is not thread-safe, so a single instance must only ever
///          be driven from one thread (see the single-threaded server below).
class state_space
{
public:
  using explorer_type = lps::explorer<false, false, lps::specification>;

  /// \brief An outgoing transition: an action label and the identifier of the
  ///        target state.
  struct transition
  {
    std::string label;
    std::size_t target;
  };

  state_space(const lps::specification& spec, const lps::explorer_options& options, const data::rewriter& rewr)
    : m_rewr(rewr),
      m_explorer(spec, options, rewr)
  {}

  /// \brief Computes the initial state, indexes it and returns its identifier.
  std::size_t initial_state()
  {
    data::mutable_indexed_substitution<> sigma;
    lps::state s0;
    const data::data_expression_list& init = m_explorer.initial_state();
    lps::make_state(s0,
        init.begin(),
        init.size(),
        [&](data::data_expression& result, const data::data_expression& x) { m_rewr(result, x, sigma); });
    return index_of(s0);
  }

  /// \brief Returns the outgoing transitions of the state with the given
  ///        identifier, registering any newly discovered target states.
  std::vector<transition> outgoing(std::size_t id)
  {
    std::vector<transition> result;
    if (id >= m_states.size())
    {
      throw mcrl2::runtime_error("unknown state identifier " + std::to_string(id));
    }
    const lps::state source = m_states.at(id);
    for (const auto& edge : m_explorer.out_edges(source))
    {
      result.push_back(transition{lps::pp(edge.action), index_of(edge.state)});
    }
    return result;
  }

  /// \brief Returns a human-readable label for a state identifier.
  std::string label(std::size_t id) const
  {
    return id < m_states.size() ? lps::pp(m_states.at(id)) : std::string();
  }

  /// \brief The number of states discovered so far.
  std::size_t discovered() const { return m_states.size(); }

private:
  std::size_t index_of(const lps::state& s) { return m_states.insert(s).first; }

  data::rewriter m_rewr;
  explorer_type m_explorer;
  atermpp::indexed_set<lps::state> m_states;
};

/// \brief Dispatches a single request and produces a JSON response.
///
/// Supported commands:
///   {"command":"initial"}                -> {"id":N,"label":"..."}
///   {"command":"expand","state":N}       -> {"id":N,"transitions":[{"label":"a","target":K},...]}
///   {"command":"info"}                   -> {"states_discovered":N}
std::string handle_request(state_space& space, const std::string& request)
{
  try
  {
    const json::object& obj = json::parse(request).as_object();
    const std::string command = json::value_to<std::string>(obj.at("command"));

    json::object response;
    if (command == "initial")
    {
      const std::size_t id = space.initial_state();
      response["id"] = id;
      response["label"] = space.label(id);
    }
    else if (command == "expand")
    {
      if (!obj.contains("state"))
      {
        return json::serialize(json::object{{"error", "missing 'state' field"}});
      }
      const std::size_t id = json::value_to<std::size_t>(obj.at("state"));
      json::array transitions;
      for (const auto& t : space.outgoing(id))
      {
        transitions.push_back(json::object{{"label", t.label}, {"target", t.target}});
      }
      response["id"] = id;
      response["transitions"] = std::move(transitions);
    }
    else if (command == "info")
    {
      response["states_discovered"] = space.discovered();
    }
    else
    {
      response["error"] = "unknown command '" + command + "'";
    }
    return json::serialize(response);
  }
  catch (const std::exception& e)
  {
    return json::serialize(json::object{{"error", e.what()}});
  }
}

/// \brief A minimal synchronous, single-connection websocket server.
///
/// \details The server accepts one connection at a time and serves requests in
///          a blocking loop. This matches the single-threaded explorer; an
///          asynchronous, multi-session design is left as future work.
class websocket_server
{
public:
  websocket_server(unsigned short port, state_space& space)
    : m_port(port),
      m_state_space(space)
  {}

  void run()
  {
    net::io_context ioc{1};
    tcp::acceptor acceptor{ioc, tcp::endpoint{tcp::v4(), m_port}};
    mCRL2log(log::info) << "mcrl2server listening on ws://0.0.0.0:" << m_port << std::endl;

    for (;;)
    {
      tcp::socket socket{ioc};
      acceptor.accept(socket);
      try
      {
        session(std::move(socket));
      }
      catch (const beast::system_error& e)
      {
        if (e.code() != websocket::error::closed)
        {
          mCRL2log(log::warning) << "session error: " << e.what() << std::endl;
        }
      }
      catch (const std::exception& e)
      {
        mCRL2log(log::warning) << "session error: " << e.what() << std::endl;
      }
    }
  }

private:
  void session(tcp::socket socket)
  {
    websocket::stream<tcp::socket> ws{std::move(socket)};
    ws.accept();
    for (;;)
    {
      beast::flat_buffer buffer;
      ws.read(buffer); // throws websocket::error::closed when the peer disconnects
      const std::string request = beast::buffers_to_string(buffer.data());
      const std::string response = handle_request(m_state_space, request);
      ws.text(true);
      ws.write(net::buffer(response));
    }
  }

  unsigned short m_port;
  state_space& m_state_space;
};

} // namespace

class mcrl2server_tool : public rewriter_tool<input_tool>
{
protected:
  using super = rewriter_tool<input_tool>;

  unsigned short m_port = 8080;

  void add_options(interface_description& desc) override
  {
    super::add_options(desc);
    desc.add_option("port", make_mandatory_argument("PORT"), "listen for websocket connections on TCP port PORT (default 8080)", 'p');
  }

  void parse_options(const command_line_parser& parser) override
  {
    super::parse_options(parser);
    if (parser.options.count("port") > 0)
    {
      m_port = parser.option_argument_as<unsigned short>("port");
    }
  }

public:
  mcrl2server_tool()
    : super(
        "mcrl2server",
        "Maurice Laveaux",
        "expose the state space of an LPS over a websocket connection",
        "Read the LPS in INFILE and start a websocket server that answers queries about its "
        "state space. Clients request the initial state and the outgoing transitions of "
        "previously discovered states. If INFILE is not present, stdin is used.")
  {}

  bool run() override
  {
    lps::specification spec;
    lps::load_lps(spec, m_input_filename);

    data::rewriter rewr(spec.data(), rewrite_strategy());

    lps::explorer_options options;
    options.number_of_threads = 1;

    state_space space(spec, options, rewr);
    websocket_server server(m_port, space);
    server.run();

    return true;
  }
};

int main(int argc, char** argv)
{
  return mcrl2server_tool().execute(argc, argv);
}
