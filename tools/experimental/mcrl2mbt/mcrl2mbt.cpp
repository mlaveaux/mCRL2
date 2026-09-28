// Author(s): Maurice Laveaux
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/// \file mcrl2mbt.cpp
/// \brief Model-Based Testing (MBT) tool that connects to an adapter over
///        WebSocket and exercises an LPS under test using IOCO semantics.

#include <boost/crc.hpp>
#include <boost/json/src.hpp>

#include <string>

#include "mcrl2/data/rewriter.h"
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/lps/explorer_options.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/utilities/input_tool.h"

#include "io_classifier.h"
#include "mbt_session.h"

using namespace mcrl2;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;
using mcrl2::data::tools::rewriter_tool;

class mcrl2mbt_tool : public rewriter_tool<input_tool>
{
  using super = rewriter_tool<input_tool>;

  std::string m_adapter_url;
  std::string m_io_file;
  lps::explorer_options m_options;

protected:
  void add_options(interface_description& desc) override
  {
    super::add_options(desc);
    desc.add_option("adapter",
      make_mandatory_argument("URL"),
      "connect to adapter at WebSocket URL (e.g. ws://localhost:8080/mbt)",
      'a');
    desc.add_option("io-file",
      make_optional_argument("FILE", ""),
      "action classification file with 'input'/'output' sections",
      'f');
    desc.add_option("cached", "use enumeration caching techniques to speed up state exploration. ");
    desc.add_hidden_option("global-cache", "use a global cache instead of a cache per summand");
    desc.add_hidden_option("no-one-point-rule-rewrite", "do not apply the one point rule rewriter");
    desc.add_hidden_option("no-replace-constants-by-variables", "do not move constant expressions to a substitution");
  }

  void parse_options(const command_line_parser& parser) override
  {
    super::parse_options(parser);
    if (parser.options.count("adapter") > 0)
    {
      m_adapter_url = parser.option_argument("adapter");
    }
    if (parser.options.count("io-file") > 0)
    {
      m_io_file = parser.option_argument("io-file");
    }

    m_options.number_of_threads = 1;
    m_options.rewrite_strategy = rewrite_strategy();
    m_options.remove_unused_rewrite_rules = true;
    m_options.cached = parser.has_option("cached");
    m_options.global_cache = parser.has_option("global-cache");
    m_options.one_point_rule_rewrite = !parser.has_option("no-one-point-rule-rewrite");
    m_options.replace_constants_by_variables = !parser.has_option("no-replace-constants-by-variables");
  }

  bool run() override
  {
    if (m_adapter_url.empty())
    {
      throw mcrl2::runtime_error("--adapter URL is required");
    }

    lps::specification spec;
    lps::load_lps(spec, m_input_filename);

    const data::rewriter rewr
      = lps::construct_rewriter(spec, m_options.rewrite_strategy, m_options.remove_unused_rewrite_rules);

    io_classifier classifier;
    if (!m_io_file.empty())
    {
      classifier.load(m_io_file);
    }

    // Use the input filename as the LPS identifier and its CRC32 as a cheap integrity hash.
    boost::crc_32_type result;
    result.process_bytes(m_input_filename.data(), m_input_filename.size());

    mbt_session session(spec, m_options, rewr, classifier, m_input_filename, std::to_string(result.checksum()));
    session.connect(m_adapter_url);
    session.run();
    return !session.had_error();
  }

public:
  mcrl2mbt_tool()
    : super("mcrl2mbt",
        "Maurice Laveaux",
        "MBT client connecting to an adapter over WebSocket",
        "Connect to an MBT adapter at URL, load LPS from INFILE, and perform "
        "IOCO-based testing according to the MBT protocol specification.")
  {}
};

int main(int argc, char** argv)
{
  return mcrl2mbt_tool().execute(argc, argv);
}
