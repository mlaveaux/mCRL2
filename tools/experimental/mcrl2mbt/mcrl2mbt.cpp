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

#include <fstream>
#include <iomanip>
#include <sstream>

#include "mcrl2/data/rewriter.h"
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/lps/explorer_options.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/utilities/input_tool.h"
#include "mcrl2/utilities/logger.h"

#include "io_classifier.h"
#include "mbt_session.h"

using namespace mcrl2;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;
using mcrl2::data::tools::rewriter_tool;

namespace
{

std::string crc32_hex(const std::string& filename)
{
  std::ifstream f(filename, std::ios::binary);
  if (!f)
  {
    return "";
  }
  boost::crc_32_type crc;
  char buf[4096];
  while (f.read(buf, sizeof(buf)))
  {
    crc.process_bytes(buf, static_cast<std::size_t>(f.gcount()));
  }
  if (f.gcount() > 0)
  {
    crc.process_bytes(buf, static_cast<std::size_t>(f.gcount()));
  }
  std::ostringstream oss;
  oss << std::hex << std::setfill('0') << std::setw(8) << crc.checksum();
  return oss.str();
}

} // namespace

class mcrl2mbt_tool : public rewriter_tool<input_tool>
{
  using super = rewriter_tool<input_tool>;

  std::string m_adapter_url;
  std::string m_io_file;

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
      'i');
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
  }

  bool run() override
  {
    if (m_adapter_url.empty())
    {
      throw mcrl2::runtime_error("--adapter URL is required");
    }

    lps::specification spec;
    lps::load_lps(spec, m_input_filename);

    const data::rewriter rewr(spec.data(), rewrite_strategy());

    io_classifier classifier;
    if (!m_io_file.empty())
    {
      classifier.load(m_io_file);
    }

    lps::explorer_options options;
    options.number_of_threads = 1;

    // Use the input filename as the LPS identifier and its CRC32 as a cheap integrity hash.
    mbt_session session(spec, options, rewr, classifier, m_input_filename, crc32_hex(m_input_filename));
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
