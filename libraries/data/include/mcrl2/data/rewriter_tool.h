// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/rewriter_tool.h
/// \brief Base class for tools that use a data rewriter.

#ifndef MCRL2_DATA_REWRITER_TOOL_H
#define MCRL2_DATA_REWRITER_TOOL_H

#include "mcrl2/data/detail/enumerator_iteration_limit.h"
#include "mcrl2/data/rewriter.h"
#include "mcrl2/utilities/command_line_interface.h"

namespace mcrl2
{

namespace data
{

namespace tools
{

/// \brief Base class for tools that use a rewriter.
template <typename Tool>
class rewriter_tool: public Tool
{
  protected:
    /// The data rewriter strategy
    data::rewrite_strategy m_rewrite_strategy;

    /// \brief Add options to an interface description. Also includes
    /// rewriter options.
    /// \param desc An interface description.
    /// \param suppress_jittyp Boolean that if true will prevent showing that jittyp is an option for rewriting.
    void add_options(utilities::interface_description& desc, bool suppress_jittyp)
    {
      Tool::add_options(desc);

      utilities::interface_description::enum_argument<data::rewrite_strategy> rewriter_option("NAME");
      rewriter_option.add_value(data::jitty, true);
      rewriter_option.add_value(data::innermost);
#ifdef MCRL2_JITTYC_AVAILABLE
      rewriter_option.add_value(data::jitty_compiling);
#endif
      if (!suppress_jittyp)
      {
        rewriter_option.add_value(data::jitty_prover);
      }

      desc.add_option(
        "rewriter",
        rewriter_option,
        "use rewrite strategy NAME:"
        ,'r'
      );

      desc.add_option(
        "qlimit",
        utilities::make_mandatory_argument("NUM"),
        "limit enumeration of quantifiers to NUM iterations. (Default NUM=1000, NUM=0 for unlimited).",
        'Q'
      );
    }

    /// \brief Add options to an interface description. Also includes
    /// rewriter options.
    /// \param desc An interface description
    void add_options(utilities::interface_description& desc) override
    {
      add_options(desc, false);
    }

    /// \brief Parse non-standard options
    /// \param parser A command line parser
    void parse_options(const utilities::command_line_parser& parser) override
    {
      Tool::parse_options(parser);
      m_rewrite_strategy = parser.option_argument_as< data::rewrite_strategy >("rewriter");

      if(parser.options.count("qlimit"))
      {
        //Set enumerator limit for quantifier enumeration
        std::size_t qlimit = parser.option_argument_as< std::size_t >("qlimit");
        data::detail::set_enumerator_iteration_limit(qlimit == 0 ? std::numeric_limits<std::size_t>::max() : qlimit);
      }
    }

  public:

    /// \brief Constructor.
    rewriter_tool(const std::string& name,
                  const std::string& author,
                  const std::string& what_is,
                  const std::string& tool_description,
                  std::string known_issues = ""
                 )
      : Tool(name, author, what_is, tool_description, known_issues),
        m_rewrite_strategy(mcrl2::data::jitty)
    {}

    /// \brief Returns the rewrite strategy
    /// \return The rewrite strategy
    data::rewrite_strategy rewrite_strategy() const
    {
      return m_rewrite_strategy;
    }

    /// \brief Creates a data rewriter as specified on the command line.
    /// \param data_spec A data specification
    /// \return A data rewriter
    data::rewriter create_rewriter(const data::data_specification& data_spec = data::data_specification())
    {
      return data::rewriter(data_spec, rewrite_strategy());
    }
};

} // namespace tools

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_REWRITER_TOOL_H
