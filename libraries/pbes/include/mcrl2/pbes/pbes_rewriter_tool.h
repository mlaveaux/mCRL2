// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbes_rewriter_tool.h
/// \brief Base class for tools that use a pbes rewriter.

#ifndef MCRL2_PBES_PBES_REWRITER_TOOL_H
#define MCRL2_PBES_PBES_REWRITER_TOOL_H

#include <set>
#include "mcrl2/utilities/command_line_interface.h"
#include "mcrl2/pbes/pbes_rewriter_type.h"

namespace mcrl2::pbes_system::tools
{

/// \brief Base class for filter tools that use a pbes rewriter.
template <typename Tool>
class pbes_rewriter_tool: public Tool
{
  protected:
    /// \brief The type of the pbes rewriter
    pbes_system::pbes_rewriter_type m_pbes_rewriter_type;

    /// \brief Returns the types of rewriters that are available for this tool.
    /// Override this method to change the standard behavior.
    /// \return The set { simplify, quantifier_all, quantifier_finite }
    virtual std::set<pbes_system::pbes_rewriter_type> available_rewriters() const
    {
      std::set<pbes_system::pbes_rewriter_type> result;
      result.insert(pbes_system::pbes_rewriter_type::simplify);
      result.insert(pbes_system::pbes_rewriter_type::quantifier_all);
      result.insert(pbes_system::pbes_rewriter_type::quantifier_finite);
      result.insert(pbes_system::pbes_rewriter_type::quantifier_inside);
      result.insert(pbes_system::pbes_rewriter_type::quantifier_one_point);
      result.insert(pbes_system::pbes_rewriter_type::pfnf);
      result.insert(pbes_system::pbes_rewriter_type::ppg);
      result.insert(pbes_system::pbes_rewriter_type::srf);
      result.insert(pbes_system::pbes_rewriter_type::pre_srf);
      result.insert(pbes_system::pbes_rewriter_type::bqnf_quantifier);
      return result;
    }

    /// \brief Returns the default pbes rewriter.
    /// Override this method to change the standard behavior.
    /// \return The string "simplify"
    virtual pbes_system::pbes_rewriter_type default_rewriter() const
    {
      return pbes_system::pbes_rewriter_type::simplify;
    }

    /// \brief Add options to an interface description. Also includes
    /// rewriter options.
    /// \param desc An interface description
    void add_options(utilities::interface_description& desc) override
    {
      Tool::add_options(desc);

      utilities::interface_description::enum_argument<pbes_system::pbes_rewriter_type> arg(utilities::make_enum_argument<pbes_system::pbes_rewriter_type>("NAME"));

      // Compute the available rewriters, and add the approriate arguments
      std::set<pbes_system::pbes_rewriter_type> types = available_rewriters();
      for (auto type : types)
      {
        arg.add_value(type, type==default_rewriter());
      }

      desc.add_option(
        "pbes-rewriter",
        arg,
        "use pbes rewrite strategy NAME:",
        'p'
      );
    }

    /// \brief Parse non-standard options
    /// \param parser A command line parser
    void parse_options(const utilities::command_line_parser& parser) override
    {
      Tool::parse_options(parser);
      m_pbes_rewriter_type = parser.option_argument_as<pbes_system::pbes_rewriter_type>("pbes-rewriter");
    }


  public:
    /// \brief Constructor.
    /// \param name The name of the tool
    /// \param author The author(s) of the tool
    /// \param what_is One-line "what is" description of the tool
    /// \param tool_description The description of the tool
    /// \param known_issues Known issues with the tool
    pbes_rewriter_tool(const std::string& name,
                       const std::string& author,
                       const std::string& what_is,
                       const std::string& tool_description,
                       std::string known_issues = ""
                      )
      : Tool(name, author, what_is, tool_description, known_issues)
    {}

    /// \brief Destructor.
    ~pbes_rewriter_tool() override = default;

    /// \brief Returns the rewriter type
    /// \return The rewriter type
    pbes_system::pbes_rewriter_type rewriter_type() const
    {
      return m_pbes_rewriter_type;
    }
};

} // namespace mcrl2::pbes_system::tools





#endif // MCRL2_PBES_PBES_REWRITER_TOOL_H
