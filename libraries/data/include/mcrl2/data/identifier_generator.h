// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/identifier_generator.h
/// \brief The class identifier_generator.

#ifndef MCRL2_DATA_IDENTIFIER_GENERATOR_H
#define MCRL2_DATA_IDENTIFIER_GENERATOR_H

#include <set>

#include "mcrl2/core/identifier_string.h"
#include "mcrl2/utilities/number_postfix_generator.h"

namespace mcrl2::data
{

/// \brief Abstract base class for identifier generators.
/// Identifier generators generate fresh names that do not appear in a
/// given context.
/// A context is maintained containing already used identifiers.
/// Using the operator()() and operator()(std::string) fresh
/// identifiers are generated that do not appear in the context.
template <typename Generator = utilities::number_postfix_generator>
class identifier_generator
{
  protected:
    Generator m_generator;

  public:
    /// \brief Constructor.
    identifier_generator() = default;

    /// \brief Destructor.
    virtual ~identifier_generator() = default;

    /// \brief Clears the context.
    virtual void clear_context() = 0;

    /// \brief Adds the identifier s to the context.
    /// \param s An identifier.
    virtual void add_identifier(const core::identifier_string& s) = 0;

    /// \brief Removes the identifier s from the context.
    /// \param s An identifier.
    virtual void remove_identifier(const core::identifier_string& s) = 0;

    /// \brief Add a set of identifiers to the context.
    void add_identifiers(const std::set<core::identifier_string>& ids)
    {
      for (const core::identifier_string& id: ids)
      {
        add_identifier(id);
      }
    }

    /// \brief Remove a set of identifiers from the context.
    void remove_identifiers(const std::set<core::identifier_string>& ids)
    {
      for (const core::identifier_string& id: ids)
      {
        remove_identifier(id);
      }
    }

    /// \brief Returns true if the identifier s appears in the context.
    /// \param s An identifier.
    /// \return True if the identifier appears in the context.
    virtual bool has_identifier(const core::identifier_string& s) const = 0;

    /// \brief Returns a fresh identifier, with the given hint as prefix.
    /// The returned identifier is added to the context.
    /// \param hint A string
    /// \param add_to_context If true, the freshly generated identifier is added
    /// to the context to make sure no duplicates are generated.
    /// \return A fresh identifier.
    virtual core::identifier_string operator()(const std::string& hint, bool add_to_context = true)
    {
      core::identifier_string id(add_to_context?hint:m_generator(hint));
      while (has_identifier(id))
      {
        id = core::identifier_string(m_generator(hint));
      }
      if (add_to_context)
      {
        add_identifier(id);
      }
      return id;
    }
};

} // namespace mcrl2::data



#endif // MCRL2_DATA_IDENTIFIER_GENERATOR_H
