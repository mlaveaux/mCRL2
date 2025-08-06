// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/number_postfix_generator.h
/// \brief add your file description here.

#ifndef MCRL2_UTILITIES_NUMBER_POSTFIX_GENERATOR_H
#define MCRL2_UTILITIES_NUMBER_POSTFIX_GENERATOR_H

#include <cassert>
#include <cctype>
#include <map>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/text_utility.cxx"
#else
  import utilities;
#endif

namespace mcrl2::utilities {

/// \brief Identifier generator that generates names consisting of a prefix followed by a number.
/// For each prefix an index is maintained, that is incremented after each call to operator()(prefix).
class number_postfix_generator
{
  protected:
    /// \brief A map that maintains the highest index for each prefix.
    mutable std::map<std::string, std::size_t> m_index;

    /// \brief The default hint.
    std::string m_hint;

  public:

    /// \brief Constructor.
    number_postfix_generator(std::string hint = "FRESH_VAR")
      : m_hint(hint)
    {}

    /// \brief Adds the strings in the range [first, last) to the context.
    /// \param id A string
    void add_identifier(const std::string& id)
    {
      std::string::size_type i = id.find_last_not_of("0123456789");
      std::size_t new_index = 0;
      std::string name;
      if (i == std::string::npos || id.size() == i + 1) // string does not end with a number
      {
        name = id;
      }
      else
      {
        name = id.substr(0, i + 1);
        std::string num = id.substr(i + 1);
        new_index = atoi(num.c_str());
      }
      std::size_t old_index = m_index.find(name) == m_index.end() ? 0 : m_index[name];
      m_index[name] = (std::max)(old_index, new_index); // Windows requires brackets around std::max. 
    }

    /// \brief Adds the strings in the range [first, last) to the context.
    /// \param first
    /// \param last
    /// [first, last) is a sequence of strings that is used as context.
    template <typename Iter>
    void add_identifiers(Iter first, Iter last)
    {
      for (Iter i = first; i != last; ++i)
      {
        add_identifier(*i);
      }
    }

    /// \brief Constructor.
    /// \param first
    /// \param last
    /// \param hint String hint to use as prefix for generated strings.
    /// [first, last) is a sequence of strings that is used as context.
    template <typename Iter>
    number_postfix_generator(Iter first, Iter last, std::string hint = "FRESH_VAR")
      : m_hint(hint)
    {
      add_identifiers(first, last);
    }

    /// \brief Generates a fresh identifier that doesn't appear in the context.
    /// \return A fresh identifier.
    std::string operator()(std::string hint, bool add_to_context = true) const
    {
      // remove digits at the end of hint
      if (std::isdigit(hint[hint.size() - 1]))
      {
        std::string::size_type i = hint.find_last_not_of("0123456789");
        hint = hint.substr(0, i + 1);
      }

      auto j = m_index.find(hint);
      if (j == m_index.end())
      {
        if (add_to_context)
        {
          m_index[hint] = 0;
        }
        return hint;
      }
      return hint + utilities::number2string(add_to_context ? ++(j->second) : j->second + 1);
    }

    /// \brief Generates a fresh identifier that doesn't appear in the context.
    /// \return A fresh identifier.
    std::string operator()() const
    {
      return (*this)(m_hint, true);
    }

    /// \brief Returns the default hint.
    const std::string& hint() const
    {
      return m_hint;
    }

    /// \brief Returns the default hint.
    std::string& hint()
    {
      return m_hint;
    }

    /// \brief Clear the context of the generator
    void clear()
    {
      m_index.clear();
    }
};

} // namespace mcrl2::utilities



#endif // MCRL2_UTILITIES_NUMBER_POSTFIX_GENERATOR_H
