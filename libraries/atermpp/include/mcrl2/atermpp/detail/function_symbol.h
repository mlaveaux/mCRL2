// Author(s): Wieger Wesselink, Jan Friso Groote. Based on the aterm library by Paul Klint and others.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/detail/function_symbol.h
/// \brief A foundational class for function symbols. 

#ifndef DETAIL_FUNCTION_SYMBOL_H
#define DETAIL_FUNCTION_SYMBOL_H

#include <string>
#include <cstddef>

#include "mcrl2/utilities/shared_reference.h"
#include "mcrl2/atermpp/detail/aterm_configuration.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/noncopyable.cxx"
#else
  import utilities;
#endif

namespace atermpp::detail
{

/// \brief Stores the data for a function symbol (name, arity) pair.
class _function_symbol : public mcrl2::utilities::shared_reference_counted<_function_symbol, mcrl2::utilities::detail::GlobalThreadSafe>, private mcrl2::utilities::noncopyable
{
public:
  /// \brief A shared reference for this object.
  using ref = mcrl2::utilities::shared_reference<const _function_symbol>;

  _function_symbol(const std::string& name, std::size_t arity) :
     m_arity(arity),
     m_name(name)
  {}

  const std::string& name() const noexcept
  {
    return m_name;
  }

  std::size_t arity() const noexcept
  {
    return m_arity;
  }

  bool operator==(const _function_symbol& f) const noexcept
  {
    return m_arity == f.m_arity && m_name == f.m_name;
  }

private:
  const std::size_t m_arity;
  const std::string m_name;
};

} // namespace atermpp::detail


#endif // DETAIL_FUNCTION_SYMBOL_H

