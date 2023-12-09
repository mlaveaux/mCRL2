// Author(s): Wieger Wesselink, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef MCRL2_ATERMPP_FUNCTION_SYMBOL_H
#define MCRL2_ATERMPP_FUNCTION_SYMBOL_H

#include <utility>
#include "mcrl2/atermpp/detail/function_symbol.h"

namespace atermpp
{
namespace detail
{
class function_symbol_pool;
class _aterm;
}

class function_symbol
{
  friend class function_symbol_generator;
  friend class detail::function_symbol_pool;
  friend class detail::_aterm;

  friend struct std::hash<function_symbol>;

public:
  function_symbol() = default;

  /// \brief Defines a function symbol from a name and arity combination.
  function_symbol(const detail::_function_symbol* symbol)
   : m_function_symbol(symbol)
  {}

  const detail::_function_symbol* address() const
  {
    return m_function_symbol.get();
  }

  /// \brief Defines a function symbol from a name and arity combination.
  function_symbol(const std::string& name, const std::size_t arity_)
   : function_symbol(name, arity_, true)
  {}

  /// \brief Defines a function symbol from a name and arity combination.
  function_symbol(std::string&& name, const std::size_t arity_)
   : function_symbol(std::forward<std::string>(name), arity_, true)
  {}

  /// This class has non-trivial destructor so declare default copy and move operators.
  function_symbol(const function_symbol& other) noexcept = default;
  function_symbol& operator=(const function_symbol& other) noexcept = default;
  function_symbol(function_symbol&& other) noexcept = default;
  function_symbol& operator=(function_symbol&& other) noexcept = default;

  bool defined() const
  {
    return m_function_symbol.defined();
  }

  /// \brief Return the name of the function_symbol.
  /// \return The name of the function symbol.
  const std::string& name() const
  {
    return m_function_symbol->name();
  }

  /// \brief Return the arity (number of arguments) of the function symbol (function_symbol).
  /// \return The arity of the function symbol.
  std::size_t arity() const
  {
    return m_function_symbol->arity();
  }

  /// \brief Equality test.
  /// \details This operator compares the indices of the function symbols. This means
  ///         that this operation takes constant time.
  /// \returns True iff the function symbols are the same.
  bool operator ==(const function_symbol& f) const
  {
    return m_function_symbol == f.m_function_symbol;
  }

  /// \brief Inequality test.
  /// \details This operator takes constant time.
  /// \returns True iff the function symbols are not equal.
  bool operator !=(const function_symbol& f) const
  {
    return m_function_symbol != f.m_function_symbol;
  }

  /// \brief Comparison operation.
  /// \details This operator takes constant time.
  /// \returns True iff this function has a lower index than the argument.
  bool operator <(const function_symbol& f) const
  {
    return m_function_symbol < f.m_function_symbol;
  }

  /// \brief Comparison operation.
  /// \details This operator takes constant time.
  /// \returns True iff this function has a higher index than the argument.
  bool operator >(const function_symbol& f) const
  {
    return m_function_symbol > f.m_function_symbol;
  }

  /// \brief Comparison operation.
  /// \details This operator takes constant time.
  /// \returns True iff this function has a lower or equal index than the argument.
  bool operator <=(const function_symbol& f) const
  {
    return m_function_symbol <= f.m_function_symbol;
  }

  /// \brief Comparison operation.
  /// \details This operator takes constant time.
  /// \returns True iff this function has a larger or equal index than the argument.
  bool operator >=(const function_symbol& f) const
  {
    return m_function_symbol >= f.m_function_symbol;
  }

  /// \brief Swap this function with its argument.
  /// \details More efficient than assigning twice.
  /// \param f The function symbol with which the swap takes place.
  void swap(function_symbol& f)
  {
    using std::swap;
    swap(f.m_function_symbol, m_function_symbol);
  }

private:
  /// \brief Constructor for internal use only.
  function_symbol(detail::_function_symbol::ref&& f)
   : m_function_symbol(std::forward<detail::_function_symbol::ref>(f))
  {}

  /// \brief Constructor for internal use only
  function_symbol(const std::string& name, const std::size_t arity, const bool check_for_registered_functions);

  /// \brief Constructor for internal use only
  function_symbol(std::string&& name, const std::size_t arity, const bool check_for_registered_functions);

  /// \brief Calls the function symbol pool to free our used memory.
  void destroy();

  /// \brief The shared reference to the underlying function symbol.
  detail::_function_symbol::ref m_function_symbol;
};

class global_function_symbol : public function_symbol
{
public:
  /// \brief Defines a function symbol from a name and arity combination.
  /// \details This constructor should be used by global function symbols.
  global_function_symbol(const std::string& name, const std::size_t arity);
};

namespace detail
{
  /// \brief These function symbols are used to indicate integer, list and empty list terms.
  /// \details They are copied from the function_symbol_pool so that type_is_{int|list|appl} can be defined in the header.
  extern function_symbol g_as_int;
  extern function_symbol g_as_list;
  extern function_symbol g_as_empty_list;
}

} // namespace atermpp

#endif // MCRL2_ATERMPP_FUNCTION_SYMBOL_H
