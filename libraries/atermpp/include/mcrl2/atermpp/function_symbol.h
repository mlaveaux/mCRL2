// Author(s): Wieger Wesselink, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef MCRL2_ATERMPP_FUNCTION_SYMBOL_H
#define MCRL2_ATERMPP_FUNCTION_SYMBOL_H

#include <mcrl3_ffi.h>

#include <string>
#include <string_view>
#include <utility>

namespace atermpp
{

class function_symbol
{
  friend class function_symbol_generator;
  friend struct std::hash<function_symbol>;
  
public:
  function_symbol() = default;

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

  [[nodiscard]] bool defined() const
  {
    return m_function_symbol.ptr != nullptr;
  }

  /// \brief Return the name of the function_symbol.
  /// \return The name of the function symbol.
  [[nodiscard]] std::string_view name() const
  {
    return mcrl3::ffi::function_symbol_get_name(&m_function_symbol);
  }

  /// \brief Return the arity (number of arguments) of the function symbol (function_symbol).
  /// \return The arity of the function symbol.
  [[nodiscard]] std::size_t arity() const
  {
    return mcrl3::ffi::function_symbol_get_arity(&m_function_symbol);
  }

  /// \brief Equality test.
  /// \details This operator compares the indices of the function symbols. This means
  ///         that this operation takes constant time.
  /// \returns True iff the function symbols are the same.
  bool operator ==(const function_symbol& f) const
  {
    return m_function_symbol.ptr == f.m_function_symbol.ptr;
  }

  /// \brief Comparison operation.
  /// \details This operator takes constant time.
  /// \returns True iff this function has a lower index than the argument.
  std::weak_ordering operator <=>(const function_symbol& f) const
  {
    return m_function_symbol.ptr <=> f.m_function_symbol.ptr;
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
  /// \brief Constructor for internal use only
  inline
  function_symbol(const std::string& name, const std::size_t arity, const bool check_for_registered_functions)
  {
    mcrl3::ffi::function_symbol_create(name.c_str(), name.length(), arity, check_for_registered_functions);
  }

  /// \brief The shared reference to the underlying function symbol.
  mcrl3::ffi::function_symbol_t m_function_symbol;
};

class global_function_symbol : public function_symbol
{
public:
  /// \brief Defines a function symbol from a name and arity combination.
  /// \details This constructor should be used by global function symbols.
  inline
  global_function_symbol(const std::string& name, const std::size_t arity)
  {

  }
};

} // namespace atermpp

namespace std 
{

  /// \brief Specialisation of the standard hash function for function_symbol.
template<>
struct hash<atermpp::function_symbol>
{
  std::size_t operator()(const atermpp::function_symbol& f) const
  {
    // Function symbols take 48 bytes in memory, so when they are packed there
    // are at least 32 bits that do not distinguish two function symbols. As
    // such these can be removed.
    return reinterpret_cast<std::uint64_t>(f.m_function_symbol.ptr) >> 5;
  }
};

} // namespace std

#endif // MCRL2_ATERMPP_FUNCTION_SYMBOL_H
