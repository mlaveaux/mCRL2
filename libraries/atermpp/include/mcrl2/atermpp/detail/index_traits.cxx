// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_INDEX_TRAITS_H
#define MCRL2_ATERMPP_INDEX_TRAITS_H
MCRL2_MODULE;

#include <cstdint>
#include <mutex>
#include <stack>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/aterm_configuration.cxx"
  #include "mcrl2/atermpp/standard_containers/unordered_map.cxx"
#else
  export module atermpp:detail.index_traits;

  import :standard_containers.unordered_map;
  import :detail.aterm_configuration;
  import :detail.aterm_int_data;

  import utilities;
#endif

MCRL2_MODULE_EXPORT namespace atermpp::detail {

template <typename Variable, typename KeyType>
atermpp::unordered_map<KeyType, std::size_t>& variable_index_map()
{
  static atermpp::unordered_map<KeyType, std::size_t> m;
  return m;
}

template <typename Variable, typename KeyType>
std::stack<std::size_t>& variable_map_free_numbers()
{
  static std::stack<std::size_t> s;
  return s;
}

template <typename Variable, typename KeyType>
std::mutex& variable_mutex()
{
  static std::mutex m;
  return m;
}


template <typename Variable, typename KeyType>
std::size_t& variable_map_max_index()
{
  static std::size_t s;
  return s;
}

/// \brief For several variable types in mCRL2 an implicit mapping of these variables
/// to integers is available. This is done for efficiency reasons. Examples are:
///
/// data::variable, process::process_identifier
///
/// The class index_traits is used to implement this mapping. A traits class was chosen to
/// prevent pollution of the public interface of the classes that represent these variables.
///
/// N is the position of the index in the aterm.
template <typename Variable, typename KeyType, const int N>
struct index_traits
{

public:
  /// \brief Returns the index of the variable.
  static inline
  std::size_t index(const Variable& x)
  {
    const _aterm_int* i = reinterpret_cast<const _aterm_int*>(address(x[N]));
    return i->value();
  }

  /// \brief Returns an upper bound for the largest index of a variable that is currently in use.
  static inline
  std::size_t max_index()
  {
    return variable_map_max_index<Variable, KeyType>();
  }

  /// \brief Note: intended for internal use only!
  /// Returns the index of the variable. If the variable was not already in the map, it is added.
  static inline
  std::size_t insert(const KeyType& x)
  {
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe) { variable_mutex<Variable, KeyType>().lock(); }

    auto& m = variable_index_map<Variable, KeyType>();
    auto i = m.find(x);
    std::size_t value;
    if (i == m.end())
    {
      auto& s = variable_map_free_numbers<Variable, KeyType>();
      if (s.empty())
      {
        value = m.size();
        variable_map_max_index<Variable, KeyType>() = value;
      }
      else
      {
        value = s.top();
        s.pop();
      }
      m[x] = value;
    }
    else 
    {
      value = i->second;
    }
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe) { variable_mutex<Variable, KeyType>().unlock(); }
    return value;
  }

  /// \brief Note: intended for internal use only!
  /// Removes the variable from the index map.
  static inline
  void erase(const KeyType& x)
  {
    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe) { variable_mutex<Variable, KeyType>().lock(); }

    auto& m = variable_index_map<Variable, KeyType>();
    auto& s = variable_map_free_numbers<Variable, KeyType>();
    auto i = m.find(x);
    assert(i != m.end());
    s.push(i->second);
    m.erase(i);

    if constexpr (mcrl2::utilities::detail::GlobalThreadSafe) { variable_mutex<Variable, KeyType>().unlock(); }
  }

  /// \brief Note: intended for internal use only!
  /// Provides the size of the variable index map.
  static inline
  std::size_t size()
  {
    auto& m = variable_index_map<Variable, KeyType>();
    return m.size();
  }
};

} // namespace atermpp::detail



#endif // MCRL2_ATERMPP_INDEX_TRAITS_H
