// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef DETAIL_FUNCTION_SYMBOL_POOL_H
#define DETAIL_FUNCTION_SYMBOL_POOL_H
MCRL2_MODULE;

#include <memory>
#include <map>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/function_symbol_hash.cxx"
  #include "mcrl2/utilities/unordered_set.cxx"
  #include "mcrl2/utilities/mutex.cxx"
  #include "mcrl2/utilities/cache_metric.cxx"
#else
  export module atermpp:detail.function_symbol_pool;

  import :function_symbol;
  import :detail.aterm_configuration;
  import :detail.function_symbol_hash;
  
  import utilities;
#endif

MCRL2_MODULE_EXPORT namespace atermpp::detail
{

/// \brief This class stores a set of function symbols.
class function_symbol_pool : private mcrl2::utilities::noncopyable
{
public:
  function_symbol_pool();

  /// \brief Creates a function symbol pair (name, arity), returns a pointer to an existing element
  ///        if this pair is already defined.
  /// \param check_for_registered_functions Check whether there is a registered prefix p such that
  ///           name equal pn where n is a number. In that case prevent that pn will be generated
  ///           as a fresh function name.
  /// \threadsafe
  function_symbol create(const std::string& name, const std::size_t arity, const bool check_for_registered_functions = false);
  function_symbol create(std::string&& name, const std::size_t arity, const bool check_for_registered_functions = false);

  /// \brief Restore the index back to index before registering this prefix.
  /// \threadsafe
  void deregister(const std::string& prefix);

  /// \returns An index that is always a safe index for the given prefix.
  /// \todo These functions are all used by the function_symbol_generator and should probably not
  ///       be public.
  /// \threadsafe
  std::shared_ptr<std::size_t> register_prefix(const std::string& prefix);

  /// \brief Resize the function symbol pool if necessary.
  void resize_if_needed();

  /// \brief Get an index such that no function symbol with name prefix + returned value
  ///        and any value above it already exists.
  std::size_t get_sufficiently_large_postfix_index(const std::string& prefix) const;

  /// \brief Collect all garbage function symbols.
  void sweep();

  /// \returns The function symbol used by integral terms.
  const function_symbol& as_int() noexcept { return m_as_int; }

  /// \returns The function symbol used by the list constructor.
  const function_symbol& as_list() noexcept { return m_as_list; }

  /// \returns The function symbol used by the term indicating the empty list.
  const function_symbol& as_empty_list() noexcept { return m_as_empty_list; }

  /// \returns The number of function symbols stored in this pool.
  std::size_t size() const noexcept { return m_symbol_set.size(); }
  
  /// \returns The maximum number of function symbols stored in this pool.
  std::size_t capacity() const noexcept { return m_symbol_set.capacity(); }

private:
  using unordered_set = mcrl2::utilities::unordered_set<
    _function_symbol,
    function_symbol_hasher,
    function_symbol_equals,
    typename std::conditional_t<EnableBlockAllocator, 
      mcrl2::utilities::block_allocator<_function_symbol, 1024, mcrl2::utilities::detail::GlobalThreadSafe>, 
      std::allocator<_function_symbol>>,
    mcrl2::utilities::detail::GlobalThreadSafe,
    false>;

  /// \brief Stores the underlying function symbols.
  unordered_set m_symbol_set;

  /// \brief A map that records a function for each prefix that must be called to set the
  ///        postfix number to a sufficiently high number if a function symbol with the same
  ///        prefix string is registered.
  std::map<std::string, std::shared_ptr<std::size_t>> m_prefix_to_register_function_map;

  mutable mcrl2::utilities::mutex m_mutex; // Mutex for m_prefix_to_register_function_map.

  // Default function symbols.
  function_symbol m_as_int;
  function_symbol m_as_list;
  function_symbol m_as_empty_list;

  // Various performance metrics.
  mcrl2::utilities::cache_metric m_function_symbol_metrics; ///< Track the number of function symbols found in or added to the set.

  // Create helper function. 
  void create_helper(const std::string& name);
};


} // namespace atermpp::detail


#endif // DETAIL_FUNCTION_SYMBOL_POOL_H

