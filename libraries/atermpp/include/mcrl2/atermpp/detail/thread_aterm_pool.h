// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef ATERMPP_DETAIL_THREAD_ATERM_POOL_H
#define ATERMPP_DETAIL_THREAD_ATERM_POOL_H

#include "mcrl2/atermpp/detail/aterm_pool.h"
#include "mcrl2/atermpp/detail/aterm_container.h"
#include "mcrl2/utilities/hashtable.h"

#include <atomic>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/shared_mutex.cxx"
#else
  import utilities;
#endif

namespace atermpp::detail
{

/// \brief This is a thread's specific access to the global aterm pool which ensures that
///        garbage collection and hash table resizing can proceed.
class thread_aterm_pool final : public mcrl2::utilities::noncopyable
{
public:
  thread_aterm_pool(aterm_pool& global_pool)
      : m_pool(global_pool),
        m_shared_mutex(global_pool.shared_mutex()),
        m_variables(new mcrl2::utilities::hashtable<aterm_core*>()),
        m_containers(new mcrl2::utilities::hashtable<detail::aterm_container*>()),
        m_thread_interface(
            global_pool,
            [this] { mark(); },
            [this] { print_local_performance_statistics(); },
            [this] { return protection_set_size(); })
  {
    /// Identify the first constructor call as the main thread.
    static bool is_main_thread = true;
    if (is_main_thread)
    {
      m_is_main_thread = true;
      is_main_thread = false;
    }
  }

  ~thread_aterm_pool() 
  {
      // We leak values for the global aterm pool since they contain global variables (for which initialisation order is undefined).
    if (!m_is_main_thread)
    {
      // We need to prematurely unregister this thread pool since we are going to delete reference variables.
      m_thread_interface.unregister();
      delete m_variables;
      delete m_containers;
    }
  }

  /// \details threadsafe
  inline function_symbol create_function_symbol(const std::string& name, const std::size_t arity, const bool check_for_registered_functions = false);

  /// \details threadsafe
  inline function_symbol create_function_symbol(std::string&& name, const std::size_t arity, const bool check_for_registered_functions = false);

  /// \details threadsafe
  inline void create_int(aterm& term, std::size_t val);

  /// \details threadsafe
  inline void create_term(aterm& term, const function_symbol& sym);

  /// \details threadsafe
  template<class ...Terms>
  inline void create_appl(aterm& term, const function_symbol& sym, const Terms&... arguments);

  /// \details threadsafe
  template<class Term, class INDEX_TYPE, class ...Terms>
  inline void create_appl_index(aterm& term, const function_symbol& sym, const Terms&... arguments);

  /// \details threadsafe
  template<typename ForwardIterator>
  inline void create_appl_dynamic(aterm& term,
      const function_symbol& sym,
      ForwardIterator begin,
      ForwardIterator end);

  /// \threadsafe
  template<typename InputIterator, typename ATermConverter>
  inline void create_appl_dynamic(aterm& term,
      const function_symbol& sym,
      ATermConverter convert_to_aterm,
      InputIterator begin,
      InputIterator end);

  /// \brief Consider the given variable when marking underlying terms.
  inline void register_variable(aterm_core* variable);

  /// \brief Removes the given variable from the active variables.
  inline void deregister_variable(aterm_core* variable);

  /// \brief Consider the given container when marking underlying terms.
  inline void register_container(aterm_container* variable);

  /// \brief Removes the given container from the active variables.
  inline void deregister_container(aterm_container* variable);

  // Implementation of thread_aterm_pool_interface
  inline void mark();
  inline void print_local_performance_statistics() const;
  inline std::size_t protection_set_size() const;

  /// Acquire a shared lock on this thread aterm pool.
  inline mcrl2::utilities::shared_guard lock_shared() { return m_shared_mutex.lock_shared(); }

  /// Acquire an exclusive lock
  inline mcrl2::utilities::lock_guard lock() { return m_shared_mutex.lock(); }

  /// Returns true iff we are in a shared section.
  inline bool is_shared_locked() { return m_shared_mutex.is_shared_locked(); }

  /// Triggers a global garbage collection
  inline void collect() { m_pool.collect(m_shared_mutex); }

private:
  aterm_pool& m_pool;

  /// Keeps track of pointers to all existing aterm variables and containers.
  mcrl2::utilities::shared_mutex m_shared_mutex;  
  mcrl2::utilities::hashtable<aterm_core*>* m_variables;
  mcrl2::utilities::hashtable<detail::aterm_container*>* m_containers;

  std::size_t m_variable_insertions = 0;
  std::size_t m_container_insertions = 0;
  std::stack<std::reference_wrapper<_aterm>> m_todo; ///< A reusable todo stack.

  bool m_is_main_thread = false;

  /// The registered thread aterm pool.
  thread_aterm_pool_interface m_thread_interface; 
};

/// \brief A reference to the thread local term pool storage
thread_aterm_pool& g_thread_term_pool();

} // namespace atermpp::detail


#include "thread_aterm_pool_implementation.h"

#endif // ATERMPP_DETAIL_ATERM_POOL_H
