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
#include "mcrl2/utilities/shared_mutex.h"

#include <atomic>


namespace atermpp::detail
{

// Forward declaration needed by the extern below.
class thread_aterm_pool;

/// \brief Pointer to the main thread's pool; set on its first access and used
///        by ~thread_aterm_pool to transfer orphaned protection-set entries
///        when a worker thread exits.
extern thread_aterm_pool* g_main_thread_pool;

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
    // We leak values for the global aterm pool since they contain global
    // variables (for which initialisation order is undefined).
    if (!m_is_main_thread)
    {
      // Transfer any remaining variables and containers to the main-thread
      // pool BEFORE unregistering.  This keeps them reachable during GC,
      // which is essential for static-local aterm_core objects (e.g. those
      // in sort/data header files) that were first initialised on this worker
      // thread and will be destroyed at program exit on the main thread.
      // Without the transfer those objects become invisible to the garbage
      // collector the moment this pool is unregistered, so the next GC cycle
      // would reclaim their underlying _aterm nodes.
      if (g_main_thread_pool != nullptr)
      {
        for (const aterm_core* v : *m_variables)
        {
          if (v != nullptr)
          {
            g_main_thread_pool->register_variable(const_cast<aterm_core*>(v));
          }
        }
        for (const aterm_container* c : *m_containers)
        {
          if (c != nullptr)
          {
            g_main_thread_pool->register_container(
                const_cast<aterm_container*>(c));
          }
        }
      }
      // We need to prematurely unregister this thread pool since we are going
      // to delete the reference-variable hashtables.
      m_thread_interface.unregister();
      delete m_variables;
      delete m_containers;
    }
  }

  /// \details threadsafe
  inline 
  function_symbol create_function_symbol(const std::string& name, std::size_t arity, bool check_for_registered_functions = false);

  /// \details threadsafe
  inline 
  function_symbol create_function_symbol(std::string&& name, std::size_t arity, bool check_for_registered_functions = false);

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
  [[nodiscard]] inline mcrl2::utilities::shared_guard lock_shared() { return mcrl2::utilities::shared_guard(m_shared_mutex); }

  /// Acquire an exclusive lock
  [[nodiscard]] inline mcrl2::utilities::lock_guard lock() { return mcrl2::utilities::lock_guard(m_shared_mutex); }

  /// Returns true iff we are in a shared section.
  inline bool is_shared_locked() { return m_shared_mutex.is_shared_locked(); }

  /// Provides access to the underlying shared mutex.
  inline mcrl2::utilities::shared_mutex& shared_mutex() { return m_shared_mutex; }

  /// Triggers a global garbage collection
  inline void collect() { m_pool.collect(m_shared_mutex); }

  /// Resizes the global hash tables when needed.
  inline void resize() { m_pool.resize_if_needed(m_shared_mutex); }

private:
  aterm_pool& m_pool;

  /// Keeps track of pointers to all existing aterm variables and containers.
  mcrl2::utilities::shared_mutex m_shared_mutex;  
  mcrl2::utilities::hashtable<aterm_core*>* m_variables;
  mcrl2::utilities::hashtable<detail::aterm_container*>* m_containers;

  std::size_t m_variable_insertions = 0;
  std::size_t m_container_insertions = 0;
  std::stack<std::reference_wrapper<_aterm>> m_todo; ///< A reusable todo stack.

  long m_count_until_check; // Counter used to check whether the data structures need a resize or recollect
                            // to avoid checking too often, and incrementing global counters too frequently.

  bool m_is_main_thread = false;

  /// The registered thread aterm pool.
  thread_aterm_pool_interface m_thread_interface; 
};

/// \brief A reference to the thread local term pool storage
thread_aterm_pool& g_thread_term_pool();

} // namespace atermpp::detail


#include "thread_aterm_pool_implementation.h"

#endif // ATERMPP_DETAIL_ATERM_POOL_H
