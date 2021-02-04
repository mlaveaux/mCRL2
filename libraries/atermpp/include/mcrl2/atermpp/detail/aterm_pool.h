// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef ATERMPP_DETAIL_ATERM_POOL_H
#define ATERMPP_DETAIL_ATERM_POOL_H

#include "mcrl2/atermpp/detail/aterm_pool_storage.h"
#include "mcrl2/atermpp/detail/function_symbol_pool.h"

#include <condition_variable>

namespace atermpp
{
namespace detail
{

/// Define several specializations of the term pool storage objects.
using integer_term_storage = aterm_pool_storage<_aterm_int, aterm_int_hasher, aterm_int_equals, 0>;
using term_storage = aterm_pool_storage<_aterm, aterm_hasher_finite<0>, aterm_equals_finite<0>, 0>;
using arbitrary_function_application_storage = aterm_pool_storage<_aterm_appl<1>,
  aterm_hasher<DynamicNumberOfArguments>,
  aterm_equals<DynamicNumberOfArguments>,
  DynamicNumberOfArguments>;

template<std::size_t N>
using function_application_storage = aterm_pool_storage<_aterm_appl<N>, aterm_hasher_finite<N>, aterm_equals_finite<N>, N>;

/// \brief A thread specific aterm pool that provides a local interface to the global term pool.
///        Ensures that terms created by this thread are protected during garbage collection.
class thread_aterm_pool_interface
{
public:
  virtual ~thread_aterm_pool_interface() {}

  /// \brief Mark the terms created by this thread to prevent them being garbage collected.
  virtual void mark() = 0;

  /// \brief Print performance statistics for data stored for this thread.
  virtual void print_local_performance_statistics() const = 0;

  /// \brief Indicates that thread aterm pool is not inserting elements.
  virtual bool busy() const = 0;
};

class thread_aterm_pool;

/// \brief The interface for the term library. Provides the storage of
///        of all classes of terms.
/// \details Internally uses different storage objects to store specific
///          classes of terms. For a given term creation it can decide what
///          storage to use at run-time using its function symbol.
class aterm_pool : public mcrl2::utilities::noncopyable
{
public:
  inline aterm_pool();
  inline ~aterm_pool();

  friend class thread_aterm_pool;

  /// \brief Creates a function symbol pair (name, arity).
  /// \see function_symbol_pool.
  inline function_symbol create_function_symbol(const std::string& name, const std::size_t arity, const bool check_for_registered_functions = false);

  /// \brief Force garbage collection on all storages.
  /// \threadsafe
  inline void collect();

  /// \brief Register a thread specific aterm pool.
  /// \threadsafe
  inline void register_thread_aterm_pool(thread_aterm_pool_interface& pool);

  /// \brief Remove thread specific aterm pool.
  /// \threadsafe
  inline void remove_thread_aterm_pool(thread_aterm_pool_interface& pool);

  /// \brief The number of terms that can be stored without resizing.
  inline std::size_t capacity() const noexcept;

  /// \returns The total number of terms residing in the pool.
  inline std::size_t size() const;

  /// \brief Prints various performance statistics for the term pool.
  inline void print_performance_statistics() const;

  /// \returns A global term that indicates the empty list.
  aterm& empty_list() noexcept { return m_empty_list; }

  /// \returns The function symbol used by integral terms.
  const function_symbol& as_int() noexcept { return m_function_symbol_pool.as_int(); }

  /// \returns The function symbol used by the list constructor.
  const function_symbol& as_list() noexcept { return m_function_symbol_pool.as_list(); }

  /// \returns The function symbol used by the term indicating the empty list.
  const function_symbol& as_empty_list() noexcept { return m_function_symbol_pool.as_empty_list(); }

  /// \brief Add a callback that is triggered whenever a term with the given function symbol is created.
  inline void add_creation_hook(function_symbol sym, term_callback callback);

  /// \brief Add a callback that is triggered whenever a term with the given function symbol is destroyed.
  inline void add_deletion_hook(function_symbol sym, term_callback callback);

  /// \brief Enable garbage collection when passing true and disable otherwise.
  inline void enable_garbage_collection(bool enable) { m_enable_garbage_collection = enable; };

  inline function_symbol_pool& get_symbol_pool() { return m_function_symbol_pool; }

  // These functions of the aterm pool should be called through a thread_aterm_pool.
private:

  /// \brief Triggers garbage collection and resizing when conditions are met.
  /// \param allow_collect Actually perform the garbage collection instead of only updating the counters.
  /// \threadsafe
  inline void created_term(bool allow_collect);

  /// \brief Collect garbage on all storages.
  /// \threadsafe
  inline void collect_impl();

  /// \brief Creates a integral term with the given value.
  inline bool create_int(aterm& term, std::size_t val);

  /// \brief Creates a term with the given function symbol.
  inline bool create_term(aterm& term, const function_symbol& sym);

  /// \brief Creates a function application with the given function symbol and arguments.
  template<class ...Terms>
  inline bool create_appl(aterm& term, const function_symbol& sym, const Terms&... arguments);

  /// \brief Creates a function application with the given function symbol and the arguments
  ///       as provided by the given iterator. This function assumes that the arity of the
  ///       function symbol is equal to the number of elements in the iterator.
  template<typename ForwardIterator>
  bool create_appl_dynamic(aterm& term,
      const function_symbol& sym,
      ForwardIterator begin,
      ForwardIterator end);

  /// \brief Creates a function application with the given function symbol and the arguments
  ///       as provided by the given iterator. This function assumes that the arity of the
  ///       function symbol is equal to the number of elements in the iterator.
  template<typename InputIterator, typename ATermConverter>
  bool create_appl_dynamic(aterm& term,
      const function_symbol& sym,
      ATermConverter convert_to_aterm,
      InputIterator begin,
      InputIterator end);

  /// \brief Resizes all storages if necessary.
  /// \threadsafe.
  inline void resize_if_needed();

  /// \returns True iff the thread aterm pool should call wait().
  inline bool should_wait();

  /// \brief Wait for the term pool to finish garbage collection.
  inline void wait();

  /// \brief Prevent any thread aterm pool from creating or retrieving terms.
  inline void halt();

  /// \brief Allow thread pools to resume their threads.
  inline void resume();

  /// \brief The set of local aterm pools.
  std::vector<thread_aterm_pool_interface*> m_thread_pools;

  /// \brief Mutex for adding/removing local pools in m_thread_pools.
  std::mutex m_mutex;

  /// \brief Storage for the function symbols.
  function_symbol_pool m_function_symbol_pool;

  /// \brief Storage for integral terms.
  integer_term_storage m_int_storage;

  /// Storage for function applications with a fixed number of arguments.
  std::tuple<
    term_storage,
    function_application_storage<1>,
    function_application_storage<2>,
    function_application_storage<3>,
    function_application_storage<4>,
    function_application_storage<5>,
    function_application_storage<6>,
    function_application_storage<7>
  > m_appl_storage;

  /// Storage for term_appl with a dynamic number of arguments larger than 7.
  arbitrary_function_application_storage m_appl_dynamic_storage;

  /// Track the number of terms destroyed and reduce the freelist.
  std::atomic<long> m_count_until_collection = 0;
  std::atomic<long> m_count_until_resize = 0;

  std::atomic<bool> m_guard = false; /// Instructs local thread pools to wait.
  std::atomic<bool> m_enable_garbage_collection = EnableGarbageCollection; /// Garbage collection is enabled.

  /// A conditional variable to avoid busy waiting on the guard.
  std::condition_variable m_wait_variable;

  /// Represents an empty list.
  aterm m_empty_list;
};

} // namespace detail
} // namespace atermpp

#include "aterm_pool_implementation.h"
#include "aterm_pool_storage_implementation.h"

#endif // ATERMPP_DETAIL_ATERM_POOL_H
