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

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/shared_mutex.h"
#else
  import utilities;
#endif

namespace atermpp::detail
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

// There are some annoying circular dependencies between a aterm_pool and the contained thread_aterm_pool_interfaces
class aterm_pool;

/// \brief A thread specific aterm pool that provides a local interface to the global term pool.
///        Ensures that terms created by this thread are protected during garbage collection.
class thread_aterm_pool_interface final
{
public:
  thread_aterm_pool_interface(aterm_pool& pool, std::function<void()> mark_function, std::function<void()> print_function, std::function<std::size_t()> protection_set_size_function);
  ~thread_aterm_pool_interface();

  /// \brief Mark the terms created by this thread to prevent them being garbage collected.
  void mark()
  {
    m_mark_function();
  }

  /// \brief Print performance statistics for data stored for this thread.
  void print_local_performance_statistics() const
  {
    m_print_function();
  }

  /// \returns The total number of terms residing in the pool.
  std::size_t protection_set_size() const 
  {
    return m_protection_set_size_function();
  }

  // Prematurely unregister the thread_aterm_pool_interface.
  void unregister();

private:
  aterm_pool& m_pool;
  std::function<void()> m_mark_function;
  std::function<void()> m_print_function;
  std::function<std::size_t()> m_protection_set_size_function;
  bool m_registered = true;
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

  /// \brief Creates a function symbol pair (name, arity).
  /// \see function_symbol_pool.
  inline function_symbol create_function_symbol(std::string&& name, const std::size_t arity, const bool check_for_registered_functions = false);

  /// \brief Register a thread specific aterm pool.
  /// \details threadsafe
  inline void register_thread_aterm_pool(thread_aterm_pool_interface& pool);

  /// \brief Remove thread specific aterm pool.
  /// \details threadsafe
  inline void remove_thread_aterm_pool(thread_aterm_pool_interface& pool);

  /// \brief The number of terms that can be stored without resizing.
  inline std::size_t capacity() const noexcept;

  /// \returns The total number of terms residing in the pool.
  inline std::size_t size() const;

  /// \brief Prints various performance statistics for the term pool.
  inline void print_performance_statistics() const;

  /// \returns A global term that indicates the empty list.
  aterm& empty_list() noexcept { return reinterpret_cast<aterm&>(m_empty_list); }  // TODO remove this reinterpret cast by letting m_empty_list become an aterm.

  /// \returns The function symbol used by integral terms.
  const function_symbol& as_int() noexcept { return m_function_symbol_pool.as_int(); }

  /// \returns The function symbol used by the list constructor.
  const function_symbol& as_list() noexcept { return m_function_symbol_pool.as_list(); }

  /// \returns The function symbol used by the term indicating the empty list.
  const function_symbol& as_empty_list() noexcept { return m_function_symbol_pool.as_empty_list(); }

  /// \brief Add a callback that is triggered whenever a term with the given function symbol is destroyed.
  inline void add_deletion_hook(function_symbol sym, term_callback callback);

  /// \brief Enable garbage collection when passing true and disable otherwise.
  inline void enable_garbage_collection(bool enable) { m_enable_garbage_collection = enable; };

  inline function_symbol_pool& get_symbol_pool() { return m_function_symbol_pool; }

  // These functions of the aterm pool should be called through a thread_aterm_pool.
private:
  /// \brief Force garbage collection on all storages.
  /// \details threadsafe
  inline void collect(mcrl2::utilities::shared_mutex& mutex);

  /// \brief Triggers garbage collection and resizing when conditions are met.
  /// \param allow_collect Actually perform the garbage collection instead of only updating the counters.
  /// \param mutex The shared mutex that should be used for locking if necessary.
  /// \details threadsafe
  inline void created_term(bool allow_collect, mcrl2::utilities::shared_mutex& mutex);

  /// \brief Collect garbage on all storages.
  /// \details threadsafe
  inline void collect_impl(mcrl2::utilities::shared_mutex& mutex);

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
  /// \details threadsafe
  inline void resize_if_needed(mcrl2::utilities::shared_mutex& shared);

  mcrl2::utilities::shared_mutex& shared_mutex() { return m_shared_mutex; }

  /// \returns The total number of term variables residing in the protection sets.
  inline std::size_t protection_set_size() const;

  /// \brief The set of local aterm pools.
  std::vector<thread_aterm_pool_interface* > m_thread_pools;

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

  std::atomic<bool> m_enable_garbage_collection = EnableGarbageCollection; /// Garbage collection is enabled.

  /// All the shared mutexes.
  mcrl2::utilities::shared_mutex m_shared_mutex;

  /// Represents an empty list.
  aterm_core m_empty_list;
};

inline
thread_aterm_pool_interface::thread_aterm_pool_interface(aterm_pool& pool, std::function<void()> mark_function, std::function<void()> print_function, std::function<std::size_t()> protection_set_size_function)
  : m_pool(pool), m_mark_function(mark_function), m_print_function(print_function), m_protection_set_size_function(protection_set_size_function)
{
  m_pool.register_thread_aterm_pool(*this);
}

inline
void thread_aterm_pool_interface::unregister()
{
  if (m_registered)
  {
    m_pool.remove_thread_aterm_pool(*this);
    m_registered = false;
  }
}

inline
thread_aterm_pool_interface::~thread_aterm_pool_interface()
{
  unregister();
}

} // namespace atermpp::detail


#include "aterm_pool_implementation.h"
#include "aterm_pool_storage_implementation.h"

#endif // ATERMPP_DETAIL_ATERM_POOL_H
