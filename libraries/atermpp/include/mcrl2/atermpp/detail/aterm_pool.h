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

/// \brief The interface for the term library. Provides the storage of
///        of all classes of terms.
/// \details Internally uses different storage objects to store specific
///          classes of terms. For a given term creation it can decide what
///          storage to use at run-time using its function symbol.
class aterm_pool : public mcrl2::utilities::noncopyable
{
public:

  /// \brief Should be able to call mark() from any storage.
  /// \todo Make mark() private and change enable this friend class.
  template<typename Element, typename Hash, typename Equals, std::size_t N>
  friend class aterm_pool_storage;

  inline aterm_pool();
  inline ~aterm_pool();

  /// \brief Add a callback that is triggered whenever a term with the given function symbol is created.
  inline void add_creation_hook(function_symbol sym, term_callback callback);

  /// \brief Add a callback that is triggered whenever a term with the given function symbol is destroyed.
  inline void add_deletion_hook(function_symbol sym, term_callback callback);

  /// \brief The number of terms that can be stored without resizing.
  inline std::size_t capacity() const noexcept;

  /// \brief Triggers garbage collection on all storages.
  inline void collect();

  /// \brief Enable garbage collection when passing true and disable otherwise.
  inline void enable_garbage_collection(bool enable);

  /// \brief Creates a integral term with the given value.
  inline void create_int(aterm& term, std::size_t val);

  /// \brief Creates a term with the given function symbol.
  inline void create_term(aterm& term, const function_symbol& sym);

  /// \brief Creates a function application with the given function symbol and arguments.
  template<class ...Terms>
  inline void create_appl(aterm& term, const function_symbol& sym, const Terms&... arguments);

  /// \brief Creates a function application with the given function symbol and the arguments
  ///       as provided by the given iterator. This function assumes that the arity of the
  ///       function symbol is equal to the number of elements in the iterator.
  template<typename ForwardIterator>
  inline void create_appl_dynamic(aterm& term,
                              const function_symbol& sym,
                              ForwardIterator begin,
                              ForwardIterator end);

  /// \brief Creates a function application with the given function symbol and the arguments
  ///       as provided by the given iterator. This function assumes that the arity of the
  ///       function symbol is equal to the number of elements in the iterator.
  template<typename InputIterator, typename ATermConverter>
  inline void create_appl_dynamic(aterm& term,
                              const function_symbol& sym,
                              ATermConverter convert_to_aterm,
                              InputIterator begin,
                              InputIterator end);

  /// \brief Prints various performance statistics for the term pool.
  inline void print_performance_statistics() const;

  /// \returns The total number of terms residing in the pool.
  inline std::size_t size() const;

  /// \returns A global term that indicates the empty list.
  aterm& empty_list() noexcept { return m_empty_list; }

  /// \returns The function symbol used by integral terms.
  const function_symbol& as_int() noexcept { return m_function_symbol_pool.as_int(); }

  /// \returns The function symbol used by the list constructor.
  const function_symbol& as_list() noexcept { return m_function_symbol_pool.as_list(); }

  /// \returns The function symbol used by the term indicating the empty list.
  const function_symbol& as_empty_list() noexcept { return m_function_symbol_pool.as_empty_list(); }

  /// \returns The pool of function symbols.
  function_symbol_pool& get_symbol_pool() { return m_function_symbol_pool; }

private:
  /// \brief Triggers garbage collection when certain conditions are met.
  inline void trigger_collection();

  /// \brief Resizes the hash tables if necessary.
  inline void resize_if_needed();

  /// Storage for the function symbols.
  function_symbol_pool m_function_symbol_pool;

  /// Storage for integral terms.
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

  /// Counter for number of terms added before garbage collection is triggered.
  std::size_t m_count_until_collection = 0;
  std::size_t m_count_until_resize = 0;

  /// It can happen that during create_appl with converter the converter generates new terms.
  /// As such these terms might only be protected after the term_appl was actually created.
  std::size_t m_creation_depth = 0;

  /// Defer garbage collection until the creation depth is equal to zero again.
  bool m_deferred_garbage_collection = false;

  /// Enable automatically triggered garbage collection.
  bool m_enable_garbage_collection = EnableGarbageCollection;

  /// Represents an empty list.
  aterm m_empty_list;
};

} // namespace detail
} // namespace atermpp

#include "aterm_pool_implementation.h"
#include "aterm_pool_storage_implementation.h"

#endif // ATERMPP_DETAIL_ATERM_POOL_H
