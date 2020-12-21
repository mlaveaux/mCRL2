// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef ATERMPP_DETAIL_ATERM_POOL_STORAGE_H
#define ATERMPP_DETAIL_ATERM_POOL_STORAGE_H

#include "mcrl2/atermpp/detail/aterm_hash.h"
#include "mcrl2/utilities/cache_metric.h"
#include "mcrl2/utilities/unordered_set.h"

#include <stack>
#include <utility>

namespace atermpp
{
namespace detail
{

// Forward declaration
class aterm_pool;

/// \brief Marks a term and recursively all arguments that are not reachable.
inline void mark_term(const _aterm& root, std::stack<std::reference_wrapper<_aterm>>& todo);

/// \brief This class provides for all types of term storage. It also
///       provides garbage collection via its mark and sweep functions.
/// \details Internally a hash set is used to ensure that the created terms are unique.
template<typename Element,
         typename Hash = aterm_hasher<>,
         typename Equals = aterm_equals<>,
         std::size_t N = DynamicNumberOfArguments>
class aterm_pool_storage : private mcrl2::utilities::noncopyable
{
public:
  using unordered_set = mcrl2::utilities::unordered_set<
    Element,
    Hash,
    Equals,
    typename std::conditional<N == DynamicNumberOfArguments,
      atermpp::detail::_aterm_appl_allocator<>,
      mcrl2::utilities::block_allocator<Element, 1024, GlobalThreadSafe>>::type,
    GlobalThreadSafe,
    false>;
  using iterator = typename unordered_set::iterator;
  using const_iterator = typename unordered_set::const_iterator;

  /// \brief The local pool is a friend class so it can mark terms.
  friend class aterm_pool;

  aterm_pool_storage(aterm_pool& pool);

  /// \brief Add a callback that is triggered whenever a term with the given function symbol is created.
  void add_creation_hook(function_symbol sym, term_callback callback);

  /// \brief Add a callback that is triggered whenever a term with the given function symbol is destroyed.
  void add_deletion_hook(function_symbol sym, term_callback callback);

  /// \returns The total number of terms that can be stored without resizing.
  std::size_t capacity() const noexcept { return m_term_set.capacity(); }

  /// \brief Creates a integral term with the given value.
  bool create_int(aterm& term, std::size_t value);

  /// \brief Creates a term with the given function symbol.
  bool create_term(aterm& term, const function_symbol& sym);

  /// \brief Creates a function application with the given function symbol and arguments.
  template<class ...Terms>
  bool create_appl(aterm& term, const function_symbol& sym, const Terms&... arguments);

  /// \brief Creates a function application with the given function symbol and the arguments
  ///       as provided by the given iterator. This function assumes that the arity of the
  ///       function symbol is equal to N and that the iterator has exactly N elements.
  template<typename ForwardIterator>
  bool create_appl_iterator(aterm& term,
                            const function_symbol& sym,
                            ForwardIterator begin,
                            ForwardIterator end);

  /// \brief Creates a function application with the given function symbol and the arguments
  ///       as provided by the given iterator. This function assumes that the arity of the
  ///       function symbol is equal to N and that the iterator has exactly N elements.
  template<typename InputIterator,
           typename TermConverter>
  bool create_appl_iterator(aterm& term,
                            const function_symbol& sym,
                            TermConverter converter,
                            InputIterator begin,
                            InputIterator end);

  /// \brief Creates a function application with the given function symbol and the arguments
  ///        as provided by the given iterator.
  template<typename ForwardIterator>
  bool create_appl_dynamic(aterm& term,
                           const function_symbol& sym,
                           ForwardIterator begin,
                           ForwardIterator end);

  /// \brief Creates a function application with the given function symbol and the arguments
  ///        as provided by the given iterator, but the converter is applied to each argument.
  template<typename InputIterator,
           typename TermConverter>
  bool create_appl_dynamic(aterm& term,
                           const function_symbol& sym,
                           TermConverter converter,
                           InputIterator begin,
                           InputIterator end);

  /// \brief Prints various performance statistics for this storage.
  /// \param identifier A string to identify the printed message for this storage.
  void print_performance_stats(const char* identifier) const;

  /// \brief Marks all terms that are reachable and should not be destroyed.
  void mark();

  /// \brief sweep Destroys all terms that are not reachable. Requires that
  ///        mark() was called first.
  void sweep();

  /// \brief Resizes the hash table if necessary.
  void resize_if_needed();

  /// \returns The number of terms stored in this storage.
  std::size_t size() const { return m_term_set.size(); }

  /// \brief Check that all arguments of a term application are marked properly.
  bool verify_mark();

  /// \brief Check that all arguments of a term application are marked properly.
  bool verify_sweep();

private:
  using callback_pair = std::pair<function_symbol, term_callback>;

  /// \brief Calls the creation hook attached to the function symbol of this term.
  void call_creation_hook(unprotected_aterm term);

  /// \brief Calls the deletion hook attached to the function symbol of this term.
  void call_deletion_hook(unprotected_aterm term);

  /// \brief Removes an element from the unordered set and deallocates it.
  iterator destroy(iterator it);

  /// \brief Inserts a term constructed by the given arguments, checks for existing term.
  template<typename ...Args>
  bool emplace(aterm& term, Args&&... args);

  /// \returns True if and only if this term storage can store term applications with a dynamic
  ///          number of arguments.
  constexpr bool is_dynamic_storage() const;

  /// \brief Verify that the given term was constructed properly.
  template<std::size_t Arity = N>
  bool verify_term(const _aterm& term);

  /// The pool that this storage belongs to.
  aterm_pool& m_pool;

  /// This is the set of term pointers to keep the terms unique.
  unordered_set m_term_set;

  /// This array stores creation, resp deletion, hooks for function symbols.
  std::vector<callback_pair> m_creation_hooks;
  std::vector<callback_pair> m_deletion_hooks;

  /// A reusable todo stack.
  std::stack<std::reference_wrapper<_aterm>> todo;

  // Various performance statistics.

  mcrl2::utilities::cache_metric m_term_metric; ///< Count the number of times a term has been found in or is added to the set.
  std::size_t m_erased_blocks = 0; /// The number of blocks that have been erased in the block allocator.
};

} // namespace detail
} // namespace atermpp


#endif // ATERMPP_DETAIL_ATERM_POOL_STORAGE_H
