// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_ATERM_POOL_STORAGE_IMPLEMENTION_H
#define MCRL2_ATERMPP_ATERM_POOL_STORAGE_IMPLEMENTION_H
#pragma once

#include "mcrl2/atermpp/detail/aterm_pool.h"

#include "mcrl2/utilities/stack_array.h"

#include <cstring>

namespace atermpp
{
namespace detail
{

/// \brief Construct the proxy where its arguments are given by applying the converter
///        to each element in the iterator.
template<std::size_t N,
         typename InputIterator,
         typename TermConverter,
         typename std::enable_if<mcrl2::utilities::is_iterator<InputIterator>::value, void>::type* = nullptr>
inline std::array<unprotected_aterm, N> construct_arguments(InputIterator it, InputIterator end, TermConverter converter)
{
  // The end is only used for debugging to ensure that the arity and std::distance(it, end) match.
  mcrl2::utilities::mcrl2_unused(end);

  // Copy the arguments into this array. Doesn't change any reference count, because they are unprotected terms.
  std::array<unprotected_aterm, N> arguments;
  for (size_t i = 0; i < N; ++i)
  {
    assert(it != end);
    arguments[i] = converter(*it);
    ++it;
  }
  assert(it == end);

  return arguments;
}

#define ATERM_POOL_STORAGE_TEMPLATES template<typename Element, typename Hash, typename Equals, std::size_t N, bool ThreadSafe>
#define ATERM_POOL_STORAGE aterm_pool_storage<Element, Hash, Equals, N, ThreadSafe>

ATERM_POOL_STORAGE_TEMPLATES
ATERM_POOL_STORAGE::aterm_pool_storage(aterm_pool& pool) :
  m_pool(pool),
  m_term_set(1 << 14)
{}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::add_creation_hook(function_symbol sym, term_callback callback)
{
  // The code handling the hooks is currently assuming that there is only one
  // hook per function symbol. If more hooks are allowed, then this code
  // should be changed.
  for (const auto& hook : m_creation_hooks)
  {
    mcrl2::utilities::mcrl2_unused(hook);
    assert(hook.first != sym);
  }
  m_creation_hooks.emplace_back(sym, callback);
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::add_deletion_hook(function_symbol sym, term_callback callback)
{
  // See the comments at add_creation_hook.
  for (const auto& hook : m_deletion_hooks)
  {
    mcrl2::utilities::mcrl2_unused(hook);
    assert(hook.first != sym);
  }
  m_deletion_hooks.emplace_back(sym, callback);
}

ATERM_POOL_STORAGE_TEMPLATES
aterm ATERM_POOL_STORAGE::create_int(std::size_t value)
{
  return emplace(value);
}

ATERM_POOL_STORAGE_TEMPLATES
aterm ATERM_POOL_STORAGE::create_term(const function_symbol& symbol)
{
  assert(symbol.arity() == 0);
  return emplace(symbol);
}

ATERM_POOL_STORAGE_TEMPLATES
template<class ...Terms>
aterm ATERM_POOL_STORAGE::create_appl(const function_symbol& symbol, const Terms&... arguments)
{
  assert(symbol.arity() == sizeof...(arguments));
  return emplace(symbol, arguments...);
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename ForwardIterator>
aterm ATERM_POOL_STORAGE::create_appl_iterator(const function_symbol& symbol,
                                        ForwardIterator begin,
                                        ForwardIterator end)
{
  return emplace(symbol, begin, end);
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename InputIterator, typename TermConverter>
aterm ATERM_POOL_STORAGE::create_appl_iterator(const function_symbol& symbol,
                                        TermConverter converter,
                                        InputIterator begin,
                                        InputIterator end)
{
  std::array<unprotected_aterm, N> arguments = construct_arguments<N>(begin, end, converter);
  return emplace(symbol, arguments);
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename ForwardIterator>
aterm ATERM_POOL_STORAGE::create_appl_dynamic(const function_symbol& symbol,
                                        ForwardIterator begin,
                                        ForwardIterator end)
{
  return emplace(symbol, begin, end);
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename InputIterator,
          typename TermConverter>
aterm ATERM_POOL_STORAGE::create_appl_dynamic(const function_symbol& symbol,
                                        TermConverter converter,
                                        InputIterator it,
                                        InputIterator end)
{
  // The end is only used for debugging to ensure that the arity and std::distance(it, end) match.
  mcrl2::utilities::mcrl2_unused(end);

  MCRL2_DECLARE_STACK_ARRAY(arguments, unprotected_aterm, symbol.arity());
  for (std::size_t i = 0; i < symbol.arity(); ++i)
  {
    assert(it != end);
    arguments[i] = converter(*it);
    ++it;
  }
  assert(it == end);

  // Find or create a new term and return it.
  return emplace(symbol, arguments.begin(), arguments.end());
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::print_performance_stats(const char* identifier) const
{
  if constexpr (EnableHashtableMetrics)
  {
    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool_storage(" << identifier << ") hashtable:\n";
    print_performance_statistics(m_term_set);
  }

  if (EnableGarbageCollectionMetrics && m_erased_blocks > 0)
  {
    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool_storage(" << identifier << "): Consolidate removed " << m_erased_blocks << " blocks.\n";
  }

  if constexpr (EnableCreationMetrics)
  {
    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool_storage(" << identifier << "): emplace() " << m_term_metric.message() << ".\n";
  }
}


ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::mark()
{
  for (const Element& term : m_term_set)
  {
    // If a term is marked its arguments have been marked as well.
    if (term.is_reachable() && !term.is_marked())
    {
      // Mark all terms (and their subterms) that are reachable, i.e the root set.
      mark_term(term);
    }
  }
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::sweep()
{
  // Iterate over all terms and removes the ones that are marked.
  for (auto it = m_term_set.begin(); it != m_term_set.end(); )
  {
    const Element& term = *it;

    if (!term.is_reachable())
    {
      it = destroy(it);
    }
    else
    {
      // Reset terms that have been marked.
      if (term.is_marked())
      {
        term.reset();
      }
      ++it;
    }
  }

  // Clean up unnecessary blocks.
  m_erased_blocks = m_term_set.get_allocator().consolidate();
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::resize_if_needed()
{
  m_term_set.rehash_if_needed();
}

/// PRIVATE FUNCTIONS

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::call_creation_hook(unprotected_aterm term)
{
  for (const auto& [symbol, callback] : m_creation_hooks)
  {
    if (symbol == term.function())
    {
      callback(static_cast<const aterm&>(term));
    }
  }
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::call_deletion_hook(unprotected_aterm term)
{
  for (const auto& [symbol, callback] : m_deletion_hooks)
  {
    if (symbol == term.function())
    {
      callback(static_cast<const aterm&>(term));
    }
  }
}


ATERM_POOL_STORAGE_TEMPLATES
bool ATERM_POOL_STORAGE::verify_mark()
{
  // Check for consistency that if a term is reachable its arguments are as well.
  for (const Element& term : m_term_set)
  {
    if (term.is_reachable() && term.function().arity() > 0)
    {
       const _term_appl& ta = static_cast<_term_appl&>(const_cast<Element&>(term));
       for (std::size_t i = 0; i < ta.function().arity(); ++i)
       {
         assert(detail::address(ta.arg(i))->is_reachable());
       }
    }
  }
  return true;
}

ATERM_POOL_STORAGE_TEMPLATES
bool ATERM_POOL_STORAGE::verify_sweep()
{
  // Check that no argument was removed from a reachable term.
  for (const Element& term : m_term_set)
  {
    (void)term;
    assert(verify_term(term));
  }
  return true;
}

/// Private definitions

ATERM_POOL_STORAGE_TEMPLATES
typename ATERM_POOL_STORAGE::iterator ATERM_POOL_STORAGE::destroy(iterator it)
{
  // Store the term temporarily to be able to deallocate it after removing it from the set.
  const Element& term = *it;
  assert(!term.is_reachable());

  // Trigger the deletion hook before the term is actually destroyed.
  call_deletion_hook(&term);

  // Remove them from the hash table, will also destroy terms with fixed arity.
  return m_term_set.erase(it);
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename ...Args>
aterm ATERM_POOL_STORAGE::emplace(Args&&... args)
{
  auto [it, added] = m_term_set.emplace(std::forward<Args>(args)...);

  aterm term(&(*it));
  if (added)
  {
    // A new term was created
    if constexpr (EnableCreationMetrics) { m_term_metric.miss(); }
    m_pool.trigger_collection();
    call_creation_hook(term);
  }
  else if constexpr (EnableCreationMetrics)
  {
    // A term was already found in the set.
    m_term_metric.hit();
  }

  return term;
}

ATERM_POOL_STORAGE_TEMPLATES
constexpr bool ATERM_POOL_STORAGE::is_dynamic_storage() const
{
  return N == DynamicNumberOfArguments;
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::mark_term(const _aterm& root)
{
  // Do not use the stack, because this might run out of stack memory for large lists.
  todo.push(const_cast<_aterm&>(root));

  // Mark the term depth-first to reduce the maximum todo size required.
  while (!todo.empty())
  {
    _aterm& term = todo.top();
    todo.pop();

    // (Statically) determine the arity of the function application.
    const std::size_t arity = term.function().arity();
    _term_appl& term_appl = static_cast<_term_appl&>(term);

    for (std::size_t i = 0; i < arity; ++i)
    {
      // Marks all arguments that are not already (marked as) reachable, because the current
      // term is reachable and as such its arguments are reachable as well.
      _aterm& argument = *detail::address(term_appl.arg(i));
      if (!argument.is_reachable())
      {
        // Only mark this term if it current is unreachable, because the marking is applied to
        // the reference counter.
        argument.mark();

        // Add the argument to be explored as well.
        todo.push(argument);
      }
    }

    // Each argument should be reachable.
    assert(term.is_reachable());
  }
}

ATERM_POOL_STORAGE_TEMPLATES
template<std::size_t Arity>
bool ATERM_POOL_STORAGE::verify_term(const _aterm& term)
{
  // Check that a valid function symbol was used and that its arity belongs to this storage.
  assert(term.function().defined());
  assert(Arity == DynamicNumberOfArguments ? true : term.function().arity() == N);

  // Check that all of its arguments are defined.
  if (term.function().arity() > 0)
  {
    const _term_appl& ta = static_cast<const _term_appl&>(term);
    for (std::size_t i = 0; i < ta.function().arity(); ++i)
    {
      assert(ta.arg(i).defined());
      verify_term<DynamicNumberOfArguments>(*detail::address(ta.arg(0)));
    }
  }
  return true;
}

#undef ATERM_POOL_STORAGE_TEMPLATES
#undef ATERM_POOL_STORAGE

} // namespace detail
} // namespace atermpp

#endif
