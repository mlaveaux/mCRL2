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

#include "mcrl2/utilities/stack_array.h"
#include "mcrl2/atermpp/detail/aterm_pool.h"

#include "mcrl2/utilities/logger.h"

#include <type_traits>
#include <cstring>

namespace atermpp::detail
{

/// \brief Construct the proxy where its arguments are given by applying the converter
///        to each element in the iterator.
template <std::size_t N,
    typename InputIterator,
    typename TermConverter,
    std::enable_if_t<mcrl2::utilities::is_iterator<InputIterator>::value, void>* = nullptr>
inline std::array<unprotected_aterm_core, N>
construct_arguments(InputIterator it, [[maybe_unused]] InputIterator end, TermConverter converter)
  requires std::is_convertible_v<std::invoke_result_t<TermConverter, typename InputIterator::value_type>, aterm>
{
  // Copy the arguments into this array.
  std::array<unprotected_aterm_core, N> arguments;
  for (size_t i = 0; i < N; ++i)
  {
    assert(it != end);
    arguments[i] = converter(*it);
 
    ++it;
  }
  assert(it == end);
 
  return arguments;
}

/// \brief Construct the proxy where its arguments are given by applying the converter
///        to each element in the iterator.
template <std::size_t N,
    typename InputIterator,
    typename TermConverter,
    std::enable_if_t<mcrl2::utilities::is_iterator<InputIterator>::value, void>* = nullptr>
inline std::array<unprotected_aterm_core, N>
construct_arguments(InputIterator it, [[maybe_unused]] InputIterator end, TermConverter converter)
  requires std::is_same_v<
      std::invoke_result_t<TermConverter, typename InputIterator::value_type&, typename InputIterator::value_type>,
      void>
{
  // Copy the arguments into this array. Doesn't change any reference count, because they are unprotected terms.
  std::array<unprotected_aterm_core, N> arguments;
  for (size_t i = 0; i < N; ++i)
  {
    assert(it != end);
    typename InputIterator::value_type& t= *reinterpret_cast<typename InputIterator::value_type*>(&arguments[i]);
    converter(t, *it);
    ++it;
  }
  assert(it == end);
 
  return arguments;
}

void mark_term(const _aterm& root, std::stack<std::reference_wrapper<_aterm>>& todo)
{
  if (!root.is_marked())
  {
    // Do not use the stack, because this might run out of stack memory for large lists.
    todo.emplace(const_cast<_aterm&>(root));

    // Mark the term depth-first to reduce the maximum todo size required.
    while (!todo.empty())
    {
      _aterm& term = todo.top();
      todo.pop();

      // Each term should be marked.
      term.mark();
      // Determine the arity of the function application.
      const std::size_t arity = term.function().arity();
      _term_appl& term_appl = static_cast<_term_appl&>(term);

      for (std::size_t i = 0; i < arity; ++i)
      {
        // Marks all arguments that are not already (marked as) reachable, because the current
        // term is reachable and as such its arguments are reachable as well.
        _aterm& argument = *detail::address(term_appl.arg(i));
        if (!argument.is_marked())
        {
          argument.mark();

          // Add the argument to be explored as well.
          todo.emplace(argument);
        }
      }

    }
  }
}

#define ATERM_POOL_STORAGE_TEMPLATES template<typename Element, typename Hash, typename Equals, std::size_t N>
#define ATERM_POOL_STORAGE aterm_pool_storage<Element, Hash, Equals, N>

ATERM_POOL_STORAGE_TEMPLATES
ATERM_POOL_STORAGE::aterm_pool_storage(aterm_pool& pool) :
  m_pool(pool),
  m_term_set(1 << 14)
{}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::add_deletion_hook(function_symbol sym, term_callback callback)
{
  for ([[maybe_unused]] 
    const auto& hook : m_deletion_hooks)
  {
    assert(hook.first != sym);
  }
  m_deletion_hooks.emplace_back(sym, callback);
}

ATERM_POOL_STORAGE_TEMPLATES
bool ATERM_POOL_STORAGE::create_int(aterm& term, std::size_t value)
{
  return emplace(reinterpret_cast<aterm_core&>(term), value);   // TODO: remove reinterpret_cast
}

ATERM_POOL_STORAGE_TEMPLATES
bool ATERM_POOL_STORAGE::create_term(aterm& term, const function_symbol& symbol)
{
  assert(symbol.arity() == 0);
  return emplace(reinterpret_cast<aterm_core&>(term), symbol);  // TODO: remove reinterpret_cast
}

template <std::size_t N>
void store_in_argument_array_(std::size_t , std::array<unprotected_aterm_core, N>& )
{}

template <std::size_t N, class FUNCTION_OR_TERM_TYPE, typename... Args>
void store_in_argument_array_(std::size_t i, 
                              std::array<unprotected_aterm_core, N>& argument_array, 
                              FUNCTION_OR_TERM_TYPE& function_or_term, 
                              const Args&... args)
{
  if constexpr (std::is_convertible_v<FUNCTION_OR_TERM_TYPE, atermpp::aterm>)
  {
    argument_array[i]=function_or_term;
  }
  // check whether the function_or_term invoked on an empty argument yields an aterm.
  else if constexpr (mcrl2::utilities::is_applicable< FUNCTION_OR_TERM_TYPE, void>::value)
  {
    argument_array[i]=function_or_term();
  }
  // Otherwise function_or_term is supposed to  have type void(term& result), putting the term in result. 
  else
  {
    // function_or_term(static_cast<Term&>(argument_array[i]));

    using traits = mcrl2::utilities::function_traits<decltype(&FUNCTION_OR_TERM_TYPE::operator())>;
    function_or_term(static_cast<typename traits::template arg<0>::type&>(argument_array[i]));
  }
  store_in_argument_array_(i+1, argument_array, args...);
}

template <std::size_t N, typename... Args>
void store_in_argument_array(std::array<unprotected_aterm_core, N>& argument_array,
                             const Args&... args)
{
  store_in_argument_array_(0, argument_array, args...);
}



ATERM_POOL_STORAGE_TEMPLATES
template<class ...Terms>
bool ATERM_POOL_STORAGE::create_appl(aterm& term, const function_symbol& symbol, const Terms&... arguments)
{
  assert(symbol.arity() == sizeof...(arguments));
  if constexpr (detail::are_terms<Terms...>::value)
  {
    return emplace(reinterpret_cast<aterm_core&>(term), symbol, arguments...); // TODO remove reinterpret_cast.
  }
  else 
  {
    std::array<unprotected_aterm_core, N> argument_array;
    
    // Evaluate the functions or terms and put the result in "argument_array".
    store_in_argument_array(argument_array, arguments...);

    return emplace(term, symbol, argument_array.begin(), argument_array.end());
  }
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename ForwardIterator>
bool ATERM_POOL_STORAGE::create_appl_iterator(aterm& term,
                                        const function_symbol& symbol,
                                        ForwardIterator begin,
                                        ForwardIterator end)
{
  return emplace(term, symbol, begin, end);
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename InputIterator, typename TermConverter>
bool ATERM_POOL_STORAGE::create_appl_iterator(aterm& term,
                                        const function_symbol& symbol,
                                        TermConverter converter,
                                        InputIterator begin,
                                        InputIterator end)
{
  std::array<unprotected_aterm_core, N> arguments = construct_arguments<N>(begin, end, converter);
  return emplace(term, symbol, arguments);
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename ForwardIterator>
bool ATERM_POOL_STORAGE::create_appl_dynamic(aterm& term,
                                        const function_symbol& symbol,
                                        ForwardIterator begin,
                                        ForwardIterator end)
{
  return emplace(term, symbol, begin, end);
}

ATERM_POOL_STORAGE_TEMPLATES
template <typename InputIterator,
    typename TermConverter>
bool ATERM_POOL_STORAGE::create_appl_dynamic(aterm& term,
    const function_symbol& symbol,
    TermConverter converter,
    InputIterator it,
    [[maybe_unused]] // The end is only used for debugging to ensure that the arity and std::distance(it, end) match.
    InputIterator end)
  requires std::is_convertible_v<std::invoke_result_t<TermConverter, typename InputIterator::value_type>, aterm>
{
  MCRL2_DECLARE_STACK_ARRAY(arguments, unprotected_aterm_core, symbol.arity());
  for (std::size_t i = 0; i < symbol.arity(); ++i)
  {
    assert(it != end);
    arguments[i] = converter(*it);
    ++it;
  }
  assert(it == end);
 
  // Find or create a new term and return it.
  return emplace(term, symbol, arguments.begin(), arguments.end());
}

ATERM_POOL_STORAGE_TEMPLATES
template <typename InputIterator,
    typename TermConverter>
bool ATERM_POOL_STORAGE::create_appl_dynamic(aterm& term,
    const function_symbol& symbol,
    TermConverter converter,
    InputIterator it,
    [[maybe_unused]] // The end is only used for debugging to ensure that the arity and std::distance(it, end) match.
    InputIterator end)
  requires std::is_same_v<
      std::invoke_result_t<TermConverter, typename InputIterator::value_type&, typename InputIterator::value_type>,
      void>
{
  MCRL2_DECLARE_STACK_ARRAY(arguments, unprotected_aterm_core, symbol.arity());
  for (std::size_t i = 0; i < symbol.arity(); ++i)
  {
    assert(it != end);
    // The construction below is possible as with protection sets the new term simply
    // overwrites the object at arguments[i].
    converter(reinterpret_cast<typename InputIterator::value_type&>(arguments[i]), *it);
    ++it;
  }
  assert(it == end);
 
  // Find or create a new term and return it.
  return emplace(term, symbol, arguments.begin(), arguments.end());
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::print_performance_stats(const char* identifier) const
{
  if (EnableHashtableMetrics)
  {
    mCRL2log(mcrl2::log::info) << "g_term_pool(" << identifier << ") hashtable:\n";
    print_performance_statistics(m_term_set);
  }

  if (EnableGarbageCollectionMetrics && m_erasedBlocks > 0)
  {
    mCRL2log(mcrl2::log::info) << "g_term_pool(" << identifier << "): Consolidate removed " << m_erasedBlocks << " blocks.\n";
  }

  if (EnableCreationMetrics)
  {
    mCRL2log(mcrl2::log::info) << "g_term_pool(" << identifier << "): emplace() " << m_term_metric.message() << ".\n";
  }
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::sweep()
{
  // Iterate over all terms and removes the ones that are marked.
  for (auto it = m_term_set.begin(); it != m_term_set.end(); )
  {
    const Element& term = *it;

    if (!term.is_marked())
    {
      // For constants, i.e., arity zero and integer terms we do not mark, but use their reachability directly. 
      it = destroy(it);
    }
    else
    {
      // Reset terms that have been marked.
      term.unmark();
      ++it;
    }
  }

  if constexpr (EnableBlockAllocator)
  {
    // Clean up unnecessary blocks.
    m_erasedBlocks = m_term_set.get_allocator().consolidate();
  }
}

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::resize_if_needed()
{
  m_term_set.rehash_if_needed();
}

/// PRIVATE FUNCTIONS

ATERM_POOL_STORAGE_TEMPLATES
void ATERM_POOL_STORAGE::call_deletion_hook(unprotected_aterm_core term)
{
  for (const auto& [symbol, callback] : m_deletion_hooks)
  {
    if (symbol == term.function())
    {
      assert(verify_term(*detail::address(term)));
      callback(reinterpret_cast<const aterm&>(term));   // TODO Check whether this cast is OK. It was a static cast. 
    }
  }
}


ATERM_POOL_STORAGE_TEMPLATES
bool ATERM_POOL_STORAGE::verify_mark()
{
  // Check for consistency that if a term is marked then its arguments are as well.
  for (const Element& term : m_term_set)
  {
    if (term.is_marked() && term.function().arity() > 0)
    {
       const _term_appl& ta = static_cast<_term_appl&>(const_cast<Element&>(term));
       for (std::size_t i = 0; i < ta.function().arity(); ++i)
       {
         assert(detail::address(ta.arg(i))->is_marked());
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

  // Trigger the deletion hook before the term is actually destroyed.
  call_deletion_hook(&term);

  // Remove them from the hash table, will also destroy terms with fixed arity.
  return m_term_set.erase(it);
}

ATERM_POOL_STORAGE_TEMPLATES
template<typename ...Args>
bool ATERM_POOL_STORAGE::emplace(aterm_core& term, Args&&... args)
{
  auto [it, added] = m_term_set.emplace(std::forward<Args>(args)...);
  new (&term) atermpp::unprotected_aterm_core(&*it); 

  if (added)
  {
    if (EnableCreationMetrics) { m_term_metric.miss(); }
  }
  else if (EnableCreationMetrics)
  {
    // A term was already found in the set.
    m_term_metric.hit();
  }

  return added;
}

ATERM_POOL_STORAGE_TEMPLATES
constexpr bool ATERM_POOL_STORAGE::is_dynamic_storage() const
{
  return N == DynamicNumberOfArguments;
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
    assert(ta.function().defined());
    for (std::size_t i = 0; i < ta.function().arity(); ++i)
    {
      assert(ta.arg(i).defined());
      assert(&ta!=detail::address(ta.arg(i)));
      verify_term<DynamicNumberOfArguments>(*detail::address(ta.arg(i)));
    }
  }
  return true;
}

#undef ATERM_POOL_STORAGE_TEMPLATES
#undef ATERM_POOL_STORAGE

} // namespace atermpp::detail


#endif
