// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef ATERMPP_DETAIL_THREAD_ATERM_POOL_IMPLEMENTATION_H
#define ATERMPP_DETAIL_THREAD_ATERM_POOL_IMPLEMENTATION_H
#pragma once

#include "thread_aterm_pool.h"

#include <chrono>

namespace atermpp
{
namespace detail
{

function_symbol thread_aterm_pool::create_function_symbol(const std::string& name, const std::size_t arity, const bool check_for_registered_functions)
{
  enter();
  function_symbol symbol = m_pool.create_function_symbol(name, arity, check_for_registered_functions);
  leave();
  return symbol;
}

void thread_aterm_pool::create_int(aterm& term, size_t val)
{
  enter();
  bool added = m_pool.create_int(term, val);
  leave();
  if (added) { m_pool.trigger_collection(); }
}

void thread_aterm_pool::create_term(aterm& term, const atermpp::function_symbol& sym)
{
  enter();
  bool added = m_pool.create_term(term, sym);
  leave();
  if (added) { m_pool.trigger_collection(); }
}

template<class ...Terms>
void thread_aterm_pool::create_appl(aterm& term, const function_symbol& sym, const Terms&... arguments)
{
  enter();
  bool added = m_pool.create_appl(term, sym, arguments...);
  leave();
  if (added) { m_pool.trigger_collection(); }
}

template<typename InputIterator>
void thread_aterm_pool::create_appl_dynamic(aterm& term,
                            const function_symbol& sym,
                            InputIterator begin,
                            InputIterator end)
{
  enter();
  bool added = m_pool.create_appl_dynamic(term, sym, begin, end);
  leave();
  if (added) { m_pool.trigger_collection(); }
}

template<typename InputIterator, typename ATermConverter>
void thread_aterm_pool::create_appl_dynamic(aterm& term,
                            const function_symbol& sym,
                            ATermConverter convert_to_aterm,
                            InputIterator begin,
                            InputIterator end)
{
  if constexpr (GlobalThreadSafe) { throw std::runtime_error("This function is not thread safe"); }

  enter();
  ++m_creation_depth;
  bool added = m_pool.create_appl_dynamic(term, sym, convert_to_aterm, begin, end);
  --m_creation_depth;
  leave();

  if (added) { m_pool.trigger_collection(); }
}

void thread_aterm_pool::register_variable(aterm* variable)
{
  if constexpr (EnableVariableRegistrationMetrics) { m_variable_cache.hit(); }

  auto [it, inserted] = m_variables.emplace(variable);

  // Apply the tag for the inserted variable.
  mcrl2::utilities::mcrl2_unused(inserted);
  it->tag();
}

void thread_aterm_pool::remove_variable(aterm* variable)
{
  // This variable was registered at some point.
  assert(m_variables.find(mcrl2::utilities::tagged_pointer<aterm>(variable)) != m_variables.end());

  if constexpr (EnableVariableRegistrationMetrics) { m_variable_cache.miss(); }

  auto it = m_variables.find(mcrl2::utilities::tagged_pointer<aterm>(variable));
  assert(it != m_variables.end());
  it->untag();
}

void thread_aterm_pool::mark()
{
  // Marks all terms that are reachable from any tagged variable. Furthermore, remove variables that are not tagged.
  for (auto it = m_variables.begin(); it != m_variables.end();)
  {
    auto variable = *it;
    if (variable.tagged())
    {
      // Mark all terms (and their subterms) that are reachable, i.e the root set.
      assert(variable != nullptr);
      _aterm* term = detail::address(*variable);
      if (variable->defined() && !term->is_marked())
      {
        // This variable is not a default term and that term has not been marked.
        mark_term(*term, todo);
      }

      ++it;
    }
    else
    {
      it = m_variables.erase(it);
    }
  }
}

void thread_aterm_pool::print_local_performance_statistics() const
{
  if constexpr (EnableVariableRegistrationMetrics)
  {
    mCRL2log(mcrl2::log::info, "Performance") << "thread_aterm_pool:\n";
    mCRL2log(mcrl2::log::info, "Performance") << m_variables.size() << " variables in root set (" << m_variable_cache.message() << ".\n";
  }
}

void thread_aterm_pool::enter()
{
  if constexpr (GlobalThreadSafe) // && m_creation_depth == 0)
  {
    m_busy_flag.store(true, std::memory_order::memory_order_release);

    if (m_pool.should_wait())
    {
      m_waiting_flag = true;

      // Wait for the guard to become false.
      while (m_pool.should_wait()) {};
    }
  }
}

void thread_aterm_pool::leave()
{
  if constexpr (GlobalThreadSafe) // && m_creation_depth == 0)
  {
    m_busy_flag.store(false, std::memory_order::memory_order_release);
  }
}


} // namespace detail
} // namespace atermpp

#endif // ATERMPP_DETAIL_ATERM_POOL_IMPLEMENTATION_H
