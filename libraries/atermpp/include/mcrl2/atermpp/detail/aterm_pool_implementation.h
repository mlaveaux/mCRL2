// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef ATERMPP_DETAIL_ATERM_POOL_IMPLEMENTATION_H
#define ATERMPP_DETAIL_ATERM_POOL_IMPLEMENTATION_H
#pragma once

#include "aterm_pool.h"

#include <chrono>

namespace atermpp
{
namespace detail
{

aterm_pool::aterm_pool() :
  m_int_storage(*this),
  m_appl_storage(
    *this,
    *this,
    *this,
    *this,
    *this,
    *this,
    *this,
    *this
  ),
  m_appl_dynamic_storage(*this)
{
  m_count_until_collection = capacity();
  m_count_until_resize = m_int_storage.capacity();
  
  // Initialize the empty list.
  create_appl(m_empty_list, m_function_symbol_pool.as_empty_list());
}

aterm_pool::~aterm_pool()
{
  print_performance_statistics();
}

void aterm_pool::add_creation_hook(function_symbol sym, term_callback callback)
{
  const std::size_t arity = sym.arity();

  switch (arity)
  {
  case 0:
  {
    if (sym == get_symbol_pool().as_int())
    {
      m_int_storage.add_creation_hook(sym, callback);
    }
    else
    {
      return std::get<0>(m_appl_storage).add_creation_hook(sym, callback);
    }
    break;
  }
  case 1:
    std::get<1>(m_appl_storage).add_creation_hook(sym, callback);
    break;
  case 2:
    std::get<2>(m_appl_storage).add_creation_hook(sym, callback);
    break;
  case 3:
    std::get<3>(m_appl_storage).add_creation_hook(sym, callback);
    break;
  case 4:
    std::get<4>(m_appl_storage).add_creation_hook(sym, callback);
    break;
  case 5:
    std::get<5>(m_appl_storage).add_creation_hook(sym, callback);
    break;
  case 6:
    std::get<6>(m_appl_storage).add_creation_hook(sym, callback);
    break;
  case 7:
    std::get<7>(m_appl_storage).add_creation_hook(sym, callback);
    break;
  default:
    m_appl_dynamic_storage.add_creation_hook(sym, callback);
  }
}

void aterm_pool::add_deletion_hook(function_symbol sym, term_callback callback)
{
  const std::size_t arity = sym.arity();

  switch (arity)
  {
  case 0:
  {
    if (sym == get_symbol_pool().as_int())
    {
      m_int_storage.add_deletion_hook(sym, callback);
    }
    else
    {
      std::get<0>(m_appl_storage).add_deletion_hook(sym, callback);
    }
  }
    break;
  case 1:
    std::get<1>(m_appl_storage).add_deletion_hook(sym, callback);
    break;
  case 2:
    std::get<2>(m_appl_storage).add_deletion_hook(sym, callback);
    break;
  case 3:
    std::get<3>(m_appl_storage).add_deletion_hook(sym, callback);
    break;
  case 4:
    std::get<4>(m_appl_storage).add_deletion_hook(sym, callback);
    break;
  case 5:
    std::get<5>(m_appl_storage).add_deletion_hook(sym, callback);
    break;
  case 6:
    std::get<6>(m_appl_storage).add_deletion_hook(sym, callback);
    break;
  case 7:
    std::get<7>(m_appl_storage).add_deletion_hook(sym, callback);
    break;
  default:
    m_appl_dynamic_storage.add_deletion_hook(sym, callback);
  }
}

void aterm_pool::collect()
{
  m_count_until_collection = 0;
  collect_impl();
}

void aterm_pool::register_thread_aterm_pool(thread_aterm_pool_interface &pool)
{
  if constexpr (GlobalThreadSafe) { m_mutex.lock(); }

  mCRL2log(mcrl2::log::debug) << "Registered thread_local aterm pool\n";
  m_thread_pools.insert(m_thread_pools.end(), &pool);

  if constexpr (GlobalThreadSafe) { m_mutex.unlock(); }
}

void aterm_pool::remove_thread_aterm_pool(thread_aterm_pool_interface& pool)
{
  if constexpr (GlobalThreadSafe) { m_mutex.lock(); }

  mCRL2log(mcrl2::log::debug) << "Removed thread_local aterm pool\n";
  auto it = std::find(m_thread_pools.begin(), m_thread_pools.end(), &pool);

  if (it != m_thread_pools.end())
  {
    m_thread_pools.erase(it);
  }

  if constexpr (GlobalThreadSafe) { m_mutex.unlock(); }
}

void aterm_pool::print_performance_statistics() const
{
  m_int_storage.print_performance_stats("integral_storage");
  std::get<0>(m_appl_storage).print_performance_stats("term_storage");
  std::get<1>(m_appl_storage).print_performance_stats("function_application_storage_1");
  std::get<2>(m_appl_storage).print_performance_stats("function_application_storage_2");
  std::get<3>(m_appl_storage).print_performance_stats("function_application_storage_3");
  std::get<4>(m_appl_storage).print_performance_stats("function_application_storage_4");
  std::get<5>(m_appl_storage).print_performance_stats("function_application_storage_5");
  std::get<6>(m_appl_storage).print_performance_stats("function_application_storage_6");
  std::get<7>(m_appl_storage).print_performance_stats("function_application_storage_7");

  m_appl_dynamic_storage.print_performance_stats("arbitrary_function_application_storage");

#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  if (mcrl2::utilities::EnableReferenceCountMetrics)
  {
    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool: all reference counts changed " << _aterm::reference_count_changes() << " times.\n";
  }
#endif
  // Print information for the local aterm pools.
  for (thread_aterm_pool_interface* local : m_thread_pools)
  {
    local->print_local_performance_statistics();
  }
}

std::size_t aterm_pool::capacity() const noexcept
{
  // Determine the total number of terms in any storage.
  return m_int_storage.capacity()
    + std::get<0>(m_appl_storage).capacity()
    + std::get<1>(m_appl_storage).capacity()
    + std::get<2>(m_appl_storage).capacity()
    + std::get<3>(m_appl_storage).capacity()
    + std::get<4>(m_appl_storage).capacity()
    + std::get<5>(m_appl_storage).capacity()
    + std::get<6>(m_appl_storage).capacity()
    + std::get<7>(m_appl_storage).capacity()
    + m_appl_dynamic_storage.capacity();
}

std::size_t aterm_pool::size() const
{
  // Determine the total number of terms in any storage.
  return m_int_storage.size()
    + std::get<0>(m_appl_storage).size()
    + std::get<1>(m_appl_storage).size()
    + std::get<2>(m_appl_storage).size()
    + std::get<3>(m_appl_storage).size()
    + std::get<4>(m_appl_storage).size()
    + std::get<5>(m_appl_storage).size()
    + std::get<6>(m_appl_storage).size()
    + std::get<7>(m_appl_storage).size()
    + m_appl_dynamic_storage.size();
}

// private

void aterm_pool::trigger_collection()
{
  // Defer garbage collection when it happens too often.
  if (m_count_until_collection.fetch_sub(1, std::memory_order_relaxed) == 0)
  {
    --m_count_until_collection;
  }
  else
  {
    collect_impl();
  }

  if (m_count_until_resize.fetch_sub(1, std::memory_order_relaxed) == 0)
  {
    resize_if_needed();
  }
}

void aterm_pool::collect_impl()
{
  if (!m_enable_garbage_collection) { return; }

  halt();
  if (m_count_until_collection > 0)
  {
    // Another thread has performed garbage collection, so we can ignore it.
    resume();
    return;
  }

  for (const auto& pool : m_thread_pools)
  {
    if (!pool->allow_collect())
    {
      // Some thread is still creating a term recursively.
      resume();
      return;
    }
  }

  auto timestamp = std::chrono::system_clock::now();
  std::size_t old_size = size();

#ifdef MCRL2_ATERMPP_REFERENCE_COUNTED
  // Marks all terms that are reachable via any reachable term to
  // not be garbage collected.
  // For integer and terms without arguments the marking is not needed, because
  // they do not have arguments that might have to be marked.
  std::get<1>(m_appl_storage).mark();
  std::get<2>(m_appl_storage).mark();
  std::get<3>(m_appl_storage).mark();
  std::get<4>(m_appl_storage).mark();
  std::get<5>(m_appl_storage).mark();
  std::get<6>(m_appl_storage).mark();
  std::get<7>(m_appl_storage).mark();
  m_appl_dynamic_storage.mark();
#else
  // Instruct each thread_local to mark its terms to be protected.
  for (const auto& pool : m_thread_pools)
  {
    pool->mark();
  }
#endif

  assert(std::get<0>(m_appl_storage).verify_mark());
  assert(std::get<1>(m_appl_storage).verify_mark());
  assert(std::get<2>(m_appl_storage).verify_mark());
  assert(std::get<3>(m_appl_storage).verify_mark());
  assert(std::get<4>(m_appl_storage).verify_mark());
  assert(std::get<5>(m_appl_storage).verify_mark());
  assert(std::get<6>(m_appl_storage).verify_mark());
  assert(std::get<7>(m_appl_storage).verify_mark());
  assert(m_appl_dynamic_storage.verify_mark());

  // Keep track of the duration for marking and reset for sweep.
  auto mark_duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - timestamp).count();
  timestamp = std::chrono::system_clock::now();

  // Garbage collect terms that are not reachable.
  m_int_storage.sweep();
  std::get<0>(m_appl_storage).sweep();
  std::get<1>(m_appl_storage).sweep();
  std::get<2>(m_appl_storage).sweep();
  std::get<3>(m_appl_storage).sweep();
  std::get<4>(m_appl_storage).sweep();
  std::get<5>(m_appl_storage).sweep();
  std::get<6>(m_appl_storage).sweep();
  std::get<7>(m_appl_storage).sweep();
  m_appl_dynamic_storage.sweep();

  // Check that after sweeping the terms are consistent.
  assert(m_int_storage.verify_sweep());
  assert(std::get<0>(m_appl_storage).verify_sweep());
  assert(std::get<1>(m_appl_storage).verify_sweep());
  assert(std::get<2>(m_appl_storage).verify_sweep());
  assert(std::get<3>(m_appl_storage).verify_sweep());
  assert(std::get<4>(m_appl_storage).verify_sweep());
  assert(std::get<5>(m_appl_storage).verify_sweep());
  assert(std::get<6>(m_appl_storage).verify_sweep());
  assert(std::get<7>(m_appl_storage).verify_sweep());
  assert(m_appl_dynamic_storage.verify_sweep());

  // Print some statistics.
  if constexpr (EnableGarbageCollectionMetrics)
  {
    auto sweep_duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - timestamp).count();

    // Print the relevant information.
    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool: Garbage collected " << old_size - size() << " terms, " << size() << " terms remaining in "
      << mark_duration + sweep_duration << " ms (marking " << mark_duration << " ms + sweep " << sweep_duration << " ms).\n";
  }

  // Garbage collect function symbols.
  m_function_symbol_pool.sweep();

  print_performance_statistics();

  // Use some heuristics to determine when the next collect should be called automatically.
  m_count_until_collection = size();

  resume();
}

function_symbol aterm_pool::create_function_symbol(const std::string& name, const std::size_t arity, const bool check_for_registered_functions)
{
  return m_function_symbol_pool.create(name, arity, check_for_registered_functions);
}

bool aterm_pool::create_int(aterm& term, size_t val)
{
  return m_int_storage.create_int(term, val);
}

bool aterm_pool::create_term(aterm& term, const atermpp::function_symbol& sym)
{
  return std::get<0>(m_appl_storage).create_term(term, sym);
}

template<class ...Terms>
bool aterm_pool::create_appl(aterm& term, const function_symbol& sym, const Terms&... arguments)
{
  return std::get<sizeof...(Terms)>(m_appl_storage).create_appl(term, sym, arguments...);
}

template<typename ForwardIterator>
bool aterm_pool::create_appl_dynamic(aterm& term,
                            const function_symbol& sym,
                            ForwardIterator begin,
                            ForwardIterator end)
{
  const std::size_t arity = sym.arity();

  switch(arity)
  {
  case 0:
    return std::get<0>(m_appl_storage).create_term(term, sym);
  case 1:
    return std::get<1>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
  case 2:
    return std::get<2>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
  case 3:
    return std::get<3>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
  case 4:
    return std::get<4>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
  case 5:
    return std::get<5>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
  case 6:
    return std::get<6>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
  case 7:
    return std::get<7>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
  default:
    return m_appl_dynamic_storage.create_appl_dynamic(term, sym, begin, end);
  }
}

template<typename InputIterator, typename ATermConverter>
bool aterm_pool::create_appl_dynamic(aterm& term,
                            const function_symbol& sym,
                            ATermConverter converter,
                            InputIterator begin,
                            InputIterator end)
{
  const std::size_t arity = sym.arity();
  switch(arity)
  {
  case 0:
    return std::get<0>(m_appl_storage).create_term(term, sym);
  case 1:
    return std::get<1>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
  case 2:
    return std::get<2>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
  case 3:
    return std::get<3>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
  case 4:
    return std::get<4>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
  case 5:
    return std::get<5>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
  case 6:
    return std::get<6>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
  case 7:
    return std::get<7>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
  default:
    return m_appl_dynamic_storage.create_appl_dynamic(term, sym, converter, begin, end);
  }
}

bool aterm_pool::should_wait()
{
  return m_guard.load(std::memory_order::memory_order_acquire);
}

void aterm_pool::wait()
{
  std::unique_lock lock(m_mutex);
  m_wait_variable.wait(lock, [this]() { return !m_guard; });
  lock.unlock();
}

void aterm_pool::resize_if_needed()
{
  halt();
  if (m_count_until_resize > 0)
  {
    // Another thread has resized the tables, so we can ignore it.
    resume();
    return;
  }

  auto timestamp = std::chrono::system_clock::now();
  std::size_t old_capacity = capacity();

  // Attempt to resize all storages.
  m_function_symbol_pool.resize_if_needed();

  m_int_storage.resize_if_needed();
  std::get<0>(m_appl_storage).resize_if_needed();
  std::get<1>(m_appl_storage).resize_if_needed();
  std::get<2>(m_appl_storage).resize_if_needed();
  std::get<3>(m_appl_storage).resize_if_needed();
  std::get<4>(m_appl_storage).resize_if_needed();
  std::get<5>(m_appl_storage).resize_if_needed();
  std::get<6>(m_appl_storage).resize_if_needed();
  std::get<7>(m_appl_storage).resize_if_needed();
  m_appl_dynamic_storage.resize_if_needed();

  // Attempt to resize ever so often.
  m_count_until_resize = 10000;

  if (EnableGarbageCollectionMetrics && old_capacity != capacity())
  {
    // Only print if a resize actually took place.
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - timestamp).count();

    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool: Resized hash tables from " << old_capacity << " to " << capacity() << " capacity in "
                                              << duration << " ms.\n";
  }

  resume();
}

void aterm_pool::halt()
{
  if constexpr (!GlobalThreadSafe) { return; }

  // Only one thread can halt everything.
  m_mutex.lock();

  // Indicate that threads must wait.
  m_guard = true;

  // Wait for all pools to indicate that they are not busy.
  bool all_finished = true;
  do
  {
    all_finished = true;
    for (const auto& pool : m_thread_pools)
    {
      all_finished = all_finished && !pool->busy();
    }
  }
  while(!all_finished);
}

void aterm_pool::resume()
{
  if constexpr (!GlobalThreadSafe) { return; }

  m_guard = false;
  m_mutex.unlock();

  // Notify waiting threads that they can proceed (guard will be false).
  m_wait_variable.notify_all();
}

} // namespace detail
} // namespace atermpp

#endif // ATERMPP_DETAIL_ATERM_POOL_IMPLEMENTATION_H
