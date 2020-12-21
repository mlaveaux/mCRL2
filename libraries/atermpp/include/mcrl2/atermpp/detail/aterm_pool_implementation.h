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

void aterm_pool::trigger_collection()
{
  if (m_count_until_collection > 0)
  {
    --m_count_until_collection;
  }
  else
  {
    if (m_enable_garbage_collection)
    {
      collect();
    }
  }

  if (m_count_until_resize > 0)
  {
    --m_count_until_resize;
  }
  else
  {
    resize_if_needed();
  }
}

void aterm_pool::collect()
{
  if (m_creation_depth > 0)
  {
    m_deferred_garbage_collection = true;
    return;
  }

  auto timestamp = std::chrono::system_clock::now();

  m_deferred_garbage_collection = false;
  std::size_t old_size = size();

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

  // Collect all terms that are not reachable or marked.
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
  if (EnableGarbageCollectionMetrics)
  {
    // Update the times
    auto sweep_duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - timestamp).count();

    // Print the relevant information.
    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool: Garbage collected " << old_size - size() << " terms, " << size() << " terms remaining in "
      << mark_duration + sweep_duration << " ms (marking " << mark_duration << " ms + sweep " << sweep_duration << " ms).\n";
  }

  print_performance_statistics();

  // Also garbage collect the symbol pool.
  m_function_symbol_pool.sweep();

  // Use some heuristics to determine when the next collection is called.
  m_count_until_collection = size();
}

void aterm_pool::enable_garbage_collection(bool enable)
{
  m_enable_garbage_collection = enable;
}

void aterm_pool::create_int(aterm& term, size_t val)
{
  if (m_int_storage.create_int(term, val))
  {
    trigger_collection();
  }
}

void aterm_pool::create_term(aterm& term, const atermpp::function_symbol& sym)
{
  if (std::get<0>(m_appl_storage).create_term(term, sym))
  {
    trigger_collection();
  }
}

template<class ...Terms>
void aterm_pool::create_appl(aterm& term, const function_symbol& sym, const Terms&... arguments)
{
  if (std::get<sizeof...(Terms)>(m_appl_storage).create_appl(term, sym, arguments...))
  {
    trigger_collection();
  }
}

template<typename ForwardIterator>
void aterm_pool::create_appl_dynamic(aterm& term,
                            const function_symbol& sym,
                            ForwardIterator begin,
                            ForwardIterator end)
{
  const std::size_t arity = sym.arity();

  bool added = false;
  switch(arity)
  {
  case 0:
    added = std::get<0>(m_appl_storage).create_term(term, sym);
    break;
  case 1:
    added = std::get<1>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
    break;
  case 2:
    added = std::get<2>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
    break;
  case 3:
    added = std::get<3>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
    break;
  case 4:
    added = std::get<4>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
    break;
  case 5:
    added = std::get<5>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
    break;
  case 6:
    added = std::get<6>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
    break;
  case 7:
    added = std::get<7>(m_appl_storage).template create_appl_iterator<ForwardIterator>(term, sym, begin, end);
    break;
  default:
    added = m_appl_dynamic_storage.create_appl_dynamic(term, sym, begin, end);
  }

  if (added)
  {
    trigger_collection();
  }
}

template<typename InputIterator, typename ATermConverter>
void aterm_pool::create_appl_dynamic(aterm& term,
                            const function_symbol& sym,
                            ATermConverter converter,
                            InputIterator begin,
                            InputIterator end)
{
  ++m_creation_depth;

  const std::size_t arity = sym.arity();

  bool added = false;
  switch(arity)
  {
  case 0:
    added = std::get<0>(m_appl_storage).create_term(term, sym);
    break;
  case 1:
    added = std::get<1>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
    break;
  case 2:
    added = std::get<2>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
    break;
  case 3:
    added = std::get<3>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
    break;
  case 4:
    added = std::get<4>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
    break;
  case 5:
    added = std::get<5>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
    break;
  case 6:
    added = std::get<6>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
    break;
  case 7:
    added = std::get<7>(m_appl_storage).template create_appl_iterator<InputIterator, ATermConverter>(term, sym, converter, begin, end);
    break;
  default:
    added = m_appl_dynamic_storage.create_appl_dynamic(term, sym, converter, begin, end);
  }

  if (added)
  {
    trigger_collection();
  }

  --m_creation_depth;

  // Trigger a deferred garbage collection when it was requested and the term has been protected.
  if (m_creation_depth == 0 && m_deferred_garbage_collection)
  {
    if (EnableGarbageCollectionMetrics)
    {
      mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool: Deferred garbage collection.\n";
    }
    collect();
  }
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

  if (mcrl2::utilities::EnableReferenceCountMetrics)
  {
    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool: all reference counts changed " << _aterm::reference_count_changes() << " times.\n";
  }
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

// private functions

void aterm_pool::resize_if_needed()
{
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

  // Find the hash table with the least amount of free buckets.
  m_count_until_resize = std::min(m_int_storage.capacity() - m_int_storage.size(),
                         std::min(std::get<0>(m_appl_storage).capacity() - std::get<0>(m_appl_storage).size(),
                         std::min(std::get<1>(m_appl_storage).capacity() - std::get<1>(m_appl_storage).size(),
                         std::min(std::get<2>(m_appl_storage).capacity() - std::get<2>(m_appl_storage).size(),
                         std::min(std::get<3>(m_appl_storage).capacity() - std::get<3>(m_appl_storage).size(),
                         std::min(std::get<4>(m_appl_storage).capacity() - std::get<4>(m_appl_storage).size(),
                         std::min(std::get<5>(m_appl_storage).capacity() - std::get<5>(m_appl_storage).size(),
                         std::min(std::get<6>(m_appl_storage).capacity() - std::get<6>(m_appl_storage).size(),
                         std::min(std::get<7>(m_appl_storage).capacity() - std::get<7>(m_appl_storage).size(),
                                  m_appl_dynamic_storage.capacity() - m_appl_dynamic_storage.size())))))))));

  if (EnableGarbageCollectionMetrics && old_capacity != capacity())
  {
    // Only print if a resize actually took place.
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - timestamp).count();

    mCRL2log(mcrl2::log::info, "Performance") << "aterm_pool: Resized hash tables from " << old_capacity << " to " << capacity() << " capacity in "
                                              << duration << " ms.\n";

    print_performance_statistics();
  }
}

} // namespace detail
} // namespace atermpp

#endif // ATERMPP_DETAIL_ATERM_POOL_IMPLEMENTATION_H
