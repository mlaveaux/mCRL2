// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
MCRL2_MODULE;

#include <mutex>
#include <chrono>

#include "mcrl2/utilities/logger.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/function_symbol_pool.cxx"
#else
  module atermpp;

  import :detail.function_symbol_pool;

  import utilities;
#endif

using namespace atermpp;
using namespace atermpp::detail;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::detail;

function_symbol_pool::function_symbol_pool()
{
  // Initialize the default function symbols.
  m_as_int = create("<aterm_int>", 0);
  m_as_list = create("<list_constructor>", 2);
  m_as_empty_list = create("<empty_list>", 0);

  // Initialize the global copies of these default function symbols.
  g_as_int = m_as_int;
  g_as_list = m_as_list;
  g_as_empty_list = m_as_empty_list;
}

void function_symbol_pool::create_helper(const std::string& name)
{
  if constexpr (GlobalThreadSafe) { m_mutex.lock(); }

  // Check whether there is a registered prefix p such that name equal pn where n is a number.
  // In that case prevent that pn will be generated as a fresh function name.
  std::size_t start_of_index = name.find_last_not_of("0123456789") + 1;

  if (start_of_index < name.size()) // Otherwise there is no trailing number.
  {
    std::string potential_number = name.substr(start_of_index); // Get the trailing string after prefix_ of function_name.
    std::string prefix = name.substr(0, start_of_index);
    auto prefix_it = m_prefix_to_register_function_map.find(prefix);
    if (prefix_it != m_prefix_to_register_function_map.end())  // points to the prefix.
    {
      try
      {
        std::size_t number = std::stoul(potential_number);
        *prefix_it->second = std::max(*prefix_it->second, number + 1); // Set the index belonging to the found prefix to at least a safe number+1.
      }
      catch (std::exception&)
      {
        // Can be std::invalid_argument or an out_of_range exception.
        // In both cases nothing needs to be done, and the exception can be ignored.
      }
    }
  }

  if  constexpr (GlobalThreadSafe) { m_mutex.unlock(); }
}

function_symbol function_symbol_pool::create(std::string&& name, const std::size_t arity, const bool check_for_registered_functions)
{
  auto it = m_symbol_set.find(name, arity);
  if (it != m_symbol_set.end())
  {
    if constexpr (EnableCreationMetrics) { m_function_symbol_metrics.hit(); }

    // The element already exists so return it.
    return function_symbol(_function_symbol::ref(&(*it)));
  }
  else
  {
    if constexpr (EnableCreationMetrics) { m_function_symbol_metrics.miss(); }

    const _function_symbol& symbol = *m_symbol_set.emplace(std::forward<std::string>(name), arity).first;
    if (check_for_registered_functions)
    {
      create_helper(symbol.name());
    }

    return function_symbol(_function_symbol::ref(&symbol));
  }
}

function_symbol function_symbol_pool::create(const std::string& name, const std::size_t arity, const bool check_for_registered_functions)
{
  auto it = m_symbol_set.find(name, arity);
  if (it != m_symbol_set.end())
  {
    if constexpr (EnableCreationMetrics) { m_function_symbol_metrics.hit(); }

    // The element already exists so return it.
    return function_symbol(_function_symbol::ref(&(*it)));
  }
  else
  {
    if constexpr (EnableCreationMetrics) { m_function_symbol_metrics.miss(); }

    const _function_symbol& symbol = *m_symbol_set.emplace(name, arity).first;
    if (check_for_registered_functions)
    {
      create_helper(name);
    }

    return function_symbol(_function_symbol::ref(&symbol));
  }
}

void function_symbol_pool::deregister(const std::string& prefix)
{
  m_mutex.lock();
  m_prefix_to_register_function_map.erase(prefix);
  m_mutex.unlock();
}

std::shared_ptr<std::size_t> function_symbol_pool::register_prefix(const std::string& prefix)
{
  m_mutex.lock();

  auto it = m_prefix_to_register_function_map.find(prefix);
  if (it != m_prefix_to_register_function_map.end())
  {
    auto result = it->second;
    m_mutex.unlock();
    return result;
  }
  else
  {
    std::size_t index = get_sufficiently_large_postfix_index(prefix);
    std::shared_ptr<std::size_t> shared_index = std::make_shared<std::size_t>(index);
    m_prefix_to_register_function_map[prefix] = shared_index;

    m_mutex.unlock();
    return shared_index;
  }
}

std::size_t function_symbol_pool::get_sufficiently_large_postfix_index(const std::string& prefix) const
{
  std::size_t index = 0;
  for (const auto& f : m_symbol_set)
  {
    const std::string& function_name = f.name();

    if (function_name.starts_with(prefix)) // The function name starts with the prefix
    {
      std::string potential_number = function_name.substr(prefix.size()); // Get the trailing string after prefix_ of function_name.
      std::size_t end_of_number;
      try
      {
        std::size_t number = std::stoul(potential_number, &end_of_number);
        if (end_of_number == potential_number.size()) // A proper number was read.
        {
          if (number >= index)
          {
            index = number + 1;
          }
        }
      }
      catch (std::exception&)
      {
        // Can be std::invalid_argument or an out_of_range exception.
        // In both cases nothing needs to be done, and the exception can be ignored.
      }
    }
  }

  return index;
}

void function_symbol_pool::sweep()
{
  // Prevents changes to the symbol_set
  std::unique_lock lock(m_mutex);

  auto timestamp = std::chrono::system_clock::now();
  std::size_t old_size = size();

  // Clean up function symbols with a zero reference count.
  for (auto it = m_symbol_set.begin(); it != m_symbol_set.end(); )
  {
    if (it->reference_count() == 0)
    {
      it = m_symbol_set.erase(it);
    }
    else
    {
      ++it;
    }
  }

  std::size_t erased_blocks = 0; //m_symbol_set.get_allocator().consolidate();

  if constexpr (EnableGarbageCollectionMetrics)
  {
    auto sweep_duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - timestamp).count();

    // Print the relevant information.
    mCRL2log(mcrl2::log::info) << "function_symbol_pool: Garbage collected " << old_size - size() << " function symbols, " << size() << " function symbols remaining in "
      << sweep_duration << " ms.\n";

    mCRL2log(mcrl2::log::info) << "function_symbol_pool: Consolidate removed " << erased_blocks << " blocks.\n";
  }

  if constexpr (EnableHashtableMetrics)
  {
    print_performance_statistics(m_symbol_set);
  }

  if constexpr (EnableCreationMetrics)
  {
    mCRL2log(mcrl2::log::info) << "g_function_symbol_pool: Stores " << size() << " function symbols. create() " << m_function_symbol_metrics.message() << ".\n";
  }

  if constexpr (EnableReferenceCountMetrics)
  {
    mCRL2log(mcrl2::log::info) << "g_function_symbol_pool: all reference counts changed " << _function_symbol::reference_count_changes() << " times.\n";
  }
}
void function_symbol_pool::resize_if_needed()
{
  m_symbol_set.rehash_if_needed();
}
