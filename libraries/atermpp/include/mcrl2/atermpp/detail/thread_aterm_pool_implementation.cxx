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

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/detail/index_traits.cxx"
  #include "mcrl2/atermpp/detail/thread_aterm_pool.cxx"
#endif

namespace atermpp::detail
{

function_symbol thread_aterm_pool::create_function_symbol(std::string&& name, const std::size_t arity, const bool check_for_registered_functions)
{
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  function_symbol symbol = m_pool.create_function_symbol(std::move(name), arity, check_for_registered_functions);
  return symbol;
}

function_symbol thread_aterm_pool::create_function_symbol(const std::string& name, const std::size_t arity, const bool check_for_registered_functions)
{
  std::string name_copy = name;
  return create_function_symbol(std::move(name_copy), arity, check_for_registered_functions);
}

void thread_aterm_pool::create_int(aterm& term, size_t val)
{
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  bool added = m_pool.create_int(term, val);
  guard.unlock_shared();
   
  if (added) { m_pool.created_term(!m_shared_mutex.is_shared_locked(), m_shared_mutex); }
}

void thread_aterm_pool::create_term(aterm& term, const atermpp::function_symbol& sym)
{
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  bool added = m_pool.create_term(term, sym);
  guard.unlock_shared();

  if (added) { m_pool.created_term(!m_shared_mutex.is_shared_locked(), m_shared_mutex); }
}

template<class ...Terms>
void thread_aterm_pool::create_appl(aterm& term, const function_symbol& sym, const Terms&... arguments)
{
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  bool added = m_pool.create_appl(term, sym, arguments...);
  guard.unlock_shared();

  if (added) { m_pool.created_term(!m_shared_mutex.is_shared_locked(), m_shared_mutex); }
}

template<class Term, class INDEX_TYPE, class ...Terms>
void thread_aterm_pool::create_appl_index(aterm& term, const function_symbol& sym, const Terms&... arguments)
{
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  std::array<unprotected_aterm_core, sizeof...(arguments)> argument_array;
  store_in_argument_array(argument_array, arguments...);

  bool added;
  if constexpr (sizeof...(arguments)==1)
  {
    /* Code below is more elegant than succeeding code, but it unnecessarily copies and protects a term.
        m_pool.create_int(term, atermpp::detail::index_traits<Term, INDEX_TYPE, 1>::
            insert(static_cast<INDEX_TYPE>(static_cast<aterm>(address(argument_array[0]))))); */
    m_pool.create_int(term, 
                      atermpp::detail::index_traits<Term, INDEX_TYPE, 1>::
                                insert(*reinterpret_cast<INDEX_TYPE*>(&(argument_array[0]))));
    added = m_pool.create_appl(term, sym, argument_array[0], term);
  }
  else
  {
    /* Code below is more elegant than the actual succeeding code, but it unnecessarily creates an extra pair, which is costly. 
        m_pool.create_int(
          term,
          atermpp::detail::index_traits<Term, INDEX_TYPE, 2>::
            insert(std::make_pair(static_cast<typename INDEX_TYPE::first_type>(static_cast<aterm>(address(argument_array[0]))),
                                  static_cast<typename INDEX_TYPE::second_type>(static_cast<aterm>(address(argument_array[1])))))); */
    m_pool.create_int(term,
                      atermpp::detail::index_traits<Term, INDEX_TYPE, 2>::
                                insert(*reinterpret_cast<INDEX_TYPE*>(&argument_array[0])));
    added = m_pool.create_appl(term, sym, argument_array[0], argument_array[1], term);
  }
  guard.unlock_shared();

  if (added) { m_pool.created_term(!m_shared_mutex.is_shared_locked(), m_shared_mutex); }
}

template<typename InputIterator>
void thread_aterm_pool::create_appl_dynamic(aterm& term,
                            const function_symbol& sym,
                            InputIterator begin,
                            InputIterator end)
{
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  bool added = m_pool.create_appl_dynamic(term, sym, begin, end);
  guard.unlock_shared();
    
  if (added) { m_pool.created_term(!m_shared_mutex.is_shared_locked(), m_shared_mutex); }
}

template<typename InputIterator, typename ATermConverter>
void thread_aterm_pool::create_appl_dynamic(aterm& term,
                            const function_symbol& sym,
                            ATermConverter convert_to_aterm,
                            InputIterator begin,
                            InputIterator end)
{  
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  bool added = m_pool.create_appl_dynamic(term, sym, convert_to_aterm, begin, end);
  guard.unlock_shared();

  if (added) { m_pool.created_term(!m_shared_mutex.is_shared_locked(), m_shared_mutex); }
}

void thread_aterm_pool::register_variable(aterm_core* variable)
{
  if constexpr (EnableVariableRegistrationMetrics) { ++m_variable_insertions; }

  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
      
  // Resizing of the protection set should not interfere with garbage collection and rehashing 
  if (m_variables->must_resize())
  {
    m_variables->resize();
  }

[[maybe_unused]]
  auto [it, inserted] = m_variables->insert(variable);

  // The variable must be inserted.
  assert(inserted);
}

void thread_aterm_pool::deregister_variable(aterm_core* variable)
{
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  m_variables->erase(variable);
}

void thread_aterm_pool::register_container(aterm_container* container)
{
  if constexpr (EnableVariableRegistrationMetrics) { ++m_container_insertions; }

  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  if (m_containers->must_resize())
  {
    m_containers->resize();
  }

  [[maybe_unused]]
  auto [it, inserted] = m_containers->insert(container);
  
  // The container must be inserted.
  assert(inserted);
}

void thread_aterm_pool::deregister_container(aterm_container* container)
{
  mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
  m_containers->erase(container);
}

void thread_aterm_pool::mark()
{
  for (const aterm_core* variable : *m_variables) 
  {
    if (variable != nullptr)
    {
      // Mark all terms (and their subterms) that are reachable, i.e the root set.
      _aterm* term = detail::address(*variable);
      if (term != nullptr && !term->is_marked()) 
      {
        // This variable is not a default term and that term has not been marked.
        mark_term(*term, m_todo);
      }
    }
  }

  for (const aterm_container* container : *m_containers)
  {
    if (container != nullptr)
    {
      // The container marks the contained terms itself.
      container->mark(m_todo);
    }
  }
}

void thread_aterm_pool::print_local_performance_statistics() const
{
  if constexpr (EnableVariableRegistrationMetrics)
  {
    mCRL2log(mcrl2::log::info) << "thread_aterm_pool: " << m_variables->size() << " variables in root set (" << m_variable_insertions << " total insertions)"
                               << " and " << m_containers->size() << " containers in root set (" << m_container_insertions << " total insertions).\n";
  }
}

std::size_t thread_aterm_pool::protection_set_size() const
{
  std::size_t result = m_variables->size();

  for (const auto& container : *m_containers)
  {
    if (container != nullptr)
    {
      result += container->size();
    }
  }

  return result;
}

} // namespace atermpp::detail


#endif // ATERMPP_DETAIL_ATERM_POOL_IMPLEMENTATION_H
