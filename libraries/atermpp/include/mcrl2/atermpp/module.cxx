// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/execution_timer.h
/// \brief Class to obtain running times of code.

#ifndef MCRL2_UTILITIES_UTILTIES_H
#define MCRL2_UTILITIES_UTILTIES_H

#ifdef MCRL2_ENABLE_MODULES
  export module atermpp;

  export import :aterm;
  export import :aterm_core;
  export import :aterm_list;
  export import :concepts;
  export import :function_symbol;
  export import :function_symbol_generator;
  export import :type_traits;
  export import :detail.aterm_appl_iterator;  
  export import :detail.aterm_configuration;
  export import :detail.aterm_container;
  export import :detail.aterm_core_data;
  export import :detail.aterm_data;
  export import :detail.aterm_hash;
  export import :detail.aterm_int_data;
  export import :detail.aterm_list_data;
  export import :detail.aterm_list_iterator;
  export import :detail.aterm_pool_storage;  
  export import :detail.aterm_pool;  
  export import :detail.function_symbol_data;
  export import :detail.function_symbol_hash;
  export import :detail.function_symbol_pool;  
  export import :detail.global_aterm_pool;  
  export import :detail.thread_aterm_pool;  
  export import :detail.index_traits; 
  export import :standard_containers.unordered_map; 
#endif

#endif // MCRL2_UTILITIES_UTILTIES_H