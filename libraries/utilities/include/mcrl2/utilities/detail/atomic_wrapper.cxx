// Author(s): Jan Friso
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_UTILITIES_DETAIL_ATOMIC_WRAPPER_H_
#define MCRL2_UTILITIES_DETAIL_ATOMIC_WRAPPER_H_
MCRL2_MODULE;

#include <atomic>

#ifdef MCRL2_ENABLE_MODULES
  export module utilities_detail:atomic_wrapper;
#endif

MCRL2_MODULE_EXPORT namespace mcrl2::utilities::detail
{

// The purpose of this wrapper is to allow copying of atomics.
template<typename T>
struct atomic_wrapper: public std::atomic<T>
{
  atomic_wrapper()
   : std::atomic<T>(0)
  {}

  atomic_wrapper(const T t)
   : std::atomic<T>(t)
  {}

  atomic_wrapper(const atomic_wrapper<T>& other)
    : std::atomic<T>(other.load())
  {}

  atomic_wrapper& operator=(const atomic_wrapper<T>& other)
  {
    this->store(other.load());
    return *this;
  }
};

} // namespace mcrl2::utilities::detail

#endif // MCRL2_UTILITIES_DETAIL_ATOMIC_WRAPPER_H_