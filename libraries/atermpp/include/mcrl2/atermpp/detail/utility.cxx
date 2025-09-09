// Author(s): Maurice Laveaux.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef ATERMPP_DETAIL_UTILITY_H
#define ATERMPP_DETAIL_UTILITY_H
MCRL2_MODULE;

#include <array>
#include <cstdint>

#ifdef MCRL2_ENABLE_MODULES
    export module atermpp:detail.utility;

    import :aterm_core;

    import utilities;
#endif

MCRL2_MODULE_EXPORT namespace atermpp::detail
{  
  template <std::size_t N>
  void store_in_argument_array_(std::size_t , std::array<unprotected_aterm_core, N>& )
  {}

  template <std::size_t N, class FUNCTION_OR_TERM_TYPE, typename... Args>
  void store_in_argument_array_(std::size_t i, 
                                std::array<unprotected_aterm_core, N>& argument_array, 
                                FUNCTION_OR_TERM_TYPE& function_or_term, 
                                const Args&... args)
  {
    if constexpr (std::is_convertible_v<FUNCTION_OR_TERM_TYPE, unprotected_aterm_core>)
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
}

#endif // ATERMPP_DETAIL_UTILITY_H