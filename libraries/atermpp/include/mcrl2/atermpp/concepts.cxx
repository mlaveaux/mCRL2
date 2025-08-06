// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRl2_ATERMPP_CONCEPTS_H
#define MCRl2_ATERMPP_CONCEPTS_H
MCRL2_MODULE;

#include <type_traits>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm_core.cxx"
#else
  export module atermpp:concepts;

  import :aterm_core;
#endif

MCRL2_MODULE_EXPORT namespace atermpp {

/// Concept that can be used to indicate that T is of an ATerm type.
template<typename T>
concept IsATerm = requires(T t)
{
    /// Term must be derived from the base aterm.
    requires std::is_base_of_v<atermpp::unprotected_aterm_core, std::remove_reference_t<T>>;

    /// aterm cast cannot be applied types derived from aterms where extra fields are added.
    requires sizeof(std::remove_reference_t<T>) == sizeof(atermpp::unprotected_aterm_core);

    /// A pointer to a standard-layout class may be converted (with reinterpret_cast) to a pointer to its first non-static data member and vice versa. 
    requires std::is_standard_layout_v<T>;
};

/// Concept that can be used to indicate that T is a function that can convert an aterm to another aterm.
template<typename T>
concept IsTermConverter = requires(T t)
{
    /// Calling the function with an aterm must yield an aterm.
    std::is_convertible_v<std::remove_reference_t<T>, atermpp::unprotected_aterm_core>;
    std::is_invocable_v<std::remove_reference_t<T>, atermpp::unprotected_aterm_core>;
};

static_assert(atermpp::IsATerm<unprotected_aterm_core>, "unprotected_aterm_core must be an aterm");
static_assert(atermpp::IsATerm<aterm_core>, "aterm_core must be an aterm");

} // namespace atermpp

#endif // MCRl2_ATERMPP_CONCEPTS_H