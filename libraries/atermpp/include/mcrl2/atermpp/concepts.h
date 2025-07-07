// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_ATERMPP_CONCEPTS_H
#define MCRL2_ATERMPP_CONCEPTS_H

#include "mcrl2/atermpp/unprotected_aterm.h"

#include <type_traits>

namespace atermpp {

/// Concept that can be used to indicate that T is of an ATerm type.
template<typename T>
concept IsATerm = requires(T t)
{
    /// Term must be derived from the base aterm.
    requires std::is_base_of<atermpp::unprotected_aterm, std::remove_reference_t<T>>::value;

    /// aterm can only contain the unprotected_aterm, and optionally the root index.
    requires sizeof(std::remove_reference_t<T>) <= sizeof(atermpp::unprotected_aterm) + sizeof(std::size_t);
};

/// Concept that can be used to indicate that T is a function that can convert an aterm to another aterm.
template<typename T>
concept IsTermConverter = requires(T t)
{
    /// Calling the function with an aterm must yield an aterm.
    std::is_convertible<std::remove_reference_t<T>, atermpp::unprotected_aterm>::value;
    std::is_invocable<std::remove_reference_t<T>, atermpp::unprotected_aterm>::value;
};

static_assert(atermpp::IsATerm<unprotected_aterm>, "unprotected_aterm_core must be an aterm");

} // namespace atermpp

#endif // MCRL2_ATERMPP_CONCEPTS_H