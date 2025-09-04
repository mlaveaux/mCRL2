// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#pragma once

#include <ranges>
#include <iostream>

namespace mcrl2::utilities 
{

/// Concept that can be used to identify C as a container.
/// TODO: Is term_list a container, and if not, why?
template<typename C>
concept IsContainer = requires(C c)
{
    /// Container must be a range.
    requires std::ranges::range<std::remove_reference_t<C>>;

    /// Container must have a value type.
    typename std::remove_reference_t<C>::value_type;
};

template<typename T>
concept IsPrintable = requires(T t) {
    { std::cout << t } -> std::same_as<std::ostream&>;
};

} // namespace mcrl2::utilities