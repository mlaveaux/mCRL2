// Author(s): Maurice Laveaux and Menno Bartels
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <ranges>

namespace mcrl2::pbes_system::detail
{

/// Fold is only available in C++23 so we provide a simple implementation here.
template<typename T, typename BinaryOperation>
  requires std::is_invocable_r_v<T, BinaryOperation, T, T>
inline T fold_left(const std::ranges::range auto& range, BinaryOperation op)
{
  auto it = std::ranges::begin(range);
  auto end = std::ranges::end(range);
  
  if (it == end)
  {
    throw std::invalid_argument("fold_left: input range is empty");
  }
  
  T result = *it;
  ++it;
  
  for (; it != end; ++it)
  {
    result = op(std::move(result), *it);
  }
  
  return result;
}

} // namespace mcrl2::pbes_system::detail