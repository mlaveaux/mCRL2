// Author(s): Jan Friso Groote. Based on bes_equation_limit.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/detail/res_equation_limit.h
/// \brief A global variable for counting the number of RES equations in presinst
/// and parity_game_generator. If the number of equations exceeds a limit, an
/// exception is thrown. A static template class variable is used, that can be
/// set from everywhere.

#ifndef MCRL2_PRES_DETAIL_RES_EQUATION_LIMIT_H
#define MCRL2_PRES_DETAIL_RES_EQUATION_LIMIT_H

#include <limits>
#include <stdexcept>

namespace mcrl2::pres_system::detail
{

template <class T> // note, T is only a dummy
struct res_equation_limit
{
  static std::size_t max_res_equations;
};

template <class T>
std::size_t res_equation_limit<T>::max_res_equations = (std::numeric_limits<std::size_t>::max)();

inline
void set_res_equation_limit(std::size_t size)
{
  res_equation_limit<std::size_t>::max_res_equations = size;
}

inline
void check_res_equation_limit(std::size_t size)
{
  if (size >= res_equation_limit<std::size_t>::max_res_equations)
  {
    throw std::out_of_range("Error: number of RES equations has exceeded the limit");
  }
}

} // namespace mcrl2::pres_system::detail

#endif // MCRL2_PRES_DETAIL_RES_EQUATION_LIMIT_H
