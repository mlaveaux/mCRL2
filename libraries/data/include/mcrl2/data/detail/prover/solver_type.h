// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/detail/prover/solver_type.h
/// \brief Enumeration for the types of solvers

#ifndef MCRL2_DATA_DETAIL_PROVER_SOLVER_TYPE_H
#define MCRL2_DATA_DETAIL_PROVER_SOLVER_TYPE_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::data::detail
{

/// \brief The enumeration type smt_solver_type enumerates all available SMT solvers.
enum smt_solver_type
{
  solver_type_cvc,
  solver_type_z3
};

/// \brief standard conversion from string to solver type
inline
smt_solver_type parse_solver_type(const std::string& s)
{
  if (s == "cvc")
  {
    return solver_type_cvc;
  }
  else if (s == "z3")
  {
    return solver_type_z3;
  }
  else
  {
    throw mcrl2::runtime_error("unknown solver type " + s);
  }
}

/// \brief standard conversion from stream to solver type
inline
std::istream& operator>>(std::istream& is, smt_solver_type& s)
{
  try
  {
    std::string str;
    is >> str;
    s = parse_solver_type(str);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

/// \brief standard conversion from solver type to string
inline
std::string print_solver_type(const smt_solver_type s)
{
  switch(s)
  {
    case solver_type_cvc: return "cvc";
    case solver_type_z3: return "z3";
    default: throw mcrl2::runtime_error("unknown solver type");
  }
}

/// \brief standard conversion from solvert type to stream
inline std::ostream& operator<<(std::ostream& os,smt_solver_type s)
{
  os << print_solver_type(s);
  return os;
}

/// \brief description of solver type
inline
std::string description(const smt_solver_type s)
{
  switch(s)
  {
    case solver_type_cvc: return "the SMT solver CVC3";
    case solver_type_z3: return "the SMT solver Z3";
    default: throw mcrl2::runtime_error("unknown solver type");
  }
}

} // namespace mcrl2::data::detail

#endif // MCRL2_DATA_DETAIL_PROVER_SOLVER_TYPE_H
