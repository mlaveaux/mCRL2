// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbes_rewriter_type.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PBES_REWRITER_TYPE_H
#define MCRL2_PBES_PBES_REWRITER_TYPE_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system
{

/// \brief An enumerated type for the available pbes rewriters
enum class pbes_rewriter_type
{
  simplify,
  quantifier_all,
  quantifier_finite,
  quantifier_inside,
  quantifier_one_point,
  prover,
  pfnf,
  ppg,
  srf,
  pre_srf,
  bqnf_quantifier
};

/// \brief Parses a pbes rewriter type
inline pbes_rewriter_type parse_pbes_rewriter_type(const std::string& type)
{
  if (type == "simplify")
  {
    return pbes_rewriter_type::simplify;
  }
  if (type == "quantifier-all")
  {
    return pbes_rewriter_type::quantifier_all;
  }
  if (type == "quantifier-finite")
  {
    return pbes_rewriter_type::quantifier_finite;
  }
  if (type == "quantifier-inside")
  {
    return pbes_rewriter_type::quantifier_inside;
  }
  if (type == "quantifier-one-point")
  {
    return pbes_rewriter_type::quantifier_one_point;
  }
  if (type == "prover")
  {
    return pbes_rewriter_type::prover;
  }
  if (type == "pfnf")
  {
    return pbes_rewriter_type::pfnf;
  }
  if (type == "ppg")
  {
    return pbes_rewriter_type::ppg;
  }
  if (type == "srf")
  {
    return pbes_rewriter_type::srf;
  }
  if (type == "pre-srf")
  {
    return pbes_rewriter_type::pre_srf;
  }

  if (type == "bqnf-quantifier")
  {
    return pbes_rewriter_type::bqnf_quantifier;
  }
  throw mcrl2::runtime_error("unknown pbes rewriter option " + type);
}

/// \brief Prints a pbes rewriter type
inline std::string print_pbes_rewriter_type(const pbes_rewriter_type type)
{
  switch (type)
  {
  case pbes_rewriter_type::simplify:
    return "simplify";
  case pbes_rewriter_type::quantifier_all:
    return "quantifier-all";
  case pbes_rewriter_type::quantifier_finite:
    return "quantifier-finite";
  case pbes_rewriter_type::quantifier_inside:
    return "quantifier-inside";
  case pbes_rewriter_type::quantifier_one_point:
    return "quantifier-one-point";
  case pbes_rewriter_type::prover:
    return "prover";
  case pbes_rewriter_type::pfnf:
    return "pfnf";
  case pbes_rewriter_type::ppg:
    return "ppg";
  case pbes_rewriter_type::bqnf_quantifier:
    return "bqnf-quantifier";
  case pbes_rewriter_type::srf:
    return "srf";
  case pbes_rewriter_type::pre_srf:
    return "pre-srf";    
  }

  return "unknown pbes rewriter";
}

/// \brief Returns a description of a pbes rewriter
inline std::string description(const pbes_rewriter_type type)
{
  switch (type)
  {
  case pbes_rewriter_type::simplify:
    return "for simplification";
  case pbes_rewriter_type::quantifier_all:
    return "for eliminating all quantifiers";
  case pbes_rewriter_type::quantifier_finite:
    return "for eliminating finite quantifier variables";
  case pbes_rewriter_type::quantifier_inside:
    return "for pushing quantifiers inside";
  case pbes_rewriter_type::quantifier_one_point:
    return "for one point rule quantifier elimination";
  case pbes_rewriter_type::prover:
    return "for rewriting using a prover";
  case pbes_rewriter_type::pfnf:
    return "for rewriting into PFNF normal form";
  case pbes_rewriter_type::ppg:
    return "for rewriting into Parameterised Parity Game form";
  case pbes_rewriter_type::srf:
    return "for rewriting into SRF normal form";
  case pbes_rewriter_type::pre_srf:
    return "for rewriting into pre-SRF normal form";
  case pbes_rewriter_type::bqnf_quantifier:
    return "for rewriting quantifiers over conjuncts to conjuncts of quantifiers (experimental)";
  }

  throw mcrl2::runtime_error("unknown pbes rewriter");
}

/// \brief Stream operator for rewriter type
/// \param is An input stream
/// \param t A rewriter type
/// \return The input stream
inline std::istream& operator>>(std::istream& is, pbes_rewriter_type& t)
{
  std::string s;
  is >> s;
  try
  {
    t = parse_pbes_rewriter_type(s);
  }
  catch (const mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

inline std::ostream& operator<<(std::ostream& os, const pbes_rewriter_type t)
{
  os << print_pbes_rewriter_type(t);
  return os;
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_PBES_REWRITER_TYPE_H
