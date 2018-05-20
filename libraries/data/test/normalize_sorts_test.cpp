// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file normalize_sorts_test.cpp
/// \brief Test for normalizing sorts.

#include "mcrl2/data/data_specification.h"
#include "mcrl2/data/normalize_sorts.h"
#include "mcrl2/data/parse.h"
#include <algorithm>
#include <boost/test/minimal.hpp>
#include <iterator>
#include <set>
#include <vector>

using namespace mcrl2;
using namespace mcrl2::data;

void test_normalize_sorts()
{
  std::string DATASPEC =
    "sort Bit = struct e0 | e1;      \n"
    "     AbsBit = struct arbitrary; \n"
    "                                \n"
    "map  inv: Bit -> Bit;           \n"
    "     h: Bit -> AbsBit;          \n"
    "                                \n"
    "eqn  inv(e0)  =  e1;            \n"
    "     inv(e1)  =  e0;            \n"
    ;

  data_specification dataspec = parse_data_specification(DATASPEC);

  data::function_symbol f;
  f = parse_function_symbol("abseq : AbsBit # AbsBit -> Set(Bool)", DATASPEC);
  dataspec.add_mapping(f);
  f = parse_function_symbol("absinv : AbsBit -> Set(AbsBit)", DATASPEC);
  dataspec.add_mapping(f);

  data_equation_vector equations = dataspec.user_defined_equations();
  data::normalize_sorts(equations, dataspec);
}

// The test below checks whether the calculation of a confluent and terminating rewrite system for types using
// Knuth-Bendix completion is efficient. Aleksi Peltonen showed in the spring of 2018 that Knuth-Bendix completion
// was exponential, causing the example below not to terminate within reasonable time. 
void test_apply_knuth_bendix_completion_on_sorts()
{
  std::string DATASPEC =
    "sort A_t = Nat; B_t = Nat; C_t = Nat; D_t = Nat; E_t = Nat; F_t = Nat; G_t = Nat; H_t = Nat; I_t = Nat; J_t=Nat; K_t=Nat; \n"
    "     L_t = Nat; M_t = Nat; N_t = Nat; O_t = Nat;\n"
    "     S_t = struct s( A:A_t, B:B_t, C:C_t, D:D_t, E:E_t, F:F_t, G:G_t, H:H_t, I:I_t, J:J_t, K:K_t, L:L_t, M:M_t, N:N_t, O:O_t); \n";

  data_specification dataspec = parse_data_specification(DATASPEC);

  data_equation_vector equations = dataspec.user_defined_equations();
  data::normalize_sorts(equations, dataspec);
}

int test_main(int argc, char* argv[])
{
  test_normalize_sorts();
  test_apply_knuth_bendix_completion_on_sorts();

  return 0;
}
