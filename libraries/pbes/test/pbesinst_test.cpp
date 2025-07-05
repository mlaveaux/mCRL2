// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbesinst_test.cpp
/// \brief Add your file description here.

#define BOOST_TEST_MODULE pbesinst_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/lps/detail/test_input.h"
#include "mcrl2/modal_formula/detail/test_input.h"
#include "mcrl2/modal_formula/parse.h"
#include "mcrl2/pbes/is_bes.h"
#include "mcrl2/pbes/lps2pbes.h"
#include "mcrl2/pbes/pbesinst_finite_algorithm.h"
#include "mcrl2/pbes/pbesinst_symbolic.h"
#include "mcrl2/pbes/txt2pbes.h"

using namespace mcrl2;
using namespace mcrl2::pbes_system;

inline
pbes pbesinst_lazy(const pbes& p)
{
  pbes q = p;
  pbesinst_algorithm algorithm(q.data());
  algorithm.run(q);
  return algorithm.get_result();
}

inline
pbes pbesinst_finite(const pbes& p)
{
  pbes q = p;
  pbesinst_finite_algorithm algorithm(data::jitty);
  algorithm.run(q);
  return q;
}

std::string test1 =
  "pbes                                                                              \n"
  "                                                                                  \n"
  "nu X(b:Bool, n:Nat) = (val(b) => X(!b, n)) && (val(!b) => X(!b, n+1));            \n"
  "                                                                                  \n"
  "init X(true,0);                                                                   \n"
  ;

std::string test2 =
  "pbes                                                                              \n"
  "                                                                                  \n"
  "nu X(b:Bool, n:Nat) = forall c:Bool. X(c,n);                                      \n"
  "                                                                                  \n"
  "init X(true,0);                                                                   \n"
  ;

std::string test3 =
  "pbes                                                                              \n"
  "                                                                                  \n"
  "nu X(b:Bool, n:Nat) = exists c:Bool. X(c,n+1);                                    \n"
  "                                                                                  \n"
  "init X(true,0);                                                                   \n"
  ;

std::string test4 =
  "pbes                                                                              \n"
  "                                                                                  \n"
  "nu X(b:Bool, n:Nat) = val(b && n < 10) => X(!b,n+1);                              \n"
  "                                                                                  \n"
  "init X(true,0);                                                                   \n"
  ;

std::string test5 =
  "sort D = struct d1 | d2;                                                          \n"
  "                                                                                  \n"
  "pbes                                                                              \n"
  "                                                                                  \n"
  "nu X(d:D, n:Nat) = val(d == d1 && n < 10) => X(d2,n+1);                           \n"
  "                                                                                  \n"
  "init X(d1,0);                                                                     \n"
  ;

std::string test6 =
  "pbes                                                                              \n"
  "nu X(b:Bool) = forall c:Bool. X(if (c,!c,c));                                     \n"
  "                                                                                  \n"
  "init X(true);                                                                     \n"
  ;

std::string test7 =
  "sort Enum2 = struct e1_5 | e0_5;                                                  \n"
  "                                                                                  \n"
  "map                                                                               \n"
  "                                                                                  \n"
  "     C5_fun2: Enum2 # Enum2 # Enum2  -> Enum2;                                    \n"
  "     C5_fun1: Enum2 # Nat # Nat  -> Nat;                                          \n"
  "                                                                                  \n"
  "var  y23,y22,y21,x5,y14,y13,y12,y11,y10,x2,e3,e2,e1: Enum2;                       \n"
  "     y20,y19,y18,x4,y9,y8,y7,y6,y5,x1: Nat;                                       \n"
  "     y17,y16,y15,x3,y4,y3,y2,y1,y,x: Bool;                                        \n"
  "eqn                                                                               \n"
  "     C5_fun2(e0_5, y14, y13)  =  y14;                                             \n"
  "     C5_fun2(e1_5, y14, y13)  =  y13;                                             \n"
  "     C5_fun2(e3, x2, x2)  =  x2;                                                  \n"
  "     C5_fun1(e0_5,  y6, y5)  =  y5;                                               \n"
  "     C5_fun1(e1_5,  y6, y5)  =  y6;                                               \n"
  "     C5_fun1(e2,  x1, x1)  =  x1;                                                 \n"
  "                                                                                  \n"
  "pbes nu X(s3_P: Enum2,  n_P: Nat) =                                               \n"
  "                                                                                  \n"
  "(forall e: Enum2.  X(C5_fun2(e, e, e1_5), C5_fun1(e, 0, n_P))                     \n"
  "                                                                                  \n"
  ")                                                                                 \n"
  "                                                                                  \n"
  "                                                                                  \n"
  ";                                                                                 \n"
  "                                                                                  \n"
  "init X(e1_5,  0);                                                                 \n"
  ;

std::string test8 =
  "pbes                                                                              \n"
  "                                                                                  \n"
  "nu X(b:Bool) = val(b) && Y(!b);                                                   \n"
  "                                                                                  \n"
  "mu Y(c:Bool) = forall d:Bool. X(d && c) || Y(d);                                  \n"
  "                                                                                  \n"
  "init X(true);                                                                     \n"
  ;

// This pbes triggered a garbage collection problem, that has been solved.
std::string random1 =
  "pbes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        \n"
  "nu X0(c:Bool, n:Nat) = (forall n:Nat.((val(n < 3)) && (((val(n < 3)) || (exists m:Nat.((val(m < 3)) || (X3(m + 1, m > 0))))) && ((forall m:Nat.((val(m < 3)) && (!X2(m + 1, 1)))) => ((val(c)) || (val(n < 3))))))) || ((val(false)) || (X0(false, n + 1)));                                                                                                                                                                                                                                                \n"
  "nu X1(b:Bool) = (!(!(forall k:Nat.((val(k < 3)) && ((forall k:Nat.((val(k < 3)) && ((X2(1, k + 1)) && (val(false))))) || ((X4(k > 0, k + 1)) && (X1(k > 1)))))))) && (!(forall m:Nat.((val(m < 3)) && (((val(m < 2)) && (val(m > 0))) && (val(true))))));                                                                                                                                                                                                                                                   \n"
  "mu X2(m:Nat, n:Nat) = (((val(m < 2)) && (X4(m == n, n + 1))) || ((val(false)) || ((val(true)) => (X0(n == m, 0))))) || (forall k:Nat.((val(k < 3)) && (exists m:Nat.((val(m < 3)) || ((val(n < 2)) && (X2(m + 1, m + 1)))))));                                                                                                                                                                                                                                                                              \n"
  "nu X3(n:Nat, c:Bool) = ((forall k:Nat.((val(k < 3)) && (!((forall m:Nat.((val(m < 3)) && (val(n > 0)))) => (val(c)))))) && ((X3(0, n < 3)) && (exists m:Nat.((val(m < 3)) || ((!(exists n:Nat.((val(n < 3)) || (val(m < 2))))) && (X2(n + 1, m + 1))))))) || ((!(!X1(n > 0))) || (val(false)));                                                                                                                                                                                                             \n"
  "nu X4(c:Bool, n:Nat) = (((exists m:Nat.((val(m < 3)) || (val(m > 0)))) && ((!(!X0(n < 3, 0))) || (!((val(n > 0)) => (!X2(0, 0)))))) => (forall k:Nat.((val(k < 3)) && (!(forall n:Nat.((val(n < 3)) && (val(n < 2)))))))) => (!(forall m:Nat.((val(m < 3)) && ((val(c)) && (forall m:Nat.((val(m < 3)) && (!X3(n + 1, false))))))));                                                                                                                                                                        \n"
  "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            \n"
  "init X0(true, 0);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \n"
  ;

// This pbes triggered a garbage collection problem, that has been solved.
std::string random2 =
  "pbes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        \n"
  "mu X0(m:Nat, b:Bool) = (X0(m + 1, m > 0)) && (((forall m:Nat.((val(m < 3)) && (forall k:Nat.((val(k < 3)) && (!((val(k > 1)) && (val(false)))))))) && (exists n:Nat.((val(n < 3)) || (!(!(val(n > 1))))))) || (exists n:Nat.((val(n < 3)) || (forall n:Nat.((val(n < 3)) && (exists n:Nat.((val(n < 3)) || (((forall m:Nat.((val(m < 3)) && (val(m == n)))) => (X3(0, n > 1))) && (!(!X2(n == m)))))))))));                                                                                                 \n"
  "mu X1(c:Bool, b:Bool) = (((!X1(c, true)) || ((!X2(b)) || (val(b)))) && (forall k:Nat.((val(k < 3)) && (((val(c)) => (X4(c, 1))) => (val(true)))))) => (val(true));                                                                                                                                                                                                                                                                                                                                          \n"
  "nu X2(b:Bool) = (!(!(((exists m:Nat.((val(m < 3)) || (val(m < 3)))) => ((X0(1, true)) => (!(val(false))))) && ((!((!X1(b, true)) && (!X3(0, false)))) => (val(false)))))) => (exists m:Nat.((val(m < 3)) || (forall k:Nat.((val(k < 3)) && (exists n:Nat.((val(n < 3)) || (forall n:Nat.((val(n < 3)) && (exists n:Nat.((val(n < 3)) || (forall n:Nat.((val(n < 3)) && (val(false))))))))))))));                                                                                                            \n"
  "mu X3(m:Nat, c:Bool) = (exists m:Nat.((val(m < 3)) || ((!(val(c))) && (forall k:Nat.((val(k < 3)) && (exists m:Nat.((val(m < 3)) || (forall k:Nat.((val(k < 3)) && (val(false))))))))))) || (exists m:Nat.((val(m < 3)) || (exists n:Nat.((val(n < 3)) || ((!((exists k:Nat.((val(k < 3)) || (!X4(false, k + 1)))) && ((!X1(m > 0, n < 2)) || ((val(n > 0)) || (val(c)))))) && (exists k:Nat.((val(k < 3)) || (exists n:Nat.((val(n < 3)) || (exists k:Nat.((val(k < 3)) || (X3(k + 1, k == n))))))))))))); \n"
  "nu X4(b:Bool, n:Nat) = (exists m:Nat.((val(m < 3)) || ((val(b)) || ((val(n > 1)) && ((X1(m > 1, n == m)) || (val(m > 0))))))) && ((exists k:Nat.((val(k < 3)) || (!(!((!X3(1, n > 0)) && (!(val(k == n)))))))) => (forall n:Nat.((val(n < 3)) && (exists k:Nat.((val(k < 3)) || (!(!(X4(n > 1, 0)))))))));                                                                                                                                                                                                  \n"
  "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            \n"
  "init X0(0, true);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \n"
  ;

// This pbes triggered an error with pbesinst finite
std::string random3 =
  "pbes                                                                                                                                                                                                                                                                                                                                                              \n"
  "mu X0(n:Nat, c:Bool) = ((!(((val(n < 2)) && (!X1)) && (exists k:Nat.((val(k < 3)) || (val(c)))))) && (((forall n:Nat.((val(n < 3)) && (!X4(n > 1)))) || (!(val(n > 1)))) => (val(true)))) && (exists n:Nat.((val(n < 3)) || (exists k:Nat.((val(k < 3)) || (forall m:Nat.((val(m < 3)) && (X0(k + 1, n < 3))))))));                                               \n"
  "mu X1 = (((!((val(true)) => (X2))) => (!(!(!(!X1))))) && (((val(false)) || (X3(1))) && (forall n:Nat.((val(n < 3)) && (val(true)))))) || (forall k:Nat.((val(k < 3)) && (!(forall k:Nat.((val(k < 3)) && (exists m:Nat.((val(m < 3)) || (val(k < 2)))))))));                                                                                                      \n"
  "mu X2 = ((exists k:Nat.((val(k < 3)) || (exists k:Nat.((val(k < 3)) || (exists m:Nat.((val(m < 3)) || ((val(m == k)) => (X1)))))))) => (((val(true)) => (!X2)) || (exists k:Nat.((val(k < 3)) || (val(false)))))) => (!((val(false)) => (!(X3(0)))));                                                                                                             \n"
  "mu X3(n:Nat) = ((exists n:Nat.((val(n < 3)) || (exists n:Nat.((val(n < 3)) || (forall n:Nat.((val(n < 3)) && (val(false)))))))) && (forall m:Nat.((val(m < 3)) && (((forall m:Nat.((val(m < 3)) && (exists m:Nat.((val(m < 3)) || (exists m:Nat.((val(m < 3)) || (!(val(n > 0))))))))) && ((!X4(m > 1)) || (!X2))) || (val(n < 2)))))) => ((val(n < 2)) => (X1)); \n"
  "nu X4(b:Bool) = ((val(true)) => ((forall m:Nat.((val(m < 3)) && (val(false)))) => ((X1) || ((!X3(0)) => (forall n:Nat.((val(n < 3)) && (val(n > 1)))))))) && (forall m:Nat.((val(m < 3)) && (forall n:Nat.((val(n < 3)) && ((X2) && (val(false)))))));                                                                                                            \n"
  "                                                                                                                                                                                                                                                                                                                                                                  \n"
  "init X0(0, true);                                                                                                                                                                                                                                                                                                                                                 \n"
  ;

void test_pbes(const std::string& pbes_spec, bool test_finite, bool test_lazy)
{
  pbes p = txt2pbes(pbes_spec);
  std::cout << "------------------------------\n" << pbes_system::pp(p) << std::endl;
  if (!p.is_closed())
  {
    std::cout << "ERROR: the pbes is not closed!" << std::endl;
    return;
  }

  if (test_finite)
  {
    std::cout << "FINITE" << std::endl;
    try
    {
      using namespace pbes_system;
      pbes q1 = ::pbesinst_finite(p);
      std::cout << pbes_system::pp(q1) << std::endl;
    }
    catch (mcrl2::runtime_error& e)
    {
      std::cout << "pbesinst failed: " << e.what() << std::endl;
    }
  }

  if (test_lazy)
  {
    std::cout << "LAZY" << std::endl;
    try
    {
      using namespace pbes_system;
      pbes q1 = pbesinst_lazy(p);
      std::cout << pbes_system::pp(q1) << std::endl;
    }
    catch (mcrl2::runtime_error& e)
    {
      std::cout << "pbesinst failed: " << e.what() << std::endl;
    }
  }
}

BOOST_AUTO_TEST_CASE(test_pbesinst)
{
  test_pbes(test1, true, false);
  test_pbes(test2, true, true);
  test_pbes(test3, true, false);
  test_pbes(test4, true, true);
  test_pbes(test5, true, true);
  test_pbes(test6, true, true);
  test_pbes(test7, true, true);
  test_pbes(test8, true, true);
  test_pbes(random3, false, true);
}

BOOST_AUTO_TEST_CASE(test_pbesinst_finite)
{
  pbes p = txt2pbes(random3);
  pbes q = pbesinst_finite(p);
  std::cerr << pbes_system::pp(q) << std::endl;

  std::string text =
    "sort D = struct d1 | d2;                                        \n"
    "                                                                \n"
    "pbes                                                            \n"
    "                                                                \n"
    "nu X(d:D) = (val(d == d1) && X(d2)) || (val(d == d2) && X(d1)); \n"
    "                                                                \n"
    "init X(d1);                                                     \n"
    ;
  pbes p1 = txt2pbes(text);
  pbesinst_finite_algorithm algorithm(data::jitty);
  pbesinst_variable_map variable_map = mcrl2::pbes_system::detail::parse_pbes_parameter_map(p1, "X(*:D)");
  algorithm.run(p1, variable_map);
}

BOOST_AUTO_TEST_CASE(test_abp_no_deadlock)
{
  lps::specification spec=remove_stochastic_operators(lps::linearise(lps::detail::ABP_SPECIFICATION()));
  state_formulas::state_formula formula = state_formulas::parse_state_formula(lps::detail::NO_DEADLOCK(), spec, false);
  bool timed = false;
  pbes p = lps2pbes(spec, formula, timed);
  data::rewriter::strategy rewriter_strategy = data::jitty;
  bool print_equations = true;
  pbes_system::pbesinst_algorithm algorithm(p.data(), rewriter_strategy, print_equations);
  algorithm.run(p);
  pbes q = algorithm.get_result();
  std::cout << "--- ABP ---" << std::endl;
  std::cout << pbes_system::pp(q) << std::endl;
  BOOST_CHECK(is_bes(q));
}

// Example supplied by Tim Willemse, 23-05-2011
BOOST_AUTO_TEST_CASE(test_functions)
{
  std::string text =
    "sort D = struct one | two;           \n"
    "                                     \n"
    "map  f: D -> D;                      \n"
    "                                     \n"
    "eqn  f  =  lambda x: D. one;         \n"
    "                                     \n"
    "pbes nu X(d: D, g: D -> D) =         \n"
    "       forall e: D. X(e, g[e -> e]); \n"
    "                                     \n"
    "init X(one, f);                      \n"
    ;
  pbes p = txt2pbes(text);
  data::rewriter::strategy rewrite_strategy = data::jitty;
  pbesinst_finite_algorithm algorithm(rewrite_strategy);
  mcrl2::pbes_system::detail::pbes_parameter_map parameter_map = mcrl2::pbes_system::detail::parse_pbes_parameter_map(p, "X(*:D)");
  algorithm.run(p, parameter_map);
}

void test_pbesinst_symbolic(const std::string& text)
{
  pbes p;
  p = txt2pbes(text);
  pbesinst_symbolic_algorithm algorithm(p);
  algorithm.run();
}

BOOST_AUTO_TEST_CASE(test_pbesinst_symbolic1)
{
  test_pbesinst_symbolic(test2);
  test_pbesinst_symbolic(test4);
  test_pbesinst_symbolic(test5);
  test_pbesinst_symbolic(test6);
}

#ifdef MCRL2_EXTENDED_TESTS
BOOST_AUTO_TEST_CASE(test_pbesinst_slow)
{
  test_pbes(random1, false, true);
  test_pbes(random2, false, true);
}
#endif // MCRL2_EXTENDED_TESTS
