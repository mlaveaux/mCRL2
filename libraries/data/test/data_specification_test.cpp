// Author(s): Jeroen Keiren, Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file data_specification_test.cpp
/// \brief Basic regression test for data specifications.

#define BOOST_TEST_MODULE data_specification_test
#include <boost/test/included/unit_test.hpp>

#include "mcrl2/data/data_io.h"
#include "mcrl2/data/merge_data_specifications.h"
#include "mcrl2/data/parse.h"
#include "mcrl2/data/print.h"

using namespace mcrl2;
using namespace mcrl2::data;

template <typename Container1, typename Container2>
bool equal_content(Container1 const& c1, Container2 const& c2)
{
  std::set<typename Container1::value_type> s1(c1.begin(), c1.end());
  std::set<typename Container2::value_type> s2(c2.begin(), c2.end());

  if (s1 != s2)
  {
    std::clog << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
                      << "Detailed comparison:" << std::endl;
    std::clog << data::pp(c1) << std::endl;
    std::clog << data::pp(c2) << std::endl;
  }
  return s1 == s2;
}

bool detailed_compare_for_equality(data_specification const& left, data_specification const& right)
{
  std::set<sort_expression> left_sorts(left.sorts().begin(), left.sorts().end());
  std::set<sort_expression> right_sorts(right.sorts().begin(), right.sorts().end());
  BOOST_CHECK(left_sorts == right_sorts);

  std::set<data::function_symbol> left_constructors(left.constructors().begin(), left.constructors().end());
  std::set<data::function_symbol> right_constructors(right.constructors().begin(), right.constructors().end());
  BOOST_CHECK(left_constructors == right_constructors);

  std::set<data::function_symbol> left_mappings(left.mappings().begin(), left.mappings().end());
  std::set<data::function_symbol> right_mappings(right.mappings().begin(), right.mappings().end());
  BOOST_CHECK(left_mappings == right_mappings);

  BOOST_CHECK(left.equations() == right.equations());

  if (/*(left_aliases != right_aliases)*/
       (left_sorts != right_sorts)
      || (left_constructors != right_constructors)
      || (left_mappings != right_mappings)
      || (left.equations() != right.equations()))
  {
    std::clog << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
                  << "Specification detailed comparison:" << std::endl;

    if (left_sorts != right_sorts)
    {
      std::clog << "Sorts (left)  " << data::pp(left.sorts()) << std::endl;
      std::clog << "Sorts (right) " << data::pp(right.sorts()) << std::endl;
    }
    if (left_constructors != right_constructors)
    {
      std::clog << "Constructors (left)  " << data::pp(left.constructors()) << std::endl;
      std::clog << "Constructors (right) " << data::pp(right.constructors()) << std::endl;
    }
    if (left_mappings != right_mappings)
    {
      std::clog << "Mappings (left)  " << data::pp(left.mappings()) << std::endl;
      std::clog << "Mappings (right) " << data::pp(right.mappings()) << std::endl;
    }
    if (left.equations() != right.equations())
    {
      std::clog << "Equations (left)  " << data::pp(left.equations()) << std::endl;
      std::clog << "Equations (right) " << data::pp(right.equations()) << std::endl;
    }

    std::clog << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    return false;
  }
  return true;
}

bool compare_for_equality(data_specification const& left, data_specification const& right)
{
  return detailed_compare_for_equality(left, right);
}

void test_sorts()
{
  std::clog << "test_sorts" << std::endl;

  basic_sort s("S");
  basic_sort s0("S0");
  alias s1(s,basic_sort("S1"));

  std::vector< sort_expression > sort_list;
  sort_list.push_back(basic_sort("S1"));
  sort_list.push_back(s0);

  data_specification spec;
  spec.add_sort(s);
  spec.add_sort(s0);
  spec.add_alias(s1);
  data_specification spec1;
  spec1.add_alias(s1);
  spec1.add_sort(s0);
  spec1.add_sort(s);

  BOOST_CHECK(compare_for_equality(spec, spec1));

  basic_sort s2("S2");
  sort_expression_vector s2l {s2};
  spec.add_context_sort(s2);
  spec1.add_context_sorts(s2l);
  BOOST_CHECK(compare_for_equality(spec, spec1));

  std::for_each(s2l.begin(), s2l.end(), [&spec](const sort_expression& s){ spec.remove_sort(s); });
  spec1.remove_sort(s2);
  compare_for_equality(spec, spec1);
}

void test_aliases()
{
  std::clog << "test_aliases" << std::endl;

  basic_sort s("S");
  basic_sort t("T");
  alias      s1(basic_sort("S1"), s);
  alias      s2(basic_sort("S2"), s);

  data_specification spec;

  std::set<basic_sort> sorts = { s, t };
  std::for_each(sorts.begin(), sorts.end(), [&spec](const basic_sort& s){ spec.add_sort(s); });

}

void test_constructors()
{
  std::clog << "test_constructors" << std::endl;

  basic_sort s("S");
  basic_sort s0("S0");
  function_sort s0s(sort_expression_list({s0}),s);
  data::function_symbol f("f", s);
  data::function_symbol g("g", s0s);
  data::function_symbol h("h", s0);
  function_symbol_vector fgl {f,g};
  function_symbol_vector hl {h};
  function_symbol_vector fghl {f,g,h};

  data_specification spec;
  spec.add_sort(s);
  spec.add_sort(s0);

  spec.add_constructor(f);
  spec.add_constructor(g);
  spec.add_constructor(h);

  data_specification spec1(spec);
  std::for_each(fghl.begin(), fghl.end(), [&spec1](const function_symbol& f){ spec1.add_constructor(f); });

  function_symbol_vector constructors(spec.constructors());
  BOOST_CHECK(spec.constructors(s) == fgl);
std::cerr << "#constructors " << spec.constructors().size() << "\n";
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  BOOST_CHECK(constructors.size() == 9); // f,g,h, true, false.
#else
  BOOST_CHECK(constructors.size() == 7); // f,g,h, true, false.
#endif
  BOOST_CHECK(std::find(constructors.begin(), constructors.end(), f) != constructors.end());
  BOOST_CHECK(std::find(constructors.begin(), constructors.end(), g) != constructors.end());
  BOOST_CHECK(std::find(constructors.begin(), constructors.end(), h) != constructors.end());

  BOOST_CHECK(compare_for_equality(spec, spec1));
  BOOST_CHECK(spec.constructors() == spec1.constructors());
  BOOST_CHECK(spec.constructors(s) == fgl);
  BOOST_CHECK(spec.constructors(s0) == hl);
  BOOST_CHECK(spec1.constructors(s) == fgl);
  BOOST_CHECK(spec1.constructors(s0) == hl);
  spec.add_constructor(data::function_symbol("i", s0));
  data::function_symbol i("i", s0);
  spec.remove_constructor(i);
  BOOST_CHECK(compare_for_equality(spec, spec1));

  spec.add_constructor(i);
  function_symbol_vector il {i};
  std::for_each(il.begin(), il.end(), [&spec1](const function_symbol& f){ spec1.add_constructor(f); });

  BOOST_CHECK(compare_for_equality(spec, spec1));

  spec.remove_constructor(i);
  std::for_each(il.begin(), il.end(), [&spec1](const function_symbol& f){ spec1.remove_constructor(f); });

  BOOST_CHECK(compare_for_equality(spec, spec1));
}

void test_functions()
{
  std::clog << "test_functions" << std::endl;

  basic_sort s("S");
  basic_sort s0("S0");
  function_sort s0s(sort_expression_list({s0}), s);
  data::function_symbol f("f", s);
  data::function_symbol g("g", s0s);
  data::function_symbol h("h", s0);

  function_symbol_vector fgl {f,g};
  function_symbol_vector hl {h};
  function_symbol_vector fghl {f,g,h};

  data_specification spec;
  spec.add_sort(s);
  spec.add_sort(s0);
  spec.add_mapping(f);
  spec.add_mapping(g);
  spec.add_mapping(h);

  data_specification spec1(spec);
  std::for_each(fghl.begin(), fghl.end(), [&spec1](const function_symbol& f){ spec1.add_mapping(f); } );

std::cerr << "#mappings " << spec.mappings().size() << "\n";
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  BOOST_CHECK(spec.mappings().size() == 107);
#else
  BOOST_CHECK(spec.mappings().size() == 54);
#endif

  function_symbol_vector mappings(spec.mappings());
  BOOST_CHECK(std::find(mappings.begin(), mappings.end(), f) != mappings.end());
  BOOST_CHECK(std::find(mappings.begin(), mappings.end(), g) != mappings.end());
  BOOST_CHECK(std::find(mappings.begin(), mappings.end(), h) != mappings.end());

  BOOST_CHECK(compare_for_equality(spec, spec1));
  BOOST_CHECK(spec.mappings(s).size() == 3);
  BOOST_CHECK(std::find(spec.mappings(s).begin(), spec.mappings(s).end(), f) != spec.mappings(s).end());
  BOOST_CHECK(std::find(spec.mappings(s).begin(), spec.mappings(s).end(), g) != spec.mappings(s).end());
  BOOST_CHECK(std::find(spec.mappings(s0).begin(), spec.mappings(s0).end(), h) != spec.mappings(s0).end());
  BOOST_CHECK(spec1.mappings(s).size() == 3);
  BOOST_CHECK(std::find(spec1.mappings(s).begin(), spec1.mappings(s).end(), f) != spec1.mappings(s).end());
  BOOST_CHECK(std::find(spec1.mappings(s).begin(), spec1.mappings(s).end(), g) != spec1.mappings(s).end());
  BOOST_CHECK(std::find(spec1.mappings(s0).begin(), spec1.mappings(s0).end(), h) != spec1.mappings(s0).end());

  data::function_symbol i("i", s0);
  spec.add_mapping(i);
  function_symbol_vector il {i};
  std::for_each(il.begin(), il.end(), [&spec1](const function_symbol& f){ spec1.add_mapping(f); });
  compare_for_equality(spec, spec1);

  std::for_each(il.begin(), il.end(), [&spec](const function_symbol& f){ spec.remove_mapping(f); });
  spec1.remove_mapping(i);
  compare_for_equality(spec, spec1);
}

void test_equations()
{
  std::clog << "test_equations" << std::endl;
  basic_sort s("S");
  basic_sort s0("S0");
  function_sort s0s(sort_expression_list({s0}), s);
  data::function_symbol f("f", s0s);
  variable x("x", s0);
  data_expression_vector xel {x};
  application fx(f, xel.begin(), xel.end());
  variable_vector xl {x};
  data_equation fxx(xl, x, fx, x);

  data_specification spec;
  data_specification spec1;
  spec.add_sort(s);
  spec.add_sort(s0);

  spec1 = spec;
  BOOST_CHECK(compare_for_equality(spec, spec1));
  spec.add_equation(fxx);
  data_equation_vector fxxl {fxx};
  std::for_each(fxxl.begin(), fxxl.end(), [&spec1](const data_equation& eqn){ spec1.add_equation(eqn); });

  BOOST_CHECK(compare_for_equality(spec, spec1));

  data_equation fxf(xl, x, fx, f);
  data_equation_vector fxfl {fxf};
  spec.add_equation(fxf);
  std::for_each(fxfl.begin(), fxfl.end(), [&spec1](const data_equation& eqn){ spec1.add_equation(eqn); });

  BOOST_CHECK(compare_for_equality(spec, spec1));

  data_equation_vector result = find_equations(spec, f);
  BOOST_CHECK(result.size() == 2);
  BOOST_CHECK(std::find(result.begin(), result.end(), fxf) != result.end());
  BOOST_CHECK(std::find(result.begin(), result.end(), fxx) != result.end());
  std::for_each(fxfl.begin(), fxfl.end(), [&spec](const data_equation& eqn){ spec.remove_equation(eqn); });
  spec1.remove_equation(fxf);
  BOOST_CHECK(compare_for_equality(spec, spec1));
}

void test_is_certainly_finite()
{
  std::clog << "test_is_certainly_finite" << std::endl;
  basic_sort s("S");
  basic_sort s0("S0");
  function_sort s0s0(sort_expression_list({s0}), s0);
  data::function_symbol f("f", s);
  data::function_symbol g("g", s0s0);
  variable x("x", s0);
  application gx(g, x);
  data_specification spec;
  spec.add_sort(s);
  spec.add_sort(s0);
  spec.add_constructor(f);
  spec.add_constructor(g);

  BOOST_CHECK(spec.is_certainly_finite(s));
  BOOST_CHECK(!spec.is_certainly_finite(s0));

  spec.add_context_sort(sort_real::real_());
  BOOST_CHECK(spec.is_certainly_finite(sort_bool::bool_()));
  BOOST_CHECK(!spec.is_certainly_finite(sort_pos::pos()));
  BOOST_CHECK(!spec.is_certainly_finite(sort_nat::nat()));
  BOOST_CHECK(!spec.is_certainly_finite(sort_int::int_()));
  BOOST_CHECK(!spec.is_certainly_finite(sort_real::real_()));

  basic_sort s1("S1");
  basic_sort s2("S2");
  spec.add_constructor(data::function_symbol("h", make_function_sort_(s1, s2)));
  spec.add_constructor(data::function_symbol("i", make_function_sort_(s2, s1)));

  spec.add_alias(alias(basic_sort("a"), s));
  spec.add_alias(alias(basic_sort("a0"), s0));
  spec.add_alias(alias(basic_sort("a1"), s1));

  BOOST_CHECK(spec.is_certainly_finite(normalize_sorts(basic_sort("a"),spec)));
  BOOST_CHECK(!spec.is_certainly_finite(normalize_sorts(basic_sort("a0"),spec)));
  BOOST_CHECK(!spec.is_certainly_finite(normalize_sorts(basic_sort("a1"),spec)));

  using namespace sort_list;

  BOOST_CHECK(!spec.is_certainly_finite(list(s)));
  BOOST_CHECK(!spec.is_certainly_finite(list(s0)));

  using namespace sort_set;

  BOOST_CHECK(spec.is_certainly_finite(set_(s)));
  BOOST_CHECK(!spec.is_certainly_finite(set_(s0)));
  BOOST_CHECK(!spec.is_certainly_finite(normalize_sorts(set_(s0),spec)));

  using namespace sort_bag;

  BOOST_CHECK(!spec.is_certainly_finite(bag(s)));
  BOOST_CHECK(!spec.is_certainly_finite(bag(s0)));
  BOOST_CHECK(spec.is_certainly_finite(make_function_sort_(s,s)));
  BOOST_CHECK(!spec.is_certainly_finite(make_function_sort_(s,s0)));
  BOOST_CHECK(!spec.is_certainly_finite(make_function_sort_(s0,s)));

  // structured sort
  /* This test should be reconsidered, as it does not work in this way.
 *   Having a test for structured sort is good.

  std::vector< data::structured_sort_constructor_argument > arguments;

  arguments.push_back(data::structured_sort_constructor_argument(s));
  arguments.push_back(data::structured_sort_constructor_argument(s0));
  arguments.push_back(data::structured_sort_constructor_argument(s1));

  std::vector< data::structured_sort_constructor > constructors;
  constructors.push_back(data::structured_sort_constructor("a", structured_sort_constructor_argument_list(arguments.begin(), arguments.begin() + 1)));
  constructors.push_back(data::structured_sort_constructor("b", structured_sort_constructor_argument_list(arguments.begin() + 1, arguments.begin() + 2)));
  constructors.push_back(data::structured_sort_constructor("b", structured_sort_constructor_argument_list(arguments.begin() + 2, arguments.begin() + 3)));

  structured_sort struct1(structured_sort_constructor_list(constructors.begin(), constructors.begin() + 1));
  structured_sort struct2(structured_sort_constructor_list(constructors.begin() + 1, constructors.begin() + 2));
  structured_sort struct3(structured_sort_constructor_list(constructors.begin() + 2, constructors.begin() + 3));
  structured_sort struct4(structured_sort_constructor_list(constructors.begin() + 0, constructors.begin() + 3));
  spec.add_sort(struct1);   add_sort does not work on complex sorts, since 2/6/2015.
  spec.add_sort(struct2);
  spec.add_sort(struct3);
  spec.add_sort(struct4);
  BOOST_CHECK(spec.is_certainly_finite(struct1));
  BOOST_CHECK(!spec.is_certainly_finite(struct2));
  BOOST_CHECK(!spec.is_certainly_finite(struct3));
  BOOST_CHECK(!spec.is_certainly_finite(struct4)); */
}

void test_constructor()
{
  std::clog << "test_constructor" << std::endl;
  std::string SPEC =
    "sort                      \n"
    "  D     = struct d1 | d2; \n"
    "  Error = struct e;       \n"
    "                          \n"
    ;
  data_specification data = parse_data_specification(SPEC);
  atermpp::aterm a = data::detail::data_specification_to_aterm(data);
  data_specification spec1(a);
}

template < typename Container, typename Expression >
bool search(Container const& container, Expression const& expression)
{
  return std::find(container.begin(), container.end(), expression) != container.end();
}

bool search_alias(const alias_vector& v, const sort_expression& s)
{
  for(const alias& a: v)
  {
    if(a.name() == s)
    {
      return true;
    }
  }
  return false;
}

BOOST_AUTO_TEST_CASE(test_system_defined)
{
  std::clog << "test_system_defined" << std::endl;

  data_specification specification;

  BOOST_CHECK(!specification.constructors(sort_bool::bool_()).empty());

  specification = parse_data_specification(
                    "sort S;"
                    "map f: Set(S);");

  BOOST_CHECK(search(specification.sorts(), sort_set::set_(basic_sort("S"))));
  // BOOST_CHECK(search(specification.sorts(), sort_fset::fset(basic_sort("S")))); MUST BE CHECKED ALSO?

  specification = parse_data_specification(
                    "sort D = Set(Nat);"
                    "sort E = D;"
                    "sort F = E;");
  BOOST_CHECK(specification.constructors(::sort_set::set_(sort_nat::nat())) == specification.constructors(basic_sort("D")));
  BOOST_CHECK(specification.constructors(normalize_sorts(basic_sort("D"),specification)) ==
              specification.constructors(normalize_sorts(basic_sort("E"),specification)));
  BOOST_CHECK(specification.mappings(normalize_sorts(basic_sort("D"),specification)) ==
              specification.mappings(normalize_sorts(basic_sort("E"),specification)));
  BOOST_CHECK(specification.constructors(normalize_sorts(basic_sort("D"),specification)) ==
              specification.constructors(normalize_sorts(basic_sort("F"),specification)));

  data_specification copy = specification;

  // A data specification that is constructed using data_specification_to_aterm is assumed not
  // not be type checked. This must be indicated explicitly.
  data_specification specification1=data_specification(data::detail::data_specification_to_aterm(copy));
  BOOST_CHECK(compare_for_equality(specification1,specification));

  specification = parse_data_specification(
                    "sort D = struct d(getBool : Bool)?is_d;"
                    "sort E = D;"
                    "sort F = E;");

  alias_vector aliases(specification.user_defined_aliases());
  BOOST_CHECK(search_alias(aliases, basic_sort("D")));
  BOOST_CHECK(search_alias(aliases, basic_sort("E")));
  BOOST_CHECK(search_alias(aliases, basic_sort("F")));
  BOOST_CHECK(specification.constructors(basic_sort("D")).size() == 1);

  BOOST_CHECK(specification.constructors(basic_sort("D")) == specification.constructors(basic_sort("E")));
  BOOST_CHECK(specification.mappings(basic_sort("D")) == specification.mappings(basic_sort("E")));
  BOOST_CHECK(specification.constructors(basic_sort("D")) == specification.constructors(basic_sort("F")));

  // Check for the non presence of function sort
  BOOST_CHECK(specification.mappings(make_function_sort_(basic_sort("D"), sort_bool::bool_())).empty());

  specification.add_mapping(data::function_symbol("f", make_function_sort_(sort_bool::bool_(), sort_bool::bool_(), sort_nat::nat())));

  BOOST_CHECK(!specification.mappings(make_function_sort_(sort_bool::bool_(), sort_bool::bool_(), sort_nat::nat())).empty());

  // Manually structured sort
  std::vector< data::structured_sort_constructor_argument > arguments;

  // sort that references itself by a name
  arguments.emplace_back(basic_sort("Q"));

  std::vector< data::structured_sort_constructor > constructors;
  constructors.emplace_back("q", structured_sort_constructor_argument_list(arguments.begin(), arguments.begin() + 1));

  specification.add_alias(alias(basic_sort("Q"), data::structured_sort(constructors)));
}

BOOST_AUTO_TEST_CASE(test_utility_functionality)
{
  data_specification spec;

  basic_sort s("S");
  basic_sort s0("S0");
  basic_sort a("a");
  function_sort s0s(sort_expression_list({s0}), s);
  data::function_symbol f("f", s);

  data::function_symbol g("g", s0s);
  data::function_symbol h("h", s0);

  {
    const std::set<sort_expression>& sorts(spec.sorts());
    BOOST_CHECK(std::find(sorts.begin(), sorts.end(), s0) == sorts.end());
    BOOST_CHECK(std::find(sorts.begin(), sorts.end(), s) == sorts.end());
    function_symbol_vector constructors(spec.constructors());
    BOOST_CHECK(std::find(constructors.begin(), constructors.end(), f) == constructors.end());
    function_symbol_vector mappings(spec.mappings());
    BOOST_CHECK(std::find(mappings.begin(), mappings.end(), f) == mappings.end());
    BOOST_CHECK(std::find(mappings.begin(), mappings.end(), g) == mappings.end());
  }

  spec.add_sort(s0);
  spec.add_constructor(f);
  spec.add_mapping(g);

  {
    const std::set<sort_expression>& sorts(spec.sorts());
    BOOST_CHECK(std::find(sorts.begin(), sorts.end(), s0) != sorts.end());
    BOOST_CHECK(std::find(sorts.begin(), sorts.end(), s) != sorts.end()); // Automatically added!
    function_symbol_vector constructors(spec.constructors());
    BOOST_CHECK(std::find(constructors.begin(), constructors.end(), f) != constructors.end());
    function_symbol_vector mappings(spec.mappings());
    BOOST_CHECK(std::find(mappings.begin(), mappings.end(), f) == mappings.end());
    BOOST_CHECK(std::find(mappings.begin(), mappings.end(), g) != mappings.end());
    BOOST_CHECK(std::find(mappings.begin(), mappings.end(), h) == mappings.end());
  }
  spec.add_mapping(h);

  spec.add_sort(s);
  spec.add_alias(alias(basic_sort("a"),s));

  const std::set<sort_expression> sorts(spec.sorts());
  function_symbol_vector constructors(spec.constructors());
  function_symbol_vector mappings(spec.mappings());


  BOOST_CHECK(std::find(sorts.begin(), sorts.end(), s0) != sorts.end());
  BOOST_CHECK(sorts.find(normalize_sorts(a, spec)) != sorts.end());
  BOOST_CHECK(std::find(sorts.begin(), sorts.end(), s) != sorts.end());

  BOOST_CHECK(std::find(mappings.begin(), mappings.end(), f) == mappings.end());
  BOOST_CHECK(std::find(constructors.begin(), constructors.end(), f) != constructors.end());

  BOOST_CHECK(std::find(constructors.begin(), constructors.end(), g) == constructors.end());
  BOOST_CHECK(std::find(mappings.begin(), mappings.end(), g) != mappings.end());
  BOOST_CHECK(std::find(constructors.begin(), constructors.end(), g) == constructors.end());
  BOOST_CHECK(std::find(constructors.begin(), constructors.end(), h) == constructors.end());
  BOOST_CHECK(std::find(mappings.begin(), mappings.end(), h) != mappings.end());
  BOOST_CHECK(spec.constructors(a) == spec.constructors(s));
  BOOST_CHECK(spec.mappings(a) == spec.mappings(s));
}

BOOST_AUTO_TEST_CASE(test_normalisation)
{
  using namespace mcrl2::data;
  using namespace mcrl2::data::sort_list;
  using namespace mcrl2::data::sort_set;
  using namespace mcrl2::data::sort_bag;

  data_specification specification;

  basic_sort A("A");

  specification.add_sort(A);

  basic_sort L("L");
  basic_sort S("S");
  basic_sort B("B");

  specification.add_alias(alias(L, list(A)));
  specification.add_alias(alias(S, set_(A)));
  specification.add_alias(alias(B, bag(A)));
  BOOST_CHECK(normalize_sorts(L,specification) == normalize_sorts(list(A),specification));
  BOOST_CHECK(normalize_sorts(list(L),specification) == normalize_sorts(list(list(A)),specification));
  BOOST_CHECK(normalize_sorts(S,specification) == normalize_sorts(set_(A),specification));
  BOOST_CHECK(normalize_sorts(list(S),specification) == normalize_sorts(list(set_(A)),specification));
  BOOST_CHECK(normalize_sorts(B,specification) == normalize_sorts(bag(A),specification));
  BOOST_CHECK(normalize_sorts(list(B),specification) == normalize_sorts(list(bag(A)),specification));

  specification = parse_data_specification(
                    "sort A = struct a(B);"
                    "sort B = struct b(A)|c;");

  std::vector< structured_sort_constructor_argument > arguments;

  arguments.emplace_back(basic_sort("B"));
  arguments.emplace_back(basic_sort("A"));

  std::vector< structured_sort_constructor > constructors;

  constructors.emplace_back("a", structured_sort_constructor_argument_list(arguments.begin(), arguments.begin() + 1));
  constructors.emplace_back("b", structured_sort_constructor_argument_list(arguments.begin() + 1, arguments.end()));
  constructors.emplace_back("c");

  structured_sort sA(data::structured_sort(structured_sort_constructor_list(constructors.begin(), constructors.begin() + 1)));
  structured_sort sB(data::structured_sort(structured_sort_constructor_list(constructors.begin() + 1, constructors.end())));

  const std::set<sort_expression> sorts(specification.sorts());
  BOOST_CHECK(sorts.find(normalize_sorts(sA, specification)) != sorts.end());
  BOOST_CHECK(sorts.find(normalize_sorts(sB, specification)) != sorts.end());

  BOOST_CHECK(normalize_sorts(sA,specification) == normalize_sorts(normalize_sorts(sA,specification),specification));
  BOOST_CHECK(normalize_sorts(sB,specification) == normalize_sorts(normalize_sorts(sB,specification),specification));

  // Check whether the sort expression struct f(struct f(A)|g) |g normalises to A
  // in the context of the specification sort A = struct f(A) | g;
  specification = parse_data_specification(
                    "sort A = struct f(A) | g;");

  std::vector< structured_sort_constructor_argument > arguments1;
  arguments1.emplace_back(basic_sort("A"));

  std::vector< structured_sort_constructor > constructors1;
  constructors1.emplace_back("f", arguments1);
  constructors1.emplace_back("g");

  sort_expression s1=structured_sort(constructors1); // s1 has the shape struct f(A)|g
  std::vector< structured_sort_constructor_argument > arguments2;
  arguments2.emplace_back(s1);

  std::vector< structured_sort_constructor > constructors2;
  constructors2.emplace_back("f", arguments2);
  constructors2.emplace_back("g");

  sort_expression s2=structured_sort(constructors2); // s2 has the shape f(struct f(A)|g) |g
  BOOST_CHECK(normalize_sorts(s2,specification)==basic_sort("A"));
}

BOOST_AUTO_TEST_CASE(test_copy)
{
  std::clog << "test_copy" << std::endl;

  data_specification specification = parse_data_specification(
                                       "sort D = struct d(bla : Bool)?is_d;"
                                       "sort A = S;"
                                       "sort S;"
                                       "map f: Set(S);");

  function_symbol_vector constructors(specification.constructors());
  BOOST_CHECK(std::find(constructors.begin(), constructors.end(), sort_bool::true_()) != constructors.end());

  data_specification other;
  other = specification;

  BOOST_CHECK(other == specification);

  specification = data_specification();


  BOOST_CHECK(normalize_sorts(basic_sort("A"),other) == normalize_sorts(basic_sort("S"),other));

  const std::set<sort_expression> sorts(specification.sorts());
  BOOST_CHECK(std::find(sorts.begin(), sorts.end(), basic_sort("A")) == sorts.end());
}

BOOST_AUTO_TEST_CASE(test_specification)
{
  data_specification spec = parse_data_specification("sort D = struct d1|d2;");
  BOOST_CHECK(spec.constructors(basic_sort("D")).size() == 2);
}

BOOST_AUTO_TEST_CASE(test_bke)
{
  std::cout << "test_bke" << std::endl;

  std::string BKE =
    "% This model is translated from the mCRL model used for analysing the Bilateral                \n"
    "% Key Exchange (BKE) protocol. The analysis is described in a paper with the                   \n"
    "% name 'Analysing the BKE-security protocol with muCRL', by Jan Friso Groote,                  \n"
    "% Sjouke Mauw and Alexander Serebrenik.                                                        \n"
    "%                                                                                              \n"
    "% The translation of the existing mCRL model into this mCRL2 model has been                    \n"
    "% performed manually.  The purpose was making use of the additional language                   \n"
    "% features of mCRL2 with respect to mCRL.                                                      \n"
    "%                                                                                              \n"
    "% The behaviour of this model should be bisimular with the original mcrl model                 \n"
    "% after renaming actions. Though this is not verified for all system                           \n"
    "% configurations below.                                                                        \n"
    "%                                                                                              \n"
    "% Eindhoven, June 11, 2008, Jeroen van der Wulp                                                \n"
    "                                                                                               \n"
    "% Agents. There are exactly three agents - A, B, E An order E < A < B is                       \n"
    "% imposed on agents to reduce the size of the state space.                                     \n"
    "sort Agent = struct A | B | E;                                                                 \n"
    "map  less : Agent # Agent -> Bool;                                                             \n"
    "var  a: Agent;                                                                                 \n"
    "eqn  less(A,a) = (a == B);                                                                     \n"
    "     less(B,a) = false;                                                                        \n"
    "     less(E,a) = (a != E);                                                                     \n"
    "                                                                                               \n"
    "sort Address = struct address(agent : Agent);                                                  \n"
    "map  bad_address : Address;                                                                    \n"
    "                                                                                               \n"
    "% A nonce is a random, unpredictable value which is used to make the                           \n"
    "% exchanged messages unique and thus helps to counter replay attacks.                          \n"
    "sort Nonce = struct nonce(value : Nat);                                                        \n"
    "                                                                                               \n"
    "% There are two kinds of keys used in the protocol: symmetric and                              \n"
    "% asymmetric ones (functional keys).                                                           \n"
    "% Symmetric keys have form K(n) where n is a natural number.                                   \n"
    "sort SymmetricKey = struct symmetric_key(value : Nat);                                         \n"
    "                                                                                               \n"
    "% Sort for representing asymmetric keys                                                        \n"
    "sort AsymmetricKey = struct public_key(Agent)?is_public |                                      \n"
    "                            secret_key(Agent)?is_secret |                                      \n"
    "                            hash(value : Nonce)?is_hash;                                       \n"
    "map  has_complementary_key: AsymmetricKey -> Bool;                                             \n"
    "     complementary_key    : AsymmetricKey -> AsymmetricKey;                                    \n"
    "var  a : Agent;                                                                                \n"
    "     n : Nonce;                                                                                \n"
    "eqn  has_complementary_key(public_key(a)) = true;                                              \n"
    "     has_complementary_key(secret_key(a)) = true;                                              \n"
    "     has_complementary_key(hash(n)) = false;                                                   \n"
    "     complementary_key(public_key(a)) = secret_key(a);                                         \n"
    "     complementary_key(secret_key(a)) = public_key(a);                                         \n"
    "                                                                                               \n"
    "sort Key = struct key(SymmetricKey)?is_symmetric | key(AsymmetricKey)?is_asymmetric;           \n"
    "map  has_complementary_key: Key -> Bool;                                                       \n"
    "     complementary_key : Key -> Key; % gets the complementary key if key is asymmetric         \n"
    "var  a,a1   : Agent;                                                                           \n"
    "     n,n1   : Nat;                                                                             \n"
    "     k,k1   : Key;                                                                             \n"
    "     ak,ak1 : AsymmetricKey;                                                                   \n"
    "eqn                                                                                            \n"
    "     % gets the complementary key if key is asymmetric                                         \n"
    "     complementary_key(key(ak)) = key(complementary_key(ak));                                  \n"
    "     has_complementary_key(key(ak)) = has_complementary_key(ak);                               \n"
    "                                                                                               \n"
    "sort Message = struct                                                                          \n"
    "        encrypt(Nonce, Address, AsymmetricKey)?is_message_1 |                                  \n"
    "        encrypt(AsymmetricKey, Nonce, SymmetricKey, AsymmetricKey)?is_message_2 |              \n"
    "        encrypt(AsymmetricKey, SymmetricKey)?is_message_3;                                     \n"
    "map  valid_message_1 : Message # AsymmetricKey -> Bool;                                        \n"
    "     valid_message_2 : Message # AsymmetricKey -> Bool;                                        \n"
    "     valid_message_3 : Message # SymmetricKey  -> Bool;                                        \n"
    "     used_key : Message -> Key;                  % key used to encrypt                         \n"
    "var  sk, sk1 : SymmetricKey;                                                                   \n"
    "     ak, ak1, ak2 : AsymmetricKey;                                                             \n"
    "     n, n1  : Nonce;                                                                           \n"
    "     m, m1   : Message;                                                                        \n"
    "     a, a1   : Address;                                                                        \n"
    "eqn  used_key(encrypt(n,a,ak)) = key(ak);                                                      \n"
    "     used_key(encrypt(ak,n1,sk,ak1)) = key(ak1);                                               \n"
    "     used_key(encrypt(ak,sk)) = key(sk);                                                       \n"
    "     valid_message_1(m, ak) = is_message_1(m) && (used_key(m) == key(ak));                     \n"
    "     valid_message_2(m, ak) = is_message_2(m) && (used_key(m) == key(ak));                     \n"
    "     valid_message_3(m, sk) = is_message_3(m) && (used_key(m) == key(sk));                     \n"
    "                                                                                               \n"
    "% Type for message sets; currently cannot use Set() because set iteration is not possible      \n"
    "sort MessageSet = List(Message);                                                               \n"
    "map  insert                : Message # MessageSet -> MessageSet;                               \n"
    "     select_crypted_by     : Key # MessageSet -> MessageSet;                                   \n"
    "     select_not_crypted_by : Key # MessageSet -> MessageSet;                                   \n"
    "     select                : (Message -> Bool) # MessageSet -> MessageSet;                     \n"
    "var  k,k1 : Key;                                                                               \n"
    "     m,m1 : Message;                                                                           \n"
    "     ms   : MessageSet;                                                                        \n"
    "     c    : Message -> Bool;                                                                   \n"
    "eqn                                                                                            \n"
    "     % inserts a message m, if it is not in the list                                           \n"
    "     insert(m,[]) = [m];                                                                       \n"
    "     m < m1  -> insert(m,m1|>ms) = m|>m1|>ms;                                                  \n"
    "     m == m1 -> insert(m,m1|>ms) = m1|>ms;                                                     \n"
    "     m1 < m  -> insert(m,m1|>ms) = m1|>insert(m,ms);                                           \n"
    "                                                                                               \n"
    "     % the set (as ordered list) of messages in m that are signed by sk                        \n"
    "     select_crypted_by(k,ms) = select(lambda x : Message.k == used_key(x),ms);                 \n"
    "     select_not_crypted_by(k,ms) = select(lambda x : Message.k != used_key(x),ms);             \n"
    "                                                                                               \n"
    "     select(c,[]) = [];                                                                        \n"
    "     select(c,m|>ms) = if(c(m),m|>r,r) whr r = select(c, ms) end;                              \n"
    "                                                                                               \n"
    "% The eavesdropper's knowledge consists of:                                                    \n"
    "%  * a list of addresses                                                                       \n"
    "%  * a list of nonces                                                                          \n"
    "%  * a list of keys (both symmetric and asymmetric)                                            \n"
    "%  * a list of messages of which the key is not known                                          \n"
    "sort Knowledge = struct                                                                        \n"
    "        knowledge(addresses : Set(Address),                                                    \n"
    "                  nonces : Set(Nonce),                                                         \n"
    "                  keys : Set(Key),                                                             \n"
    "                  messages : MessageSet);                                                      \n"
    "map  update_knowledge : Message # Knowledge -> Knowledge;                                      \n"
    "     propagate : MessageSet # Knowledge -> Knowledge;                                          \n"
    "     propagate : Key # Knowledge -> Knowledge;                                                 \n"
    "     add_key : Key # Knowledge -> Knowledge;                                                   \n"
    "     add_nonce : Nonce # Knowledge -> Knowledge;                                               \n"
    "     add_address : Address # Knowledge -> Knowledge;                                           \n"
    "var  m  : Message;                                                                             \n"
    "     as : Set(Address);                                                                        \n"
    "     ns : Set(Nonce);                                                                          \n"
    "     ks : Set(Key);                                                                            \n"
    "     ms : MessageSet;                                                                          \n"
    "     k  : Knowledge;                                                                           \n"
    "     sk : SymmetricKey;                                                                        \n"
    "     ak,hk : AsymmetricKey;                                                                    \n"
    "     ck : Key;                                                                                 \n"
    "     n, n1 : Nonce;                                                                            \n"
    "     a  : Address;                                                                             \n"
    "eqn                                                                                            \n"
    "     % adds keys to knowledge that are part of known messages encrypted with a new key         \n"
    "     has_complementary_key(ak) && complementary_key(key(ak)) in keys(k) ->                     \n"
    "         update_knowledge(encrypt(n,a,ak),k) =                                                 \n"
    "                propagate(key(ak), add_key(key(ak),                                            \n"
    "                        add_address(a, add_nonce(n, k))));                                     \n"
    "     has_complementary_key(ak) && complementary_key(key(ak)) in keys(k) ->                     \n"
    "         update_knowledge(encrypt(hk,n1,sk,ak),k) =                                            \n"
    "                propagate(key(sk), propagate(key(ak),                                          \n"
    "                      add_key(key(sk),add_key(key(ak),k))));                                   \n"
    "     key(sk) in keys(k) ->                                                                     \n"
    "         update_knowledge(encrypt(ak,sk),k) =                                                  \n"
    "                propagate(key(ak), add_key(key(ak),k));                                        \n"
    "                                                                                               \n"
    "     % adds a message that cannot be decrypted with any known key                              \n"
    "     ((is_symmetric(ck) && !(ck in keys(k))) ||                                                \n"
    "       (is_asymmetric(ck) && !(has_complementary_key(ck) &&                                    \n"
    "                (complementary_key(ck) in keys(k))))) whr ck = used_key(m) end ->              \n"
    "         update_knowledge(m,k) =                                                               \n"
    "                knowledge(addresses(k),nonces(k),keys(k),insert(m,messages(k)));               \n"
    "                                                                                               \n"
    "     % adds a key to knowledge                                                                 \n"
    "     add_key(ck,knowledge(as,ns,ks,ms)) = knowledge(as,ns,ks + {ck},ms);                       \n"
    "     add_nonce(n,knowledge(as,ns,ks,ms)) = knowledge(as,ns + {n},ks,ms);                       \n"
    "     add_address(a,knowledge(as,ns,ks,ms)) = knowledge(as + {a},ns,ks,ms);                     \n"
    "                                                                                               \n"
    "     % adds keys to knowledge that are part of messages encrypted with a key k                 \n"
    "     propagate([],k) = k;                                                                      \n"
    "     propagate(m|>ms,k) = propagate(ms, update_knowledge(m, k));                               \n"
    "     propagate(ck,knowledge(as,ns,ks,ms)) =                                                    \n"
    "        propagate(select_crypted_by(ck,ms),                                                    \n"
    "                knowledge(as,ns,ks,select_not_crypted_by(ck,ms)));                             \n"
    ;

  data_specification data_spec = parse_data_specification(BKE);
  const alias_vector& aliases = data_spec.user_defined_aliases();
  for (const alias& a: aliases)
  {
    std::cout << "alias " << a << std::endl;
    const sort_expression& s = a.reference();
    if (is_structured_sort(s))
    {
      structured_sort sort(s);
      for (const structured_sort_constructor& constructor: sort.constructors())
      {
        for (const structured_sort_constructor_argument& argument: constructor.arguments())
        {
          std::cout << "argument: " << argument << " " << argument << std::endl;
          const atermpp::aterm& name = argument.name();
          if (name != core::empty_identifier_string())
          {
            std::cout << "name = " << name << std::endl;
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_abuse_of_tail)
{
  std::cout << "Test abuse of tail\n";
  const std::string spec_string =
    "map tail:Nat#List(Bool) -> List(Bool);\n"
    "var vs: List(Bool);\n"
    "n: Nat;\n"
    "eqn tail(n, vs) = if(n==0, vs, tail(Int2Nat(n-1),vs));\n";
  try
  {
    data_specification data_spec=parse_data_specification(spec_string);
    BOOST_CHECK(false); // Typechecking is supposed to fail; one cannot get here.
  }
  catch (mcrl2::runtime_error&) // NOLINT(bugprone-empty-catch)
  {
    // It is expected that a runtime error is thrown.
  }
}

BOOST_AUTO_TEST_CASE(test_merge_data_specifications)
{
  std::string DATASPEC =
    "sort D;\n"
    "cons s: D;"
  ;
  data_specification dataspec1 = parse_data_specification(DATASPEC);
  data_specification dataspec2 = parse_data_specification(DATASPEC);
  data_specification dataspec3 = merge_data_specifications(dataspec1, dataspec2);
  BOOST_CHECK(dataspec1 == dataspec3);
}

BOOST_AUTO_TEST_CASE(test_standard_sorts_mappings_functions)
{
   std::set < sort_expression > sorts;
   std::set < function_symbol > constructors;
   std::set <function_symbol > mappings;

   data_specification spec;
   spec.get_system_defined_sorts_constructors_and_mappings(sorts, constructors, mappings);

std::cerr << "#sorts " << sorts.size() << "\n";
std::cerr << "#constructorss " << constructors.size() << "\n";
std::cerr << "#mappings " << mappings.size() << "\n";

#ifdef MCRL2_ENABLE_MACHINENUMBERS
   BOOST_CHECK(sorts.size()==11);
   BOOST_CHECK(constructors.size()==19);
   BOOST_CHECK(mappings.size()==333);
#else
   BOOST_CHECK(sorts.size()==10);
   BOOST_CHECK(constructors.size()==17);
   BOOST_CHECK(mappings.size()==243);
#endif
}


// The next check tests whether elementary sorts in function symbols occur in data_spec.sorts().
// The sort natnatpair occurring in Nat64 was forgotten to be added.
BOOST_AUTO_TEST_CASE(test_constructors_and_sorts)
{
  std::cerr << "Test: test_constructors_and_sorts\n";
  data::data_specification data_spec = data::parse_data_specification(
    "sort Bit = struct b0(r:Real) | b1;\n"
    "map  invert: Bit -> Bit;\n"
    "var  r1: Real;\n"
    "eqn  invert(b1)= b0(1);\n"
    "     invert(b0(r1))= b1;\n");

  std::set < sort_expression > sorts;
  std::set < function_symbol > constructors;
  std::set <function_symbol > mappings;

  data_spec.get_system_defined_sorts_constructors_and_mappings(sorts, constructors, mappings);

  for(const function_symbol& f: data_spec.constructors())
  {
    for(const sort_expression& sort: find_sort_expressions(f))
    {
      BOOST_CHECK(is_function_sort(sort)||data_spec.sorts().count(sort)>0);
    }
  }
}
