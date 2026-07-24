// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/absinthe.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_ABSINTHE_H
#define MCRL2_PBES_ABSINTHE_H

#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/data/data_expression.h"
#include "mcrl2/data/consistency.h"
#include "mcrl2/data/detail/data_construction.h"
#include "mcrl2/pbes/builder.h"
#include "mcrl2/utilities/detail/separate_keyword_section.h"
#include "mcrl2/data/detail/print_parse_check.h"

namespace mcrl2::pbes_system
{

template <typename Term>
std::string print_term(const Term& x)
{
  return data::pp(x) + " " + data::pp(x);
  }

  template <typename Term>
  std::string print_symbol(const Term& x)
  {
    return data::pp(x) + ": " + data::pp(x.sort());
  }

namespace detail {

  inline
  data::data_specification& absinthe_data_specification()
  {
    static data::data_specification dataspec;
    return dataspec;
  }

  template <typename T>
  inline void absinthe_check_expression(const T& x)
  {
    data::detail::print_parse_check(x, absinthe_data_specification());
  }

  // Returns true if f appears as a structured sort constructor in dataspec.
  inline
  bool is_structured_sort_constructor(const data::data_specification& dataspec, const data::function_symbol& f)
  {
    for (const data::alias& a: dataspec.user_defined_aliases())
    {
      if (f.sort() != a.name())
      {
        continue;
      }
      const data::sort_expression& s = a.reference();
      if (data::is_structured_sort(s))
      {
        const auto& ss = atermpp::down_cast<data::structured_sort>(s);
        for (const data::function_symbol& g: ss.constructor_functions())
        {
          if (f.name() == g.name())
          {
            return true;
          }
        }
      }
    }
    return false;
  }

  inline
  void print_used_function_symbols(const pbes& p)
  {
    mCRL2log(log::debug) << "--- used function symbols ---" << std::endl;
    for (const data::function_symbol& f: pbes_system::find_function_symbols(p))
    {
      mCRL2log(log::debug) << print_symbol(f) << std::endl;
    }
  }

  // TODO: Is this correct if s has the shape A -> (B -> C)? Should the result be (B -> C) or C?
  inline
  data::sort_expression target_sort(const data::sort_expression& s)
  {
    if (data::is_basic_sort(s))
    {
      return s;
    }
    else if (data::is_function_sort(s))
    {
      const auto& fs = atermpp::down_cast<data::function_sort>(s);
      return fs.codomain();
    }
    else if (data::is_container_sort(s))
    {
      return s;
    }
    throw mcrl2::runtime_error("target_sort: unsupported sort " + print_term(s) + " detected!");
    return data::sort_expression();
  }

} // namespace detail

struct absinthe_algorithm
{
  using sort_expression_substitution_map = std::map<data::sort_expression, data::sort_expression>;
  using function_symbol_substitution_map = std::map<data::function_symbol, data::function_symbol>;
  using abstraction_map = std::map<data::sort_expression, data::function_symbol>;

  // Used for generating variables of sort comprehensions.
  data::set_identifier_generator m_generator;

  struct absinthe_sort_expression_builder: public sort_expression_builder<absinthe_sort_expression_builder>
  {
    using super = sort_expression_builder<absinthe_sort_expression_builder>;
    using super::apply;

    const abstraction_map& sigmaH;
    const sort_expression_substitution_map& sigmaS;
    const function_symbol_substitution_map& sigmaF;
    data::set_identifier_generator& generator;

    absinthe_sort_expression_builder(const abstraction_map& sigmaA_,
                                     const sort_expression_substitution_map& sigmaS_,
                                     const function_symbol_substitution_map& sigmaF_,
                                     data::set_identifier_generator& generator_
                                    )
      : sigmaH(sigmaA_),
        sigmaS(sigmaS_),
        sigmaF(sigmaF_),
        generator(generator_)
    {}

    template <class T>
    void apply(T& result, const data::sort_expression& x)
    {
      auto i = sigmaS.find(x);
      if (i == sigmaS.end())
      {
        super::apply(result, x);
      }
      else
      {
        result = i->second;
      }
    }

    template <class T>
    void apply(T& result, const data::function_symbol& x)
    {
      auto i = sigmaF.find(x);
      if (i != sigmaF.end())
      {
        result = i->second;
        return;
      }
      throw mcrl2::runtime_error("function symbol " + print_symbol(x) + " not present in the function symbol mapping!");
    }

    template <class T>
    void apply(T& result, const data::application& x)
    {
      if (data::is_variable(x.head()))
      {
        data::variable v = atermpp::down_cast<data::variable>(x.head());
        data::sort_expression sort;
        super::apply(sort, v.sort());
        v = data::variable(v.name(), sort);
        result = data::detail::create_finite_set(data::application(v, x.begin(), x.end()));
        return;
      }
      else if (data::is_function_symbol(x.head()))
      {
        super::apply(result, x);
        return;
      }
      else
      {
        throw mcrl2::runtime_error("don't know how to handle arbitrary expression as head: " + data::pp(x));
      }
    }

    template <class T>
    void apply(T& result, const data::lambda& x)
    {
      data::data_expression body;
      super::apply(body, x);
      data::sort_expression s = body.sort();
      generator.add_identifiers(data::find_identifiers(x));
      data::variable v(generator("v"), s);
      result = data::detail::create_set_comprehension(v, data::equal_to(v, body));
    }

    template <class T>
    void apply(T& result, const data::data_expression& x)
    {
      if (data::is_variable(x))
      {
        super::apply(result, x);
        result = data::detail::create_finite_set(result);
      }
      else
      {
        // check if it is a "ground term", i.e. it does not contain any variables
        auto i = sigmaH.find(x.sort());
        if (i != sigmaH.end() && data::find_all_variables(x).empty())
        {
          data::data_expression_list args = { x };
          result = data::detail::create_finite_set(data::application(i->second, args));
        }
        else
        {
          // first apply the sort and function symbol transformations
          super::apply(result, x);
        }
      }
    }
  };

  // Applies sigmaS to a sort expression
  struct sort_function
  {
    using argument_type = data::sort_expression;
    using result_type = data::sort_expression;

    absinthe_sort_expression_builder f;

    sort_function(const abstraction_map& sigmaH,
                  const sort_expression_substitution_map& sigmaS,
                  const function_symbol_substitution_map& sigmaF,
                  data::set_identifier_generator& generator
                 )
      : f(sigmaH, sigmaS, sigmaF, generator)
    {}

    data::sort_expression operator()(const data::sort_expression& x)
    {
      data::sort_expression sort;
      f.apply(sort, x);
      return sort;
    }
  };

  struct absinthe_data_expression_builder: public pbes_expression_builder<absinthe_data_expression_builder>
  {
    using super = pbes_expression_builder<absinthe_data_expression_builder>;
    using super::apply;
    using super::update;

    data::variable_list make_variables(const data::data_expression_list& x, const std::string& hint, sort_function sigma) const
    {
      std::vector<data::variable> result;
      std::size_t i = 0;
      for (auto j = x.begin(); j != x.end(); ++i, ++j)
      {
        result.emplace_back(hint + utilities::number2string(i), sigma(j->sort()));
      }
      return data::variable_list(result.begin(), result.end());
    }

    const abstraction_map& sigmaH;
    const sort_expression_substitution_map& sigmaS;
    const function_symbol_substitution_map& sigmaF;
    data::set_identifier_generator& generator;
    bool m_is_over_approximation;

    data::data_expression lift(const data::data_expression& x)
    {
      data::data_expression result;
      absinthe_sort_expression_builder(sigmaH, sigmaS, sigmaF, generator).apply(result, x);
      return result;
    }

    data::data_expression_list lift(const data::data_expression_list& x)
    {
      data::data_expression_list result;
      absinthe_sort_expression_builder(sigmaH, sigmaS, sigmaF, generator).apply(result, x);
      return result;
    }

    data::variable_list lift(const data::variable_list& x)
    {
      data::variable_list result;
      absinthe_sort_expression_builder(sigmaH, sigmaS, sigmaF, generator).apply(result, x);
      return result;
    }

    pbes_system::propositional_variable lift(const pbes_system::propositional_variable& x)
    {
      pbes_system::propositional_variable result;
      absinthe_sort_expression_builder(sigmaH, sigmaS, sigmaF, generator).apply(result, x);
      return result;
    }

    absinthe_data_expression_builder(const abstraction_map& sigmaA_,
                                     const sort_expression_substitution_map& sigmaS_,
                                     const function_symbol_substitution_map& sigmaF_,
                                     data::set_identifier_generator& generator_,
                                     bool is_over_approximation)
      : sigmaH(sigmaA_),
        sigmaS(sigmaS_),
        sigmaF(sigmaF_),
        generator(generator_),
        m_is_over_approximation(is_over_approximation)
    {}

    template <class T>
    void apply(T& result, const data::data_expression& x)
    {
      data::data_expression x1 = lift(x);
      if (m_is_over_approximation)
      {
        result = atermpp::down_cast<T>(data::detail::create_set_in(data::true_(), x1));
      }
      else
      {
        result =  atermpp::down_cast<T>(data::not_(data::detail::create_set_in(data::false_(), x1)));
      }
    }

    template <class T>
    void apply(T& result, const propositional_variable_instantiation& x)
    {
      data::data_expression_list e = lift(x.parameters());
      data::variable_list variables = make_variables(x.parameters(), "x", sort_function(sigmaH, sigmaS, sigmaF, generator));
      data::data_expression_list::iterator i = e.begin();
      data::variable_list::iterator j = variables.begin();
      data::data_expression_vector z;
      for (; i != e.end(); ++i, ++j)
      {
        z.push_back(data::detail::create_set_in(*j, *i));
      }
      data::data_expression q = data::lazy::join_and(z.begin(), z.end());
      if (m_is_over_approximation)
      {
        result = optimized_exists(variables, and_(atermpp::down_cast<pbes_expression>(q), 
          propositional_variable_instantiation(x.name(), data::data_expression_list(variables))));
      }
      else
      {
        result = optimized_forall(variables, imp(atermpp::down_cast<pbes_expression>(q), propositional_variable_instantiation(x.name(), data::data_expression_list(variables))));
      }
    }

    template <class T>
    void apply(T& result, const pbes_system::forall& x)
    {
      pbes_system::pbes_expression body;
      super::apply(body, x.body());
      result = optimized_forall(lift(x.variables()), body);
    }

    template <class T>
    void apply(T& result, const pbes_system::exists& x)
    {
      pbes_system::pbes_expression body;
      super::apply(body, x.body());
      result = optimized_exists(lift(x.variables()), body);
    }

    void update(pbes_system::pbes_equation& x)
    {
      x.variable() = lift(x.variable());
      pbes_system::pbes_expression result;
      super::apply(result, x.formula());
      x.formula() = result;
    }

    void update(pbes_system::pbes& x)
    {
      super::update(x.equations());
      pbes_expression kappa;
      apply(kappa, x.initial_state());
      core::identifier_string name("GeneratedZ");
      propositional_variable Z(name, data::variable_list());
      x.equations().emplace_back(fixpoint_symbol::mu(), Z, kappa);
      x.initial_state() = propositional_variable_instantiation(name, data::data_expression_list());
    }
  };

  sort_expression_substitution_map parse_sort_expression_mapping(const std::string& text, const data::data_specification& dataspec)
  {
    sort_expression_substitution_map result;

    for (const std::string& line: utilities::regex_split(text, "\\n"))
    {
      std::vector<std::string> words = utilities::regex_split(line, ":=");
      if (words.size() == 2)
      {
        data::sort_expression lhs = data::parse_sort_expression(words[0], dataspec);
        data::sort_expression rhs = data::parse_sort_expression(words[1], dataspec);
        result[lhs] = rhs;
      }
    }
    return result;
  }

  // Parse the right hand sides of the function symbol mapping in text, and add them to dataspec if needed
  void parse_right_hand_sides(const std::string& text, data::data_specification& dataspec)
  {
    std::string dataspec_text = data::pp(dataspec);
    for (const std::string& line: utilities::regex_split(text, "\\n"))
    {
      std::vector<std::string> words = utilities::regex_split(line, ":=");
      if (words.size() == 2)
      {
        data::function_symbol f = data::parse_function_symbol(words[1], dataspec_text);
        if (!pbes_system::detail::is_structured_sort_constructor(dataspec, f))
        {
          dataspec.add_mapping(f);
        }
      }
    }
  }

  function_symbol_substitution_map parse_function_symbol_mapping(const std::string& text, const data::data_specification& dataspec)
  {
    function_symbol_substitution_map result;
    std::string dataspec_text = data::pp(dataspec);

    for (const std::string& line: utilities::regex_split(text, "\\n"))
    {
      std::vector<std::string> words = utilities::regex_split(line, ":=");
      if (words.size() == 2)
      {
        data::function_symbol lhs = data::parse_function_symbol(words[0], dataspec_text);
        std::string s = words[1];
        s = utilities::regex_replace(";\\s*$", "", s);
        data::function_symbol rhs = data::parse_function_symbol(s, dataspec_text);
        result[lhs] = rhs;
      }
    }
    return result;
  }

  // text is a data_specification; extract the user defined mappings
  abstraction_map parse_abstraction_map(const std::string& text)
  {
    abstraction_map result;
    data::data_specification dataspec = data::parse_data_specification(text);
    for (const data::function_symbol& i: dataspec.user_defined_mappings())
    {
      const auto& f = atermpp::down_cast<data::function_sort>(i.sort());
      if (f.domain().size() != 1)
      {
        throw mcrl2::runtime_error("cannot abstract the function " + data::pp(i) + " since the arity of the domain is not equal to one!");
      }
      result[f.domain().front()] = i;
    }
    return result;
  }

  // creates a finite set containing one data expression
  struct make_data_expression_set
  {
    data::data_expression operator()(const data::data_expression& x) const
    {
      data::sort_expression s = x.sort();
      data::data_expression result = data::sort_fset::empty(s);
      result = data::sort_fset::cons_(s, x, result);
      return result;
    }
  };

  // transforms the sort expression s to Set(s)
  struct make_set
  {
    data::sort_expression operator()(const data::sort_expression& s) const
    {
      return data::sort_set::set_(s);
    }
  };

  // function that transforms a function symbol
  struct lift_function_symbol_1_2
  {
    std::map<std::string, std::string> unprintable;
    std::set<std::string> suffix_with_sort;

    lift_function_symbol_1_2()
    {
      unprintable["&&"] = "and";
      unprintable["||"] = "or";
      unprintable["!"] = "not";
      unprintable["#"] = "len";
      unprintable["."] = "element_at";
      unprintable["+"] = "plus";
      unprintable["-"] = "minus";
      unprintable[">"] = "greater";
      unprintable["<"] = "less";
      unprintable[">="] = "ge";
      unprintable["<="] = "le";
      unprintable["=="] = "eq";
      unprintable["!="] = "neq";
      unprintable["[]"] = "emptylist";
      unprintable["++"] = "concat";
      unprintable["<|"] = "snoc";
      unprintable["|>"] = "cons";
      unprintable["@cNat"] = "cNat";
      unprintable["@cDub"] = "cDub";
      unprintable["@c0"] = "c0";
      unprintable["@c1"] = "c1";
      unprintable["@func_update"] = "func_update";
      unprintable["@cInt"] = "cInt";
      unprintable["@cNeg"] = "cNeg";
      unprintable["@most_significant_digit"] = "most_significant_digit";
      unprintable["@succ_pos"] = "succpos";
      unprintable["@most_significant_digitNat"] = "most_significant_digitNat";
      unprintable["@concat_digit"] = "concat_digit";
      unprintable["@succ_nat"] = "succ_nat";
      

      suffix_with_sort.insert("[]");
      suffix_with_sort.insert("|>");
    }

    std::string print_cleaned(const data::sort_expression& s) const
    {
      std::string result = data::pp(s);
      result = utilities::regex_replace("\\(", "_", result);
      result = utilities::regex_replace("\\)", "_", result);
      result = utilities::regex_replace("#", "_", result);
      result = utilities::regex_replace("->", "_", result);
      result = utilities::remove_whitespace(result);
      return result;
    }


    data::function_symbol operator()(const data::function_symbol& f, sort_function sigma) const
    {
    	using utilities::detail::contains;

      std::string name = std::string(f.name());

      bool print_sort = contains(suffix_with_sort, std::string(f.name()));
      auto i = unprintable.find(name);
      if (i != unprintable.end())
      {
        name = i->second;
      }
      name = "Generated_" + name;
      if (print_sort)
      {
        name = name + print_cleaned(f.sort());
      }

      const data::sort_expression& s = f.sort();
      if (data::is_basic_sort(s))
      {
        return data::function_symbol(name, sigma(s));
      }
      else if (data::is_function_sort(s))
      {
        // Apply sigmaS recursively to s
        //   f:        tail:           List(Nat) -> List(Nat)
        //   result:   generated_tail: List(AbsNat) -> Set(List(AbsNat))
        data::function_sort fs = atermpp::down_cast<data::function_sort>(sigma(s));
        return data::function_symbol(name, data::function_sort(fs.domain(), make_set()(fs.codomain())));
      }
      else if (data::is_container_sort(s))
      {
        // Apply sigmaS recursively to s
        // Example: List(Nat) -> List(AbsNat)
        return data::function_symbol(name, sigma(s));
      }
      throw mcrl2::runtime_error("absinthe algorithm: unsupported sort " + print_term(s) + " detected!");
      return data::function_symbol();
    }
  };

  // function that lifts a function symbol
  struct lift_function_symbol_2_3
  {
    data::function_symbol operator()(const data::function_symbol& f) const
    {
      using namespace data;
      std::string name = "Lift" + utilities::trim_copy(std::string(f.name()));
      const sort_expression &s = f.sort();
      if (is_basic_sort(s))
      {
        return function_symbol(name, make_set()(s));
      }
      else if (is_function_sort(s))
      {
        const auto& fs = atermpp::down_cast<data::function_sort>(s);
        const sort_expression_list& sl = fs.domain();
        return function_symbol(name, function_sort(sort_expression_list(sl.begin(),sl.end(), make_set()), fs.codomain()));
      }
      else if (is_container_sort(s))
      {
        return data::function_symbol(name, make_set()(s));
      }
      throw mcrl2::runtime_error("absinthe algorithm (lift): unsupported sort " + print_term(s) + " detected!");
    }
  };

  // function that generates an equation from a function symbol and it's corresponding 'generated' version
  struct lift_equation_1_2
  {
    lift_equation_1_2() = default;

    std::vector<data::variable> make_variables(const data::sort_expression_list& sorts, const std::string& hint, sort_function sigma) const
    {
      std::vector<data::variable> result;
      std::size_t i = 0;
      for (auto j = sorts.begin(); j != sorts.end(); ++i, ++j)
      {
        result.emplace_back(hint + utilities::number2string(i), sigma(*j));
      }
      return result;
    }

    // sigmaH is used for checking consistency of the user input
    data::data_equation operator()(const data::function_symbol& f1, const data::function_symbol& f2, sort_function sigma, const abstraction_map& sigmaH) const
    {
      mCRL2log(log::debug) << "lift_equation_1_2 f1 = " << print_symbol(f1) << " f2 = " << print_symbol(f2) << std::endl;
      data::variable_list variables;
      const data::data_expression& condition = data::true_();
      data::data_expression lhs;
      data::data_expression rhs;

      data::sort_expression s1 = f1.sort();

      if (data::is_basic_sort(s1))
      {
        lhs = f2;
        data::function_symbol f1_sigma(f1.name(), sigma(f1.sort()));

        auto i = sigmaH.find(f1.sort());
        if (i == sigmaH.end())
        {
          rhs = f1_sigma;
        }
        else
        {
          data::function_symbol h = i->second;
          rhs = data::application(h, f1);
        }
      }
      else if (data::is_function_sort(s1))
      {
        data::function_sort fs1(f1.sort());
        data::function_sort fs2(f2.sort());

        // check validity
        if (data::is_function_update_application(f1))
        {
          // TODO: add check that the domain of the updated function does not contain abstraction sorts
        }
        else if (fs1.domain() != fs2.domain())
        {
          throw std::runtime_error("can not generalize functions with abstraction sorts in the domain: " + data::pp(f1) + ": " + data::pp(s1));
        }

        data::variable_vector x = make_variables(fs2.domain(), "x", sigma);
        variables = data::variable_list(x.begin(),x.end());
        lhs = data::application(f2, x.begin(), x.end());
        data::application f_x(f1, x.begin(), x.end());

        data::function_symbol f1_sigma(f1.name(), sigma(f1.sort()));

        auto i = sigmaH.find(detail::target_sort(f1.sort()));
        if (i == sigmaH.end())
        {
          data::application f1_sigma_x(f1_sigma, x.begin(), x.end());
          rhs = data::detail::create_finite_set(f_x);
        }
        else
        {
          data::function_symbol h = i->second;
          rhs = data::detail::create_finite_set(data::application(h, f_x));
        }
      }
      else if (data::is_container_sort(s1))
      {
        // Example:
        //   f1:   [] : List(Nat)
        //   f2:   generated_emptylist : List(AbsNat)
        //  eqn:   generated_emptylist = []                    met [] : List(AbsNat)
        // tail: List(AbsNat) -> List(Nat)
        lhs = f2;
        data::function_symbol f1_sigma(f1.name(), sigma(f1.sort()));
        rhs = f1_sigma;
      }
      else
      {
        throw mcrl2::runtime_error("absinthe algorithm (lift_equation_1_2): unsupported sort " + print_term(s1) + " detected!");
      }

      if (lhs.sort() != rhs.sort())
      {
        throw mcrl2::runtime_error("absinthe algorithm (lift_equation_1_2): lhs.sort() and rhs.sort are not equal: " + data::pp(lhs.sort()) + " <-> " + data::pp(rhs.sort()));
      }

      return data::data_equation(variables, condition, lhs, rhs);
    }
  };

  // function that generates an equation from a function symbol and it's corresponding lifted version
  struct lift_equation_2_3
  {
    lift_equation_2_3() = default;

    std::vector<data::variable> make_variables(const data::sort_expression_list& sorts, const std::string& hint, sort_function sigma) const
    {
      std::vector<data::variable> result;
      std::size_t i = 0;
      for (auto j = sorts.begin(); j != sorts.end(); ++i, ++j)
      {
        result.emplace_back(hint + utilities::number2string(i), sigma(*j));
      }
      return result;
    }

    // Let x = [x1:D1, ..., xn:Dn] and X = [X1:Set(D1), ..., Xn:Set(Dn)]. Returns the expression
    //
    // (x1 in X1 /\ ... /\ xn in Xn)
    data::data_expression enumerate_domain(const std::vector<data::variable>& x, const std::vector<data::variable>& X) const
    {
      std::vector<data::data_expression> a;
      auto i = x.begin();
      auto j = X.begin();
      for (; i != x.end(); ++i, ++j)
      {
        a.push_back(data::detail::create_set_in(*i, *j));
      }
      data::data_expression body = data::lazy::join_and(a.begin(), a.end());
      return body;
    }

    data::data_equation operator()(const data::function_symbol& f2, const data::function_symbol& f3, sort_function sigma) const
    {
      data::variable_list variables;
      const data::data_expression& condition = data::true_();
      data::data_expression lhs;
      data::data_expression rhs;

      const data::sort_expression& s2 = f2.sort();

      if (data::is_basic_sort(s2))
      {
        lhs = f3;
        rhs = data::detail::create_finite_set(f2);
      }
      else if (data::is_function_sort(s2))
      {
        data::function_sort fs2(f2.sort());
        data::function_sort fs3(f3.sort());

        if (!data::is_container_sort(fs2.codomain()))
        {
          throw mcrl2::runtime_error("The codomain " + data::pp(fs2.codomain()) + " of function " + data::pp(f2.name()) +  " should be a container sort!");
        }

        // TODO: generate these variables in a proper way
        std::vector<data::variable> x = make_variables(fs2.domain(), "x", sigma);
        std::vector<data::variable> X = make_variables(fs3.domain(), "X", sigma);

        variables = data::variable_list(X.begin(), X.end());
        lhs = data::application(f3, X.begin(), X.end());
        data::variable y("y", data::detail::get_set_sort(atermpp::down_cast<data::container_sort>(fs2.codomain())));
        data::data_expression Y = data::application(f2, x.begin(), x.end());
        data::data_expression body = data::and_(enumerate_domain(x, X), data::detail::create_set_in(y, Y));
        rhs = data::detail::create_set_comprehension(y, data::exists(x, body));
      }
      else if (data::is_container_sort(s2))
      {
        lhs = f3;
        rhs = data::detail::create_finite_set(f2);
      }
      else
      {
        throw mcrl2::runtime_error("absinthe algorithm (lift_equation_2_3): unsupported sort " + print_term(s2) + " detected!");
      }
      return data::data_equation(variables, condition, lhs, rhs);
    }
  };

  template <typename Map>
  std::string print_mapping(const Map& m)
  {
    std::ostringstream out;
    for (auto i = m.begin(); i != m.end(); ++i)
    {
      out << i->first << " -> " << i->second << std::endl;
    }
    return out.str();
  }

  void check_consistency(const data::sort_expression& s, const data::sort_expression& t, const data::function_symbol& f, sort_expression_substitution_map& sigmaS) const
  {
    sort_expression_substitution_map::const_iterator i = sigmaS.find(s);
    if (i != sigmaS.end() && i->second != t)
    {
      throw mcrl2::runtime_error("inconsistent abstraction " + data::pp(s) + " := " + data::pp(t) + " detected in the abstraction of " + print_symbol(f) + " (elsewhere it is abstracted as " + data::pp(s) + " := " + data::pp(i->second) + ").");
    }
    else
    {
      sigmaS[s] = t;
    }
  }

  // Let f1: s1 x ... sn -> s and f2: t1 x ... tn -> t
  // This function checks if the correspondence si -> ti conflicts with sigmaS.
  void check_consistency(const data::function_symbol& f1, const data::function_symbol& f2, sort_expression_substitution_map& sigmaS) const
  {
    const data::sort_expression& s1 = f1.sort();
    const data::sort_expression& s2 = f2.sort();

    if (data::is_basic_sort(s1))
    {
      check_consistency(s1, s2, f1, sigmaS);
    }
    else if (data::is_function_sort(s1))
    {
      data::function_sort fs1(s1);
      data::function_sort fs2(s2);

      const data::sort_expression_list& domain1 = fs1.domain();
      const data::sort_expression_list& domain2 = fs2.domain();

      data::sort_expression_list::iterator i1 = domain1.begin();
      data::sort_expression_list::iterator i2 = domain2.begin();

      for (; i1 != domain1.end(); ++i1, ++i2)
      {
        check_consistency(*i1, *i2, f1, sigmaS);
      }
    }
//    else if (data::is_container_sort(s1))
//    {
//    }
  }

  // add lifted mappings and equations to the data specification
  void lift_data_specification(const pbes& p, const abstraction_map& sigmaH, const sort_expression_substitution_map& sigmaS, function_symbol_substitution_map& sigmaF, data::data_specification& dataspec)
  {
    using utilities::detail::has_key;

    sort_expression_substitution_map sigmaS_consistency = sigmaS; // is only used for consistency checking
    sort_function sigma(sigmaH, sigmaS, sigmaF, m_generator);

    // add lifted versions of used function symbols that are not specified by the user to sigmaF, and adds them to the data specification as well
    std::set<data::function_symbol> used_function_symbols = pbes_system::find_function_symbols(p);

    // add List containers for user defined sorts, since they are used in the translation
    const data::basic_sort_vector& sorts = dataspec.user_defined_sorts();
    for (const data::basic_sort& sort: sorts)
    {
      data::sort_expression s = data::container_sort(data::list_container(), sort);
      dataspec.add_context_sort(s);
    }

    // add constructor functions of List containers of abstracted sorts to sigmaF
    for (const auto& i: sigmaH)
    {
      data::sort_expression s = data::container_sort(data::list_container(), i.first);
      dataspec.add_context_sort(s);
      data::function_symbol_vector list_constructors = dataspec.constructors(s);
      for (const data::function_symbol& f1: list_constructors)
      {
        data::function_symbol f2 = lift_function_symbol_1_2()(f1, sigma);
        sigmaF[f1] = f2;
        dataspec.add_mapping(f2);
        mCRL2log(log::debug) << "adding list constructor " << f1 << " to sigmaF" << std::endl;
      }
    }

    for (const data::function_symbol& f1: used_function_symbols)
    {
      mCRL2log(log::debug) << "lifting function symbol: " << f1 << std::endl;
      if (!has_key(sigmaF, f1))
      {
        data::function_symbol f2 = lift_function_symbol_1_2()(f1, sigma);
        mCRL2log(log::debug) << "lifted function symbol: " << f1 << " to " << f2 << std::endl;
        check_consistency(f1, f2, sigmaS_consistency);
        sigmaF[f1] = f2;
        dataspec.add_mapping(f2);

        data::data_equation eq = lift_equation_1_2()(f1, f2, sigma, sigmaH);
        mCRL2log(log::debug) << "adding equation: " << eq << std::endl;
        dataspec.add_equation(eq);
      }
    }

    for (auto& i : sigmaF)
    {
      data::function_symbol f2 = i.second;
      data::function_symbol f3 = lift_function_symbol_2_3()(f2);

      mCRL2log(log::debug) << "adding mapping: " << f3 << " " << f3.sort() << std::endl;
      dataspec.add_mapping(f3);

      // update sigmaF
      i.second = f3;

      // make an equation for the lifted function symbol f
      data::data_equation eq = lift_equation_2_3()(f2, f3, sigma);
      mCRL2log(log::debug) << "adding equation: " << eq << std::endl;
      dataspec.add_equation(eq);
    }
  }

  void print_fsvec(const data::function_symbol_vector& v, const std::string& msg) const
  {
    mCRL2log(log::debug) << "--- " << msg << std::endl;
    for (const data::function_symbol& f: v)
    {
      mCRL2log(log::debug) << print_symbol(f) << std::endl;
    }
  }

  void print_fsmap(const function_symbol_substitution_map& v, const std::string& msg) const
  {
    mCRL2log(log::debug) << "--- " << msg << std::endl;
    for (const auto& i: v)
    {
      mCRL2log(log::debug) << print_symbol(i.first) << "  -->  " << print_symbol(i.second) << std::endl;
    }
  }

  void enable_logging()
  {
    log::logger::set_reporting_level(log::debug);
  }

  void run(pbes& p, const std::string& abstraction_text, bool is_over_approximation)
  {
    // split the string abstraction_text into four different parts
    std::string function_symbol_mapping_text;
    std::string user_sorts_text;
    std::string user_equations_text;
    std::string abstraction_mapping_text;
    std::string pbes_sorts_text;

    std::string text = abstraction_text;
    std::vector<std::string> all_keywords = { "sort", "var", "eqn", "map", "cons", "absfunc", "absmap" };
    std::pair<std::string, std::string> q;

    q = utilities::detail::separate_keyword_section(text, "sort", all_keywords);
    user_sorts_text = q.first;
    text = q.second;

    q = utilities::detail::separate_keyword_section(text, "absmap", all_keywords);

    abstraction_mapping_text = q.first;
    text = q.second;

    // must be the last one!
    q = utilities::detail::separate_keyword_section(text, "absfunc", all_keywords);
    function_symbol_mapping_text = q.first;
    user_equations_text = q.second;

    // extract pbes sorts
    std::string ptext = data::pp(p.data());
    q = utilities::detail::separate_keyword_section(ptext, "sort", all_keywords);
    pbes_sorts_text = q.first;

    // 0) split user_dataspec_text into user_sorts_text and user_equations_text
    mCRL2log(log::debug) << "--- user sorts ---\n" << user_sorts_text << std::endl;
    mCRL2log(log::debug) << "--- user equations ---\n" << user_equations_text << std::endl;
    mCRL2log(log::debug) << "--- function mapping ---\n" << function_symbol_mapping_text << std::endl;
    mCRL2log(log::debug) << "--- abstraction mapping ---\n" << abstraction_mapping_text << std::endl;
    mCRL2log(log::debug) << "--- pbes sorts ---\n" << pbes_sorts_text << std::endl;

    if (!abstraction_mapping_text.starts_with("absmap"))
    {
      throw mcrl2::runtime_error("the abstraction mapping may not be empty!");
    }

    // 1) create the data specification dataspec, which consists of user_sorts_text, abstract_mapping_text and p.data()
    data::data_specification dataspec = data::parse_data_specification(data::pp(p.data()) + "\n" + user_sorts_text + "\n" + abstraction_mapping_text.substr(3));
    mCRL2log(log::debug) << "--- data specification 1) ---\n" << dataspec << std::endl;

    // 2) parse the right hand sides of the function symbol mapping, and add them to dataspec
    parse_right_hand_sides(function_symbol_mapping_text, dataspec);
    mCRL2log(log::debug) << "--- data specification 2) ---\n" << dataspec << std::endl;

    // 3) add user_equations_text to dataspec
    dataspec = data::parse_data_specification(data::pp(dataspec) + "\n" + user_equations_text);
    mCRL2log(log::debug) << "--- data specification 3) ---\n" << dataspec << std::endl;

    // abstraction functions (specified by the user)
    abstraction_map sigmaH = parse_abstraction_map(pbes_sorts_text + "\n" + user_sorts_text + "\n" + abstraction_mapping_text.substr(3));

    // sort expressions replacements (extracted from sigmaH)
    sort_expression_substitution_map sigmaS;
    for (auto& i: sigmaH)
    {
      data::function_symbol f = i.second;
      const data::function_sort& fs = atermpp::down_cast<data::function_sort>(f.sort());
      sigmaS[i.first] = fs.codomain();
    }
    mCRL2log(log::debug) << "\n--- sort expression mapping ---\n" << print_mapping(sigmaS) << std::endl;

    // function symbol replacements (specified by the user)
    function_symbol_substitution_map sigmaF = parse_function_symbol_mapping(function_symbol_mapping_text, dataspec);
    mCRL2log(log::debug) << "\n--- function symbol mapping ---\n" << print_mapping(sigmaF) << std::endl;

    m_generator.add_identifiers(data::function_and_mapping_identifiers(p.data()));
    m_generator.add_identifiers(data::function_and_mapping_identifiers(dataspec));

    // 4) add lifted sorts, mappings and equations to dataspec
    // before: the mapping sigmaF is f1 -> f2
    // after: the mapping sigmaF is f1 -> f3
    // after: f2 and f3 have been added to dataspec
    // after: equations for f3 have been added to dataspec
    // generate mapping f1 -> f2 for missing function symbols
    lift_data_specification(p, sigmaH, sigmaS, sigmaF, dataspec);
    mCRL2log(log::debug) << "--- data specification 4) ---\n" << dataspec << std::endl;

    mCRL2log(log::debug) << "\n--- function symbol mapping after lifting ---\n" << print_mapping(sigmaF) << std::endl;

    mCRL2log(log::debug) << "--- pbes before ---\n" << p << std::endl;

    p.data() = dataspec;

    // then transform the data expressions and the propositional variable instantiations
    absinthe_data_expression_builder(sigmaH, sigmaS, sigmaF, m_generator, is_over_approximation).update(p);

    mCRL2log(log::debug) << "--- pbes after ---\n" << p << std::endl;
  }
};

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_ABSINTHE_H
