// Author(s): Jeroen Keiren, Jeroen van der Wulp, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <iostream>

#include "mcrl2/core/load_aterm.h"
#include "mcrl2/data/data_specification.h"
#include "mcrl2/data/detail/data_utility.h"
#include "mcrl2/data/detail/io.h"
#include "mcrl2/data/replace.h"
#include "mcrl2/data/substitutions/sort_expression_assignment.h"

// Predefined datatypes
#include "mcrl2/data/function_update.h"
#include "mcrl2/data/list.h"


namespace mcrl2::data
{

class finiteness_helper
{
  protected:

    const data_specification& m_specification;
    std::set< sort_expression > m_visiting;

    bool is_finite_aux(const sort_expression& s)
    {
      const function_symbol_vector& constructors=m_specification.constructors(s);
      if (constructors.empty())
      {
        return false;
      }

      for(const function_symbol& f: constructors)
      {
        if (is_function_sort(f.sort()))
        {
          const function_sort& f_sort=atermpp::down_cast<function_sort>(f.sort());
          const sort_expression_list& l=f_sort.domain();

          for(const sort_expression& e: l)
          {
            if (!is_finite(e))
            {
              return false;
            }
          }
        }
      }
      return true;
    }

  public:

    finiteness_helper(const data_specification& specification) : m_specification(specification)
    { }

    bool is_finite(const sort_expression& s)
    {
      assert(s==normalize_sorts(s,m_specification));
      if (m_visiting.count(s)>0)
      {
        return false;
      }

      m_visiting.insert(s);

      bool result=false;
      if (is_basic_sort(s))
      {
        result=is_finite(basic_sort(s));
      }
      else if (is_container_sort(s))
      {
        result=is_finite(container_sort(s));
      }
      else if (is_function_sort(s))
      {
        result=is_finite(function_sort(s));
      }
      else if (is_structured_sort(s))
      {
        result=is_finite(structured_sort(s));
      }

      m_visiting.erase(s);
      return result;
    }

    bool is_finite(const basic_sort& s)
    {
      return is_finite_aux(s);
    }

    bool is_finite(const function_sort& s)
    {
      for(const sort_expression& sort: s.domain())
      {
        if (!is_finite(sort))
        {
          return false;
        }
      }

      return is_finite(s.codomain());
    }

    bool is_finite(const container_sort& s)
    {
      return (s.container_name() == set_container()) ? is_finite(s.element_sort()) : false;
    }

    bool is_finite(const alias&)
    {
      assert(0);
      return false;
    }

    bool is_finite(const structured_sort& s)
    {
      return is_finite_aux(s);
    }
};

/// \brief Checks whether a sort is certainly finite.
///
/// \param[in] s A sort expression
/// \return true if s can be determined to be finite,
///      false otherwise.
bool data_specification::is_certainly_finite(const sort_expression& s) const
{
  const bool result=finiteness_helper(*this).is_finite(s);
  return result;
}


// The function below checks whether there is an alias loop, e.g. aliases
// of the form A=B; B=A; or more complex A=B->C; B=Set(D); D=List(A); Loops
// through structured sorts are allowed. If a loop is detected, an exception
// is thrown.
void sort_specification::check_for_alias_loop(
  const sort_expression& s,
  std::set<sort_expression> sorts_already_seen,
  const bool toplevel) const
{
  if (is_basic_sort(s))
  {
    if (sorts_already_seen.count(s)>0)
    {
      throw mcrl2::runtime_error("Sort alias " + pp(s) + " is defined in terms of itself.");
    }
    for(const alias& a: m_user_defined_aliases)
    {
      if (a.name() == s)
      {
        sorts_already_seen.insert(s);
        check_for_alias_loop(a.reference(), sorts_already_seen, true);
        sorts_already_seen.erase(s);
        return;
      }
    }
    return;
  }

  if (is_container_sort(s))
  {
    check_for_alias_loop(container_sort(s).element_sort(),sorts_already_seen,false);
    return;
  }

  if (is_function_sort(s))
  {
    sort_expression_list s_domain(function_sort(s).domain());
    for(const sort_expression& sort: s_domain)
    {
      check_for_alias_loop(sort,sorts_already_seen,false);
    }

    check_for_alias_loop(function_sort(s).codomain(),sorts_already_seen,false);
    return;
  }

  // A sort declaration with a struct on toplevel can be recursive. Otherwise a
  // check needs to be made.
  if (is_structured_sort(s) && !toplevel)
  {
    const structured_sort ss(s);
    structured_sort_constructor_list constructors=ss.constructors();
    for(const structured_sort_constructor& constructor: constructors)
    {
      structured_sort_constructor_argument_list ssca=constructor.arguments();
      for(const structured_sort_constructor_argument& a: ssca)
      {
        check_for_alias_loop(a.sort(),sorts_already_seen,false);
      }
    }
  }

}


// This function returns the normal form of e, under the the map map1.
// This normal form is obtained by repeatedly applying map1, until this
// is not possible anymore. It is assumed that this procedure terminates. There is
// no check for loops.
static sort_expression find_normal_form(
  const sort_expression& e,
  const std::multimap< sort_expression, sort_expression >& map1,
  std::set < sort_expression > sorts_already_seen = std::set < sort_expression >())
{
  assert(sorts_already_seen.find(e)==sorts_already_seen.end()); // e has not been seen already.
  assert(!is_untyped_sort(e));
  assert(!is_untyped_possible_sorts(e));

  if (is_function_sort(e))
  {
    const function_sort fs(e);
    const sort_expression normalised_codomain=
      find_normal_form(fs.codomain(),map1,sorts_already_seen);
    const sort_expression_list& domain=fs.domain();
    sort_expression_list normalised_domain;
    for(const sort_expression& s: domain)
    {
      normalised_domain.push_front(find_normal_form(s,map1,sorts_already_seen));
    }
    return function_sort(reverse(normalised_domain),normalised_codomain);
  }

  if (is_container_sort(e))
  {
    const container_sort cs(e);
    return container_sort(cs.container_name(),find_normal_form(cs.element_sort(),map1,sorts_already_seen));
  }

  sort_expression result_sort;

  if (is_structured_sort(e))
  {
    const structured_sort ss(e);
    structured_sort_constructor_list constructors=ss.constructors();
    structured_sort_constructor_list normalised_constructors;
    for(const structured_sort_constructor& constructor: constructors)
    {
      structured_sort_constructor_argument_list normalised_ssa;
      for(const structured_sort_constructor_argument& a: constructor.arguments())
      {
        normalised_ssa.push_front(structured_sort_constructor_argument(a.name(),
                                      find_normal_form(a.sort(),map1,sorts_already_seen)));
      }

      normalised_constructors.push_front(
                                structured_sort_constructor(
                                  constructor.name(),
                                  reverse(normalised_ssa),
                                  constructor.recogniser()));

    }
    result_sort=structured_sort(reverse(normalised_constructors));
  }

  if (is_basic_sort(e))
  {
    result_sort=e;
  }


  assert(is_basic_sort(result_sort) || is_structured_sort(result_sort));
  const std::multimap< sort_expression, sort_expression >::const_iterator i1=map1.find(result_sort);
  if (i1!=map1.end()) // found
  {
#ifndef NDEBUG
    sorts_already_seen.insert(result_sort);
#endif
    return find_normal_form(i1->second,map1,sorts_already_seen);
  }
  return result_sort;
}

void sort_specification::add_predefined_basic_sorts()
{
    add_system_defined_sort(sort_bool::bool_());
    add_system_defined_sort(sort_pos::pos());
}

void sort_specification::import_system_defined_sort(const sort_expression& sort)
{
  if (is_untyped_sort(sort) || is_untyped_possible_sorts(sort))
  {
    mCRL2log(mcrl2::log::debug) << "Erroneous attempt to insert an untyped sort into the a sort specification\n";
    return;
  }
  // Add an element, and stop if it was already added.
  if (!m_sorts_in_context.insert(sort).second)
  {
    return;
  }

  sorts_are_not_necessarily_normalised_anymore();
  // add the sorts on which this sorts depends.
  if (sort == sort_real::real_())
  {
    // Int is required as the rewrite rules of Real rely on it.
    import_system_defined_sort(sort_int::int_());
  }
  else if (sort == sort_int::int_())
  {
    // See above, Int requires Nat.
    import_system_defined_sort(sort_nat::nat());
  }
  else if (sort == sort_nat::nat())
  {
#ifdef MCRL2_ENABLE_MACHINENUMBERS
    import_system_defined_sort(sort_machine_word::machine_word());
    import_system_defined_sort(sort_nat::natnatpair());
#else
    // Nat requires NatPair.
    import_system_defined_sort(sort_nat::natpair());
#endif
  }
  else if (sort == sort_pos::pos())
  {
#ifdef MCRL2_ENABLE_MACHINENUMBERS
    import_system_defined_sort(sort_machine_word::machine_word());
#endif
  }
  else if (is_function_sort(sort))
  {
    const function_sort& fsort=atermpp::down_cast<function_sort>(sort);
    import_system_defined_sorts(fsort.domain());
    import_system_defined_sort(fsort.codomain());
  }
  else if (is_container_sort(sort))
  {
    const sort_expression element_sort(container_sort(sort).element_sort());
    // Import the element sort (which may be a complex sort also).
    import_system_defined_sort(element_sort);
    if (sort_list::is_list(sort))
    {
      import_system_defined_sort(sort_nat::nat()); // Required for lists.
    }
    else if (sort_set::is_set(sort))
    {
      import_system_defined_sort(sort_fset::fset(element_sort));
      // Import the functions from element_sort->Bool.
      sort_expression_list element_sorts;
      element_sorts.push_front(element_sort);
      import_system_defined_sort(function_sort(element_sorts,sort_bool::bool_()));
    }
    else if (sort_fset::is_fset(sort))
    {
      import_system_defined_sort(sort_nat::nat());
    }
    else if (sort_bag::is_bag(sort))
    {
      // Add the sorts Nat and set_(element_sort) to the specification.
      import_system_defined_sort(sort_nat::nat()); // Required for bags.
      import_system_defined_sort(sort_set::set_(element_sort));
      import_system_defined_sort(sort_fbag::fbag(element_sort));

      // Add the function sort element_sort->Nat to the specification
      sort_expression_list element_sorts ;
      element_sorts.push_front(element_sort);
      import_system_defined_sort(function_sort(element_sorts,sort_nat::nat()));
    }
    else if (sort_fbag::is_fbag(sort))
    {
      import_system_defined_sort(sort_nat::nat()); // Required for bags.
      import_system_defined_sort(sort_fset::fset(sort));
    }
  }
  else if (is_structured_sort(sort))
  {
    structured_sort s_sort(sort);
    function_symbol_vector f(s_sort.constructor_functions(sort));
    for(const function_symbol& f: s_sort.constructor_functions(sort))
    {
      import_system_defined_sort(f.sort());
    }
  }
}

// The function below recalculates m_normalised_aliases, such that
// it forms a confluent terminating rewriting system using which
// sorts can be normalised.
// This algorithm is described in the document: algorithm-for-sort-equivalence.tex in
// the developers library of the mCRL2 toolset.
void sort_specification::reconstruct_m_normalised_aliases() const
{
  // First reset the normalised aliases and the mappings and constructors that have been
  // inherited to basic sort aliases during a previous round of sort normalisation.
  m_normalised_aliases.clear();

  // This is the first step of the algorithm.
  // Check for loops in the aliases. The type checker should already have done this,
  // but we check it again here. If there is a loop m_normalised_aliases will not be
  // built.
    for(const alias& a: m_user_defined_aliases)
    {
      std::set < sort_expression > sorts_already_seen; // Empty set.
    try
    {
      check_for_alias_loop(a.name(),sorts_already_seen,true);
    }
    catch (mcrl2::runtime_error &)
    {
      mCRL2log(log::debug) << "Encountered an alias loop in the alias for " << a.name() <<". The normalised aliases are not constructed\n";
      return;
    }
  }

  // This is the second step of the algorithm.
  // Copy m_normalised_aliases. All aliases are stored from left to right,
  // except structured sorts, which are stored from right to left. The reason is
  // that structured sorts can be recursive, and therefore, they cannot be
  // rewritten from left to right, as this can cause sorts to be infinitely rewritten.

  std::multimap< sort_expression, sort_expression > sort_aliases_to_be_investigated;
  std::multimap< sort_expression, sort_expression > resulting_normalized_sort_aliases;

  for(const alias& a: m_user_defined_aliases)
  {
    if (is_structured_sort(a.reference()))
    {
      sort_aliases_to_be_investigated.insert(std::pair<sort_expression,sort_expression>(a.reference(),a.name()));
    }
    else
    {
      resulting_normalized_sort_aliases.insert(std::pair<sort_expression,sort_expression>(a.name(),a.reference()));
    }
  }

  // Apply Knuth-Bendix completion to the rules in m_normalised_aliases.
  for(; !sort_aliases_to_be_investigated.empty() ;)
  {
    const std::multimap< sort_expression, sort_expression >::iterator it=sort_aliases_to_be_investigated.begin();
    const sort_expression lhs=it->first;
    const sort_expression rhs=it->second;
    sort_aliases_to_be_investigated.erase(it);

    for(const std::pair<const sort_expression, sort_expression >& p: resulting_normalized_sort_aliases)
    {
      const sort_expression s1=data::replace_sort_expressions(lhs,sort_expression_assignment(p.first,p.second), true);

      if (s1!=lhs)
      {
        // There is a conflict between the two sort rewrite rules.
        assert(is_basic_sort(rhs));
        // Choose the normal form on the basis of a lexicographical ordering. This guarantees
        // uniqueness of normal forms over different tools. Ordering on addresses (as used previously)
        // proved to be unstable over different tools.
        const bool rhs_to_s1 = is_basic_sort(s1) && pp(basic_sort(s1))<=pp(rhs);
        const sort_expression left_hand_side=(rhs_to_s1?rhs:s1);
        const sort_expression pre_normal_form=(rhs_to_s1?s1:rhs);
        assert(is_basic_sort(pre_normal_form));
        const sort_expression& e1=pre_normal_form;
        if (e1!=left_hand_side)
        {
          const sort_expression normalised_lhs=find_normal_form(left_hand_side,resulting_normalized_sort_aliases);
          // Check whether the inserted sort rewrite rule is already in sort_aliases_to_be_investigated.
          if (std::find_if(sort_aliases_to_be_investigated.lower_bound(normalised_lhs),
                        sort_aliases_to_be_investigated.upper_bound(normalised_lhs),
                        [&rhs](const std::pair<sort_expression,sort_expression>& x){ return x.second==rhs; })
                   == sort_aliases_to_be_investigated.upper_bound(normalised_lhs)) // Not found.
          {
            sort_aliases_to_be_investigated.insert(
                  std::pair<sort_expression,sort_expression > (normalised_lhs, e1));
          }
        }
      }
      else
      {
        const sort_expression s2 = data::replace_sort_expressions(p.first,sort_expression_assignment(lhs,rhs), true);
        if (s2!=p.first)
        {
          assert(is_basic_sort(p.second));
          // Choose the normal form on the basis of a lexicographical ordering. This guarantees
          // uniqueness of normal forms over different tools.
          const bool i_second_to_s2 = is_basic_sort(s2) && pp(basic_sort(s2))<=pp(p.second);
          const sort_expression left_hand_side=(i_second_to_s2?p.second:s2);
          const sort_expression pre_normal_form=(i_second_to_s2?s2:p.second);
          assert(is_basic_sort(pre_normal_form));
          const sort_expression& e2=pre_normal_form;
          if (e2!=left_hand_side)
          {
            const sort_expression normalised_lhs=find_normal_form(left_hand_side,resulting_normalized_sort_aliases);
            // Check whether the inserted sort rewrite rule is already in sort_aliases_to_be_investigated.
            if (std::find_if(sort_aliases_to_be_investigated.lower_bound(normalised_lhs),
                          sort_aliases_to_be_investigated.upper_bound(normalised_lhs),
                          [&rhs](const std::pair<sort_expression,sort_expression>& x){ return x.second==rhs; })
                     == sort_aliases_to_be_investigated.upper_bound(normalised_lhs)) // Not found.
            {
              sort_aliases_to_be_investigated.insert(
                    std::pair<sort_expression,sort_expression > (normalised_lhs,e2));
            }
          }
        }
      }
    }
    assert(lhs!=rhs);
    const sort_expression normalised_lhs = find_normal_form(lhs,resulting_normalized_sort_aliases);
    const sort_expression normalised_rhs = find_normal_form(rhs,resulting_normalized_sort_aliases);
    if (normalised_lhs!=normalised_rhs)
    {
      resulting_normalized_sort_aliases.insert(std::pair<sort_expression,sort_expression >(normalised_lhs,normalised_rhs));
    }
  }
  // Copy resulting_normalized_sort_aliases into m_normalised_aliases, i.e. from multimap to map.
  // If there are rules with equal left hand side, only one is arbitrarily chosen. Rewrite the
  // right hand side to normal form.

  for(const std::pair<const sort_expression,sort_expression>& p: resulting_normalized_sort_aliases)
  {
    const sort_expression normalised_rhs = find_normal_form(p.second,resulting_normalized_sort_aliases);
    m_normalised_aliases[p.first]=normalised_rhs;

    assert(p.first!=normalised_rhs);
  }
}

///\brief Adds the system defined sorts to the sets with constructors, mappings, and equations for
//        a given sort. If the boolean skip_equations is true, no equations are added.
void data_specification::find_associated_system_defined_data_types_for_a_sort(
                   const sort_expression& sort,
                   std::set < function_symbol >& constructors,
                   std::set < function_symbol >& mappings,
                   std::set < data_equation >& equations,
                   implementation_map& cpp_implemented_functions,
                   const bool skip_equations) const
{
  // add sorts, constructors, mappings and equations
  if (sort == sort_bool::bool_())
  {
    function_symbol_vector f(sort_bool::bool_generate_constructors_code());
    constructors.insert(f.begin(), f.end());
    f = sort_bool::bool_generate_functions_code();
    mappings.insert(f.begin(), f.end());
    implementation_map f1 = sort_bool::bool_cpp_implementable_mappings();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    f1 = sort_bool::bool_cpp_implementable_constructors();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    if (!skip_equations)
    {
      data_equation_vector e(sort_bool::bool_generate_equations_code());
      equations.insert(e.begin(),e.end());
    }
  }
  else if (sort == sort_real::real_())
  {
    function_symbol_vector f(sort_real::real_generate_constructors_code());
    constructors.insert(f.begin(),f.end());
    f = sort_real::real_generate_functions_code();
    mappings.insert(f.begin(),f.end());
    implementation_map f1 = sort_int::int_cpp_implementable_mappings();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    f1 = sort_int::int_cpp_implementable_constructors();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    if (!skip_equations)
    {
      data_equation_vector e(sort_real::real_generate_equations_code());
      equations.insert(e.begin(),e.end());
    }
  }
  else if (sort == sort_int::int_())
  {
    function_symbol_vector f(sort_int::int_generate_constructors_code());
    constructors.insert(f.begin(),f.end());
    f = sort_int::int_generate_functions_code();
    mappings.insert(f.begin(),f.end());
    implementation_map f1 = sort_int::int_cpp_implementable_mappings();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    f1 = sort_int::int_cpp_implementable_constructors();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    if (!skip_equations)
    {
      data_equation_vector e(sort_int::int_generate_equations_code());
      equations.insert(e.begin(),e.end());
    }
  }
  else if (sort == sort_nat::nat())
  {
    function_symbol_vector f(sort_nat::nat_generate_constructors_code());
    constructors.insert(f.begin(),f.end());
    f = sort_nat::nat_generate_functions_code();
    mappings.insert(f.begin(),f.end());
    implementation_map f1 = sort_nat::nat_cpp_implementable_mappings();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    f1 = sort_nat::nat_cpp_implementable_constructors();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    if (!skip_equations)
    {
      data_equation_vector e(sort_nat::nat_generate_equations_code());
      equations.insert(e.begin(),e.end());
    }
  }
  else if (sort == sort_pos::pos())
  {
    function_symbol_vector f(sort_pos::pos_generate_constructors_code());
    constructors.insert(f.begin(),f.end());
    f = sort_pos::pos_generate_functions_code();
    mappings.insert(f.begin(),f.end());
    implementation_map f1 = sort_pos::pos_cpp_implementable_mappings();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    f1 = sort_pos::pos_cpp_implementable_constructors();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    if (!skip_equations)
    {
      data_equation_vector e(sort_pos::pos_generate_equations_code());
      equations.insert(e.begin(),e.end());
    }
  }
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  else if (sort == sort_machine_word::machine_word())
  {
    function_symbol_vector f(sort_machine_word::machine_word_generate_constructors_code());
    constructors.insert(f.begin(),f.end());
    f = sort_machine_word::machine_word_generate_functions_code();
    mappings.insert(f.begin(),f.end());
    implementation_map f1 = sort_machine_word::machine_word_cpp_implementable_mappings();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    f1 = sort_machine_word::machine_word_cpp_implementable_constructors();
    cpp_implemented_functions.insert(f1.begin(), f1.end());
    if (!skip_equations)
    {
      data_equation_vector e(sort_machine_word::machine_word_generate_equations_code());
      equations.insert(e.begin(),e.end());
    }
  }
#endif
  else if (is_function_sort(sort))
  {
    const sort_expression& t = static_cast<const function_sort&>(sort).codomain();
    const sort_expression_list& l = static_cast<const function_sort&>(sort).domain();
    if (l.size()==1)
    {
      const function_symbol_vector f = function_update_generate_functions_code(l.front(),t);
      mappings.insert(f.begin(),f.end());
      implementation_map f1 = function_update_cpp_implementable_mappings(l.front(),t);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      f1 = function_update_cpp_implementable_constructors();
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      if (!skip_equations)
      {
        data_equation_vector e(function_update_generate_equations_code(l.front(),t));
        equations.insert(e.begin(),e.end());
      }
    }
  }
  else if (is_container_sort(sort))
  {
    sort_expression element_sort(container_sort(sort).element_sort());
    if (sort_list::is_list(sort))
    {
      function_symbol_vector f(sort_list::list_generate_constructors_code(element_sort));
      constructors.insert(f.begin(),f.end());
      f = sort_list::list_generate_functions_code(element_sort);
      mappings.insert(f.begin(),f.end());
      implementation_map f1 = sort_list::list_cpp_implementable_mappings(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      f1 = sort_list::list_cpp_implementable_constructors(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      if (!skip_equations)
      {
        data_equation_vector e(sort_list::list_generate_equations_code(element_sort));
        equations.insert(e.begin(),e.end());
      }
    }
    else if (sort_set::is_set(sort))
    {
      sort_expression_list element_sorts;
      element_sorts.push_front(element_sort);
      function_symbol_vector f(sort_set::set_generate_constructors_code(element_sort));
      constructors.insert(f.begin(),f.end());
      f = sort_set::set_generate_functions_code(element_sort);
      mappings.insert(f.begin(),f.end());
      implementation_map f1 = sort_set::set_cpp_implementable_mappings(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      f1 = sort_set::set_cpp_implementable_constructors(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      if (!skip_equations)
      {
        data_equation_vector e(sort_set::set_generate_equations_code(element_sort));
        equations.insert(e.begin(),e.end());
      }
    }
    else if (sort_fset::is_fset(sort))
    {
      function_symbol_vector f = sort_fset::fset_generate_constructors_code(element_sort);
      constructors.insert(f.begin(),f.end());
      f = sort_fset::fset_generate_functions_code(element_sort);
      mappings.insert(f.begin(),f.end());
      implementation_map f1 = sort_fset::fset_cpp_implementable_mappings(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      f1 = sort_fset::fset_cpp_implementable_constructors(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      if (!skip_equations)
      {
        data_equation_vector e = sort_fset::fset_generate_equations_code(element_sort);
        equations.insert(e.begin(),e.end());
      }
    }
    else if (sort_bag::is_bag(sort))
    {
      sort_expression_list element_sorts;
      element_sorts.push_front(element_sort);
      function_symbol_vector f(sort_bag::bag_generate_constructors_code(element_sort));
      constructors.insert(f.begin(),f.end());
      f = sort_bag::bag_generate_functions_code(element_sort);
      mappings.insert(f.begin(),f.end());
      implementation_map f1 = sort_bag::bag_cpp_implementable_mappings(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      f1 = sort_bag::bag_cpp_implementable_constructors(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      if (!skip_equations)
      {
        data_equation_vector e(sort_bag::bag_generate_equations_code(element_sort));
        equations.insert(e.begin(),e.end());
      }
    }
    else if (sort_fbag::is_fbag(sort))
    {
      function_symbol_vector f = sort_fbag::fbag_generate_constructors_code(element_sort);
      constructors.insert(f.begin(),f.end());
      f = sort_fbag::fbag_generate_functions_code(element_sort);
      mappings.insert(f.begin(),f.end());
      implementation_map f1 = sort_fbag::fbag_cpp_implementable_mappings(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      f1 = sort_fbag::fbag_cpp_implementable_constructors(element_sort);
      cpp_implemented_functions.insert(f1.begin(), f1.end());
      if (!skip_equations)
      {
        data_equation_vector e = sort_fbag::fbag_generate_equations_code(element_sort);
        equations.insert(e.begin(),e.end());
      }
    }
  }
  else if (is_structured_sort(sort))
  {
    insert_mappings_constructors_for_structured_sort(
                    atermpp::down_cast<structured_sort>(sort),
                    constructors, mappings, equations, skip_equations);
  }
  add_standard_mappings_and_equations(sort, mappings, equations, skip_equations);
}

void data_specification::get_system_defined_sorts_constructors_and_mappings(
            std::set < sort_expression >& sorts,
            std::set < function_symbol >& constructors,
            std::set <function_symbol >& mappings) const
{
  implementation_map cpp_implemented_functions;

  sorts.insert(sort_bool::bool_());
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  sorts.insert(sort_machine_word::machine_word());
#endif
  sorts.insert(sort_pos::pos());
  sorts.insert(sort_nat::nat());
  sorts.insert(sort_int::int_());
  sorts.insert(sort_real::real_());
  sorts.insert(sort_list::list(sort_pos::pos()));
  sorts.insert(sort_fset::fset(sort_pos::pos()));
  sorts.insert(sort_set::set_(sort_pos::pos()));
  sorts.insert(sort_fbag::fbag(sort_pos::pos()));
  sorts.insert(sort_bag::bag(sort_pos::pos()));

  std::set < data_equation > dummy_equations;
  for(const sort_expression& s: sorts)
  {
    find_associated_system_defined_data_types_for_a_sort(s, constructors, mappings, dummy_equations, cpp_implemented_functions, true);
  }
  assert(dummy_equations.size()==0);
}

bool data_specification::is_well_typed() const
{
  // check 1)
  if (!detail::check_data_spec_sorts(constructors(), sorts()))
  {
    std::clog << "data_specification::is_well_typed() failed: not all of the sorts appearing in the constructors "
              << data::pp(constructors()) << " are declared in " << data::pp(sorts()) << std::endl;
    return false;
  }

  // check 2)
  if (!detail::check_data_spec_sorts(mappings(), sorts()))
  {
    std::clog << "data_specification::is_well_typed() failed: not all of the sorts appearing in the mappings "
              << data::pp(mappings()) << " are declared in " << data::pp(sorts()) << std::endl;
    return false;
  }

  return true;
}

/// There are two types of representations of ATerms:
///  - the bare specification that does not contain constructor, mappings
///    and equations for system defined sorts
///  - specification that includes all system defined information (legacy)
/// The last type must eventually disappear but is unfortunately still in
/// use in a substantial amount of source code.
/// Note, all sorts with name prefix \@legacy_ are eliminated
void data_specification::build_from_aterm(const atermpp::aterm& term)
{
  assert(core::detail::check_rule_DataSpec(term));

  // Note backwards compatibility measure: alias is no longer a sort_expression
  const atermpp::term_list<atermpp::aterm> term_sorts=
                 atermpp::down_cast<atermpp::term_list<atermpp::aterm> >(term[0][0]);
  const data::function_symbol_list term_constructors=
                 atermpp::down_cast<data::function_symbol_list>(term[1][0]);
  const data::function_symbol_list term_mappings=
                 atermpp::down_cast<data::function_symbol_list>(term[2][0]);
  const data::data_equation_list term_equations=
                 atermpp::down_cast<data::data_equation_list>(term[3][0]);

  // Store the sorts and aliases.
  for(const atermpp::aterm& t: term_sorts)
  {
    if (data::is_alias(t)) // Compatibility with legacy code
    {
      add_alias(atermpp::down_cast<data::alias>(t));
    }
    else
    {
      add_sort(atermpp::down_cast<basic_sort>(t));
    }
  }

  // Store the constructors.
  for(const function_symbol& f: term_constructors)
  {
    add_constructor(f);
  }

  // Store the mappings.
  for(const function_symbol& f: term_mappings)
  {
    add_mapping(f);
  }

  // Store the equations.
  for(const data_equation& e: term_equations)
  {
    add_equation(e);
  }
}

data_specification::data_specification(const basic_sort_vector& sorts,
  const alias_vector& aliases,
  const function_symbol_vector& constructors,
  const function_symbol_vector& user_defined_mappings,
  const data_equation_vector& user_defined_equations)
  : sort_specification(sorts, aliases)
{
  // Store the constructors.
  for(const function_symbol& f: constructors)
  {
    add_constructor(f);
  }

  // Store the mappings.
  for(const function_symbol& f: user_defined_mappings)
  {
    add_mapping(f);
  }

  // Store the equations.
  for(const data_equation& e: user_defined_equations)
  {
    add_equation(e);
  }

  assert(is_well_typed());
}

} // namespace mcrl2::data

