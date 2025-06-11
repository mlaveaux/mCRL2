// Author(s): Jeroen van der Wulp, Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/print.h
/// \brief Provides utilities for pretty printing.

#ifndef MCRL2_DATA_PRINT_H
#define MCRL2_DATA_PRINT_H

#include "mcrl2/data/data_specification.h"
#include "mcrl2/data/detail/is_untyped.h"
#include "mcrl2/data/detail/print_utility.h"
#include "mcrl2/data/standard_container_utility.h"

namespace mcrl2
{

namespace data
{

namespace detail {

// Check that the argument, which must be of sort Pos, represents 1.
inline 
bool is_one(const data_expression& x)
{
  assert(x.sort() == sort_pos::pos());
#ifdef MCRL2_ENABLE_MACHINENUMBERS
  return (sort_pos::is_most_significant_digit_application(x) &&
          (( is_machine_number(sort_pos::arg(x)) && atermpp::down_cast<machine_number>(sort_pos::arg(x)).value() == 1) ||
           (  sort_machine_word::is_one_word_function_symbol(sort_pos::arg(x))))
         ) ||
         sort_pos::is_c1_function_symbol(x);
#else
  return sort_pos::is_c1_function_symbol(x);
#endif
}

inline
bool is_numeric_cast(const data_expression& x)
{
  return data::sort_nat::is_pos2nat_application(x)
         || data::sort_int::is_pos2int_application(x)
         || data::sort_real::is_pos2real_application(x)
         || data::sort_int::is_nat2int_application(x)
         || data::sort_real::is_nat2real_application(x)
         || data::sort_real::is_int2real_application(x)
#ifndef MCRL2_ENABLE_MACHINENUMBERS
         || data::sort_nat::is_cnat_application(x)
         || data::sort_int::is_cint_application(x)
#endif
         || data::sort_real::is_creal_application(x)
          ;
}

inline bool look_through_numeric_casts(const data_expression& x, std::function<bool(const data_expression&)> f)
{
  if (is_numeric_cast(x))
  {
    return look_through_numeric_casts(atermpp::down_cast<application>(x)[0], f);
  }
  return f(x);
}

// Check whether this term is an applicatino with a head of the shape "UntypedIdentifier(s)" for
// the given string s.
inline bool head_matches_undefined_symbol(const data_expression& x, const core::identifier_string& s)
{
  if (is_application(x))
  {
    const data_expression& h=atermpp::down_cast<application>(x).head();
    if (is_untyped_identifier(h))
    {
      const untyped_identifier& u=atermpp::down_cast<untyped_identifier>(h);
      return u.name()==s;
    }
  }
  return false;
}

/* inline
data_expression remove_numeric_casts(data_expression x)
{
  while (is_numeric_cast(x))
  {
    x = *atermpp::down_cast<application>(x).begin();
  }
  return x;
} */

inline
bool is_plus(const application& x)
{
  return look_through_numeric_casts(
                x,
                [](const data_expression& x)
                   { 
                     return // sort_int::is_plus_application(x)  ||
                            // sort_nat::is_plus_application(x)  ||
                            // sort_pos::is_plus_application(x)  ||
                            sort_real::is_plus_application(x) ||
                            head_matches_undefined_symbol(x, sort_pos::plus_name()); });
}

inline
bool is_minus(const application& x)
{
  return look_through_numeric_casts(
                x,
                [](const data_expression& x)
                   { return // sort_int::is_minus_application(x)  ||
                            sort_real::is_minus_application(x) ||
                            head_matches_undefined_symbol(x, sort_int::minus_name()); });
}

inline
bool is_mod(const application& x)
{
  return look_through_numeric_casts(
                x,
                [](const data_expression& x)
                   { return sort_int::is_mod_application(x) ||
                            sort_nat::is_mod_application(x); });
}

inline
bool is_div(const application& x)
{
  return look_through_numeric_casts(
                x,
                [](const data_expression& x)
                   { return sort_int::is_div_application(x) ||
                            sort_nat::is_div_application(x); });
}

#ifndef MCRL2_ENABLE_MACHINENUMBERS
inline
bool is_divmod(const application& x)
{
  return look_through_numeric_casts(
                x,
                [](const data_expression& x)
                   { return sort_nat::is_divmod_application(x); });
}
#endif

inline
bool is_divides(const application& x)
{
  return look_through_numeric_casts(
                x,
                [](const data_expression& x)
                   { return sort_real::is_divides_application(x)  ||
                            head_matches_undefined_symbol(x, sort_real::divides_name()); });
}

inline
bool is_implies(const application& x)
{
  return sort_bool::is_implies_application(x);
}

inline
bool is_set_union(const application& x)
{
  return sort_set::is_union_application(x);
}

inline
bool is_set_difference(const application& x)
{
  return sort_set::is_difference_application(x);
}

inline
bool is_bag_join(const application& x)
{
  return sort_bag::is_union_application(x);
}

inline
bool is_bag_difference(const application& x)
{
  return sort_bag::is_difference_application(x);
}

inline
bool is_and(const application& x)
{
  return sort_bool::is_and_application(x);
}

inline
bool is_or(const application& x)
{
  return sort_bool::is_or_application(x);
}

inline
bool is_equal_to(const application& x)
{
  return data::is_equal_to_application(x);
}

inline
bool is_not_equal_to(const application& x)
{
  return data::is_not_equal_to_application(x);
}

inline
bool is_less(const application& x)
{
  return data::is_less_application(x);
}

inline
bool is_less_equal(const application& x)
{
  return data::is_less_equal_application(x);
}

inline
bool is_greater(const application& x)
{
  return data::is_greater_application(x);
}

inline
bool is_greater_equal(const application& x)
{
  return data::is_greater_equal_application(x);
}

inline
bool is_in(const application& x)
{
  return sort_list::is_in_application(x);
}

inline
bool is_times(const application& x)
{
  return look_through_numeric_casts(
                x,
                [](const data_expression& x)
                   { return sort_int::is_times_application(x) ||
                            head_matches_undefined_symbol(x, sort_pos::times_name()); });
}

inline
bool is_element_at(const application& x)
{
  return sort_list::is_element_at_application(x);
}

inline
bool is_set_intersection(const application& x)
{
  return sort_set::is_intersection_application(x);
}

inline
bool is_bag_intersection(const application& x)
{
  return sort_bag::is_intersection_application(x);
}

inline
bool is_concat(const application& x)
{
  return sort_list::is_concat_application(x);
}

inline
bool is_cons_list(data_expression x)
{
  while (sort_list::is_cons_application(x))
  {
    x = sort_list::right(x);
  }
  return sort_list::is_empty_function_symbol(x);
}

inline
bool is_snoc_list(data_expression x)
{
  while (sort_list::is_snoc_application(x))
  {
    x = sort_list::left(x);
  }
  return sort_list::is_empty_function_symbol(x);
}

inline
bool is_cons(const application& x)
{
  return sort_list::is_cons_application(x) && !is_cons_list(x);
}

inline
bool is_snoc(const application& x)
{
  return sort_list::is_snoc_application(x) && !is_snoc_list(x);
}

} // namespace detail

int precedence(const data_expression& x);

inline
int precedence(const application& x)
{
  // N.B. this code should match printing of a creal
  if (sort_real::is_creal_application(x))
  {
    const data_expression& numerator = sort_real::left(x);
    const data_expression& denominator = sort_real::right(x);
    // if (sort_pos::is_c1_function_symbol(denominator))
    if (detail::is_one(denominator))
    {
      return precedence(numerator);
    }
    else
    {
      return precedence(sort_real::divides(numerator, sort_int::pos2int(denominator)));
    }
  }
  else if (detail::is_implies(x))
  {
    return 2;
  }
  else if (detail::is_or(x))
  {
    return 3;
  }
  else if (detail::is_and(x))
  {
    return 4;
  }
  else if (detail::is_equal_to(x) ||
           detail::is_not_equal_to(x)
          )
  {
    return 5;
  }
  else if (   detail::is_less(x)
              || detail::is_less_equal(x)
              || detail::is_greater(x)
              || detail::is_greater_equal(x)
              || detail::is_in(x)
          )
  {
    return 6;
  }
  else if (detail::is_cons(x))
  {
    return 7;
  }
  else if (detail::is_snoc(x))
  {
    return 8;
  }
  else if (detail::is_concat(x))
  {
    return 9;
  }
  else if (   detail::is_plus(x)
              || detail::is_minus(x)
              || detail::is_set_union(x)
              || detail::is_set_difference(x)
              || detail::is_bag_join(x)
              || detail::is_bag_difference(x)
          )
  {
    return 10;
  }
  else if (   detail::is_div(x)
              || detail::is_mod(x)
#ifndef MCRL2_ENABLE_MACHINENUMBERS
              || detail::is_divmod(x)
#endif
              || detail::is_divides(x)
          )
  {
    return 11;
  }
  else if (   detail::is_times(x)
              || detail::is_element_at(x)
              || detail::is_set_intersection(x)
              || detail::is_bag_intersection(x)
          )
  {
    return 12;
  }
  else if (is_function_update_application(x) || is_function_update_stable_application(x))
  {
    return 13;
  }
  // TODO: add function application (there seems to be no recognizer for it)
  return core::detail::max_precedence;
}

constexpr int precedence(const forall&)            { return 1; }
constexpr int precedence(const exists&)            { return 1; }
constexpr int precedence(const lambda&)            { return 1; }
constexpr int precedence(const set_comprehension&) { return core::detail::max_precedence; }
constexpr int precedence(const bag_comprehension&) { return core::detail::max_precedence; }
constexpr int precedence(const where_clause&)      { return 0; }
inline int precedence(const data_expression& x)
{
  if (data::is_application(x))            { return precedence(atermpp::down_cast<application>(x)); }
  else if (data::is_forall(x))            { return precedence(atermpp::down_cast<forall>(x)); }
  else if (data::is_exists(x))            { return precedence(atermpp::down_cast<exists>(x)); }
  else if (data::is_lambda(x))            { return precedence(atermpp::down_cast<lambda>(x)); }
  else if (data::is_set_comprehension(x)) { return precedence(atermpp::down_cast<set_comprehension>(x)); }
  else if (data::is_bag_comprehension(x)) { return precedence(atermpp::down_cast<bag_comprehension>(x)); }
  else if (data::is_where_clause(x))      { return precedence(atermpp::down_cast<where_clause>(x)); }
  return core::detail::max_precedence;
}

inline
bool is_left_associative(const data_expression& x)
{
  return !sort_bool::is_implies_application(x) && !sort_list::is_cons_application(x);
}

inline
bool is_right_associative(const data_expression& x)
{
  if (!is_application(x))
  {
    return false;
  }
  const auto& x_ = atermpp::down_cast<application>(x);
  return !detail::is_minus(x_) && !is_equal_to_application(x);
}

namespace detail
{

template <typename Derived>
struct printer: public data::add_traverser_sort_expressions<core::detail::printer, Derived>
{
  typedef data::add_traverser_sort_expressions<core::detail::printer, Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;
  using super::derived;
  using super::print_expression;
  using super::print_unary_operand;
  using super::print_list;

  void print_unary_data_operation(const application& x, const std::string& op)
  {
    derived().print(op);
    print_expression(x[0], precedence(x[0]) < precedence(x));
  }

  void print_binary_data_operation(const application& x, const data_expression& x1, const data_expression& x2, const std::string& op)
  {
    const int p = precedence(x);
    const int p1 = precedence(x1);
    const int p2 = precedence(x2);
    print_expression(x1, (p1 < p) || (p1 == p && !is_left_associative(x)));
    derived().print(op);
    print_expression(x2, (p2 < p) || (p2 == p && !is_right_associative(x)));
  }

  void print_binary_data_operation(const application& x, const std::string& op)
  {
    const data_expression& x1 = x[0];
    const data_expression& x2 = x[1];
    print_binary_data_operation(x, x1, x2, op);
  }

  bool is_infix_operation(const application& x) const
  {
    if (x.size() != 2)
    {
      return false;
    }
    core::identifier_string name;
    if (is_function_symbol(x.head()))
    {
      name = function_symbol(x.head()).name();
    }
    else if (is_untyped_identifier(x.head()))
    {
      name = untyped_identifier(x.head()).name();
    }
    else
    {
      return false;
    }
    return
      (name == data::sort_bool::implies_name())      ||
      (name == data::sort_bool::and_name())          ||
      (name == data::sort_bool::or_name())           ||
      data::detail::equal_symbol::is_symbol(name)    ||
      data::detail::not_equal_symbol::is_symbol(name)     ||
      data::detail::less_symbol::is_symbol(name)          ||
      data::detail::less_equal_symbol::is_symbol(name)    ||
      data::detail::greater_symbol::is_symbol(name)       ||
      data::detail::greater_equal_symbol::is_symbol(name) ||
      (name == data::sort_list::in_name())           ||
      (name == data::sort_list::cons_name())         ||
      (name == data::sort_list::snoc_name())         ||
      (name == data::sort_list::concat_name())       ||
      (name == data::sort_real::plus_name())         ||
      (name == data::sort_real::minus_name())        ||
      (name == data::sort_set::union_name())         ||
      (name == data::sort_fset::union_name())        ||
      (name == data::sort_set::difference_name())    ||
      (name == data::sort_fset::difference_name())   ||
      (name == data::sort_bag::union_name())         ||
      (name == data::sort_fbag::union_name())        ||
      (name == data::sort_bag::difference_name())    ||
      (name == data::sort_fbag::difference_name())   ||
      (name == data::sort_int::div_name())           ||
      (name == data::sort_int::mod_name())           ||
      (name == data::sort_real::divides_name())      ||
      (name == data::sort_int::times_name())         ||
      (name == data::sort_list::element_at_name())   ||
      (name == data::sort_set::intersection_name())  ||
      (name == data::sort_bag::intersection_name());
  }

  core::identifier_string generate_identifier(const std::string& prefix, const data_expression& context) const
  {
    data::set_identifier_generator generator;
    std::set<variable> variables = data::find_all_variables(context);
    for (const variable& v: variables)
    {
      generator.add_identifier(v.name());
    }
    return generator(prefix);
  }

  template <typename Container>
  void print_container(const Container& container,
                       int container_precedence = -1,
                       const std::string& separator = ", ",
                       const std::string& open_bracket = "(",
                       const std::string& close_bracket = ")"
                      )
  {
    for (auto i = container.begin(); i != container.end(); ++i)
    {
      if (i != container.begin())
      {
        derived().print(separator);
      }
      bool print_brackets = (container.size() > 1) && (precedence(*i) < container_precedence);
      if (print_brackets)
      {
        derived().print(open_bracket);
      }
      derived().apply(*i);
      if (print_brackets)
      {
        derived().print(close_bracket);
      }
    }
  }

  template <typename Variable>
  void print_variable(const Variable& x, bool print_sort = false)
  {
    derived().apply(x);
    if (print_sort)
    {
      derived().print(": ");
      derived().apply(x.sort());
    }
  }

  struct get_sort_default
  {
    template <typename T>
    sort_expression operator()(const T& t) const
    {
      return t.sort();
    }
  };

  template <typename Container, typename SortAccessor>
  void print_sorted_declarations(const Container& container,
                                 bool print_sorts = true,
                                 bool join_sorts = true,
                                 bool maximally_shared = false,
                                 const std::string& opener = "(",
                                 const std::string& closer = ")",
                                 const std::string& separator = ", ",
                                 SortAccessor get_sort = get_sort_default()
                                )
  {
    auto first = container.begin();
    auto last = container.end();
    if (first == last)
    {
      return;
    }

    derived().print(opener);

    if (maximally_shared)
    {
      typedef typename Container::value_type T;

      // sort_map[s] will contain all elements t of container with t.sort() == s.
      std::map<sort_expression, std::vector<T> > sort_map;

      // sorts will contain all sort expressions s that appear as a key in sort_map,
      // in the order they are encountered in container
      std::vector<sort_expression> sorts;

      for (auto i = container.begin(); i != container.end(); ++i)
      {
        if (sort_map.find(i->sort()) == sort_map.end())
        {
          sorts.push_back(i->sort());
        }
        sort_map[i->sort()].push_back(*i);
      }

      // do the actual printing
      for (auto i = sorts.begin(); i != sorts.end(); ++i)
      {
        if (i != sorts.begin())
        {
          derived().print(separator);
        }
        const std::vector<T>& v = sort_map[*i];
        print_list(v, "", "", ",");
        derived().print(": ");
        derived().apply(*i);
      }
    }
    else
    {
      while (first != last)
      {
        if (first != container.begin())
        {
          derived().print(separator);
        }

        if (print_sorts && join_sorts)
        {
          // determine a consecutive interval [first, i) with elements of the same sorts
          auto i = first;
          do
          {
            ++i;
          }

          // print the elements of the interval [first, i)
          while (i != last && i->sort() == first->sort());

          for (auto j = first; j != i; ++j)
          {
            if (j != first)
            {
              derived().print(",");
            }
            derived().apply(*j);
          }

          // print the sort
          if (print_sorts)
          {
            derived().print(": ");
            derived().apply(get_sort(*first));
          }

          // update first
          first = i;
        }
        else
        {
          derived().apply(*first);

          // print the sort
          if (print_sorts)
          {
            derived().print(": ");
            derived().apply(get_sort(*first));
          }

          // update first
          ++first;
        }
      }
    }
    derived().print(closer);
  }

  template <typename Container>
  void print_variables(const Container& container,
                       bool print_sorts = true,
                       bool join_sorts = true,
                       bool maximally_shared = false,
                       const std::string& opener = "(",
                       const std::string& closer = ")",
                       const std::string& separator = ", "
                      )
  {
    print_sorted_declarations(container, print_sorts, join_sorts, maximally_shared, opener, closer, separator, get_sort_default());
  }

  template <typename Container>
  void print_assignments(const Container& container,
                         bool print_lhs = true,
                         const std::string& opener = "",
                         const std::string& closer = "",
                         const std::string& separator = ", ",
                         const std::string& assignment_symbol = " = "
                        )
  {
    if (container.empty())
    {
      return;
    }
    derived().print(opener);
    for (typename Container::iterator i = container.begin(); i != container.end(); ++i)
    {
      if (i != container.begin())
      {
        derived().print(separator);
      }
      if (print_lhs)
      {
        derived().apply(i->lhs());
        derived().print(assignment_symbol);
      }
      derived().apply(i->rhs());
    }
    derived().print(closer);
  }

  template <typename T>
  void print_condition(const T& x, const std::string& arrow = "  ->  ")
  {
    if (!sort_bool::is_true_function_symbol(x))
    {
      print_expression(x, true);
      derived().print(arrow);
    }
  }

  template <typename Container>
  void print_sort_list(const Container& container,
                       const std::string& opener = "(",
                       const std::string& closer = ")",
                       const std::string& separator = ", "
                      )
  {
    if (container.empty())
    {
      return;
    }
    derived().print(opener);
    for (auto i = container.begin(); i != container.end(); ++i)
    {
      if (i != container.begin())
      {
        derived().print(separator);
      }
      bool print_brackets = is_function_sort(*i);
      if (print_brackets)
      {
        derived().print("(");
      }
      derived().apply(*i);
      if (print_brackets)
      {
        derived().print(")");
      }
    }
    derived().print(closer);
  }

  void print_list_enumeration(const application& x)
  {
    derived().print("[");
    print_container(x, precedence(x));
    derived().print("]");
  }

  void print_set_enumeration(const application& x)
  {
    derived().print("{ ");
    print_container(x, precedence(x));
    derived().print(" }");
  }

  void print_bag_enumeration(const application& x)
  {
    derived().print("{ ");
    application::const_iterator i = x.begin();
    while (i != x.end())
    {
      if (i != x.begin())
      {
        derived().print(", ");
      }
      derived().apply(*i++);
      derived().print(": ");
      derived().apply(*i++);
    }
    derived().print(" }");
  }

  void print_setbag_comprehension(const abstraction& x)
  {
    derived().print("{ ");
    print_variables(x.variables(), true, true, false, "", "", ", ");
    derived().print(" | ");
    derived().apply(x.body());
    derived().print(" }");
  }

  bool is_abstraction_application(const application& x) const
  {
    return is_abstraction(x.head());
  }

  bool is_cons_list(data_expression x) const
  {
    while (sort_list::is_cons_application(x))
    {
      x = sort_list::right(x);
    }
    return sort_list::is_empty_function_symbol(x);
  }

  bool is_snoc_list(data_expression x) const
  {
    while (sort_list::is_snoc_application(x))
    {
      x = sort_list::left(x);
    }
    return sort_list::is_empty_function_symbol(x);
  }

  bool is_fset_cons_list(data_expression x)
  {
    while (sort_fset::is_cons_application(x) || sort_fset::is_insert_application(x))
    {
      x = sort_fset::right(x);
    }
    return sort_fset::is_empty_function_symbol(x);
  }

  /// \brief Returns true if x is a list composed of cons, insert and cinsert applications.
  bool is_fbag_cons_list(data_expression x)
  {
    while (sort_fbag::is_cons_application(x) || sort_fbag::is_insert_application(x) || sort_fbag::is_cinsert_application(x))
    {
      x = sort_fbag::arg3(x);
    }
    return sort_fbag::is_empty_function_symbol(x);
  }

  bool is_numeric_expression(const application& x)
  {
    return    sort_pos::is_pos(x.sort())
           || sort_nat::is_nat(x.sort())
           || sort_int::is_int(x.sort())
           || sort_real::is_real(x.sort());
  }

  bool is_standard_sort(const sort_expression& x)
  {
    return    sort_pos::is_pos(x)
           || sort_bool::is_bool(x)
           || sort_nat::is_nat(x)
           || sort_int::is_int(x)
           || sort_real::is_real(x);
  }

  bool is_fset_true(const data_expression& x)
  {
    return sort_set::is_true_function_function_symbol(sort_set::left(x));
  }

  bool is_fset_false(const data_expression& x)
  {
    return sort_set::is_false_function_function_symbol(sort_set::left(x));
  }

  bool is_fset_lambda(const data_expression& x)
  {
    return is_lambda(sort_set::left(x)) && sort_fset::is_empty_function_symbol(sort_set::right(x));
  }

  bool is_fbag_zero(const data_expression& x)
  {
    return sort_bag::is_zero_function_function_symbol(sort_bag::left(x));
  }

  bool is_fbag_one(const data_expression& x)
  {
    return sort_bag::is_one_function_function_symbol(sort_bag::left(x));
  }

  bool is_fbag_lambda(const data_expression& x)
  {
    return is_lambda(sort_bag::left(x)) && sort_fbag::is_empty_function_symbol(sort_bag::right(x));
  }

  void print_cons_list(data_expression x)
  {
    data_expression_vector arguments;
    while (sort_list::is_cons_application(x))
    {
      arguments.push_back(sort_list::left(x));
      x = sort_list::right(x);
    }
    derived().print("[");
    print_container(arguments, 6);
    derived().print("]");
  }

  void print_snoc_list(data_expression x)
  {
    data_expression_vector arguments;
    while (sort_list::is_snoc_application(x))
    {
      arguments.insert(arguments.begin(), sort_list::right(x));
      x = sort_list::left(x);
    }
    derived().print("[");
    print_container(arguments, 7);
    derived().print("]");
  }

  void print_fset_cons_list(data_expression x)
  {
    data_expression_vector arguments;
    while (sort_fset::is_cons_application(x) || sort_fset::is_insert_application(x))
    {
      arguments.push_back(sort_fset::left(x));
      x = sort_fset::right(x);
    }
    derived().print("{");
    print_container(arguments, 6);
    derived().print("}");
  }

  void print_fbag_zero(const data_expression& x)
  {
    // TODO: check if this is the correct way to handle this case
    data_expression y = sort_bag::right(x);
    if (sort_fbag::is_empty_function_symbol(y))
    {
      derived().print("{:}");
    }
    else if (data::is_variable(y))
    {
      derived().print("@bagfbag(");
      derived().apply(variable(y).name());
      derived().print(")");
    }
    else
    {
      derived().apply(y);
    }
  }

  void print_fbag_one(const data_expression& x)
  {
    sort_expression s = function_sort(sort_bag::left(x).sort()).domain().front(); // the sort of the bag elements
    core::identifier_string name = generate_identifier("x", x);
    variable var(name, s);
    data_expression body = number(sort_nat::nat(), "1");
    if (!sort_fbag::is_empty_function_symbol(sort_bag::right(x)))
    {
      body = sort_nat::swap_zero(body, sort_bag::count(s, var, sort_bag::bag_fbag(s, sort_bag::right(x))));
    }
    derived().print("{ ");
    print_variable(var, true);
    derived().print(" | ");
    derived().apply(body);
    derived().print(" }");
  }

  void print_fbag_lambda(const data_expression& x)
  {
    sort_expression s = function_sort(sort_bag::left(x).sort()).domain().front(); // the sort of the bag elements
    core::identifier_string name = generate_identifier("x", x);
    variable var(name, s);
    data::lambda left(sort_bag::left(x));
    data_expression body = left.body();
    if (!sort_fbag::is_empty_function_symbol(sort_bag::right(x)))
    {
      body = sort_nat::swap_zero(body, sort_bag::count(s, var, sort_bag::bag_fbag(s, sort_bag::right(x))));
    }
    derived().print("{ ");
    print_variables(left.variables(), true, true, false, "", "", ", ");
    derived().print(" | ");
    derived().apply(body);
    derived().print(" }");
  }

  void print_fbag_default(const data_expression& x)
  {
    sort_expression s = function_sort(sort_bag::left(x).sort()).domain().front();
    core::identifier_string name = generate_identifier("x", x);
    variable var(name, s);
    data_expression body = sort_bag::left(x)(var);
    if (!sort_fbag::is_empty_function_symbol(sort_bag::right(x)))
    {
      body = sort_nat::swap_zero(body, sort_bag::count(s, var, sort_bag::bag_fbag(s, sort_bag::right(x))));
    }
    derived().print("{ ");
    print_variable(var, true);
    derived().print(" | ");
    derived().apply(body);
    derived().print(" }");
  }

  void print_fbag_cons_list(data_expression x)
  {
    std::vector<std::pair<data_expression, data_expression> > arguments;
    while (sort_fbag::is_cons_application(x) || sort_fbag::is_insert_application(x) || sort_fbag::is_cinsert_application(x))
    {
      if (sort_fbag::is_cons_application(x))
      {
        arguments.emplace_back(sort_fbag::arg1(x), sort_fbag::arg2(x));
        x = sort_fbag::arg3(x);
      }
      else if (sort_fbag::is_insert_application(x))
      {
#ifdef MCRL2_ENABLE_MACHINENUMBERS
        arguments.emplace_back(sort_fbag::arg1(x), sort_nat::pos2nat(sort_fbag::arg2(x)));
#else
        arguments.emplace_back(sort_fbag::arg1(x), sort_nat::cnat(sort_fbag::arg2(x)));
#endif
        x = sort_fbag::arg3(x);
      }
      else // if (sort_fbag::is_fbagcinsert_application(x))
      {
        arguments.emplace_back(sort_fbag::arg1(x), sort_fbag::arg2(x));
        x = sort_fbag::arg3(x);
      }
    }
    print_list(arguments, "{", "}");
  }

  void print_fset_true(const data_expression& x)
  {
    derived().print("!");
    derived().apply(sort_set::right(x));
  }

  void print_fset_false(const data_expression& x)
  {
    if (sort_fset::is_empty_function_symbol(sort_set::right(x)))
    {
      derived().print("{}");
    }
    else
    {
      derived().apply(sort_set::right(x));
    }
  }

  void print_fset_lambda(const data_expression& x)
  {
    data::lambda left(sort_set::left(x));
    derived().print("{ ");
    print_variables(left.variables(), true, true, false, "", "", ", ");
    derived().print(" | ");
    derived().apply(left.body());
    derived().print(" }");
  }

  void print_fset_set_operation(const data_expression& x, const std::string& op)
  {
    data_expression f = sort_set::arg1(x);
    data_expression g = sort_set::arg2(x);

    // print lhs
    if (sort_set::is_false_function_function_symbol(g))
    {
      derived().apply(sort_set::arg3(x));
    }
    else if (sort_set::is_true_function_function_symbol(g))
    {
      derived().print("!");
      derived().apply(sort_set::arg3(x));
    }
    else if (is_function_sort(sort_set::arg1(x).sort()))
    {
      const sort_expression sort = sort_set::arg1(x).sort();
      const sort_expression& s = atermpp::down_cast<function_sort>(sort).domain().front();
      core::identifier_string name = generate_identifier("x", x);
      variable var(name, s);
      data_expression body = sort_bool::and_(sort_bool::not_(g(var)), sort_set::in(s, var, sort_set::arg3(x)));
      derived().print("{ ");
      print_variable(var, true);
      derived().print(" | ");
      derived().apply(body);
      derived().print(" }");
    }
    else 
    {
      // In this case the term is not well formed, for instance because it contains a "Rewritten@@term" function.
      // We print the residue as an aterm. 
      derived().print(pp(atermpp::aterm(x)));
      return;
    }

    // print operator
    derived().print(op);

    // print rhs
    if (sort_set::is_false_function_function_symbol(f))
    {
      derived().apply(sort_set::arg4(x));
    }
    else if (sort_set::is_true_function_function_symbol(f))
    {
      derived().print("!");
      derived().apply(sort_set::arg4(x));
    }
    else
    {
      sort_expression s = function_sort(sort_set::arg1(x).sort()).domain().front();
      core::identifier_string name = generate_identifier("x", x);
      variable var(name, s);
      data_expression body = sort_bool::and_(sort_bool::not_(f(var)), sort_set::in(s, var, sort_set::arg4(x)));
      derived().print("{ ");
      print_variable(var, true);
      derived().print(" | ");
      derived().apply(body);
      derived().print(" }");
    }
  }

  void print_fset_default(const data_expression& x)
  {
    data_expression right = sort_set::right(x);
    // TODO: check if this is the correct way to handle this case
    if (sort_fset::is_empty_function_symbol(right))
    {
        sort_expression s = function_sort(sort_set::left(x).sort()).domain().front();
        core::identifier_string name = generate_identifier("x", x);
        variable var(name, s);
        data_expression body(sort_set::left(x)(var));
        derived().print("{ ");
        print_variable(var, true);
        derived().print(" | ");
        derived().apply(body);
        derived().print(" }");
    }
    else
    {
        sort_expression s = function_sort(sort_set::left(x).sort()).domain().front();
        core::identifier_string name = generate_identifier("x", x);
        variable var(name, s);
        data_expression lhs(sort_set::left(x)(var));
        data_expression rhs(sort_set::in(s, var, sort_set::set_fset(s, right)));
        data_expression body = not_equal_to(lhs, rhs);
        derived().print("{ ");
        print_variable(var, true);
        derived().print(" | ");
        derived().apply(body);
        derived().print(" }");
    }
  }

  template <typename Abstraction>
  void print_abstraction(const Abstraction& x, const std::string& op)
  {
    derived().enter(x);
    derived().print(op + " ");
    print_variables(x.variables(), true, true, false, "", "", ", ");
    derived().print(". ");
    derived().apply(x.body());
    derived().leave(x);
  }

  void print_function_application(const application& x)
  {
    // Add special handling of list/set/bag enumeration types. This case applies to printing
    // terms after parsing and before type checking.
    if (sort_list::is_list_enumeration_application(x))
    {
      print_list_enumeration(x);
      return;
    }
    else if (sort_set::is_set_enumeration_application(x))
    {
      print_set_enumeration(x);
      return;
    }
    else if (sort_bag::is_bag_enumeration_application(x))
    {
      print_bag_enumeration(x);
      return;
    }

    if (is_infix_operation(x))
    {
      assert(detail::is_untyped(x));
      auto i = x.begin();
      data_expression left = *i++;
      data_expression right = *i;


      print_expression(left, precedence(left)<precedence(x));
      derived().print(" ");
      derived().apply(x.head());
      derived().print(" ");
      print_expression(right, precedence(right)<precedence(x));
      return;
    }

    // print the head
    bool print_parentheses = is_abstraction(x.head());
    if (print_parentheses)
    {
      derived().print("(");
    }
    derived().apply(x.head());
    if (print_parentheses)
    {
      derived().print(")");
    }

    // print the arguments
    print_parentheses = !x.empty();
    if (is_function_symbol(x.head()) && x.size() == 1)
    {
      std::string name(function_symbol(x.head()).name());
      if (name == "!" || name == "#")
      {
        print_parentheses = precedence(*x.begin()) < core::detail::max_precedence;
      }
    }
    if (print_parentheses)
    {
      derived().print("(");
    }
    print_container(x);
    if (print_parentheses)
    {
      derived().print(")");
    }
  }

  // N.B. This is interpreted as the bag element 'x.first: x.second'
  void apply(const std::pair<data_expression, data_expression>& x)
  {
    derived().apply(x.first);
    derived().print(": ");
    derived().apply(x.second);
  }

  // TODO: this code should be generated!
  void apply(const data::container_type& x)
  {
    derived().enter(x);
    if (data::is_list_container(x))
    {
      derived().apply(data::list_container(atermpp::aterm(x)));
    }
    else if (data::is_set_container(x))
    {
      derived().apply(data::set_container(atermpp::aterm(x)));
    }
    else if (data::is_bag_container(x))
    {
      derived().apply(data::bag_container(atermpp::aterm(x)));
    }
    else if (data::is_fset_container(x))
    {
      derived().apply(data::fset_container(atermpp::aterm(x)));
    }
    else if (data::is_fbag_container(x))
    {
      derived().apply(data::fbag_container(atermpp::aterm(x)));
    }
    derived().leave(x);
  }

  void apply(const data::assignment& x)
  {
    derived().enter(x);
    derived().apply(x.lhs());
    derived().print(" = ");
    derived().apply(x.rhs());
    derived().leave(x);
  }

  // variable lists have their own notation
  void apply(const data::variable_list& x)
  {
    derived().enter(x);
    print_list(x, "", "", ", ", false);
    derived().leave(x);
  }

  void apply(const data::untyped_data_parameter& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    print_list(x.arguments(), "(", ")", ", ");
    derived().leave(x);
  }

  void apply(const data::untyped_identifier_assignment& x)
  {
    derived().enter(x);
    derived().apply(x.lhs());
    derived().print("=");
    derived().apply(x.rhs());
    derived().leave(x);
  }

  void apply(const data::untyped_set_or_bag_comprehension_binder& x)
  {
    derived().enter(x);
    derived().leave(x);
  }

  void apply(const data::set_comprehension_binder& x)
  {
    derived().enter(x);
    derived().leave(x);
  }

  void apply(const data::bag_comprehension_binder& x)
  {
    derived().enter(x);
    derived().leave(x);
  }

  void apply(const data::forall_binder& x)
  {
    derived().enter(x);
    derived().leave(x);
  }

  void apply(const data::exists_binder& x)
  {
    derived().enter(x);
    derived().leave(x);
  }

  void apply(const data::lambda_binder& x)
  {
    derived().enter(x);
    derived().leave(x);
  }

  void apply(const data::structured_sort_constructor_argument& x)
  {
    derived().enter(x);
    if (x.name() != core::empty_identifier_string())
    {
      derived().apply(x.name());
      derived().print(": ");
    }
    derived().apply(x.sort());
    derived().leave(x);
  }

  void apply(const data::structured_sort_constructor& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    print_list(x.arguments(), "(", ")", ", ");
    if (x.recogniser() != core::empty_identifier_string())
    {
      derived().print("?");
      derived().apply(x.recogniser());
    }
    derived().leave(x);
  }

  void apply(const data::alias& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    derived().print(" = ");
    derived().apply(x.reference());
    derived().leave(x);
  }

  void apply(const data::list_container& x)
  {
    derived().enter(x);
    derived().print("List");
    derived().leave(x);
  }

  void apply(const data::set_container& x)
  {
    derived().enter(x);
    derived().print("Set");
    derived().leave(x);
  }

  void apply(const data::bag_container& x)
  {
    derived().enter(x);
    derived().print("Bag");
    derived().leave(x);
  }

  void apply(const data::fset_container& x)
  {
    derived().enter(x);
    derived().print("FSet");
    derived().leave(x);
  }

  void apply(const data::fbag_container& x)
  {
    derived().enter(x);
    derived().print("FBag");
    derived().leave(x);
  }

  void apply(const data::basic_sort& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    derived().leave(x);
  }

  void apply(const data::container_sort& x)
  {
    derived().enter(x);
    derived().apply(x.container_name());
    derived().print("(");
    derived().apply(x.element_sort());
    derived().print(")");
    derived().leave(x);
  }

  void apply(const data::structured_sort& x)
  {
    derived().enter(x);
    print_list(x.constructors(), "struct ", "", " | ");
    derived().leave(x);
  }

  void apply(const data::function_sort& x)
  {
    derived().enter(x);
    print_sort_list(x.domain(), "", " -> ", " # ");
    derived().apply(x.codomain());
    derived().leave(x);
  }

  void apply(const data::untyped_sort& x)
  {
    derived().enter(x);
    derived().print("untyped_sort");
    derived().leave(x);
  }

  void apply(const data::untyped_possible_sorts& x)
  {
    derived().enter(x);
    derived().print("@untyped_possible_sorts[");
    derived().apply(x.sorts());
    derived().print("]");
    derived().leave(x);
  }

  void apply(const data::untyped_sort_variable& x)
  {
    derived().enter(x);
    derived().print("@s");
    derived().apply(x.value());
    derived().leave(x);
  }

  void apply(const data::untyped_identifier& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    derived().leave(x);
  }

  void apply(const data::variable& x)
  {
    derived().enter(x);
    derived().apply(x.name());
    derived().leave(x);
  }

  void apply(const data::function_symbol& x)
  {
    derived().enter(x);
    if (sort_nat::is_c0_function_symbol(x))
    {
      derived().print("0");
    }
    else if (sort_pos::is_c1_function_symbol(x))
    {
      derived().print("1");
    }
    else if (sort_fbag::is_empty_function_symbol(x))
    {
      derived().print("{:}");
    }
    else if (sort_fset::is_empty_function_symbol(x))
    {
      derived().print("{}");
    }
    else
    {
      derived().print(x.name());
    }
    derived().leave(x);
  }

  void apply(const data::application& x)
  {
    derived().enter(x);

    //-------------------------------------------------------------------//
    //                            bool
    //-------------------------------------------------------------------//

    if (sort_bool::is_not_application(x))
    {
      print_unary_data_operation(x, "!");
    }
    else if (sort_bool::is_and_application(x))
    {
      print_binary_data_operation(x, " && ");
    }
    else if (sort_bool::is_or_application(x))
    {
      print_binary_data_operation(x, " || ");
    }
    else if (sort_bool::is_implies_application(x))
    {
      print_binary_data_operation(x, " => ");
    }

    //-------------------------------------------------------------------//
    //                            data
    //-------------------------------------------------------------------//

    else if (data::is_equal_to_application(x))
    {
      print_binary_data_operation(x, " == ");
    }
    else if (data::is_not_equal_to_application(x))
    {
      print_binary_data_operation(x, " != ");
    }
    else if (data::is_if_application(x))
    {
      // TODO: is this correct?
      print_function_application(x);
    }
    else if (data::is_less_application(x))
    {
      print_binary_data_operation(x, " < ");
    }
    else if (data::is_less_equal_application(x))
    {
      print_binary_data_operation(x, " <= ");
    }
    else if (data::is_greater_application(x))
    {
      print_binary_data_operation(x, " > ");
    }
    else if (data::is_greater_equal_application(x))
    {
      print_binary_data_operation(x, " >= ");
    }

#ifdef MCRL2_ENABLE_MACHINENUMBERS    
    //-------------------------------------------------------------------//
    //                            machine_word
    //-------------------------------------------------------------------//

    else if (sort_machine_word::is_zero_word_function_symbol(x))
    {
      derived().print(std::to_string(sort_machine_word::detail::zero_word().value()));
    }
    else if (sort_machine_word::is_one_word_function_symbol(x))
    {
      derived().print(std::to_string(sort_machine_word::detail::one_word().value()));
    }
    else if (sort_machine_word::is_max_word_function_symbol(x))
    {
      derived().print(std::to_string(sort_machine_word::detail::max_word().value()));
    }
    else if (sort_machine_word::is_succ_word_application(x))
    {
      derived().print("((");
      derived().apply(sort_machine_word::arg(x));
      derived().print(" + 1) mod ");
      derived().print(max_machine_number_string());
      derived().print(" )");
    }
    else if (sort_machine_word::is_add_word_application(x))
    {
      derived().print("((");
      derived().apply(sort_machine_word::left(x));
      derived().print(" + ");
      derived().apply(sort_machine_word::right(x));
      derived().print(") mod ");
      derived().print(max_machine_number_string());
      derived().print(" )");
    }
    else if (sort_machine_word::is_add_overflow_word_application(x))
    {
      derived().print("((");
      derived().apply(sort_machine_word::left(x));
      derived().print(" + ");
      derived().apply(sort_machine_word::right(x));
      derived().print(") div ");
      derived().print(max_machine_number_string());
      derived().print(" )");
    }
    else if (sort_machine_word::is_times_word_application(x))
    {
      derived().print("((");
      derived().apply(sort_machine_word::left(x));
      derived().print(" * ");
      derived().apply(sort_machine_word::right(x));
      derived().print(") mod ");
      derived().print(max_machine_number_string());
      derived().print(" )");
    }
    else if (sort_machine_word::is_times_overflow_word_application(x))
    {
      derived().print("((");
      derived().apply(sort_machine_word::left(x));
      derived().print(" * ");
      derived().apply(sort_machine_word::right(x));
      derived().print(") div ");
      derived().print(max_machine_number_string());
      derived().print(" )");
    }
    /* else if (sort_machine_word::is_minus_word_application(x))
    {
      derived().print("(if(");
      derived().print_binary_operation(x, " >= ");
      derived().print(", ");
      derived().print_binary_operation(x, " - ");
      derived().print(", ");
      derived().print(max_machine_number_string());
      derived().print(" + ");
      derived().print_binary_operation(x, " - ");
      derived().print(")");
    } 
    else if (sort_machine_word::is_div_word_application(x))
    {
      derived().print_binary_operation(x, " div ");
    }
    else if (sort_machine_word::is_mod_word_application(x))
    {
      derived().print_binary_operation(x, " mod ");
    } */

/* TODO: Handle the following cases. 
 
     @div_doubleword <"div_doubleword">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> -> @word                                                  
     @div_double_doubleword <"div_double_doubleword">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> # @word <"arg4"> -> @word                   
     @div_triple_doubleword <"div_triple_doubleword">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> # @word <"arg4"> # @word <"arg5"> -> @word  
     @mod_doubleword <"mod_doubleword">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> -> @word                                                  
     @sqrt_doubleword <"sqrt_doubleword">: @word <"arg1"> # @word <"arg2"> -> @word                                                                 
     @sqrt_tripleword <"sqrt_tripleword">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> -> @word                                                
     @sqrt_tripleword_overflow <"sqrt_tripleword_overflow">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> -> @word                              
     @sqrt_quadrupleword <"sqrt_quadrupleword">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> # @word <"arg4"> -> @word                         
     @sqrt_quadrupleword_overflow <"sqrt_quadrupleword_overflow">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> # @word <"arg4"> -> @word       
     @pred_word <"pred_word">: @word <"arg"> ->@word                                                                                                
*/


#endif

    //-------------------------------------------------------------------//
    //                            pos
    //-------------------------------------------------------------------//

#ifdef MCRL2_ENABLE_MACHINENUMBERS
    else if (data::sort_pos::is_positive_constant(x))
    {
      derived().print(data::sort_pos::positive_constant_as_string(x));
    }
    else if (sort_pos::is_most_significant_digit_application(x))
    {
      derived().apply(sort_pos::arg(x));
    }
    else if (sort_pos::is_concat_digit_application(x))
    {
      derived().print(max_machine_number_string() + "* (");
      derived().apply(sort_pos::arg1(x));
      derived().print(") + ");
      derived().apply(sort_pos::arg2(x));
    }
#else
    else if (sort_pos::is_cdub_application(x))
    {
      if (data::sort_pos::is_positive_constant(x))
      {
        derived().print(data::sort_pos::positive_constant_as_string(x));
      }
      else
      {
        std::vector<char> number = data::detail::string_to_vector_number("1");
        derived().apply(detail::reconstruct_pos_mult(x, number));
      }
    }
#endif
    // TODO: handle @pospred
    else if (sort_pos::is_plus_application(x))
    {
      print_binary_data_operation(x, " + ");
    }
    else if (sort_pos::is_add_with_carry_application(x))
    {
#ifdef MCRL2_ENABLE_MACHINENUMBERS
      const data_expression& x1 = sort_pos::left(x);
      const data_expression& x2 = sort_pos::right(x);
      derived().apply(sort_pos::plus(x1, x2));
#else
      const data_expression& b = sort_pos::arg1(x);
      const data_expression& x1 = sort_pos::arg2(x);
      const data_expression& x2 = sort_pos::arg3(x);
      if (b == data::sort_bool::true_())
      {
        derived().apply(sort_pos::succ(sort_pos::plus(x1, x2)));
      }
      else if (b == sort_bool::false_())
      {
        derived().apply(sort_pos::plus(x1, x2));
      }
      else
      {
        derived().apply(if_(b, x1, x2));
      }
#endif
    }
    else if (sort_pos::is_times_application(x))
    {
      print_binary_data_operation(x, " * ");
    }
    // TODO: handle @powerlog2

    //-------------------------------------------------------------------//
    //                            natpair
    //-------------------------------------------------------------------//

#ifdef MCRL2_ENABLE_MACHINENUMBERS
    else if (sort_nat::is_first_application(x))
    {
      // TODO: verify if this is the correct way of dealing with first/divmod
      data_expression y = sort_nat::arg(x);
      if (!sort_nat::is_divmod_aux_application(y))
      {
        print_function_application(x);
      }
      else
      {
        print_expression(sort_nat::left(y), true);     // Removed  left_precedence(y) and replaced it with true, which may be undesireable. Happens 3x more below. 
        derived().print(" div ");
        print_expression(sort_nat::right(y), true);
      } 
    }
    else if (sort_nat::is_last_application(x))
    {
      // TODO: verify if this is the correct way of dealing with last/divmod
      data_expression y = sort_nat::arg(x);
      if (!sort_nat::is_divmod_aux_application(y))
      {
        print_function_application(x);
      }
      else
      {
        print_expression(sort_nat::left(y), true);
        derived().print(" mod ");
        print_expression(sort_nat::right(y), true);
      }
    }
#else
    else if (sort_nat::is_first_application(x))
    {
        // TODO: verify if this is the correct way of dealing with first/divmod
        const auto& y = atermpp::down_cast<application>(sort_nat::arg(x));
      if (!sort_nat::is_divmod_application(y))
      {
        print_function_application(x);
      }
      else
      {
        print_binary_data_operation(y, " div ");
      }
    }
    else if (sort_nat::is_last_application(x))
    {
      // TODO: verify if this is the correct way of dealing with last/divmod
      const auto& y = atermpp::down_cast<application>(sort_nat::arg(x));
      if (!sort_nat::is_divmod_application(y))
      {
        print_function_application(x);
      }         
      else      
      {            
        print_binary_data_operation(y, " mod ");
      }                     
    }                       
#endif

    //-------------------------------------------------------------------//
    //                            nat
    //-------------------------------------------------------------------//

#ifdef MCRL2_ENABLE_MACHINENUMBERS
    else if (data::sort_nat::is_natural_constant(x))
    {
      derived().print(data::sort_nat::natural_constant_as_string(x));
    }
    else if (sort_nat::is_most_significant_digit_nat_application(x))
    {
      derived().apply(sort_nat::arg(x));
    }   
    else if (sort_nat::is_concat_digit_application(x))
    { 
      derived().print("(" + max_machine_number_string() + "*");
      derived().apply(sort_nat::arg1(x));
      derived().print(") + ");
      derived().apply(sort_nat::arg2(x));
    } 
#else
    else if (sort_nat::is_cnat_application(x))
    {
      derived().apply(sort_nat::arg(x));
    }
#endif
    else if (sort_nat::is_pos2nat_application(x))
    {
      derived().apply(*x.begin());
    }
    // TODO: handle @dub
    // TODO: handle @dubsucc
    else if (sort_nat::is_plus_application(x))
    {
      print_binary_data_operation(x, " + ");
    }
    // TODO: handle @gtesubtb
    else if (sort_nat::is_times_application(x))
    {
      print_binary_data_operation(x, " * ");
    }
    else if (sort_nat::is_div_application(x))
    {
      print_binary_data_operation(x, sort_nat::left(x), sort_nat::right(x), " div ");
    }
    else if (sort_nat::is_mod_application(x))
    {
      print_binary_data_operation(x, sort_nat::left(x), sort_nat::right(x), " mod ");
    }
    // TODO: handle @monus
    // TODO: handle @swap_zero*
    // TODO: handle @sqrt_nat

    //-------------------------------------------------------------------//
    //                            int
    //-------------------------------------------------------------------//

    else if (sort_int::is_cint_application(x))
    {
      derived().apply(sort_int::arg(x));
    }
    else if (sort_int::is_cneg_application(x))
    {
      derived().apply(sort_int::negate(sort_int::arg(x)));
    }
    else if (sort_int::is_nat2int_application(x))
    {
      derived().apply(*x.begin());
    }
    else if (sort_int::is_pos2int_application(x))
    {
      derived().apply(sort_int::arg(x));
    }
    else if (sort_int::is_negate_application(x))
    {
      print_unary_data_operation(x, "-");
    }
    else if (sort_int::is_plus_application(x))
    {
      print_binary_data_operation(x, " + ");
    }
    else if (sort_int::is_minus_application(x))
    {
      print_binary_data_operation(x, " - ");
    }
    else if (sort_int::is_times_application(x))
    {
      print_binary_data_operation(x, " * ");
    }
    else if (sort_int::is_div_application(x))
    {
      print_binary_data_operation(x, sort_int::left(x), sort_int::right(x), " div ");
    }
    else if (sort_int::is_mod_application(x))
    {
      print_binary_data_operation(x, sort_int::left(x), sort_int::right(x), " mod ");
    }

    //-------------------------------------------------------------------//
    //                            real
    //-------------------------------------------------------------------//

    else if (sort_real::is_creal_application(x))
    {
      data_expression numerator = sort_real::left(x);
      const data_expression& denominator = sort_real::right(x);
      if (detail::is_one(denominator))
      {
        derived().apply(numerator);
      }
      else
      {
        print_binary_data_operation(x, " / ");
      }
    }
    else if (sort_real::is_pos2real_application(x))
    {
      derived().apply(*x.begin());
    }
    else if (sort_real::is_nat2real_application(x))
    {
      derived().apply(*x.begin());
    }
    else if (sort_real::is_int2real_application(x))
    {
      derived().apply(*x.begin());
    }
    else if (sort_real::is_negate_application(x))
    {
      print_unary_data_operation(x, "-");
    }
    else if (sort_real::is_plus_application(x))
    {
      print_binary_data_operation(x, " + ");
    }
    else if (sort_real::is_minus_application(x))
    {
      print_binary_data_operation(x, " - ");
    }
    else if (sort_real::is_times_application(x))
    {
      print_binary_data_operation(x, " * ");
    }
    else if (sort_real::is_divides_application(x))
    {
      print_binary_data_operation(x, " / ");
    }
    else if (sort_real::is_reduce_fraction_application(x))
    {
      derived().apply(sort_real::divides(sort_real::left(x),sort_real::right(x)));
    }
    else if (sort_real::is_reduce_fraction_where_application(x))
    {
      derived().apply(sort_real::plus(sort_real::int2real(sort_real::arg2(x)), sort_real::divides(sort_real::arg3(x), sort_nat::pos2nat(sort_real::arg1(x)))));
    }
    // TODO: handle @redfrachlp

    //-------------------------------------------------------------------//
    //                            list
    //-------------------------------------------------------------------//

    else if (sort_list::is_list_enumeration_application(x))
    {
      print_list_enumeration(x);
    }
    else if (sort_list::is_cons_application(x))
    {
      if (is_cons_list(x))
      {
        print_cons_list(x);
      }
      else
      {
        print_binary_data_operation(x, " |> ");
      }
    }
    else if (sort_list::is_in_application(x))
    {
      print_binary_data_operation(x, " in ");
    }
    else if (sort_list::is_count_application(x))
    {
      derived().print("#");
      print_unary_operand(x, sort_list::arg(x));
    }
    else if (sort_list::is_snoc_application(x))
    {
      if (is_snoc_list(x))
      {
        print_snoc_list(x);
      }
      else
      {
        print_binary_data_operation(x, " <| ");
      }
    }
    else if (sort_list::is_concat_application(x))
    {
      print_binary_data_operation(x, " ++ ");
    }
    else if (sort_list::is_element_at_application(x))
    {
      print_binary_data_operation(x, " . ");
    }

    //-------------------------------------------------------------------//
    //                            set
    //-------------------------------------------------------------------//

    else if (sort_set::is_constructor_application(x))
    {
      if (is_fset_true(x))
      {
        print_fset_true(x);
      }
      else if (is_fset_false(x))
      {
        print_fset_false(x);
      }
      else if (is_fset_lambda(x))
      {
        print_fset_lambda(x);
      }
      else
      {
        print_fset_default(x);
      }
    }
    else if (sort_set::is_set_fset_application(x))
    {
      data_expression y = sort_set::arg(x);
      if (sort_fset::is_empty_function_symbol(y))
      {
        derived().print("{}");
      }
      else if (data::is_variable(y))
      {
        derived().print("@setfset(");
        derived().apply(variable(y).name());
        derived().print(")");
      }
      else
      {
        derived().apply(y);
      }
    }
    else if (sort_set::is_set_comprehension_application(x))
    {
      sort_expression s = function_sort(sort_set::arg(x).sort()).domain().front();
      core::identifier_string name = generate_identifier("x", x);
      variable var(name, s);
      data_expression body(sort_set::arg(x)(var));
      derived().print("{ ");
      print_variable(var, true);
      derived().print(" | ");
      derived().apply(body);
      derived().print(" }");
    }
    else if (sort_set::is_in_application(x))
    {
      print_binary_data_operation(x, " in ");
    }
    else if (sort_set::is_complement_application(x))
    {
      print_unary_data_operation(x, "!");
    }
    else if (sort_set::is_union_application(x))
    {
      print_binary_data_operation(x, " + ");
    }
    else if (sort_set::is_intersection_application(x))
    {
      print_binary_data_operation(x, " * ");
    }
    else if (sort_set::is_difference_application(x))
    {
      print_binary_data_operation(x, " - ");
    }
    else if (sort_set::is_fset_union_application(x))
    {
      print_fset_set_operation(x, " + ");
    }
    else if (sort_set::is_fset_intersection_application(x))
    {
      print_fset_set_operation(x, " * ");
    }

    //-------------------------------------------------------------------//
    //                            fset
    //-------------------------------------------------------------------//

    else if (is_fset_cons_list(x))
    {
      print_fset_cons_list(x);
    }
    else if (sort_fset::is_in_application(x))
    {
      print_binary_data_operation(x, " in ");
    }
    else if (sort_fset::is_difference_application(x))
    {
      derived().apply(sort_fset::left(x));
      derived().print(" - ");
      derived().apply(sort_fset::right(x));
    }
    else if (sort_fset::is_union_application(x))
    {
      derived().apply(sort_fset::left(x));
      derived().print(" + ");
      derived().apply(sort_fset::right(x));
    }
    else if (sort_fset::is_intersection_application(x))
    {
      derived().apply(sort_fset::left(x));
      derived().print(" * ");
      derived().apply(sort_fset::right(x));
    }
    else if (sort_fset::is_count_application(x))
    {
      derived().print("#");
      derived().apply(sort_fset::arg(x));
    }

    //-------------------------------------------------------------------//
    //                            bag
    //-------------------------------------------------------------------//

    else if (sort_bag::is_constructor_application(x))
    {
      if (is_fbag_zero(x))
      {
        print_fbag_zero(x);
      }
      else if (is_fbag_one(x))
      {
        print_fbag_one(x);
      }
      else if (is_fbag_lambda(x))
      {
        print_fbag_lambda(x);
      }
      else
      {
        print_fbag_default(x);
      }
    }
    else if (sort_bag::is_bag_fbag_application(x))
    {
      data_expression y = sort_bag::arg(x);
      if (sort_fbag::is_empty_function_symbol(y))
      {
        derived().print("{:}");
      }
      else if (data::is_variable(y))
      {
        derived().print("@bagfbag(");
        derived().apply(variable(y).name());
        derived().print(")");
      }
      else
      {
        derived().apply(y);
      }
    }
    else if (sort_bag::is_bag_comprehension_application(x))
    {
      sort_expression s = function_sort(sort_bag::arg(x).sort()).domain().front();
      core::identifier_string name = generate_identifier("x", x);
      variable var(name, s);
      data_expression body(sort_bag::arg(x)(var));
      derived().print("{ ");
      print_variable(var, true);
      derived().print(" | ");
      derived().apply(body);
      derived().print(" }");
    }
    else if (sort_bag::is_in_application(x))
    {
      print_binary_data_operation(x, " in ");
    }
    else if (sort_bag::is_union_application(x))
    {
      print_binary_data_operation(x, " + ");
    }
    else if (sort_bag::is_intersection_application(x))
    {
      print_binary_data_operation(x, " * ");
    }
    else if (sort_bag::is_difference_application(x))
    {
      print_binary_data_operation(x, " - ");
    }

    //-------------------------------------------------------------------//
    //                            fbag
    //-------------------------------------------------------------------//

    // cons / insert / cinsert
    else if (is_fbag_cons_list(x))
    {
      print_fbag_cons_list(x);
    }
    else if (sort_fbag::is_in_application(x))
    {
      print_binary_data_operation(x, " in ");
    }
    else if (sort_fbag::is_union_application(x))
    {
      print_binary_data_operation(x, " + ");
    }
    else if (sort_fbag::is_intersection_application(x))
    {
      print_binary_data_operation(x, " * ");
    }

    else if (sort_fbag::is_difference_application(x))
    {
      print_binary_data_operation(x, " - ");
    }
    else if (sort_fbag::is_count_all_application(x))
    {
      derived().print("#");
      derived().apply(sort_fbag::arg(x));
    }

    //-------------------------------------------------------------------//
    //                            function update
    //-------------------------------------------------------------------//
    else if (is_function_update_application(x) || is_function_update_stable_application(x))
    {
      data_expression x1 = data::arg1(x);
      data_expression x2 = data::arg2(x);
      data_expression x3 = data::arg3(x);
      bool print_parentheses = is_abstraction(x1);
      if (print_parentheses)
      {
        derived().print("(");
      }
      derived().apply(x1);
      if (print_parentheses)
      {
        derived().print(")");
      }
      derived().print("[");
      derived().apply(x2);
      derived().print(" -> ");
      derived().apply(x3);
      derived().print("]");
    }

    //-------------------------------------------------------------------//
    //                            abstraction
    //-------------------------------------------------------------------//
    else if (is_abstraction_application(x))
    {
      if (!x.empty()) {
        derived().print("(");
      }
      derived().apply(x.head());
      if (!x.empty())
      {
        derived().print(")(");
      }
      print_container(x);
      if (!x.empty())
      {
        derived().print(")");
      }
    }

    //-------------------------------------------------------------------//
    //                            function application
    //-------------------------------------------------------------------//
    else
    {
      print_function_application(x);
    }
    derived().leave(x);
  }

  void apply(const machine_number& x)
  {     
      derived().print(std::to_string(x.value()));
  }     

  void apply(const data::where_clause& x)
  {
    derived().enter(x);
    derived().apply(x.body());
    derived().print(" whr ");
    const assignment_expression_list& declarations = x.declarations();
    for (assignment_expression_list::const_iterator i = declarations.begin(); i != declarations.end(); ++i)
    {
      if (i != declarations.begin())
      {
        derived().print(", ");
      }
      derived().apply(*i);
    }
    derived().print(" end");
    derived().leave(x);
  }

  void apply(const data::forall& x)
  {
    print_abstraction(x, "forall");
  }

  void apply(const data::exists& x)
  {
    print_abstraction(x, "exists");
  }

  void apply(const data::lambda& x)
  {
    print_abstraction(x, "lambda");
  }

  void apply(const data::data_equation& x)
  {
    derived().enter(x);
    print_condition(x.condition());
    derived().apply(x.lhs());
    derived().print("  =  ");
    derived().apply(x.rhs());
    derived().leave(x);
  }

  // Adds variables v and function symbols f to variable_map and function_symbol_names respectively.
  void update_mappings(const data_equation& eqn,
                       std::vector<variable>& variables,
                       std::map<core::identifier_string, variable>& variable_map,
                       std::set<core::identifier_string>& function_symbol_names
                      )
  {
    for (const function_symbol& f: data::find_function_symbols(eqn))
    {
      function_symbol_names.insert(f.name());
    }
    for (const variable& v: eqn.variables())
    {
      std::pair<std::map<core::identifier_string, variable>::iterator, bool> k = variable_map.insert(std::make_pair(v.name(), v));
      if (k.second) // new variable encountered
      {
        variables.push_back(v);
      }
    }
  }

  bool has_conflict(const data_equation& eqn,
                    const std::map<core::identifier_string, variable>& variable_map
                   )
  {
    for (const variable& v: eqn.variables())
    {
      auto j = variable_map.find(v.name());
      if (j != variable_map.end() && v != j->second)
      {
        return true;
      }
    }
    return false;
  }

  /// \brief Searches in the range of equations [first, last) for the first equation
  /// that conflicts with one of the previous equations. We say that equation e1 conflicts
  /// with another equation e2 if their declared variables contain a variable with the same
  /// name and a different sort, or if a declared variable in e1 has the same name as a
  /// function symbol appearing in equation e2.
  template <typename Iter>
  Iter find_conflicting_equation(Iter first, Iter last, std::vector<variable>& variables)
  {
    std::map<core::identifier_string, variable> variable_map;
    std::set<core::identifier_string> function_symbol_names;
    for (Iter i = first; i != last; ++i)
    {
      if (has_conflict(*i, variable_map))
      {
        return i;
      }
      else
      {
        update_mappings(*i, variables, variable_map, function_symbol_names);
      }
    }
    return last; // no conflict found
  }

  template <typename AliasContainer, typename SortContainer>
  void print_sort_declarations(const AliasContainer& aliases,
                               const SortContainer& sorts,
                               const std::string& opener = "(",
                               const std::string& closer = ")",
                               const std::string& separator = ", "
                              )
  {
    if (aliases.empty() && sorts.empty())
    {
      return;
    }
    bool first_element = true;
    derived().print(opener);

    // print sorts
    for (auto i = sorts.begin(); i != sorts.end(); ++i)
    {
      if (!first_element)
      {
        derived().print(separator);
      }
      derived().apply(*i);
      first_element = false;
    }

    // print aliases
    for (auto i = aliases.begin(); i != aliases.end(); ++i)
    {
      if (!first_element)
      {
        derived().print(separator);
      }
      derived().apply(*i);
      first_element = false;
    }

    derived().print(closer);
  }

  template <typename Container>
  void print_equations(const Container& equations,
                       const data_specification& data_spec,
                       const std::string& opener = "(",
                       const std::string& closer = ")",
                       const std::string& separator = ", "
                      )
  {
    auto first = equations.begin();
    auto last = equations.end();

    Container normalized_equations = equations;
    data::normalize_sorts(normalized_equations, data_spec);

    while (first != last)
    {
      std::vector<variable> variables;
      auto i = find_conflicting_equation(first, last, variables);
      print_variables(variables, true, true, true, "var  ", ";\n", ";\n     ");

      // N.B. We print normalized equations instead of user defined equations.
      // print_list(std::vector<data_equation>(first, i), opener, closer, separator);
      auto first1 = normalized_equations.begin() + (first - equations.begin());
      auto i1 = normalized_equations.begin() + (i - equations.begin());
      print_list(std::vector<data_equation>(first1, i1), opener, closer, separator);

      first = i;
    }
  }

  void apply(const data::data_specification& x)
  {
    derived().enter(x);
    print_sort_declarations(x.user_defined_aliases(), x.user_defined_sorts(), "sort ", ";\n\n", ";\n     ");
    print_sorted_declarations(x.user_defined_constructors(), true, true, false, "cons ",";\n\n", ";\n     ", get_sort_default());
    print_sorted_declarations(x.user_defined_mappings(), true, true, false, "map  ",";\n\n", ";\n     ", get_sort_default());
    print_equations(x.user_defined_equations(), x, "eqn  ", ";\n\n", ";\n     ");
    derived().leave(x);
  }

  // Override, because there are set/bag/setbag comprehension classes that exist after parsing and before type checking.
  void apply(const data::abstraction& x)
  {
    derived().enter(x);
    if (data::is_forall(x))
    {
      derived().apply(atermpp::down_cast<data::forall>(x));
    }
    else if (data::is_exists(x))
    {
      derived().apply(atermpp::down_cast<data::exists>(x));
    }
    else if (data::is_lambda(x))
    {
      derived().apply(atermpp::down_cast<data::lambda>(x));
    }
    else if (data::is_set_comprehension(x))
    {
      print_setbag_comprehension(x);
    }
    else if (data::is_bag_comprehension(x))
    {
      print_setbag_comprehension(x);
    }
    else if (data::is_untyped_set_or_bag_comprehension(x))
    {
      print_setbag_comprehension(x);
    }
    derived().leave(x);
  }
};

} // namespace detail

/// \brief Prints the object x to a stream.
struct stream_printer
{
  template <typename T>
  void operator()(const T& x, std::ostream& out, bool precendence_aware)
  {
    core::detail::apply_printer<data::detail::printer> printer(out, precendence_aware);
    printer.apply(x);
  }
};

/// \brief Returns a string representation of the object x.
template <typename T>
std::string pp(const T& x, bool precendence_aware = true)
{
  std::ostringstream out;
  stream_printer()(x, out, precendence_aware);
  return out.str();
}

} // namespace data

} // namespace mcrl2

#endif
