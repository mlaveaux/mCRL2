// Author(s): Jeroen Keiren, Jeroen van der Wulp
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/structured_sort.h
/// \brief The class structured_sort.

#ifndef MCRL2_DATA_STRUCTURED_SORT_H
#define MCRL2_DATA_STRUCTURED_SORT_H

#include "mcrl2/data/set_identifier_generator.h"
#include "mcrl2/data/standard_numbers_utility.h"
#include "mcrl2/data/structured_sort_constructor.h"

namespace mcrl2
{

namespace data
{

/// \cond INTERNAL_DOCS
// declare for friendship
namespace sort_fset
{
function_symbol_vector fset_generate_constructors_code(const sort_expression&);
function_symbol_vector fset_generate_functions_code(const sort_expression&);
data_equation_vector fset_generate_equations_code(const sort_expression&);
} // namespace sort_fset

// declare for friendship
namespace sort_fbag
{
function_symbol_vector fbag_generate_constructors_code(const sort_expression&);
function_symbol_vector fbag_generate_functions_code(const sort_expression&);
data_equation_vector fbag_generate_equations_code(const sort_expression&);
} // namespace sort_fbag
/// \endcond

/// \brief structured sort.
///
/// A structured sort is a sort with the following structure:
///  struct c1(pr1,1:S1,1, ..., pr1,k1:S1,k1)?is_c1
///       | ...
///       | cn(prn,1:Sn,1, ..., prn,kn:Sn,kn)?is_cn
/// in this:
/// * c1, ..., cn are called constructors.
/// * pri,j are the projection functions (or constructor arguments).
/// * is_ci are the recognisers.
//--- start generated class structured_sort ---//
/// \\brief A structured sort
class structured_sort: public sort_expression
{
  public:
    /// \\brief Default constructor X3.
    structured_sort()
      : sort_expression(core::detail::default_values::SortStruct)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit structured_sort(const atermpp::aterm& term)
      : sort_expression(term)
    {
      assert(core::detail::check_term_SortStruct(*this));
    }

    /// \\brief Constructor Z14.
    explicit structured_sort(const structured_sort_constructor_list& constructors)
      : sort_expression(atermpp::aterm(core::detail::function_symbol_SortStruct(), constructors))
    {}

    /// \\brief Constructor Z2.
    template <typename Container>
    structured_sort(const Container& constructors, typename atermpp::enable_if_container<Container, structured_sort_constructor>::type* = nullptr)
      : sort_expression(atermpp::aterm(core::detail::function_symbol_SortStruct(), structured_sort_constructor_list(constructors.begin(), constructors.end())))
    {}

    /// Move semantics
    structured_sort(const structured_sort&) noexcept = default;
    structured_sort(structured_sort&&) noexcept = default;
    structured_sort& operator=(const structured_sort&) noexcept = default;
    structured_sort& operator=(structured_sort&&) noexcept = default;

    const structured_sort_constructor_list& constructors() const
    {
      return atermpp::down_cast<structured_sort_constructor_list>((*this)[0]);
    }
//--- start user section structured_sort ---//

    friend class data_specification;
    friend class sort_specification;

  private:

    friend function_symbol_vector sort_fset::fset_generate_constructors_code(const sort_expression&);
    friend function_symbol_vector sort_fset::fset_generate_functions_code(const sort_expression&);
    friend data_equation_vector sort_fset::fset_generate_equations_code(const sort_expression&);
    friend function_symbol_vector sort_fbag::fbag_generate_constructors_code(const sort_expression&);
    friend function_symbol_vector sort_fbag::fbag_generate_functions_code(const sort_expression&);
    friend data_equation_vector sort_fbag::fbag_generate_equations_code(const sort_expression&);

    static bool has_recogniser(structured_sort_constructor const& s)
    {
      return !s.recogniser().empty();
    }

    function_symbol to_pos_function(const sort_expression& s) const
    {
      function_symbol to_pos_function_("@to_pos", make_function_sort_(s, sort_pos::pos()));
      return to_pos_function_;
    }

    function_symbol equal_arguments_function(const sort_expression& s) const
    {
      function_symbol equal_arguments_function_("@equal_arguments", make_function_sort_(s, s, sort_bool::bool_()));
      return equal_arguments_function_;
    }

    function_symbol smaller_arguments_function(const sort_expression& s) const
    {
      function_symbol smaller_arguments_function_("@less_arguments", make_function_sort_(s, s, sort_bool::bool_()));
      return smaller_arguments_function_;
    }

    function_symbol smaller_equal_arguments_function(const sort_expression& s) const
    {
      function_symbol smaller_equal_arguments_function_("@less_equal_arguments", make_function_sort_(s, s, sort_bool::bool_()));
      return smaller_equal_arguments_function_;
    }

    function_symbol_vector constructor_functions(const sort_expression& s) const
    {
      function_symbol_vector result;
      for (const auto & i : constructors())
      {
        result.push_back(i.constructor_function(s));
      }
      return result;
    }

    function_symbol_vector comparison_functions(const sort_expression& s) const
    {
      function_symbol_vector result;
      result.push_back(to_pos_function(s));
      result.push_back(equal_arguments_function(s));
      result.push_back(smaller_arguments_function(s));
      result.push_back(smaller_equal_arguments_function(s));
      return result;
    }

    function_symbol_vector projection_functions(const sort_expression& s) const
    {
      function_symbol_vector result;
      for (const structured_sort_constructor& i: constructors())
      {
        function_symbol_vector projections(i.projection_functions(s));

        for (function_symbol_vector::const_iterator j = projections.begin(); j != projections.end(); ++j)
        {
          result.push_back(*j);
        }
      }
      return result;
    }

    function_symbol_vector recogniser_functions(const sort_expression& s) const
    {
      function_symbol_vector result;
      for (const structured_sort_constructor& i: constructors())
      {
        if (i.recogniser() != core::empty_identifier_string())
        {
          result.push_back(i.recogniser_function(s));
        }
      }
      return result;
    }


    data_equation_vector comparison_equations(const sort_expression& s) const
    {
      data_equation_vector result;

      // First add an equation "equal_arguments(v,v)=true", as it is sometimes useful.
      variable v("v",s);
      result.push_back(data_equation({ v }, sort_bool::true_(),application(equal_arguments_function(s),v,v),sort_bool::true_()));

      // give every constructor an index.
      std::size_t index = 1;
      for(structured_sort_constructor_list::const_iterator i = constructors().begin(); i != constructors().end(); ++i, ++index)
      {
        // constructor does not take arguments
        if (i->arguments().empty())
        {
          result.push_back(data_equation(to_pos_function(s)(i->constructor_function(s)), sort_pos::pos(index)));
          result.push_back(data_equation(equal_arguments_function(s)(i->constructor_function(s),i->constructor_function(s)), sort_bool::true_()));
          result.push_back(data_equation(smaller_arguments_function(s)(i->constructor_function(s),i->constructor_function(s)), sort_bool::false_()));
          result.push_back(data_equation(smaller_equal_arguments_function(s)(i->constructor_function(s),i->constructor_function(s)), sort_bool::true_()));
        }
        else
        {
          set_identifier_generator generator;

          std::vector< variable > variables1;
          std::vector< variable > variables2;
          for (const structured_sort_constructor_argument& k: i->arguments())
          {
            variables1.push_back(variable(generator("v"), k.sort()));
            variables2.push_back(variable(generator("w"), k.sort()));
          }
          application instance1(i->constructor_function(s), variables1);
          application instance2(i->constructor_function(s), variables2);

          result.push_back(data_equation(variables1, sort_bool::true_(), to_pos_function(s)(instance1), sort_pos::pos(index)));

          // constructors are the same, generate right hand sides of equal_arguments_function, etc.
          variable_vector::const_reverse_iterator end(variables1.rend());
          variable_vector::const_reverse_iterator k(variables1.rbegin());
          variable_vector::const_reverse_iterator l(variables2.rbegin());

          data_expression right_equal         = equal_to(*k, *l);
          data_expression right_smaller       = less(*k, *l);
          data_expression right_smaller_equal = less_equal(*k, *l);

          for (++l, ++k; k != end; ++l, ++k)
          {
            // Constructors have one or more arguments:
            // - rhs for c(x0,...,xn) == c(y0,..,yn):
            //     x0 == y0 && ... && xn == yn
            right_equal         = sort_bool::and_(equal_to(*k, *l), right_equal);
            // - rhs for c(x0,...,xn) < c(y0,..,yn):
            //     x0 < y0                                                     , when n = 0
            //     x0 < y0 || (x0 == y0 && x1 < y1)                            , when n = 1
            //     x0 < y0 || (x0 == y0 && (x1 < y1 || (x1 == y1 && x2 < y2))) , when n = 2
            //     etcetera
            right_smaller       = sort_bool::or_(less(*k, *l),
                                                 sort_bool::and_(equal_to(*k, *l), right_smaller));
            // - rhs for c(x0,...,xn) <= c(y0,..,yn):
            //     x0 <= y0                                                    , when n = 0
            //     x0 < y0 || (x0 == y0 && x1 <= y1)                           , when n = 1
            //     x0 < y0 || (x0 == y0 && (x1 < y1 || (x1 == y1 && x2 <= y2))), when n = 2
            //     etcetera
            right_smaller_equal = sort_bool::or_(less(*k, *l),
                                                 sort_bool::and_(equal_to(*k, *l), right_smaller_equal));
          }

          application left_equal = application(equal_arguments_function(s), instance1, instance2);
          application left_smaller = application(smaller_arguments_function(s), instance1, instance2);
          application left_smaller_equal = application(smaller_equal_arguments_function(s), instance1, instance2);
          variables1.insert(variables1.end(),variables2.begin(),variables2.end());
          result.push_back(data_equation(variables1, sort_bool::true_(),left_equal, right_equal));
          result.push_back(data_equation(variables1, sort_bool::true_(),left_smaller, right_smaller));
          result.push_back(data_equation(variables1, sort_bool::true_(),left_smaller_equal, right_smaller_equal));
        }
      }
      return result;
    }

    data_equation_vector constructor_equations(const sort_expression& s) const
    {
      data_equation_vector result;

      variable x("x", s);
      variable y("y", s);
      variable_list xy = { x, y };
      application to_pos_x = application(to_pos_function(s), x);
      application to_pos_y = application(to_pos_function(s), y);
      application equal_arguments_xy         = application(equal_arguments_function(s), x, y);
      application smaller_arguments_xy       = application(smaller_arguments_function(s), x, y);
      application smaller_equal_arguments_xy = application(smaller_equal_arguments_function(s), x, y);
      result.push_back(data_equation(xy, equal_to(to_pos_x, to_pos_y),     equal_to(x,y), equal_arguments_xy));
      result.push_back(data_equation(xy, not_equal_to(to_pos_x, to_pos_y), equal_to(x,y), sort_bool::false_()));
      result.push_back(data_equation(xy, less(to_pos_x, to_pos_y),         less(x,y), sort_bool::true_()));
      result.push_back(data_equation(xy, equal_to(to_pos_x, to_pos_y),     less(x,y), smaller_arguments_xy));
      result.push_back(data_equation(xy, greater(to_pos_x, to_pos_y),      less(x,y), sort_bool::false_()));
      result.push_back(data_equation(xy, less(to_pos_x, to_pos_y),         less_equal(x,y), sort_bool::true_()));
      result.push_back(data_equation(xy, equal_to(to_pos_x, to_pos_y),     less_equal(x,y), smaller_equal_arguments_xy));
      result.push_back(data_equation(xy, greater(to_pos_x, to_pos_y),      less_equal(x,y), sort_bool::false_()));
      // (JK) The following encoding of the equations would be more desirable; however,
      // rewriting fails if we use this.
/*
      result.push_back(data_equation(xy, sort_bool::true_(), equal_to(x,y), sort_bool::and_(equal_to(to_pos_x, to_pos_y), equal_arguments_xy)));
      result.push_back(data_equation(xy, sort_bool::true_(), less(x,y), sort_bool::or_(less(to_pos_x, to_pos_y), sort_bool::and_(equal_to(to_pos_x, to_pos_y), smaller_arguments_xy))));
      result.push_back(data_equation(xy, sort_bool::true_(), less_equal(x,y), sort_bool::or_(less(to_pos_x, to_pos_y), sort_bool::and_(equal_to(to_pos_x, to_pos_y), smaller_equal_arguments_xy))));
*/
      return result;
    }

    data_equation_vector projection_equations(const sort_expression& s) const
    {
      data_equation_vector result;

      for (const auto & i : constructors())
      {
        if (!i.arguments().empty())
        {
          structured_sort_constructor_argument_list arguments(i.arguments());

          set_identifier_generator generator;

          std::vector< variable > variables;

          // Create variables for equation
          for (const auto & argument : arguments)
          {
            variables.push_back(variable(generator("v"), argument.sort()));
          }

          std::vector< variable >::const_iterator v = variables.begin();

          for (structured_sort_constructor_argument_list::const_iterator j(arguments.begin()); j != arguments.end(); ++j, ++v)
          {
            if (j->name() != core::empty_identifier_string())
            {
              application lhs(function_symbol(j->name(), make_function_sort_(s, j->sort()))
                              (application(i.constructor_function(s), variables)));

              result.push_back(data_equation(variables, lhs, *v));
            }
          }
        }
      }

      return result;
    }

    data_equation_vector recogniser_equations(const sort_expression& s) const
    {
      data_equation_vector result;

      structured_sort_constructor_list cl(constructors());

      for (structured_sort_constructor_list::const_iterator i = cl.begin(); i != cl.end(); ++i)
      {
        for (const auto & j : cl)
        {
          if (j.recogniser() != core::empty_identifier_string())
          {
            data_expression right = (*i == j) ? sort_bool::true_() : sort_bool::false_();

            if (i->arguments().empty())
            {
              result.push_back(data_equation(j.recogniser_function(s)(i->constructor_function(s)), right));
            }
            else
            {

              set_identifier_generator generator;

              structured_sort_constructor_argument_list arguments(i->arguments());

              // Create variables for equation
              variable_vector variables;

              for (const auto & argument : arguments)
              {
                variables.push_back(variable(generator("v"), argument.sort()));
              }

              result.push_back(data_equation(variables, j.recogniser_function(s)(
                                               application(i->constructor_function(s), variables)), right));
            }
          }
        }
      }

      return result;
    }

public:
    /// \brief Returns the constructor functions of this sort, such that the
    ///        result can be used by the rewriter
    function_symbol_vector constructor_functions() const
    {
      return constructor_functions(*this);
    }

    /// \brief Returns the additional functions of this sort, used to implement
    ///        its comparison operators
    function_symbol_vector comparison_functions() const
    {
      return comparison_functions(*this);
    }

    /// \brief Returns the projection functions of this sort, such that the
    ///        result can be used by the rewriter
    function_symbol_vector projection_functions() const
    {
      return projection_functions(*this);
    }


    /// \brief Returns the recogniser functions of this sort, such that the
    ///        result can be used by the rewriter
    function_symbol_vector recogniser_functions() const
    {
      return recogniser_functions(*this);
    }

    /// \brief Returns the equations for ==, < and <= for this sort, such that the
    ///        result can be used by the rewriter internally
    data_equation_vector constructor_equations() const
    {
      return constructor_equations(*this);
    }

    /// \brief Returns the equations for the functions used to implement comparison
    ///        operators on this sort.
    ///        Needed to do proper rewriting
    data_equation_vector comparison_equations() const
    {
      return comparison_equations(*this);
    }

    /// \brief Generate equations for the projection functions of this sort
    /// \return A vector of equations for the projection functions of this sort.
    data_equation_vector projection_equations() const
    {
      return projection_equations(*this);
    }

    /// \brief Generate equations for the recognisers of this sort, assuming
    ///        that this sort is referred to with s.
    /// \return A vector of equations for the recognisers of this sort.
    data_equation_vector recogniser_equations() const
    {
      return recogniser_equations(*this);
    }
//--- end user section structured_sort ---//
};

/// \\brief Make_structured_sort constructs a new term into a given address.
/// \\ \param t The reference into which the new structured_sort is constructed. 
template <class... ARGUMENTS>
inline void make_structured_sort(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_SortStruct(), args...);
}

/// \\brief list of structured_sorts
typedef atermpp::term_list<structured_sort> structured_sort_list;

/// \\brief vector of structured_sorts
typedef std::vector<structured_sort>    structured_sort_vector;

// prototype declaration
std::string pp(const structured_sort& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const structured_sort& x)
{
  return out << data::pp(x);
}

/// \\brief swap overload
inline void swap(structured_sort& t1, structured_sort& t2)
{
  t1.swap(t2);
}
//--- end generated class structured_sort ---//

} // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_STRUCTURED_SORT_H

