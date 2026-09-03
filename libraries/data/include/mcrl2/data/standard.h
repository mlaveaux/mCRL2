// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/standard.h
/// \brief Standard functions that are available for all sorts.

#ifndef MCRL2_DATA_STANDARD_H
#define MCRL2_DATA_STANDARD_H

#include "mcrl2/core/detail/construction_utility.h"
#include "mcrl2/data/abstraction.h"
#include "mcrl2/data/data_equation.h"

namespace mcrl2::data
{

// predeclare
namespace sort_bool
{
const basic_sort& bool_();
const function_symbol& false_();
const function_symbol& true_();
application and_(const data_expression&,const data_expression&);
application not_(const data_expression&);
application implies(const data_expression&,const data_expression&);
bool is_bool(const sort_expression&);
} // namespace sort_bool

/// \cond INTERNAL_DOCS
namespace detail
{

/// \ Component to facilitate code generation
template < typename Derived >
struct symbol : public core::detail::singleton_identifier< Derived >
{
  static bool is_symbol(const core::identifier_string& e)
  {
    return e == core::detail::singleton_identifier< Derived >::instance();
  }

  static bool is_application(const data_expression& e)
  {
    return data::is_application(e) ? is_application(atermpp::down_cast<application>(e)) : false;
  }

  static bool is_application(const application& e)
  {
    return is_function_symbol(e.head());
  }

  static bool is_function_symbol(const data_expression& e)
  {
    return data::is_function_symbol(e) ? is_function_symbol(atermpp::down_cast<function_symbol>(e)) : false;
  }

  static bool is_function_symbol(const function_symbol& e)
  {
    return is_symbol(e.name());
  }
};

struct equal_symbol : public symbol< equal_symbol >
{
  static char const* initialise()
  {
    return "==";
  }
};
struct not_equal_symbol : public symbol< not_equal_symbol >
{
  static char const* initialise()
  {
    return "!=";
  }
};
struct if_symbol : public symbol< if_symbol >
{
  static char const* initialise()
  {
    return "if";
  }
};
struct less_symbol : public symbol< less_symbol >
{
  static char const* initialise()
  {
    return "<";
  }
};
struct less_equal_symbol : public symbol< less_equal_symbol >
{
  static char const* initialise()
  {
    return "<=";
  }
};
struct greater_symbol : public symbol< greater_symbol >
{
  static char const* initialise()
  {
    return ">";
  }
};
struct greater_equal_symbol : public symbol< greater_equal_symbol >
{
  static char const* initialise()
  {
    return ">=";
  }
};
struct less_total_symbol : public symbol< less_total_symbol >
{
  static char const* initialise()
  {
    return "less_total";
  }
};
} // namespace detail
/// \endcond

/// \brief Constructor for function symbol ==
/// \param[in] s A sort expression
/// \return function symbol equal_to
inline function_symbol equal_to(const sort_expression& s)
{
  return function_symbol(detail::equal_symbol::instance(), make_function_sort_(s, s, sort_bool::bool_()));
}

/// \brief Recogniser for function ==
/// \param[in] e A data expression
/// \return true iff e is the function symbol matching ==
template < typename DataExpression >
inline bool is_equal_to_function_symbol(const DataExpression& e)
{
  return detail::equal_symbol::is_function_symbol(e);
}

/// \brief Application of function symbol ==
/// \param[in] arg0 A data expression
/// \param[in] arg1 A data expression
/// \return Application of == to a number of arguments
inline application equal_to(const data_expression& arg0, const data_expression& arg1)
{
  assert(arg0.sort() == arg1.sort());
  return equal_to(arg0.sort())(arg0, arg1);
}

/// \brief Recogniser for application of ==
/// \param[in] e A data expression
/// \return true iff e is an application of function symbol equal_to to a
///     number of arguments
template < typename DataExpression >
inline bool is_equal_to_application(const DataExpression& e)
{
  return detail::equal_symbol::is_application(e);
}

/// \brief Constructor for function symbol !=
/// \param[in] s A sort expression
/// \return function symbol not_equal_to
inline function_symbol not_equal_to(const sort_expression& s)
{
  return function_symbol(detail::not_equal_symbol::instance(), make_function_sort_(s, s, sort_bool::bool_()));
}

/// \brief Recogniser for function !=
/// \param[in] e A data expression
/// \return true iff e is the function symbol matching !=
template < typename DataExpression >
inline bool is_not_equal_to_function_symbol(const DataExpression& e)
{
  return detail::not_equal_symbol::is_function_symbol(e);
}

/// \brief Application of function symbol !=
/// \param[in] arg0 A data expression
/// \param[in] arg1 A data expression
/// \return Application of != to a number of arguments
inline application not_equal_to(const data_expression& arg0, const data_expression& arg1)
{
  assert(arg0.sort() == arg1.sort());
  return not_equal_to(arg0.sort())(arg0, arg1);
}

/// \brief Recogniser for application of !=
/// \param[in] e A data expression
/// \return true iff e is an application of function symbol not_equal_to to a
///     number of arguments
template < typename DataExpression >
inline bool is_not_equal_to_application(const DataExpression& e)
{
  return detail::not_equal_symbol::is_application(e);
}

/// \brief Constructor for function symbol if
/// \param[in] s A sort expression
/// \return function symbol if_
inline function_symbol if_(const sort_expression& s)
{
  return function_symbol(detail::if_symbol::instance(), make_function_sort_(sort_bool::bool_(), s, s, s));
}

/// \brief Recogniser for function if
/// \param[in] e A data expression
/// \return true iff e is the function symbol matching if_
template < typename DataExpression >
inline bool is_if_function_symbol(const DataExpression& e)
{
  return detail::if_symbol::is_function_symbol(e);
}

/// \brief Application of function symbol if
/// \param[in] arg0 A data expression
/// \param[in] arg1 A data expression
/// \param[in] arg2 A data expression
/// \return Application of if to a number of arguments
inline application if_(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
{
  assert(sort_bool::is_bool(arg0.sort()));
  assert(arg1.sort() == arg2.sort());

  return if_(arg1.sort())(arg0, arg1, arg2);
}

/// \brief Recogniser for application of if
/// \param[in] e A data expression
/// \return true iff e is an application of function symbol if_ to a
///     number of arguments
template < typename DataExpression >
inline bool is_if_application(const DataExpression& e)
{
  return detail::if_symbol::is_application(e);
}

/// \brief Constructor for function symbol <
/// \param[in] s A sort expression
/// \return function symbol less
inline function_symbol less(const sort_expression& s)
{
  return function_symbol(detail::less_symbol::instance(), make_function_sort_(s, s, sort_bool::bool_()));
}

/// \brief Recogniser for function <
/// \param[in] e A data expression
/// \return true iff e is the function symbol matching <
template < typename DataExpression >
inline bool is_less_function_symbol(const DataExpression& e)
{
  return detail::less_symbol::is_function_symbol(e);
}

/// \brief Application of function symbol <
/// \param[in] arg0 A data expression
/// \param[in] arg1 A data expression
/// \return Application of < to a number of arguments
inline application less(const data_expression& arg0, const data_expression& arg1)
{
  assert(arg0.sort() == arg1.sort());
  return less(arg0.sort())(arg0, arg1);
}

/// \brief Recogniser for application of <
/// \param[in] e A data expression
/// \return true iff e is an application of function symbol less to a
///     number of arguments
template < typename DataExpression >
inline bool is_less_application(const DataExpression& e)
{
  return detail::less_symbol::is_application(e);
}

/// \brief Constructor for function symbol <=
/// \param[in] s A sort expression
/// \return function symbol less_equal
inline function_symbol less_equal(const sort_expression& s)
{
  return function_symbol(detail::less_equal_symbol::instance(), make_function_sort_(s, s, sort_bool::bool_()));
}

/// \brief Recogniser for function <=
/// \param[in] e A data expression
/// \return true iff e is the function symbol matching <=
template < typename DataExpression >
inline bool is_less_equal_function_symbol(const DataExpression& e)
{
  return detail::less_equal_symbol::is_function_symbol(e);
}

/// \brief Application of function symbol <=
/// \param[in] arg0 A data expression
/// \param[in] arg1 A data expression
/// \return Application of <= to a number of arguments
inline application less_equal(const data_expression& arg0, const data_expression& arg1)
{
  assert(arg0.sort() == arg1.sort());
  return less_equal(arg0.sort())(arg0, arg1);
}

/// \brief Recogniser for application of <=
/// \param[in] e A data expression
/// \return true iff e is an application of function symbol less_equal to a
///     number of arguments
template < typename DataExpression >
inline bool is_less_equal_application(const DataExpression& e)
{
  return detail::less_equal_symbol::is_application(e);
}

/// \brief Constructor for function symbol >
/// \param[in] s A sort expression
/// \return function symbol greater
inline function_symbol greater(const sort_expression& s)
{
  return function_symbol(detail::greater_symbol::instance(), make_function_sort_(s, s, sort_bool::bool_()));
}

/// \brief Recogniser for function >
/// \param[in] e A data expression
/// \return true iff e is the function symbol matching >
template < typename DataExpression >
inline bool is_greater_function_symbol(const DataExpression& e)
{
  return detail::greater_symbol::is_function_symbol(e);
}

/// \brief Application of function symbol >
/// \param[in] arg0 A data expression
/// \param[in] arg1 A data expression
/// \return Application of > to a number of arguments
inline application greater(const data_expression& arg0, const data_expression& arg1)
{
  assert(arg0.sort() == arg1.sort());
  return greater(arg0.sort())(arg0, arg1);
}

/// \brief Recogniser for application of >
/// \param[in] e A data expression
/// \return true iff e is an application of function symbol greater to a
///     number of arguments
template < typename DataExpression >
inline bool is_greater_application(const DataExpression& e)
{
  return detail::greater_symbol::is_application(e);
}

/// \brief Constructor for function symbol >=
/// \param[in] s A sort expression
/// \return function symbol greater_equal
inline function_symbol greater_equal(const sort_expression& s)
{
  return function_symbol(detail::greater_equal_symbol::instance(), make_function_sort_(s, s, sort_bool::bool_()));
}

/// \brief Recogniser for function >=
/// \param[in] e A data expression
/// \return true iff e is the function symbol matching >=
template < typename DataExpression >
inline bool is_greater_equal_function_symbol(const DataExpression& e)
{
  return detail::greater_equal_symbol::is_function_symbol(e);
}

/// \brief Application of function symbol >=
/// \param[in] arg0 A data expression
/// \param[in] arg1 A data expression
/// \return Application of >= to a number of arguments
inline application greater_equal(const data_expression& arg0, const data_expression& arg1)
{
  assert(arg0.sort() == arg1.sort());
  return greater_equal(arg0.sort())(arg0, arg1);
}

/// \brief Recogniser for application of >=
/// \param[in] e A data expression
/// \return true iff e is an application of function symbol greater_equal to a
///     number of arguments
template < typename DataExpression >
inline bool is_greater_equal_application(const DataExpression& e)
{
  return detail::greater_equal_symbol::is_application(e);
}

/// \brief Constructor for function symbol less_total
/// \param[in] s A sort expression
/// \return function symbol less_total
inline function_symbol less_total(const sort_expression& s)
{
  return function_symbol(detail::less_total_symbol::instance(), make_function_sort_(s, s, sort_bool::bool_()));
}

/// \brief Recogniser for function less_total
/// \param[in] e A data expression
/// \return true iff e is the function symbol matching less_total
template < typename DataExpression >
inline bool is_less_total_function_symbol(const DataExpression& e)
{
  return detail::less_total_symbol::is_function_symbol(e);
}

/// \brief Application of function symbol less_total
/// \param[in] arg0 A data expression
/// \param[in] arg1 A data expression
/// \return Application of less_total to a number of arguments
inline application less_total(const data_expression& arg0, const data_expression& arg1)
{
  assert(arg0.sort() == arg1.sort());
  return less_total(arg0.sort())(arg0, arg1);
}

/// \brief Recogniser for application of less_total
/// \param[in] e A data expression
/// \return true iff e is an application of function symbol less_total to a
///     number of arguments
template < typename DataExpression >
inline bool is_less_total_application(const DataExpression& e)
{
  return detail::less_total_symbol::is_application(e);
}

/// \brief Give all standard system defined functions for sort s
/// \param[in] s A sort expression
/// \return All standard system defined functions for sort s
inline function_symbol_vector standard_generate_functions_code(const sort_expression& s)
{
  function_symbol_vector result;
  result.push_back(equal_to(s));
  result.push_back(not_equal_to(s));
  result.push_back(if_(s));
  result.push_back(less(s));
  result.push_back(less_equal(s));
  result.push_back(greater_equal(s));
  result.push_back(greater(s));
  result.push_back(less_total(s));

  return result;
}

/// \brief Give all standard system defined equations for sort s
/// \param[in] s A sort expression
/// \return All standard system defined equations for sort s
inline data_equation_vector standard_generate_equations_code(const sort_expression& s)
{
  data_equation_vector result;
  variable b("b", sort_bool::bool_());
  variable x("x", s);
  variable y("y", s);
  result.emplace_back(variable_list({x}), equal_to(x, x), sort_bool::true_());
  result.emplace_back(variable_list({x, y}), not_equal_to(x, y), sort_bool::not_(equal_to(x, y)));
  result.emplace_back(variable_list({x, y}), if_(sort_bool::true_(), x, y), x);
  result.emplace_back(variable_list({x, y}), if_(sort_bool::false_(), x, y), y);
  result.emplace_back(variable_list({b, x}), if_(b, x, x), x);
  result.emplace_back(variable_list({x}), less(x, x), sort_bool::false_());
  result.emplace_back(variable_list({x}), less_equal(x, x), sort_bool::true_());
  result.emplace_back(variable_list({x, y}), greater_equal(x, y), less_equal(y, x));
  result.emplace_back(variable_list({x, y}), greater(x, y), less(y, x));
  result.emplace_back(variable_list({x}), less_total(x, x), sort_bool::false_());

  // For a function sort, add the equation f==g iff forall x.f(x)==g(x). This equation is not in the Specification and Analysis of Communicating Systems of 2014.
  if (is_function_sort(s))
  {
    const function_sort& fs = atermpp::down_cast<function_sort>(s);
    variable_vector xvars;
    variable_vector yvars;
    std::size_t index=0;
    for(const sort_expression& sort: fs.domain())
    {
      std::stringstream xs;
      xs << "x" << index;
      ++index;
      variable x(xs.str(),sort);
      xvars.push_back(x);
    }
    variable f("f",s);
    variable g("g",s);
    variable_list xvar_list=variable_list(xvars.begin(),xvars.end());
    result.emplace_back(variable_list({f, g}) + xvar_list,
        equal_to(f, g),
        abstraction(forall_binder(),
            xvar_list,
            equal_to(application(f, xvars.begin(), xvars.end()), application(g, xvars.begin(), xvars.end()))));

    // less_total(f,g) holds iff there is a point x, smaller (w.r.t. less_total on the domain) than every
    // point where f and g differ, such that f(x) is less_total than g(x). This defines a total order on
    // functions given a total order on the domain and codomain (see the domain/codomain totality proof).
    if (xvars.size() == 1)
    {
      const variable& x = xvars.front();
      variable y("y", x.sort());
      result.emplace_back(variable_list({f, g}),
          less_total(f, g),
          abstraction(exists_binder(),
              variable_list({x}),
              sort_bool::and_(
                  abstraction(forall_binder(),
                      variable_list({y}),
                      sort_bool::implies(less_total(y, x), equal_to(application(f, y), application(g, y)))),
                  less_total(application(f, x), application(g, x)))));
    }
  }

  return result;
}

} // namespace mcrl2::data


#endif // MCRL2_DATA_STANDARD_H
