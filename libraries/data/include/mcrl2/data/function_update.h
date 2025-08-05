// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/function_update.h
/// \brief The standard sort function_update.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/function_update.spec.

#ifndef MCRL2_DATA_FUNCTION_UPDATE_H
#define MCRL2_DATA_FUNCTION_UPDATE_H

#include "mcrl2/utilities/exception.h"
#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/bool.h"

#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

      /// \brief Give all system defined constructors for function_update.
      /// \return All system defined constructors for function_update.
      inline
      function_symbol_vector function_update_generate_constructors_code()
      {
        function_symbol_vector result;
        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for function_update.
      /// \return All system defined constructors that can be used in an mCRL2 specification for function_update.
      inline
      function_symbol_vector function_update_mCRL2_usable_constructors()
      {
        function_symbol_vector result;
        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for function_update.
      /// \return All system defined constructors that are to be implemented in C++ for function_update.
      inline
      implementation_map function_update_cpp_implementable_constructors()
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier \@func_update.
      /// \return Identifier \@func_update.
      inline
      const core::identifier_string& function_update_name()
      {
        static core::identifier_string function_update_name = core::identifier_string("@func_update");
        return function_update_name;
      }

      /// \brief Constructor for function symbol \@func_update.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \return Function symbol function_update.
      inline
      function_symbol function_update(const sort_expression& s, const sort_expression& t)
      {
        function_symbol function_update(function_update_name(), make_function_sort_(make_function_sort_(s, t), s, t, make_function_sort_(s, t)));
        return function_update;
      }

      /// \brief Recogniser for function \@func_update.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@func_update.
      inline
      bool is_function_update_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == function_update_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@func_update.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@func_update to a number of arguments.
      inline
      application function_update(const sort_expression& s, const sort_expression& t, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return function_update(s, t)(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@func_update.
      /// \param result The data expression where the \@func_update expression is put.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_function_update(data_expression& result, const sort_expression& s, const sort_expression& t, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, function_update(s, t),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@func_update.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol function_update to a
      ///     number of arguments.
      inline
      bool is_function_update_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_function_update_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@func_update_stable.
      /// \return Identifier \@func_update_stable.
      inline
      const core::identifier_string& function_update_stable_name()
      {
        static core::identifier_string function_update_stable_name = core::identifier_string("@func_update_stable");
        return function_update_stable_name;
      }

      /// \brief Constructor for function symbol \@func_update_stable.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \return Function symbol function_update_stable.
      inline
      function_symbol function_update_stable(const sort_expression& s, const sort_expression& t)
      {
        function_symbol function_update_stable(function_update_stable_name(), make_function_sort_(make_function_sort_(s, t), s, t, make_function_sort_(s, t)));
        return function_update_stable;
      }

      /// \brief Recogniser for function \@func_update_stable.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@func_update_stable.
      inline
      bool is_function_update_stable_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == function_update_stable_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@func_update_stable.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@func_update_stable to a number of arguments.
      inline
      application function_update_stable(const sort_expression& s, const sort_expression& t, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return function_update_stable(s, t)(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@func_update_stable.
      /// \param result The data expression where the \@func_update_stable expression is put.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_function_update_stable(data_expression& result, const sort_expression& s, const sort_expression& t, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, function_update_stable(s, t),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@func_update_stable.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol function_update_stable to a
      ///     number of arguments.
      inline
      bool is_function_update_stable_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_function_update_stable_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@is_not_an_update.
      /// \return Identifier \@is_not_an_update.
      inline
      const core::identifier_string& is_not_a_function_update_name()
      {
        static core::identifier_string is_not_a_function_update_name = core::identifier_string("@is_not_an_update");
        return is_not_a_function_update_name;
      }

      /// \brief Constructor for function symbol \@is_not_an_update.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \return Function symbol is_not_a_function_update.
      inline
      function_symbol is_not_a_function_update(const sort_expression& s, const sort_expression& t)
      {
        function_symbol is_not_a_function_update(is_not_a_function_update_name(), make_function_sort_(make_function_sort_(s, t), sort_bool::bool_()));
        return is_not_a_function_update;
      }

      /// \brief Recogniser for function \@is_not_an_update.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@is_not_an_update.
      inline
      bool is_is_not_a_function_update_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == is_not_a_function_update_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@is_not_an_update.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@is_not_an_update to a number of arguments.
      inline
      application is_not_a_function_update(const sort_expression& s, const sort_expression& t, const data_expression& arg0)
      {
        return is_not_a_function_update(s, t)(arg0);
      }

      /// \brief Make an application of function symbol \@is_not_an_update.
      /// \param result The data expression where the \@is_not_an_update expression is put.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_is_not_a_function_update(data_expression& result, const sort_expression& s, const sort_expression& t, const data_expression& arg0)
      {
        make_application(result, is_not_a_function_update(s, t),arg0);
      }

      /// \brief Recogniser for application of \@is_not_an_update.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol is_not_a_function_update to a
      ///     number of arguments.
      inline
      bool is_is_not_a_function_update_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_is_not_a_function_update_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@is_not_an_update.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@is_not_an_update to a number of arguments.
      inline
      void is_not_a_function_update_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void is_not_a_function_update_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==is_not_a_function_update());
        is_not_a_function_update_manual_implementation(result, a[0]);
      }


      /// \brief Generate identifier \@if_always_else.
      /// \return Identifier \@if_always_else.
      inline
      const core::identifier_string& if_always_else_name()
      {
        static core::identifier_string if_always_else_name = core::identifier_string("@if_always_else");
        return if_always_else_name;
      }

      /// \brief Constructor for function symbol \@if_always_else.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \return Function symbol if_always_else.
      inline
      function_symbol if_always_else(const sort_expression& s, const sort_expression& t)
      {
        function_symbol if_always_else(if_always_else_name(), make_function_sort_(sort_bool::bool_(), make_function_sort_(s, t), make_function_sort_(s, t), make_function_sort_(s, t)));
        return if_always_else;
      }

      /// \brief Recogniser for function \@if_always_else.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@if_always_else.
      inline
      bool is_if_always_else_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == if_always_else_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@if_always_else.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@if_always_else to a number of arguments.
      inline
      application if_always_else(const sort_expression& s, const sort_expression& t, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return if_always_else(s, t)(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@if_always_else.
      /// \param result The data expression where the \@if_always_else expression is put.
      /// \param s A sort expression.
      /// \param t A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_if_always_else(data_expression& result, const sort_expression& s, const sort_expression& t, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, if_always_else(s, t),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@if_always_else.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol if_always_else to a
      ///     number of arguments.
      inline
      bool is_if_always_else_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_if_always_else_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@if_always_else.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return The data expression corresponding to an application of \@if_always_else to a number of arguments.
      inline
      void if_always_else_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void if_always_else_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==if_always_else());
        if_always_else_manual_implementation(result, a[0], a[1], a[2]);
      }

      /// \brief Give all system defined mappings for function_update
      /// \param s A sort expression
      /// \param t A sort expression
      /// \return All system defined mappings for function_update
      inline
      function_symbol_vector function_update_generate_functions_code(const sort_expression& s, const sort_expression& t)
      {
        function_symbol_vector result;
        result.push_back(function_update(s, t));
        result.push_back(function_update_stable(s, t));
        result.push_back(is_not_a_function_update(s, t));
        result.push_back(if_always_else(s, t));
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for function_update
      /// \param s A sort expression
      /// \param t A sort expression
      /// \return All system defined mappings for function_update
      inline
      function_symbol_vector function_update_generate_constructors_and_functions_code(const sort_expression& s, const sort_expression& t)
      {
        function_symbol_vector result=function_update_generate_functions_code(s, t);
        for(const function_symbol& f: function_update_generate_constructors_code())
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for function_update
      /// \param s A sort expression
      /// \param t A sort expression
      /// \return All system defined mappings for that can be used in mCRL2 specificationis function_update
      inline
      function_symbol_vector function_update_mCRL2_usable_mappings(const sort_expression& s, const sort_expression& t)
      {
        function_symbol_vector result;
        result.push_back(function_update(s, t));
        result.push_back(function_update_stable(s, t));
        result.push_back(is_not_a_function_update(s, t));
        result.push_back(if_always_else(s, t));
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for function_update
      /// \param s A sort expression
      /// \param t A sort expression
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for function_update
      inline
      implementation_map function_update_cpp_implementable_mappings(const sort_expression& s, const sort_expression& t)
      {
        implementation_map result;
        result[is_not_a_function_update(s, t)]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(is_not_a_function_update_application,"is_not_a_function_update_manual_implementation");
        result[if_always_else(s, t)]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(if_always_else_application,"if_always_else_manual_implementation");
        return result;
      }
      ///\brief Function for projecting out argument.
      ///        arg1 from an application.
      /// \param e A data expression.
      /// \pre arg1 is defined for e.
      /// \return The argument of e that corresponds to arg1.
      inline
      const data_expression& arg1(const data_expression& e)
      {
        assert(is_function_update_application(e) || is_function_update_stable_application(e) || is_is_not_a_function_update_application(e) || is_if_always_else_application(e));
        return atermpp::down_cast<application>(e)[0];
      }

      ///\brief Function for projecting out argument.
      ///        arg2 from an application.
      /// \param e A data expression.
      /// \pre arg2 is defined for e.
      /// \return The argument of e that corresponds to arg2.
      inline
      const data_expression& arg2(const data_expression& e)
      {
        assert(is_function_update_application(e) || is_function_update_stable_application(e) || is_if_always_else_application(e));
        return atermpp::down_cast<application>(e)[1];
      }

      ///\brief Function for projecting out argument.
      ///        arg3 from an application.
      /// \param e A data expression.
      /// \pre arg3 is defined for e.
      /// \return The argument of e that corresponds to arg3.
      inline
      const data_expression& arg3(const data_expression& e)
      {
        assert(is_function_update_application(e) || is_function_update_stable_application(e) || is_if_always_else_application(e));
        return atermpp::down_cast<application>(e)[2];
      }

      /// \brief Give all system defined equations for function_update
      /// \param s A sort expression
      /// \param t A sort expression
      /// \return All system defined equations for sort function_update
      inline
      data_equation_vector function_update_generate_equations_code(const sort_expression& s, const sort_expression& t)
      {
        variable vx("x",s);
        variable vy("y",s);
        variable vv("v",t);
        variable vw("w",t);
        variable vf("f",make_function_sort_(s, t));

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vf, vv, vx}), is_not_a_function_update(s, t, vf), function_update(s, t, vf, vx, vv), if_always_else(s, t, equal_to(vf(vx), vv), vf, function_update_stable(s, t, vf, vx, vv))));
        result.push_back(data_equation(variable_list({vf, vv, vw, vx}), function_update(s, t, function_update_stable(s, t, vf, vx, vw), vx, vv), if_always_else(s, t, equal_to(vf(vx), vv), vf, function_update_stable(s, t, vf, vx, vv))));
        result.push_back(data_equation(variable_list({vf, vv, vw, vx, vy}), less(vy, vx), function_update(s, t, function_update_stable(s, t, vf, vy, vw), vx, vv), function_update_stable(s, t, function_update(s, t, vf, vx, vv), vy, vw)));
        result.push_back(data_equation(variable_list({vf, vv, vw, vx, vy}), less(vx, vy), function_update(s, t, function_update_stable(s, t, vf, vy, vw), vx, vv), if_always_else(s, t, equal_to(vf(vx), vv), function_update_stable(s, t, vf, vy, vw), function_update_stable(s, t, function_update_stable(s, t, vf, vy, vw), vx, vv))));
        result.push_back(data_equation(variable_list({vf, vv, vx, vy}), not_equal_to(vx, vy), function_update_stable(s, t, vf, vx, vv)(vy), vf(vy)));
        result.push_back(data_equation(variable_list({vf, vv, vx}), function_update_stable(s, t, vf, vx, vv)(vx), vv));
        result.push_back(data_equation(variable_list({vf, vv, vx, vy}), function_update(s, t, vf, vx, vv)(vy), if_(equal_to(vx, vy), vv, vf(vy))));
        return result;
      }

  } // namespace data

} // namespace mcrl2

#include "mcrl2/data/detail/function_update.h" // This file contains the manual implementations of rewrite functions.
#endif // MCRL2_DATA_FUNCTION_UPDATE_H
