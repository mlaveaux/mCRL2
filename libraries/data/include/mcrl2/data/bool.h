// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/bool.h
/// \brief The standard sort bool_.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/bool.spec.

#ifndef MCRL2_DATA_BOOL_H
#define MCRL2_DATA_BOOL_H

#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"

#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort bool_.
    namespace sort_bool {

      inline
      const core::identifier_string& bool_name()
      {
        static core::identifier_string bool_name = core::identifier_string("Bool");
        return bool_name;
      }

      /// \brief Constructor for sort expression Bool.
      /// \return Sort expression Bool.
      inline
      const basic_sort& bool_()
      {
        static basic_sort bool_ = basic_sort(bool_name());
        return bool_;
      }

      /// \brief Recogniser for sort expression Bool
      /// \param e A sort expression
      /// \return true iff e == bool_()
      inline
      bool is_bool(const sort_expression& e)
      {
        if (is_basic_sort(e))
        {
          return basic_sort(e) == bool_();
        }
        return false;
      }


      /// \brief Generate identifier true.
      /// \return Identifier true.
      inline
      const core::identifier_string& true_name()
      {
        static core::identifier_string true_name = core::identifier_string("true");
        return true_name;
      }

      /// \brief Constructor for function symbol true.
      
      /// \return Function symbol true_.
      inline
      const function_symbol& true_()
      {
        static function_symbol true_(true_name(), bool_());
        return true_;
      }

      /// \brief Recogniser for function true.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching true.
      inline
      bool is_true_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == true_();
        }
        return false;
      }

      /// \brief Generate identifier false.
      /// \return Identifier false.
      inline
      const core::identifier_string& false_name()
      {
        static core::identifier_string false_name = core::identifier_string("false");
        return false_name;
      }

      /// \brief Constructor for function symbol false.
      
      /// \return Function symbol false_.
      inline
      const function_symbol& false_()
      {
        static function_symbol false_(false_name(), bool_());
        return false_;
      }

      /// \brief Recogniser for function false.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching false.
      inline
      bool is_false_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == false_();
        }
        return false;
      }
      /// \brief Give all system defined constructors for bool_.
      /// \return All system defined constructors for bool_.
      inline
      function_symbol_vector bool_generate_constructors_code()
      {
        function_symbol_vector result;
        result.push_back(sort_bool::true_());
        result.push_back(sort_bool::false_());

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for bool_.
      /// \return All system defined constructors that can be used in an mCRL2 specification for bool_.
      inline
      function_symbol_vector bool_mCRL2_usable_constructors()
      {
        function_symbol_vector result;
        result.push_back(sort_bool::true_());
        result.push_back(sort_bool::false_());

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for bool_.
      /// \return All system defined constructors that are to be implemented in C++ for bool_.
      inline
      implementation_map bool_cpp_implementable_constructors()
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier !.
      /// \return Identifier !.
      inline
      const core::identifier_string& not_name()
      {
        static core::identifier_string not_name = core::identifier_string("!");
        return not_name;
      }

      /// \brief Constructor for function symbol !.
      
      /// \return Function symbol not_.
      inline
      const function_symbol& not_()
      {
        static function_symbol not_(not_name(), make_function_sort_(bool_(), bool_()));
        return not_;
      }

      /// \brief Recogniser for function !.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching !.
      inline
      bool is_not_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == not_();
        }
        return false;
      }

      /// \brief Application of function symbol !.
      
      /// \param arg0 A data expression.
      /// \return Application of ! to a number of arguments.
      inline
      application not_(const data_expression& arg0)
      {
        return sort_bool::not_()(arg0);
      }

      /// \brief Make an application of function symbol !.
      /// \param result The data expression where the ! expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_not_(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_bool::not_(),arg0);
      }

      /// \brief Recogniser for application of !.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol not_ to a
      ///     number of arguments.
      inline
      bool is_not_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_not_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier &&.
      /// \return Identifier &&.
      inline
      const core::identifier_string& and_name()
      {
        static core::identifier_string and_name = core::identifier_string("&&");
        return and_name;
      }

      /// \brief Constructor for function symbol &&.
      
      /// \return Function symbol and_.
      inline
      const function_symbol& and_()
      {
        static function_symbol and_(and_name(), make_function_sort_(bool_(), bool_(), bool_()));
        return and_;
      }

      /// \brief Recogniser for function &&.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching &&.
      inline
      bool is_and_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == and_();
        }
        return false;
      }

      /// \brief Application of function symbol &&.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of && to a number of arguments.
      inline
      application and_(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_bool::and_()(arg0, arg1);
      }

      /// \brief Make an application of function symbol &&.
      /// \param result The data expression where the && expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_and_(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bool::and_(),arg0, arg1);
      }

      /// \brief Recogniser for application of &&.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol and_ to a
      ///     number of arguments.
      inline
      bool is_and_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_and_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier ||.
      /// \return Identifier ||.
      inline
      const core::identifier_string& or_name()
      {
        static core::identifier_string or_name = core::identifier_string("||");
        return or_name;
      }

      /// \brief Constructor for function symbol ||.
      
      /// \return Function symbol or_.
      inline
      const function_symbol& or_()
      {
        static function_symbol or_(or_name(), make_function_sort_(bool_(), bool_(), bool_()));
        return or_;
      }

      /// \brief Recogniser for function ||.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching ||.
      inline
      bool is_or_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == or_();
        }
        return false;
      }

      /// \brief Application of function symbol ||.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of || to a number of arguments.
      inline
      application or_(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_bool::or_()(arg0, arg1);
      }

      /// \brief Make an application of function symbol ||.
      /// \param result The data expression where the || expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_or_(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bool::or_(),arg0, arg1);
      }

      /// \brief Recogniser for application of ||.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol or_ to a
      ///     number of arguments.
      inline
      bool is_or_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_or_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier =>.
      /// \return Identifier =>.
      inline
      const core::identifier_string& implies_name()
      {
        static core::identifier_string implies_name = core::identifier_string("=>");
        return implies_name;
      }

      /// \brief Constructor for function symbol =>.
      
      /// \return Function symbol implies.
      inline
      const function_symbol& implies()
      {
        static function_symbol implies(implies_name(), make_function_sort_(bool_(), bool_(), bool_()));
        return implies;
      }

      /// \brief Recogniser for function =>.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching =>.
      inline
      bool is_implies_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == implies();
        }
        return false;
      }

      /// \brief Application of function symbol =>.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of => to a number of arguments.
      inline
      application implies(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_bool::implies()(arg0, arg1);
      }

      /// \brief Make an application of function symbol =>.
      /// \param result The data expression where the => expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_implies(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bool::implies(),arg0, arg1);
      }

      /// \brief Recogniser for application of =>.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol implies to a
      ///     number of arguments.
      inline
      bool is_implies_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_implies_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for bool_
      /// \return All system defined mappings for bool_
      inline
      function_symbol_vector bool_generate_functions_code()
      {
        function_symbol_vector result;
        result.push_back(sort_bool::not_());
        result.push_back(sort_bool::and_());
        result.push_back(sort_bool::or_());
        result.push_back(sort_bool::implies());
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for bool_
      /// \return All system defined mappings for bool_
      inline
      function_symbol_vector bool_generate_constructors_and_functions_code()
      {
        function_symbol_vector result=bool_generate_functions_code();
        for(const function_symbol& f: bool_generate_constructors_code())
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for bool_
      /// \return All system defined mappings for that can be used in mCRL2 specificationis bool_
      inline
      function_symbol_vector bool_mCRL2_usable_mappings()
      {
        function_symbol_vector result;
        result.push_back(sort_bool::not_());
        result.push_back(sort_bool::and_());
        result.push_back(sort_bool::or_());
        result.push_back(sort_bool::implies());
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for bool_
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for bool_
      inline
      implementation_map bool_cpp_implementable_mappings()
      {
        implementation_map result;
        return result;
      }
      ///\brief Function for projecting out argument.
      ///        arg from an application.
      /// \param e A data expression.
      /// \pre arg is defined for e.
      /// \return The argument of e that corresponds to arg.
      inline
      const data_expression& arg(const data_expression& e)
      {
        assert(is_not_application(e));
        return atermpp::down_cast<application>(e)[0];
      }

      ///\brief Function for projecting out argument.
      ///        left from an application.
      /// \param e A data expression.
      /// \pre left is defined for e.
      /// \return The argument of e that corresponds to left.
      inline
      const data_expression& left(const data_expression& e)
      {
        assert(is_and_application(e) || is_or_application(e) || is_implies_application(e));
        return atermpp::down_cast<application>(e)[0];
      }

      ///\brief Function for projecting out argument.
      ///        right from an application.
      /// \param e A data expression.
      /// \pre right is defined for e.
      /// \return The argument of e that corresponds to right.
      inline
      const data_expression& right(const data_expression& e)
      {
        assert(is_and_application(e) || is_or_application(e) || is_implies_application(e));
        return atermpp::down_cast<application>(e)[1];
      }

      /// \brief Give all system defined equations for bool_
      /// \return All system defined equations for sort bool_
      inline
      data_equation_vector bool_generate_equations_code()
      {
        variable vb("b",bool_());

        data_equation_vector result;
        result.push_back(data_equation(variable_list(), not_(true_()), false_()));
        result.push_back(data_equation(variable_list(), not_(false_()), true_()));
        result.push_back(data_equation(variable_list({vb}), not_(not_(vb)), vb));
        result.push_back(data_equation(variable_list({vb}), and_(vb, true_()), vb));
        result.push_back(data_equation(variable_list({vb}), and_(vb, false_()), false_()));
        result.push_back(data_equation(variable_list({vb}), and_(true_(), vb), vb));
        result.push_back(data_equation(variable_list({vb}), and_(false_(), vb), false_()));
        result.push_back(data_equation(variable_list({vb}), or_(vb, true_()), true_()));
        result.push_back(data_equation(variable_list({vb}), or_(vb, false_()), vb));
        result.push_back(data_equation(variable_list({vb}), or_(true_(), vb), true_()));
        result.push_back(data_equation(variable_list({vb}), or_(false_(), vb), vb));
        result.push_back(data_equation(variable_list({vb}), implies(vb, true_()), true_()));
        result.push_back(data_equation(variable_list({vb}), implies(vb, false_()), not_(vb)));
        result.push_back(data_equation(variable_list({vb}), implies(true_(), vb), vb));
        result.push_back(data_equation(variable_list({vb}), implies(false_(), vb), true_()));
        result.push_back(data_equation(variable_list({vb}), equal_to(true_(), vb), vb));
        result.push_back(data_equation(variable_list({vb}), equal_to(false_(), vb), not_(vb)));
        result.push_back(data_equation(variable_list({vb}), equal_to(vb, true_()), vb));
        result.push_back(data_equation(variable_list({vb}), equal_to(vb, false_()), not_(vb)));
        result.push_back(data_equation(variable_list({vb}), less(false_(), vb), vb));
        result.push_back(data_equation(variable_list({vb}), less(true_(), vb), false_()));
        result.push_back(data_equation(variable_list({vb}), less(vb, false_()), false_()));
        result.push_back(data_equation(variable_list({vb}), less(vb, true_()), not_(vb)));
        result.push_back(data_equation(variable_list({vb}), less_equal(false_(), vb), true_()));
        result.push_back(data_equation(variable_list({vb}), less_equal(true_(), vb), vb));
        result.push_back(data_equation(variable_list({vb}), less_equal(vb, false_()), not_(vb)));
        result.push_back(data_equation(variable_list({vb}), less_equal(vb, true_()), true_()));
        return result;
      }

    } // namespace sort_bool_

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_BOOL_H
