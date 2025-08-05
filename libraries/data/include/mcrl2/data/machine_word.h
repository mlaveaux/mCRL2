// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/machine_word.h
/// \brief The standard sort machine_word.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/machine_word.spec.

#ifndef MCRL2_DATA_MACHINE_WORD_H
#define MCRL2_DATA_MACHINE_WORD_H

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

    /// \brief Namespace for system defined sort machine_word.
    namespace sort_machine_word {

      inline
      const core::identifier_string& machine_word_name()
      {
        static core::identifier_string machine_word_name = core::identifier_string("@word");
        return machine_word_name;
      }

      /// \brief Constructor for sort expression \@word.
      /// \return Sort expression \@word.
      inline
      const basic_sort& machine_word()
      {
        static basic_sort machine_word = basic_sort(machine_word_name());
        return machine_word;
      }

      /// \brief Recogniser for sort expression \@word
      /// \param e A sort expression
      /// \return true iff e == machine_word()
      inline
      bool is_machine_word(const sort_expression& e)
      {
        if (is_basic_sort(e))
        {
          return basic_sort(e) == machine_word();
        }
        return false;
      }


      /// \brief Generate identifier \@zero_word.
      /// \return Identifier \@zero_word.
      inline
      const core::identifier_string& zero_word_name()
      {
        static core::identifier_string zero_word_name = core::identifier_string("@zero_word");
        return zero_word_name;
      }

      /// \brief Constructor for function symbol \@zero_word.
      
      /// \return Function symbol zero_word.
      inline
      const function_symbol& zero_word()
      {
        static function_symbol zero_word(zero_word_name(), machine_word());
        return zero_word;
      }

      /// \brief Recogniser for function \@zero_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@zero_word.
      inline
      bool is_zero_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == zero_word();
        }
        return false;
      }

      /// \brief The data expression of an application of the constant symbol \@zero_word.
      /// \details This function is to be implemented manually.       /// \return The data expression corresponding to an application of \@zero_word to a number of arguments.
      inline
      void zero_word_manual_implementation(data_expression& result);

      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void zero_word_application(data_expression& result, [[maybe_unused]] const data_expression& a)
      {
        assert(is_function_symbol(a));
        // assert(a==zero_word());
        zero_word_manual_implementation(result);
      }


      /// \brief Generate identifier \@succ_word.
      /// \return Identifier \@succ_word.
      inline
      const core::identifier_string& succ_word_name()
      {
        static core::identifier_string succ_word_name = core::identifier_string("@succ_word");
        return succ_word_name;
      }

      /// \brief Constructor for function symbol \@succ_word.
      
      /// \return Function symbol succ_word.
      inline
      const function_symbol& succ_word()
      {
        static function_symbol succ_word(succ_word_name(), make_function_sort_(machine_word(), machine_word()));
        return succ_word;
      }

      /// \brief Recogniser for function \@succ_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@succ_word.
      inline
      bool is_succ_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == succ_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@succ_word.
      
      /// \param arg0 A data expression.
      /// \return Application of \@succ_word to a number of arguments.
      inline
      application succ_word(const data_expression& arg0)
      {
        return sort_machine_word::succ_word()(arg0);
      }

      /// \brief Make an application of function symbol \@succ_word.
      /// \param result The data expression where the \@succ_word expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_succ_word(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_machine_word::succ_word(),arg0);
      }

      /// \brief Recogniser for application of \@succ_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol succ_word to a
      ///     number of arguments.
      inline
      bool is_succ_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_succ_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@succ_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@succ_word to a number of arguments.
      inline
      void succ_word_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void succ_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==succ_word());
        succ_word_manual_implementation(result, a[0]);
      }

      /// \brief Give all system defined constructors for machine_word.
      /// \return All system defined constructors for machine_word.
      inline
      function_symbol_vector machine_word_generate_constructors_code()
      {
        function_symbol_vector result;
        result.push_back(sort_machine_word::zero_word());
        result.push_back(sort_machine_word::succ_word());

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for machine_word.
      /// \return All system defined constructors that can be used in an mCRL2 specification for machine_word.
      inline
      function_symbol_vector machine_word_mCRL2_usable_constructors()
      {
        function_symbol_vector result;
        result.push_back(sort_machine_word::zero_word());
        result.push_back(sort_machine_word::succ_word());

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for machine_word.
      /// \return All system defined constructors that are to be implemented in C++ for machine_word.
      inline
      implementation_map machine_word_cpp_implementable_constructors()
      {
        implementation_map result;
        result[sort_machine_word::zero_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::zero_word_application,"sort_machine_word::zero_word_manual_implementation");
        result[sort_machine_word::succ_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::succ_word_application,"sort_machine_word::succ_word_manual_implementation");

        return result;
      }

      /// \brief Generate identifier \@one_word.
      /// \return Identifier \@one_word.
      inline
      const core::identifier_string& one_word_name()
      {
        static core::identifier_string one_word_name = core::identifier_string("@one_word");
        return one_word_name;
      }

      /// \brief Constructor for function symbol \@one_word.
      
      /// \return Function symbol one_word.
      inline
      const function_symbol& one_word()
      {
        static function_symbol one_word(one_word_name(), machine_word());
        return one_word;
      }

      /// \brief Recogniser for function \@one_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@one_word.
      inline
      bool is_one_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == one_word();
        }
        return false;
      }

      /// \brief The data expression of an application of the constant symbol \@one_word.
      /// \details This function is to be implemented manually.       /// \return The data expression corresponding to an application of \@one_word to a number of arguments.
      inline
      void one_word_manual_implementation(data_expression& result);

      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void one_word_application(data_expression& result, [[maybe_unused]] const data_expression& a)
      {
        assert(is_function_symbol(a));
        // assert(a==one_word());
        one_word_manual_implementation(result);
      }


      /// \brief Generate identifier \@two_word.
      /// \return Identifier \@two_word.
      inline
      const core::identifier_string& two_word_name()
      {
        static core::identifier_string two_word_name = core::identifier_string("@two_word");
        return two_word_name;
      }

      /// \brief Constructor for function symbol \@two_word.
      
      /// \return Function symbol two_word.
      inline
      const function_symbol& two_word()
      {
        static function_symbol two_word(two_word_name(), machine_word());
        return two_word;
      }

      /// \brief Recogniser for function \@two_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@two_word.
      inline
      bool is_two_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == two_word();
        }
        return false;
      }

      /// \brief The data expression of an application of the constant symbol \@two_word.
      /// \details This function is to be implemented manually.       /// \return The data expression corresponding to an application of \@two_word to a number of arguments.
      inline
      void two_word_manual_implementation(data_expression& result);

      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void two_word_application(data_expression& result, [[maybe_unused]] const data_expression& a)
      {
        assert(is_function_symbol(a));
        // assert(a==two_word());
        two_word_manual_implementation(result);
      }


      /// \brief Generate identifier \@three_word.
      /// \return Identifier \@three_word.
      inline
      const core::identifier_string& three_word_name()
      {
        static core::identifier_string three_word_name = core::identifier_string("@three_word");
        return three_word_name;
      }

      /// \brief Constructor for function symbol \@three_word.
      
      /// \return Function symbol three_word.
      inline
      const function_symbol& three_word()
      {
        static function_symbol three_word(three_word_name(), machine_word());
        return three_word;
      }

      /// \brief Recogniser for function \@three_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@three_word.
      inline
      bool is_three_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == three_word();
        }
        return false;
      }

      /// \brief The data expression of an application of the constant symbol \@three_word.
      /// \details This function is to be implemented manually.       /// \return The data expression corresponding to an application of \@three_word to a number of arguments.
      inline
      void three_word_manual_implementation(data_expression& result);

      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void three_word_application(data_expression& result, [[maybe_unused]] const data_expression& a)
      {
        assert(is_function_symbol(a));
        // assert(a==three_word());
        three_word_manual_implementation(result);
      }


      /// \brief Generate identifier \@four_word.
      /// \return Identifier \@four_word.
      inline
      const core::identifier_string& four_word_name()
      {
        static core::identifier_string four_word_name = core::identifier_string("@four_word");
        return four_word_name;
      }

      /// \brief Constructor for function symbol \@four_word.
      
      /// \return Function symbol four_word.
      inline
      const function_symbol& four_word()
      {
        static function_symbol four_word(four_word_name(), machine_word());
        return four_word;
      }

      /// \brief Recogniser for function \@four_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@four_word.
      inline
      bool is_four_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == four_word();
        }
        return false;
      }

      /// \brief The data expression of an application of the constant symbol \@four_word.
      /// \details This function is to be implemented manually.       /// \return The data expression corresponding to an application of \@four_word to a number of arguments.
      inline
      void four_word_manual_implementation(data_expression& result);

      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void four_word_application(data_expression& result, [[maybe_unused]] const data_expression& a)
      {
        assert(is_function_symbol(a));
        // assert(a==four_word());
        four_word_manual_implementation(result);
      }


      /// \brief Generate identifier \@max_word.
      /// \return Identifier \@max_word.
      inline
      const core::identifier_string& max_word_name()
      {
        static core::identifier_string max_word_name = core::identifier_string("@max_word");
        return max_word_name;
      }

      /// \brief Constructor for function symbol \@max_word.
      
      /// \return Function symbol max_word.
      inline
      const function_symbol& max_word()
      {
        static function_symbol max_word(max_word_name(), machine_word());
        return max_word;
      }

      /// \brief Recogniser for function \@max_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@max_word.
      inline
      bool is_max_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == max_word();
        }
        return false;
      }

      /// \brief The data expression of an application of the constant symbol \@max_word.
      /// \details This function is to be implemented manually.       /// \return The data expression corresponding to an application of \@max_word to a number of arguments.
      inline
      void max_word_manual_implementation(data_expression& result);

      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void max_word_application(data_expression& result, [[maybe_unused]] const data_expression& a)
      {
        assert(is_function_symbol(a));
        // assert(a==max_word());
        max_word_manual_implementation(result);
      }


      /// \brief Generate identifier \@equals_zero_word.
      /// \return Identifier \@equals_zero_word.
      inline
      const core::identifier_string& equals_zero_word_name()
      {
        static core::identifier_string equals_zero_word_name = core::identifier_string("@equals_zero_word");
        return equals_zero_word_name;
      }

      /// \brief Constructor for function symbol \@equals_zero_word.
      
      /// \return Function symbol equals_zero_word.
      inline
      const function_symbol& equals_zero_word()
      {
        static function_symbol equals_zero_word(equals_zero_word_name(), make_function_sort_(machine_word(), sort_bool::bool_()));
        return equals_zero_word;
      }

      /// \brief Recogniser for function \@equals_zero_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@equals_zero_word.
      inline
      bool is_equals_zero_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == equals_zero_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@equals_zero_word.
      
      /// \param arg0 A data expression.
      /// \return Application of \@equals_zero_word to a number of arguments.
      inline
      application equals_zero_word(const data_expression& arg0)
      {
        return sort_machine_word::equals_zero_word()(arg0);
      }

      /// \brief Make an application of function symbol \@equals_zero_word.
      /// \param result The data expression where the \@equals_zero_word expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_equals_zero_word(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_machine_word::equals_zero_word(),arg0);
      }

      /// \brief Recogniser for application of \@equals_zero_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol equals_zero_word to a
      ///     number of arguments.
      inline
      bool is_equals_zero_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_equals_zero_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@equals_zero_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@equals_zero_word to a number of arguments.
      inline
      void equals_zero_word_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void equals_zero_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==equals_zero_word());
        equals_zero_word_manual_implementation(result, a[0]);
      }


      /// \brief Generate identifier \@not_equals_zero_word.
      /// \return Identifier \@not_equals_zero_word.
      inline
      const core::identifier_string& not_equals_zero_word_name()
      {
        static core::identifier_string not_equals_zero_word_name = core::identifier_string("@not_equals_zero_word");
        return not_equals_zero_word_name;
      }

      /// \brief Constructor for function symbol \@not_equals_zero_word.
      
      /// \return Function symbol not_equals_zero_word.
      inline
      const function_symbol& not_equals_zero_word()
      {
        static function_symbol not_equals_zero_word(not_equals_zero_word_name(), make_function_sort_(machine_word(), sort_bool::bool_()));
        return not_equals_zero_word;
      }

      /// \brief Recogniser for function \@not_equals_zero_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@not_equals_zero_word.
      inline
      bool is_not_equals_zero_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == not_equals_zero_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@not_equals_zero_word.
      
      /// \param arg0 A data expression.
      /// \return Application of \@not_equals_zero_word to a number of arguments.
      inline
      application not_equals_zero_word(const data_expression& arg0)
      {
        return sort_machine_word::not_equals_zero_word()(arg0);
      }

      /// \brief Make an application of function symbol \@not_equals_zero_word.
      /// \param result The data expression where the \@not_equals_zero_word expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_not_equals_zero_word(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_machine_word::not_equals_zero_word(),arg0);
      }

      /// \brief Recogniser for application of \@not_equals_zero_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol not_equals_zero_word to a
      ///     number of arguments.
      inline
      bool is_not_equals_zero_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_not_equals_zero_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@not_equals_zero_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@not_equals_zero_word to a number of arguments.
      inline
      void not_equals_zero_word_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void not_equals_zero_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==not_equals_zero_word());
        not_equals_zero_word_manual_implementation(result, a[0]);
      }


      /// \brief Generate identifier \@equals_one_word.
      /// \return Identifier \@equals_one_word.
      inline
      const core::identifier_string& equals_one_word_name()
      {
        static core::identifier_string equals_one_word_name = core::identifier_string("@equals_one_word");
        return equals_one_word_name;
      }

      /// \brief Constructor for function symbol \@equals_one_word.
      
      /// \return Function symbol equals_one_word.
      inline
      const function_symbol& equals_one_word()
      {
        static function_symbol equals_one_word(equals_one_word_name(), make_function_sort_(machine_word(), sort_bool::bool_()));
        return equals_one_word;
      }

      /// \brief Recogniser for function \@equals_one_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@equals_one_word.
      inline
      bool is_equals_one_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == equals_one_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@equals_one_word.
      
      /// \param arg0 A data expression.
      /// \return Application of \@equals_one_word to a number of arguments.
      inline
      application equals_one_word(const data_expression& arg0)
      {
        return sort_machine_word::equals_one_word()(arg0);
      }

      /// \brief Make an application of function symbol \@equals_one_word.
      /// \param result The data expression where the \@equals_one_word expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_equals_one_word(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_machine_word::equals_one_word(),arg0);
      }

      /// \brief Recogniser for application of \@equals_one_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol equals_one_word to a
      ///     number of arguments.
      inline
      bool is_equals_one_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_equals_one_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@equals_one_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@equals_one_word to a number of arguments.
      inline
      void equals_one_word_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void equals_one_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==equals_one_word());
        equals_one_word_manual_implementation(result, a[0]);
      }


      /// \brief Generate identifier \@equals_max_word.
      /// \return Identifier \@equals_max_word.
      inline
      const core::identifier_string& equals_max_word_name()
      {
        static core::identifier_string equals_max_word_name = core::identifier_string("@equals_max_word");
        return equals_max_word_name;
      }

      /// \brief Constructor for function symbol \@equals_max_word.
      
      /// \return Function symbol equals_max_word.
      inline
      const function_symbol& equals_max_word()
      {
        static function_symbol equals_max_word(equals_max_word_name(), make_function_sort_(machine_word(), sort_bool::bool_()));
        return equals_max_word;
      }

      /// \brief Recogniser for function \@equals_max_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@equals_max_word.
      inline
      bool is_equals_max_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == equals_max_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@equals_max_word.
      
      /// \param arg0 A data expression.
      /// \return Application of \@equals_max_word to a number of arguments.
      inline
      application equals_max_word(const data_expression& arg0)
      {
        return sort_machine_word::equals_max_word()(arg0);
      }

      /// \brief Make an application of function symbol \@equals_max_word.
      /// \param result The data expression where the \@equals_max_word expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_equals_max_word(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_machine_word::equals_max_word(),arg0);
      }

      /// \brief Recogniser for application of \@equals_max_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol equals_max_word to a
      ///     number of arguments.
      inline
      bool is_equals_max_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_equals_max_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@equals_max_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@equals_max_word to a number of arguments.
      inline
      void equals_max_word_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void equals_max_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==equals_max_word());
        equals_max_word_manual_implementation(result, a[0]);
      }


      /// \brief Generate identifier \@add_word.
      /// \return Identifier \@add_word.
      inline
      const core::identifier_string& add_word_name()
      {
        static core::identifier_string add_word_name = core::identifier_string("@add_word");
        return add_word_name;
      }

      /// \brief Constructor for function symbol \@add_word.
      
      /// \return Function symbol add_word.
      inline
      const function_symbol& add_word()
      {
        static function_symbol add_word(add_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return add_word;
      }

      /// \brief Recogniser for function \@add_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@add_word.
      inline
      bool is_add_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == add_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@add_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@add_word to a number of arguments.
      inline
      application add_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::add_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@add_word.
      /// \param result The data expression where the \@add_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_add_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::add_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@add_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol add_word to a
      ///     number of arguments.
      inline
      bool is_add_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_add_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@add_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@add_word to a number of arguments.
      inline
      void add_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void add_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==add_word());
        add_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@add_with_carry_word.
      /// \return Identifier \@add_with_carry_word.
      inline
      const core::identifier_string& add_with_carry_word_name()
      {
        static core::identifier_string add_with_carry_word_name = core::identifier_string("@add_with_carry_word");
        return add_with_carry_word_name;
      }

      /// \brief Constructor for function symbol \@add_with_carry_word.
      
      /// \return Function symbol add_with_carry_word.
      inline
      const function_symbol& add_with_carry_word()
      {
        static function_symbol add_with_carry_word(add_with_carry_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return add_with_carry_word;
      }

      /// \brief Recogniser for function \@add_with_carry_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@add_with_carry_word.
      inline
      bool is_add_with_carry_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == add_with_carry_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@add_with_carry_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@add_with_carry_word to a number of arguments.
      inline
      application add_with_carry_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::add_with_carry_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@add_with_carry_word.
      /// \param result The data expression where the \@add_with_carry_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_add_with_carry_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::add_with_carry_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@add_with_carry_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol add_with_carry_word to a
      ///     number of arguments.
      inline
      bool is_add_with_carry_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_add_with_carry_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@add_with_carry_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@add_with_carry_word to a number of arguments.
      inline
      void add_with_carry_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void add_with_carry_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==add_with_carry_word());
        add_with_carry_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@add_overflow_word.
      /// \return Identifier \@add_overflow_word.
      inline
      const core::identifier_string& add_overflow_word_name()
      {
        static core::identifier_string add_overflow_word_name = core::identifier_string("@add_overflow_word");
        return add_overflow_word_name;
      }

      /// \brief Constructor for function symbol \@add_overflow_word.
      
      /// \return Function symbol add_overflow_word.
      inline
      const function_symbol& add_overflow_word()
      {
        static function_symbol add_overflow_word(add_overflow_word_name(), make_function_sort_(machine_word(), machine_word(), sort_bool::bool_()));
        return add_overflow_word;
      }

      /// \brief Recogniser for function \@add_overflow_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@add_overflow_word.
      inline
      bool is_add_overflow_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == add_overflow_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@add_overflow_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@add_overflow_word to a number of arguments.
      inline
      application add_overflow_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::add_overflow_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@add_overflow_word.
      /// \param result The data expression where the \@add_overflow_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_add_overflow_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::add_overflow_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@add_overflow_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol add_overflow_word to a
      ///     number of arguments.
      inline
      bool is_add_overflow_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_add_overflow_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@add_overflow_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@add_overflow_word to a number of arguments.
      inline
      void add_overflow_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void add_overflow_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==add_overflow_word());
        add_overflow_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@add_with_carry_overflow_word.
      /// \return Identifier \@add_with_carry_overflow_word.
      inline
      const core::identifier_string& add_with_carry_overflow_word_name()
      {
        static core::identifier_string add_with_carry_overflow_word_name = core::identifier_string("@add_with_carry_overflow_word");
        return add_with_carry_overflow_word_name;
      }

      /// \brief Constructor for function symbol \@add_with_carry_overflow_word.
      
      /// \return Function symbol add_with_carry_overflow_word.
      inline
      const function_symbol& add_with_carry_overflow_word()
      {
        static function_symbol add_with_carry_overflow_word(add_with_carry_overflow_word_name(), make_function_sort_(machine_word(), machine_word(), sort_bool::bool_()));
        return add_with_carry_overflow_word;
      }

      /// \brief Recogniser for function \@add_with_carry_overflow_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@add_with_carry_overflow_word.
      inline
      bool is_add_with_carry_overflow_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == add_with_carry_overflow_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@add_with_carry_overflow_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@add_with_carry_overflow_word to a number of arguments.
      inline
      application add_with_carry_overflow_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::add_with_carry_overflow_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@add_with_carry_overflow_word.
      /// \param result The data expression where the \@add_with_carry_overflow_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_add_with_carry_overflow_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::add_with_carry_overflow_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@add_with_carry_overflow_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol add_with_carry_overflow_word to a
      ///     number of arguments.
      inline
      bool is_add_with_carry_overflow_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_add_with_carry_overflow_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@add_with_carry_overflow_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@add_with_carry_overflow_word to a number of arguments.
      inline
      void add_with_carry_overflow_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void add_with_carry_overflow_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==add_with_carry_overflow_word());
        add_with_carry_overflow_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@times_word.
      /// \return Identifier \@times_word.
      inline
      const core::identifier_string& times_word_name()
      {
        static core::identifier_string times_word_name = core::identifier_string("@times_word");
        return times_word_name;
      }

      /// \brief Constructor for function symbol \@times_word.
      
      /// \return Function symbol times_word.
      inline
      const function_symbol& times_word()
      {
        static function_symbol times_word(times_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return times_word;
      }

      /// \brief Recogniser for function \@times_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@times_word.
      inline
      bool is_times_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == times_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@times_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@times_word to a number of arguments.
      inline
      application times_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::times_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@times_word.
      /// \param result The data expression where the \@times_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::times_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@times_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol times_word to a
      ///     number of arguments.
      inline
      bool is_times_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_times_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@times_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@times_word to a number of arguments.
      inline
      void times_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void times_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==times_word());
        times_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@times_with_carry_word.
      /// \return Identifier \@times_with_carry_word.
      inline
      const core::identifier_string& times_with_carry_word_name()
      {
        static core::identifier_string times_with_carry_word_name = core::identifier_string("@times_with_carry_word");
        return times_with_carry_word_name;
      }

      /// \brief Constructor for function symbol \@times_with_carry_word.
      
      /// \return Function symbol times_with_carry_word.
      inline
      const function_symbol& times_with_carry_word()
      {
        static function_symbol times_with_carry_word(times_with_carry_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word()));
        return times_with_carry_word;
      }

      /// \brief Recogniser for function \@times_with_carry_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@times_with_carry_word.
      inline
      bool is_times_with_carry_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == times_with_carry_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@times_with_carry_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@times_with_carry_word to a number of arguments.
      inline
      application times_with_carry_word(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_machine_word::times_with_carry_word()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@times_with_carry_word.
      /// \param result The data expression where the \@times_with_carry_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_times_with_carry_word(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_machine_word::times_with_carry_word(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@times_with_carry_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol times_with_carry_word to a
      ///     number of arguments.
      inline
      bool is_times_with_carry_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_times_with_carry_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@times_with_carry_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return The data expression corresponding to an application of \@times_with_carry_word to a number of arguments.
      inline
      void times_with_carry_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void times_with_carry_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==times_with_carry_word());
        times_with_carry_word_manual_implementation(result, a[0], a[1], a[2]);
      }


      /// \brief Generate identifier \@times_overflow_word.
      /// \return Identifier \@times_overflow_word.
      inline
      const core::identifier_string& times_overflow_word_name()
      {
        static core::identifier_string times_overflow_word_name = core::identifier_string("@times_overflow_word");
        return times_overflow_word_name;
      }

      /// \brief Constructor for function symbol \@times_overflow_word.
      
      /// \return Function symbol times_overflow_word.
      inline
      const function_symbol& times_overflow_word()
      {
        static function_symbol times_overflow_word(times_overflow_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return times_overflow_word;
      }

      /// \brief Recogniser for function \@times_overflow_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@times_overflow_word.
      inline
      bool is_times_overflow_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == times_overflow_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@times_overflow_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@times_overflow_word to a number of arguments.
      inline
      application times_overflow_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::times_overflow_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@times_overflow_word.
      /// \param result The data expression where the \@times_overflow_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times_overflow_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::times_overflow_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@times_overflow_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol times_overflow_word to a
      ///     number of arguments.
      inline
      bool is_times_overflow_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_times_overflow_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@times_overflow_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@times_overflow_word to a number of arguments.
      inline
      void times_overflow_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void times_overflow_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==times_overflow_word());
        times_overflow_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@times_with_carry_overflow_word.
      /// \return Identifier \@times_with_carry_overflow_word.
      inline
      const core::identifier_string& times_with_carry_overflow_word_name()
      {
        static core::identifier_string times_with_carry_overflow_word_name = core::identifier_string("@times_with_carry_overflow_word");
        return times_with_carry_overflow_word_name;
      }

      /// \brief Constructor for function symbol \@times_with_carry_overflow_word.
      
      /// \return Function symbol times_with_carry_overflow_word.
      inline
      const function_symbol& times_with_carry_overflow_word()
      {
        static function_symbol times_with_carry_overflow_word(times_with_carry_overflow_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word()));
        return times_with_carry_overflow_word;
      }

      /// \brief Recogniser for function \@times_with_carry_overflow_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@times_with_carry_overflow_word.
      inline
      bool is_times_with_carry_overflow_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == times_with_carry_overflow_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@times_with_carry_overflow_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@times_with_carry_overflow_word to a number of arguments.
      inline
      application times_with_carry_overflow_word(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_machine_word::times_with_carry_overflow_word()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@times_with_carry_overflow_word.
      /// \param result The data expression where the \@times_with_carry_overflow_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_times_with_carry_overflow_word(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_machine_word::times_with_carry_overflow_word(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@times_with_carry_overflow_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol times_with_carry_overflow_word to a
      ///     number of arguments.
      inline
      bool is_times_with_carry_overflow_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_times_with_carry_overflow_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@times_with_carry_overflow_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return The data expression corresponding to an application of \@times_with_carry_overflow_word to a number of arguments.
      inline
      void times_with_carry_overflow_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void times_with_carry_overflow_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==times_with_carry_overflow_word());
        times_with_carry_overflow_word_manual_implementation(result, a[0], a[1], a[2]);
      }


      /// \brief Generate identifier \@minus_word.
      /// \return Identifier \@minus_word.
      inline
      const core::identifier_string& minus_word_name()
      {
        static core::identifier_string minus_word_name = core::identifier_string("@minus_word");
        return minus_word_name;
      }

      /// \brief Constructor for function symbol \@minus_word.
      
      /// \return Function symbol minus_word.
      inline
      const function_symbol& minus_word()
      {
        static function_symbol minus_word(minus_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return minus_word;
      }

      /// \brief Recogniser for function \@minus_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@minus_word.
      inline
      bool is_minus_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == minus_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@minus_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@minus_word to a number of arguments.
      inline
      application minus_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::minus_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@minus_word.
      /// \param result The data expression where the \@minus_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_minus_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::minus_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@minus_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol minus_word to a
      ///     number of arguments.
      inline
      bool is_minus_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_minus_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@minus_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@minus_word to a number of arguments.
      inline
      void minus_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void minus_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==minus_word());
        minus_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@monus_word.
      /// \return Identifier \@monus_word.
      inline
      const core::identifier_string& monus_word_name()
      {
        static core::identifier_string monus_word_name = core::identifier_string("@monus_word");
        return monus_word_name;
      }

      /// \brief Constructor for function symbol \@monus_word.
      
      /// \return Function symbol monus_word.
      inline
      const function_symbol& monus_word()
      {
        static function_symbol monus_word(monus_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return monus_word;
      }

      /// \brief Recogniser for function \@monus_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@monus_word.
      inline
      bool is_monus_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == monus_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@monus_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@monus_word to a number of arguments.
      inline
      application monus_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::monus_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@monus_word.
      /// \param result The data expression where the \@monus_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_monus_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::monus_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@monus_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol monus_word to a
      ///     number of arguments.
      inline
      bool is_monus_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_monus_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@monus_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@monus_word to a number of arguments.
      inline
      void monus_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void monus_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==monus_word());
        monus_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@div_word.
      /// \return Identifier \@div_word.
      inline
      const core::identifier_string& div_word_name()
      {
        static core::identifier_string div_word_name = core::identifier_string("@div_word");
        return div_word_name;
      }

      /// \brief Constructor for function symbol \@div_word.
      
      /// \return Function symbol div_word.
      inline
      const function_symbol& div_word()
      {
        static function_symbol div_word(div_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return div_word;
      }

      /// \brief Recogniser for function \@div_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div_word.
      inline
      bool is_div_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@div_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@div_word to a number of arguments.
      inline
      application div_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::div_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@div_word.
      /// \param result The data expression where the \@div_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_div_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::div_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@div_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div_word to a
      ///     number of arguments.
      inline
      bool is_div_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@div_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@div_word to a number of arguments.
      inline
      void div_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void div_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==div_word());
        div_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@mod_word.
      /// \return Identifier \@mod_word.
      inline
      const core::identifier_string& mod_word_name()
      {
        static core::identifier_string mod_word_name = core::identifier_string("@mod_word");
        return mod_word_name;
      }

      /// \brief Constructor for function symbol \@mod_word.
      
      /// \return Function symbol mod_word.
      inline
      const function_symbol& mod_word()
      {
        static function_symbol mod_word(mod_word_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return mod_word;
      }

      /// \brief Recogniser for function \@mod_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@mod_word.
      inline
      bool is_mod_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == mod_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@mod_word.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@mod_word to a number of arguments.
      inline
      application mod_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::mod_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@mod_word.
      /// \param result The data expression where the \@mod_word expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_mod_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::mod_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@mod_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol mod_word to a
      ///     number of arguments.
      inline
      bool is_mod_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_mod_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@mod_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@mod_word to a number of arguments.
      inline
      void mod_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void mod_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==mod_word());
        mod_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@sqrt_word.
      /// \return Identifier \@sqrt_word.
      inline
      const core::identifier_string& sqrt_word_name()
      {
        static core::identifier_string sqrt_word_name = core::identifier_string("@sqrt_word");
        return sqrt_word_name;
      }

      /// \brief Constructor for function symbol \@sqrt_word.
      
      /// \return Function symbol sqrt_word.
      inline
      const function_symbol& sqrt_word()
      {
        static function_symbol sqrt_word(sqrt_word_name(), make_function_sort_(machine_word(), machine_word()));
        return sqrt_word;
      }

      /// \brief Recogniser for function \@sqrt_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_word.
      inline
      bool is_sqrt_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_word.
      
      /// \param arg0 A data expression.
      /// \return Application of \@sqrt_word to a number of arguments.
      inline
      application sqrt_word(const data_expression& arg0)
      {
        return sort_machine_word::sqrt_word()(arg0);
      }

      /// \brief Make an application of function symbol \@sqrt_word.
      /// \param result The data expression where the \@sqrt_word expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_sqrt_word(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_machine_word::sqrt_word(),arg0);
      }

      /// \brief Recogniser for application of \@sqrt_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_word to a
      ///     number of arguments.
      inline
      bool is_sqrt_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@sqrt_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@sqrt_word to a number of arguments.
      inline
      void sqrt_word_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void sqrt_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==sqrt_word());
        sqrt_word_manual_implementation(result, a[0]);
      }


      /// \brief Generate identifier \@div_doubleword.
      /// \return Identifier \@div_doubleword.
      inline
      const core::identifier_string& div_doubleword_name()
      {
        static core::identifier_string div_doubleword_name = core::identifier_string("@div_doubleword");
        return div_doubleword_name;
      }

      /// \brief Constructor for function symbol \@div_doubleword.
      
      /// \return Function symbol div_doubleword.
      inline
      const function_symbol& div_doubleword()
      {
        static function_symbol div_doubleword(div_doubleword_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word()));
        return div_doubleword;
      }

      /// \brief Recogniser for function \@div_doubleword.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div_doubleword.
      inline
      bool is_div_doubleword_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div_doubleword();
        }
        return false;
      }

      /// \brief Application of function symbol \@div_doubleword.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@div_doubleword to a number of arguments.
      inline
      application div_doubleword(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_machine_word::div_doubleword()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@div_doubleword.
      /// \param result The data expression where the \@div_doubleword expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_div_doubleword(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_machine_word::div_doubleword(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@div_doubleword.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div_doubleword to a
      ///     number of arguments.
      inline
      bool is_div_doubleword_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_doubleword_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@div_doubleword.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return The data expression corresponding to an application of \@div_doubleword to a number of arguments.
      inline
      void div_doubleword_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void div_doubleword_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==div_doubleword());
        div_doubleword_manual_implementation(result, a[0], a[1], a[2]);
      }


      /// \brief Generate identifier \@div_double_doubleword.
      /// \return Identifier \@div_double_doubleword.
      inline
      const core::identifier_string& div_double_doubleword_name()
      {
        static core::identifier_string div_double_doubleword_name = core::identifier_string("@div_double_doubleword");
        return div_double_doubleword_name;
      }

      /// \brief Constructor for function symbol \@div_double_doubleword.
      
      /// \return Function symbol div_double_doubleword.
      inline
      const function_symbol& div_double_doubleword()
      {
        static function_symbol div_double_doubleword(div_double_doubleword_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word(), machine_word()));
        return div_double_doubleword;
      }

      /// \brief Recogniser for function \@div_double_doubleword.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div_double_doubleword.
      inline
      bool is_div_double_doubleword_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div_double_doubleword();
        }
        return false;
      }

      /// \brief Application of function symbol \@div_double_doubleword.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@div_double_doubleword to a number of arguments.
      inline
      application div_double_doubleword(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_machine_word::div_double_doubleword()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@div_double_doubleword.
      /// \param result The data expression where the \@div_double_doubleword expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_div_double_doubleword(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_machine_word::div_double_doubleword(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@div_double_doubleword.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div_double_doubleword to a
      ///     number of arguments.
      inline
      bool is_div_double_doubleword_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_double_doubleword_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@div_double_doubleword.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return The data expression corresponding to an application of \@div_double_doubleword to a number of arguments.
      inline
      void div_double_doubleword_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void div_double_doubleword_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==div_double_doubleword());
        div_double_doubleword_manual_implementation(result, a[0], a[1], a[2], a[3]);
      }


      /// \brief Generate identifier \@div_triple_doubleword.
      /// \return Identifier \@div_triple_doubleword.
      inline
      const core::identifier_string& div_triple_doubleword_name()
      {
        static core::identifier_string div_triple_doubleword_name = core::identifier_string("@div_triple_doubleword");
        return div_triple_doubleword_name;
      }

      /// \brief Constructor for function symbol \@div_triple_doubleword.
      
      /// \return Function symbol div_triple_doubleword.
      inline
      const function_symbol& div_triple_doubleword()
      {
        static function_symbol div_triple_doubleword(div_triple_doubleword_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word(), machine_word(), machine_word()));
        return div_triple_doubleword;
      }

      /// \brief Recogniser for function \@div_triple_doubleword.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div_triple_doubleword.
      inline
      bool is_div_triple_doubleword_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div_triple_doubleword();
        }
        return false;
      }

      /// \brief Application of function symbol \@div_triple_doubleword.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@div_triple_doubleword to a number of arguments.
      inline
      application div_triple_doubleword(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_machine_word::div_triple_doubleword()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@div_triple_doubleword.
      /// \param result The data expression where the \@div_triple_doubleword expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_div_triple_doubleword(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_machine_word::div_triple_doubleword(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@div_triple_doubleword.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div_triple_doubleword to a
      ///     number of arguments.
      inline
      bool is_div_triple_doubleword_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_triple_doubleword_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@div_triple_doubleword.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return The data expression corresponding to an application of \@div_triple_doubleword to a number of arguments.
      inline
      void div_triple_doubleword_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void div_triple_doubleword_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==div_triple_doubleword());
        div_triple_doubleword_manual_implementation(result, a[0], a[1], a[2], a[3], a[4]);
      }


      /// \brief Generate identifier \@mod_doubleword.
      /// \return Identifier \@mod_doubleword.
      inline
      const core::identifier_string& mod_doubleword_name()
      {
        static core::identifier_string mod_doubleword_name = core::identifier_string("@mod_doubleword");
        return mod_doubleword_name;
      }

      /// \brief Constructor for function symbol \@mod_doubleword.
      
      /// \return Function symbol mod_doubleword.
      inline
      const function_symbol& mod_doubleword()
      {
        static function_symbol mod_doubleword(mod_doubleword_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word()));
        return mod_doubleword;
      }

      /// \brief Recogniser for function \@mod_doubleword.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@mod_doubleword.
      inline
      bool is_mod_doubleword_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == mod_doubleword();
        }
        return false;
      }

      /// \brief Application of function symbol \@mod_doubleword.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@mod_doubleword to a number of arguments.
      inline
      application mod_doubleword(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_machine_word::mod_doubleword()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@mod_doubleword.
      /// \param result The data expression where the \@mod_doubleword expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_mod_doubleword(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_machine_word::mod_doubleword(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@mod_doubleword.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol mod_doubleword to a
      ///     number of arguments.
      inline
      bool is_mod_doubleword_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_mod_doubleword_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@mod_doubleword.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return The data expression corresponding to an application of \@mod_doubleword to a number of arguments.
      inline
      void mod_doubleword_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void mod_doubleword_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==mod_doubleword());
        mod_doubleword_manual_implementation(result, a[0], a[1], a[2]);
      }


      /// \brief Generate identifier \@sqrt_doubleword.
      /// \return Identifier \@sqrt_doubleword.
      inline
      const core::identifier_string& sqrt_doubleword_name()
      {
        static core::identifier_string sqrt_doubleword_name = core::identifier_string("@sqrt_doubleword");
        return sqrt_doubleword_name;
      }

      /// \brief Constructor for function symbol \@sqrt_doubleword.
      
      /// \return Function symbol sqrt_doubleword.
      inline
      const function_symbol& sqrt_doubleword()
      {
        static function_symbol sqrt_doubleword(sqrt_doubleword_name(), make_function_sort_(machine_word(), machine_word(), machine_word()));
        return sqrt_doubleword;
      }

      /// \brief Recogniser for function \@sqrt_doubleword.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_doubleword.
      inline
      bool is_sqrt_doubleword_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_doubleword();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_doubleword.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@sqrt_doubleword to a number of arguments.
      inline
      application sqrt_doubleword(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::sqrt_doubleword()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@sqrt_doubleword.
      /// \param result The data expression where the \@sqrt_doubleword expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_sqrt_doubleword(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::sqrt_doubleword(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@sqrt_doubleword.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_doubleword to a
      ///     number of arguments.
      inline
      bool is_sqrt_doubleword_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_doubleword_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@sqrt_doubleword.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@sqrt_doubleword to a number of arguments.
      inline
      void sqrt_doubleword_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void sqrt_doubleword_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==sqrt_doubleword());
        sqrt_doubleword_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@sqrt_tripleword.
      /// \return Identifier \@sqrt_tripleword.
      inline
      const core::identifier_string& sqrt_tripleword_name()
      {
        static core::identifier_string sqrt_tripleword_name = core::identifier_string("@sqrt_tripleword");
        return sqrt_tripleword_name;
      }

      /// \brief Constructor for function symbol \@sqrt_tripleword.
      
      /// \return Function symbol sqrt_tripleword.
      inline
      const function_symbol& sqrt_tripleword()
      {
        static function_symbol sqrt_tripleword(sqrt_tripleword_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word()));
        return sqrt_tripleword;
      }

      /// \brief Recogniser for function \@sqrt_tripleword.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_tripleword.
      inline
      bool is_sqrt_tripleword_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_tripleword();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_tripleword.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@sqrt_tripleword to a number of arguments.
      inline
      application sqrt_tripleword(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_machine_word::sqrt_tripleword()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@sqrt_tripleword.
      /// \param result The data expression where the \@sqrt_tripleword expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_sqrt_tripleword(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_machine_word::sqrt_tripleword(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@sqrt_tripleword.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_tripleword to a
      ///     number of arguments.
      inline
      bool is_sqrt_tripleword_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_tripleword_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@sqrt_tripleword.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return The data expression corresponding to an application of \@sqrt_tripleword to a number of arguments.
      inline
      void sqrt_tripleword_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void sqrt_tripleword_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==sqrt_tripleword());
        sqrt_tripleword_manual_implementation(result, a[0], a[1], a[2]);
      }


      /// \brief Generate identifier \@sqrt_tripleword_overflow.
      /// \return Identifier \@sqrt_tripleword_overflow.
      inline
      const core::identifier_string& sqrt_tripleword_overflow_name()
      {
        static core::identifier_string sqrt_tripleword_overflow_name = core::identifier_string("@sqrt_tripleword_overflow");
        return sqrt_tripleword_overflow_name;
      }

      /// \brief Constructor for function symbol \@sqrt_tripleword_overflow.
      
      /// \return Function symbol sqrt_tripleword_overflow.
      inline
      const function_symbol& sqrt_tripleword_overflow()
      {
        static function_symbol sqrt_tripleword_overflow(sqrt_tripleword_overflow_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word()));
        return sqrt_tripleword_overflow;
      }

      /// \brief Recogniser for function \@sqrt_tripleword_overflow.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_tripleword_overflow.
      inline
      bool is_sqrt_tripleword_overflow_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_tripleword_overflow();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_tripleword_overflow.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@sqrt_tripleword_overflow to a number of arguments.
      inline
      application sqrt_tripleword_overflow(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_machine_word::sqrt_tripleword_overflow()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@sqrt_tripleword_overflow.
      /// \param result The data expression where the \@sqrt_tripleword_overflow expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_sqrt_tripleword_overflow(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_machine_word::sqrt_tripleword_overflow(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@sqrt_tripleword_overflow.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_tripleword_overflow to a
      ///     number of arguments.
      inline
      bool is_sqrt_tripleword_overflow_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_tripleword_overflow_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@sqrt_tripleword_overflow.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return The data expression corresponding to an application of \@sqrt_tripleword_overflow to a number of arguments.
      inline
      void sqrt_tripleword_overflow_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void sqrt_tripleword_overflow_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==sqrt_tripleword_overflow());
        sqrt_tripleword_overflow_manual_implementation(result, a[0], a[1], a[2]);
      }


      /// \brief Generate identifier \@sqrt_quadrupleword.
      /// \return Identifier \@sqrt_quadrupleword.
      inline
      const core::identifier_string& sqrt_quadrupleword_name()
      {
        static core::identifier_string sqrt_quadrupleword_name = core::identifier_string("@sqrt_quadrupleword");
        return sqrt_quadrupleword_name;
      }

      /// \brief Constructor for function symbol \@sqrt_quadrupleword.
      
      /// \return Function symbol sqrt_quadrupleword.
      inline
      const function_symbol& sqrt_quadrupleword()
      {
        static function_symbol sqrt_quadrupleword(sqrt_quadrupleword_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word(), machine_word()));
        return sqrt_quadrupleword;
      }

      /// \brief Recogniser for function \@sqrt_quadrupleword.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_quadrupleword.
      inline
      bool is_sqrt_quadrupleword_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_quadrupleword();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_quadrupleword.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@sqrt_quadrupleword to a number of arguments.
      inline
      application sqrt_quadrupleword(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_machine_word::sqrt_quadrupleword()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@sqrt_quadrupleword.
      /// \param result The data expression where the \@sqrt_quadrupleword expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_sqrt_quadrupleword(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_machine_word::sqrt_quadrupleword(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@sqrt_quadrupleword.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_quadrupleword to a
      ///     number of arguments.
      inline
      bool is_sqrt_quadrupleword_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_quadrupleword_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@sqrt_quadrupleword.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return The data expression corresponding to an application of \@sqrt_quadrupleword to a number of arguments.
      inline
      void sqrt_quadrupleword_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void sqrt_quadrupleword_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==sqrt_quadrupleword());
        sqrt_quadrupleword_manual_implementation(result, a[0], a[1], a[2], a[3]);
      }


      /// \brief Generate identifier \@sqrt_quadrupleword_overflow.
      /// \return Identifier \@sqrt_quadrupleword_overflow.
      inline
      const core::identifier_string& sqrt_quadrupleword_overflow_name()
      {
        static core::identifier_string sqrt_quadrupleword_overflow_name = core::identifier_string("@sqrt_quadrupleword_overflow");
        return sqrt_quadrupleword_overflow_name;
      }

      /// \brief Constructor for function symbol \@sqrt_quadrupleword_overflow.
      
      /// \return Function symbol sqrt_quadrupleword_overflow.
      inline
      const function_symbol& sqrt_quadrupleword_overflow()
      {
        static function_symbol sqrt_quadrupleword_overflow(sqrt_quadrupleword_overflow_name(), make_function_sort_(machine_word(), machine_word(), machine_word(), machine_word(), machine_word()));
        return sqrt_quadrupleword_overflow;
      }

      /// \brief Recogniser for function \@sqrt_quadrupleword_overflow.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_quadrupleword_overflow.
      inline
      bool is_sqrt_quadrupleword_overflow_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_quadrupleword_overflow();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_quadrupleword_overflow.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@sqrt_quadrupleword_overflow to a number of arguments.
      inline
      application sqrt_quadrupleword_overflow(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_machine_word::sqrt_quadrupleword_overflow()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@sqrt_quadrupleword_overflow.
      /// \param result The data expression where the \@sqrt_quadrupleword_overflow expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_sqrt_quadrupleword_overflow(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_machine_word::sqrt_quadrupleword_overflow(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@sqrt_quadrupleword_overflow.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_quadrupleword_overflow to a
      ///     number of arguments.
      inline
      bool is_sqrt_quadrupleword_overflow_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_quadrupleword_overflow_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@sqrt_quadrupleword_overflow.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return The data expression corresponding to an application of \@sqrt_quadrupleword_overflow to a number of arguments.
      inline
      void sqrt_quadrupleword_overflow_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void sqrt_quadrupleword_overflow_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==sqrt_quadrupleword_overflow());
        sqrt_quadrupleword_overflow_manual_implementation(result, a[0], a[1], a[2], a[3]);
      }


      /// \brief Generate identifier \@pred_word.
      /// \return Identifier \@pred_word.
      inline
      const core::identifier_string& pred_word_name()
      {
        static core::identifier_string pred_word_name = core::identifier_string("@pred_word");
        return pred_word_name;
      }

      /// \brief Constructor for function symbol \@pred_word.
      
      /// \return Function symbol pred_word.
      inline
      const function_symbol& pred_word()
      {
        static function_symbol pred_word(pred_word_name(), make_function_sort_(machine_word(), machine_word()));
        return pred_word;
      }

      /// \brief Recogniser for function \@pred_word.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@pred_word.
      inline
      bool is_pred_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == pred_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@pred_word.
      
      /// \param arg0 A data expression.
      /// \return Application of \@pred_word to a number of arguments.
      inline
      application pred_word(const data_expression& arg0)
      {
        return sort_machine_word::pred_word()(arg0);
      }

      /// \brief Make an application of function symbol \@pred_word.
      /// \param result The data expression where the \@pred_word expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pred_word(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_machine_word::pred_word(),arg0);
      }

      /// \brief Recogniser for application of \@pred_word.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol pred_word to a
      ///     number of arguments.
      inline
      bool is_pred_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_pred_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@pred_word.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@pred_word to a number of arguments.
      inline
      void pred_word_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void pred_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==pred_word());
        pred_word_manual_implementation(result, a[0]);
      }


      /// \brief Generate identifier \@equal.
      /// \return Identifier \@equal.
      inline
      const core::identifier_string& equal_word_name()
      {
        static core::identifier_string equal_word_name = core::identifier_string("@equal");
        return equal_word_name;
      }

      /// \brief Constructor for function symbol \@equal.
      
      /// \return Function symbol equal_word.
      inline
      const function_symbol& equal_word()
      {
        static function_symbol equal_word(equal_word_name(), make_function_sort_(machine_word(), machine_word(), sort_bool::bool_()));
        return equal_word;
      }

      /// \brief Recogniser for function \@equal.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@equal.
      inline
      bool is_equal_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == equal_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@equal.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@equal to a number of arguments.
      inline
      application equal_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::equal_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@equal.
      /// \param result The data expression where the \@equal expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_equal_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::equal_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@equal.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol equal_word to a
      ///     number of arguments.
      inline
      bool is_equal_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_equal_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@equal.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@equal to a number of arguments.
      inline
      void equal_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void equal_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==equal_word());
        equal_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@not_equal.
      /// \return Identifier \@not_equal.
      inline
      const core::identifier_string& not_equal_word_name()
      {
        static core::identifier_string not_equal_word_name = core::identifier_string("@not_equal");
        return not_equal_word_name;
      }

      /// \brief Constructor for function symbol \@not_equal.
      
      /// \return Function symbol not_equal_word.
      inline
      const function_symbol& not_equal_word()
      {
        static function_symbol not_equal_word(not_equal_word_name(), make_function_sort_(machine_word(), machine_word(), sort_bool::bool_()));
        return not_equal_word;
      }

      /// \brief Recogniser for function \@not_equal.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@not_equal.
      inline
      bool is_not_equal_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == not_equal_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@not_equal.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@not_equal to a number of arguments.
      inline
      application not_equal_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::not_equal_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@not_equal.
      /// \param result The data expression where the \@not_equal expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_not_equal_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::not_equal_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@not_equal.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol not_equal_word to a
      ///     number of arguments.
      inline
      bool is_not_equal_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_not_equal_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@not_equal.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@not_equal to a number of arguments.
      inline
      void not_equal_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void not_equal_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==not_equal_word());
        not_equal_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@less.
      /// \return Identifier \@less.
      inline
      const core::identifier_string& less_word_name()
      {
        static core::identifier_string less_word_name = core::identifier_string("@less");
        return less_word_name;
      }

      /// \brief Constructor for function symbol \@less.
      
      /// \return Function symbol less_word.
      inline
      const function_symbol& less_word()
      {
        static function_symbol less_word(less_word_name(), make_function_sort_(machine_word(), machine_word(), sort_bool::bool_()));
        return less_word;
      }

      /// \brief Recogniser for function \@less.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@less.
      inline
      bool is_less_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == less_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@less.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@less to a number of arguments.
      inline
      application less_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::less_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@less.
      /// \param result The data expression where the \@less expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_less_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::less_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@less.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol less_word to a
      ///     number of arguments.
      inline
      bool is_less_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_less_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@less.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@less to a number of arguments.
      inline
      void less_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void less_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==less_word());
        less_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@less_equal.
      /// \return Identifier \@less_equal.
      inline
      const core::identifier_string& less_equal_word_name()
      {
        static core::identifier_string less_equal_word_name = core::identifier_string("@less_equal");
        return less_equal_word_name;
      }

      /// \brief Constructor for function symbol \@less_equal.
      
      /// \return Function symbol less_equal_word.
      inline
      const function_symbol& less_equal_word()
      {
        static function_symbol less_equal_word(less_equal_word_name(), make_function_sort_(machine_word(), machine_word(), sort_bool::bool_()));
        return less_equal_word;
      }

      /// \brief Recogniser for function \@less_equal.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@less_equal.
      inline
      bool is_less_equal_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == less_equal_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@less_equal.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@less_equal to a number of arguments.
      inline
      application less_equal_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::less_equal_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@less_equal.
      /// \param result The data expression where the \@less_equal expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_less_equal_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::less_equal_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@less_equal.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol less_equal_word to a
      ///     number of arguments.
      inline
      bool is_less_equal_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_less_equal_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@less_equal.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@less_equal to a number of arguments.
      inline
      void less_equal_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void less_equal_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==less_equal_word());
        less_equal_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@greater.
      /// \return Identifier \@greater.
      inline
      const core::identifier_string& greater_word_name()
      {
        static core::identifier_string greater_word_name = core::identifier_string("@greater");
        return greater_word_name;
      }

      /// \brief Constructor for function symbol \@greater.
      
      /// \return Function symbol greater_word.
      inline
      const function_symbol& greater_word()
      {
        static function_symbol greater_word(greater_word_name(), make_function_sort_(machine_word(), machine_word(), sort_bool::bool_()));
        return greater_word;
      }

      /// \brief Recogniser for function \@greater.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@greater.
      inline
      bool is_greater_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == greater_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@greater.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@greater to a number of arguments.
      inline
      application greater_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::greater_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@greater.
      /// \param result The data expression where the \@greater expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_greater_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::greater_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@greater.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol greater_word to a
      ///     number of arguments.
      inline
      bool is_greater_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_greater_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@greater.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@greater to a number of arguments.
      inline
      void greater_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void greater_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==greater_word());
        greater_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@greater_equal.
      /// \return Identifier \@greater_equal.
      inline
      const core::identifier_string& greater_equal_word_name()
      {
        static core::identifier_string greater_equal_word_name = core::identifier_string("@greater_equal");
        return greater_equal_word_name;
      }

      /// \brief Constructor for function symbol \@greater_equal.
      
      /// \return Function symbol greater_equal_word.
      inline
      const function_symbol& greater_equal_word()
      {
        static function_symbol greater_equal_word(greater_equal_word_name(), make_function_sort_(machine_word(), machine_word(), sort_bool::bool_()));
        return greater_equal_word;
      }

      /// \brief Recogniser for function \@greater_equal.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@greater_equal.
      inline
      bool is_greater_equal_word_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == greater_equal_word();
        }
        return false;
      }

      /// \brief Application of function symbol \@greater_equal.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@greater_equal to a number of arguments.
      inline
      application greater_equal_word(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::greater_equal_word()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@greater_equal.
      /// \param result The data expression where the \@greater_equal expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_greater_equal_word(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::greater_equal_word(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@greater_equal.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol greater_equal_word to a
      ///     number of arguments.
      inline
      bool is_greater_equal_word_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_greater_equal_word_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@greater_equal.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@greater_equal to a number of arguments.
      inline
      void greater_equal_word_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void greater_equal_word_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==greater_equal_word());
        greater_equal_word_manual_implementation(result, a[0], a[1]);
      }


      /// \brief Generate identifier \@rightmost_bit.
      /// \return Identifier \@rightmost_bit.
      inline
      const core::identifier_string& rightmost_bit_name()
      {
        static core::identifier_string rightmost_bit_name = core::identifier_string("@rightmost_bit");
        return rightmost_bit_name;
      }

      /// \brief Constructor for function symbol \@rightmost_bit.
      
      /// \return Function symbol rightmost_bit.
      inline
      const function_symbol& rightmost_bit()
      {
        static function_symbol rightmost_bit(rightmost_bit_name(), make_function_sort_(machine_word(), sort_bool::bool_()));
        return rightmost_bit;
      }

      /// \brief Recogniser for function \@rightmost_bit.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@rightmost_bit.
      inline
      bool is_rightmost_bit_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == rightmost_bit();
        }
        return false;
      }

      /// \brief Application of function symbol \@rightmost_bit.
      
      /// \param arg0 A data expression.
      /// \return Application of \@rightmost_bit to a number of arguments.
      inline
      application rightmost_bit(const data_expression& arg0)
      {
        return sort_machine_word::rightmost_bit()(arg0);
      }

      /// \brief Make an application of function symbol \@rightmost_bit.
      /// \param result The data expression where the \@rightmost_bit expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_rightmost_bit(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_machine_word::rightmost_bit(),arg0);
      }

      /// \brief Recogniser for application of \@rightmost_bit.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol rightmost_bit to a
      ///     number of arguments.
      inline
      bool is_rightmost_bit_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_rightmost_bit_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@rightmost_bit.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \return The data expression corresponding to an application of \@rightmost_bit to a number of arguments.
      inline
      void rightmost_bit_manual_implementation(data_expression& result, const data_expression& arg0);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void rightmost_bit_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==rightmost_bit());
        rightmost_bit_manual_implementation(result, a[0]);
      }


      /// \brief Generate identifier \@shift_right.
      /// \return Identifier \@shift_right.
      inline
      const core::identifier_string& shift_right_name()
      {
        static core::identifier_string shift_right_name = core::identifier_string("@shift_right");
        return shift_right_name;
      }

      /// \brief Constructor for function symbol \@shift_right.
      
      /// \return Function symbol shift_right.
      inline
      const function_symbol& shift_right()
      {
        static function_symbol shift_right(shift_right_name(), make_function_sort_(sort_bool::bool_(), machine_word(), machine_word()));
        return shift_right;
      }

      /// \brief Recogniser for function \@shift_right.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@shift_right.
      inline
      bool is_shift_right_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == shift_right();
        }
        return false;
      }

      /// \brief Application of function symbol \@shift_right.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@shift_right to a number of arguments.
      inline
      application shift_right(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_machine_word::shift_right()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@shift_right.
      /// \param result The data expression where the \@shift_right expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_shift_right(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_machine_word::shift_right(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@shift_right.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol shift_right to a
      ///     number of arguments.
      inline
      bool is_shift_right_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_shift_right_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief The data expression of an application of the function symbol \@shift_right.
      /// \details This function is to be implemented manually.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return The data expression corresponding to an application of \@shift_right to a number of arguments.
      inline
      void shift_right_manual_implementation(data_expression& result, const data_expression& arg0, const data_expression& arg1);


      /// \brief Application of a function that is user defined instead of by rewrite rules. It does not have sort parameters.
      inline
      void shift_right_application(data_expression& result, const data_expression& a1)
      {
        assert(is_application(a1));
        const application& a=atermpp::down_cast<application>(a1);
        // assert(a.head()==shift_right());
        shift_right_manual_implementation(result, a[0], a[1]);
      }

      /// \brief Give all system defined mappings for machine_word
      /// \return All system defined mappings for machine_word
      inline
      function_symbol_vector machine_word_generate_functions_code()
      {
        function_symbol_vector result;
        result.push_back(sort_machine_word::one_word());
        result.push_back(sort_machine_word::two_word());
        result.push_back(sort_machine_word::three_word());
        result.push_back(sort_machine_word::four_word());
        result.push_back(sort_machine_word::max_word());
        result.push_back(sort_machine_word::equals_zero_word());
        result.push_back(sort_machine_word::not_equals_zero_word());
        result.push_back(sort_machine_word::equals_one_word());
        result.push_back(sort_machine_word::equals_max_word());
        result.push_back(sort_machine_word::add_word());
        result.push_back(sort_machine_word::add_with_carry_word());
        result.push_back(sort_machine_word::add_overflow_word());
        result.push_back(sort_machine_word::add_with_carry_overflow_word());
        result.push_back(sort_machine_word::times_word());
        result.push_back(sort_machine_word::times_with_carry_word());
        result.push_back(sort_machine_word::times_overflow_word());
        result.push_back(sort_machine_word::times_with_carry_overflow_word());
        result.push_back(sort_machine_word::minus_word());
        result.push_back(sort_machine_word::monus_word());
        result.push_back(sort_machine_word::div_word());
        result.push_back(sort_machine_word::mod_word());
        result.push_back(sort_machine_word::sqrt_word());
        result.push_back(sort_machine_word::div_doubleword());
        result.push_back(sort_machine_word::div_double_doubleword());
        result.push_back(sort_machine_word::div_triple_doubleword());
        result.push_back(sort_machine_word::mod_doubleword());
        result.push_back(sort_machine_word::sqrt_doubleword());
        result.push_back(sort_machine_word::sqrt_tripleword());
        result.push_back(sort_machine_word::sqrt_tripleword_overflow());
        result.push_back(sort_machine_word::sqrt_quadrupleword());
        result.push_back(sort_machine_word::sqrt_quadrupleword_overflow());
        result.push_back(sort_machine_word::pred_word());
        result.push_back(sort_machine_word::equal_word());
        result.push_back(sort_machine_word::not_equal_word());
        result.push_back(sort_machine_word::less_word());
        result.push_back(sort_machine_word::less_equal_word());
        result.push_back(sort_machine_word::greater_word());
        result.push_back(sort_machine_word::greater_equal_word());
        result.push_back(sort_machine_word::rightmost_bit());
        result.push_back(sort_machine_word::shift_right());
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for machine_word
      /// \return All system defined mappings for machine_word
      inline
      function_symbol_vector machine_word_generate_constructors_and_functions_code()
      {
        function_symbol_vector result=machine_word_generate_functions_code();
        for(const function_symbol& f: machine_word_generate_constructors_code())
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for machine_word
      /// \return All system defined mappings for that can be used in mCRL2 specificationis machine_word
      inline
      function_symbol_vector machine_word_mCRL2_usable_mappings()
      {
        function_symbol_vector result;
        result.push_back(sort_machine_word::one_word());
        result.push_back(sort_machine_word::two_word());
        result.push_back(sort_machine_word::three_word());
        result.push_back(sort_machine_word::four_word());
        result.push_back(sort_machine_word::max_word());
        result.push_back(sort_machine_word::equals_zero_word());
        result.push_back(sort_machine_word::not_equals_zero_word());
        result.push_back(sort_machine_word::equals_one_word());
        result.push_back(sort_machine_word::equals_max_word());
        result.push_back(sort_machine_word::add_word());
        result.push_back(sort_machine_word::add_with_carry_word());
        result.push_back(sort_machine_word::add_overflow_word());
        result.push_back(sort_machine_word::add_with_carry_overflow_word());
        result.push_back(sort_machine_word::times_word());
        result.push_back(sort_machine_word::times_with_carry_word());
        result.push_back(sort_machine_word::times_overflow_word());
        result.push_back(sort_machine_word::times_with_carry_overflow_word());
        result.push_back(sort_machine_word::minus_word());
        result.push_back(sort_machine_word::monus_word());
        result.push_back(sort_machine_word::div_word());
        result.push_back(sort_machine_word::mod_word());
        result.push_back(sort_machine_word::sqrt_word());
        result.push_back(sort_machine_word::div_doubleword());
        result.push_back(sort_machine_word::div_double_doubleword());
        result.push_back(sort_machine_word::div_triple_doubleword());
        result.push_back(sort_machine_word::mod_doubleword());
        result.push_back(sort_machine_word::sqrt_doubleword());
        result.push_back(sort_machine_word::sqrt_tripleword());
        result.push_back(sort_machine_word::sqrt_tripleword_overflow());
        result.push_back(sort_machine_word::sqrt_quadrupleword());
        result.push_back(sort_machine_word::sqrt_quadrupleword_overflow());
        result.push_back(sort_machine_word::pred_word());
        result.push_back(sort_machine_word::equal_word());
        result.push_back(sort_machine_word::not_equal_word());
        result.push_back(sort_machine_word::less_word());
        result.push_back(sort_machine_word::less_equal_word());
        result.push_back(sort_machine_word::greater_word());
        result.push_back(sort_machine_word::greater_equal_word());
        result.push_back(sort_machine_word::rightmost_bit());
        result.push_back(sort_machine_word::shift_right());
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for machine_word
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for machine_word
      inline
      implementation_map machine_word_cpp_implementable_mappings()
      {
        implementation_map result;
        result[sort_machine_word::one_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::one_word_application,"sort_machine_word::one_word_manual_implementation");
        result[sort_machine_word::two_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::two_word_application,"sort_machine_word::two_word_manual_implementation");
        result[sort_machine_word::three_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::three_word_application,"sort_machine_word::three_word_manual_implementation");
        result[sort_machine_word::four_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::four_word_application,"sort_machine_word::four_word_manual_implementation");
        result[sort_machine_word::max_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::max_word_application,"sort_machine_word::max_word_manual_implementation");
        result[sort_machine_word::equals_zero_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::equals_zero_word_application,"sort_machine_word::equals_zero_word_manual_implementation");
        result[sort_machine_word::not_equals_zero_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::not_equals_zero_word_application,"sort_machine_word::not_equals_zero_word_manual_implementation");
        result[sort_machine_word::equals_one_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::equals_one_word_application,"sort_machine_word::equals_one_word_manual_implementation");
        result[sort_machine_word::equals_max_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::equals_max_word_application,"sort_machine_word::equals_max_word_manual_implementation");
        result[sort_machine_word::add_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::add_word_application,"sort_machine_word::add_word_manual_implementation");
        result[sort_machine_word::add_with_carry_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::add_with_carry_word_application,"sort_machine_word::add_with_carry_word_manual_implementation");
        result[sort_machine_word::add_overflow_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::add_overflow_word_application,"sort_machine_word::add_overflow_word_manual_implementation");
        result[sort_machine_word::add_with_carry_overflow_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::add_with_carry_overflow_word_application,"sort_machine_word::add_with_carry_overflow_word_manual_implementation");
        result[sort_machine_word::times_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::times_word_application,"sort_machine_word::times_word_manual_implementation");
        result[sort_machine_word::times_with_carry_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::times_with_carry_word_application,"sort_machine_word::times_with_carry_word_manual_implementation");
        result[sort_machine_word::times_overflow_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::times_overflow_word_application,"sort_machine_word::times_overflow_word_manual_implementation");
        result[sort_machine_word::times_with_carry_overflow_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::times_with_carry_overflow_word_application,"sort_machine_word::times_with_carry_overflow_word_manual_implementation");
        result[sort_machine_word::minus_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::minus_word_application,"sort_machine_word::minus_word_manual_implementation");
        result[sort_machine_word::monus_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::monus_word_application,"sort_machine_word::monus_word_manual_implementation");
        result[sort_machine_word::div_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::div_word_application,"sort_machine_word::div_word_manual_implementation");
        result[sort_machine_word::mod_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::mod_word_application,"sort_machine_word::mod_word_manual_implementation");
        result[sort_machine_word::sqrt_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::sqrt_word_application,"sort_machine_word::sqrt_word_manual_implementation");
        result[sort_machine_word::div_doubleword()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::div_doubleword_application,"sort_machine_word::div_doubleword_manual_implementation");
        result[sort_machine_word::div_double_doubleword()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::div_double_doubleword_application,"sort_machine_word::div_double_doubleword_manual_implementation");
        result[sort_machine_word::div_triple_doubleword()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::div_triple_doubleword_application,"sort_machine_word::div_triple_doubleword_manual_implementation");
        result[sort_machine_word::mod_doubleword()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::mod_doubleword_application,"sort_machine_word::mod_doubleword_manual_implementation");
        result[sort_machine_word::sqrt_doubleword()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::sqrt_doubleword_application,"sort_machine_word::sqrt_doubleword_manual_implementation");
        result[sort_machine_word::sqrt_tripleword()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::sqrt_tripleword_application,"sort_machine_word::sqrt_tripleword_manual_implementation");
        result[sort_machine_word::sqrt_tripleword_overflow()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::sqrt_tripleword_overflow_application,"sort_machine_word::sqrt_tripleword_overflow_manual_implementation");
        result[sort_machine_word::sqrt_quadrupleword()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::sqrt_quadrupleword_application,"sort_machine_word::sqrt_quadrupleword_manual_implementation");
        result[sort_machine_word::sqrt_quadrupleword_overflow()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::sqrt_quadrupleword_overflow_application,"sort_machine_word::sqrt_quadrupleword_overflow_manual_implementation");
        result[sort_machine_word::pred_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::pred_word_application,"sort_machine_word::pred_word_manual_implementation");
        result[sort_machine_word::equal_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::equal_word_application,"sort_machine_word::equal_word_manual_implementation");
        result[sort_machine_word::not_equal_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::not_equal_word_application,"sort_machine_word::not_equal_word_manual_implementation");
        result[sort_machine_word::less_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::less_word_application,"sort_machine_word::less_word_manual_implementation");
        result[sort_machine_word::less_equal_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::less_equal_word_application,"sort_machine_word::less_equal_word_manual_implementation");
        result[sort_machine_word::greater_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::greater_word_application,"sort_machine_word::greater_word_manual_implementation");
        result[sort_machine_word::greater_equal_word()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::greater_equal_word_application,"sort_machine_word::greater_equal_word_manual_implementation");
        result[sort_machine_word::rightmost_bit()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::rightmost_bit_application,"sort_machine_word::rightmost_bit_manual_implementation");
        result[sort_machine_word::shift_right()]=std::pair<std::function<void(data_expression&, const data_expression&)>, std::string>(sort_machine_word::shift_right_application,"sort_machine_word::shift_right_manual_implementation");
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
        assert(is_succ_word_application(e) || is_equals_zero_word_application(e) || is_not_equals_zero_word_application(e) || is_equals_one_word_application(e) || is_equals_max_word_application(e) || is_sqrt_word_application(e) || is_pred_word_application(e) || is_rightmost_bit_application(e));
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
        assert(is_add_word_application(e) || is_add_with_carry_word_application(e) || is_add_overflow_word_application(e) || is_add_with_carry_overflow_word_application(e) || is_times_word_application(e) || is_times_overflow_word_application(e) || is_minus_word_application(e) || is_monus_word_application(e) || is_div_word_application(e) || is_mod_word_application(e) || is_equal_word_application(e) || is_not_equal_word_application(e) || is_less_word_application(e) || is_less_equal_word_application(e) || is_greater_word_application(e) || is_greater_equal_word_application(e));
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
        assert(is_add_word_application(e) || is_add_with_carry_word_application(e) || is_add_overflow_word_application(e) || is_add_with_carry_overflow_word_application(e) || is_times_word_application(e) || is_times_overflow_word_application(e) || is_minus_word_application(e) || is_monus_word_application(e) || is_div_word_application(e) || is_mod_word_application(e) || is_equal_word_application(e) || is_not_equal_word_application(e) || is_less_word_application(e) || is_less_equal_word_application(e) || is_greater_word_application(e) || is_greater_equal_word_application(e));
        return atermpp::down_cast<application>(e)[1];
      }

      ///\brief Function for projecting out argument.
      ///        arg1 from an application.
      /// \param e A data expression.
      /// \pre arg1 is defined for e.
      /// \return The argument of e that corresponds to arg1.
      inline
      const data_expression& arg1(const data_expression& e)
      {
        assert(is_times_with_carry_word_application(e) || is_times_with_carry_overflow_word_application(e) || is_div_doubleword_application(e) || is_div_double_doubleword_application(e) || is_div_triple_doubleword_application(e) || is_mod_doubleword_application(e) || is_sqrt_doubleword_application(e) || is_sqrt_tripleword_application(e) || is_sqrt_tripleword_overflow_application(e) || is_sqrt_quadrupleword_application(e) || is_sqrt_quadrupleword_overflow_application(e) || is_shift_right_application(e));
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
        assert(is_times_with_carry_word_application(e) || is_times_with_carry_overflow_word_application(e) || is_div_doubleword_application(e) || is_div_double_doubleword_application(e) || is_div_triple_doubleword_application(e) || is_mod_doubleword_application(e) || is_sqrt_doubleword_application(e) || is_sqrt_tripleword_application(e) || is_sqrt_tripleword_overflow_application(e) || is_sqrt_quadrupleword_application(e) || is_sqrt_quadrupleword_overflow_application(e) || is_shift_right_application(e));
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
        assert(is_times_with_carry_word_application(e) || is_times_with_carry_overflow_word_application(e) || is_div_doubleword_application(e) || is_div_double_doubleword_application(e) || is_div_triple_doubleword_application(e) || is_mod_doubleword_application(e) || is_sqrt_tripleword_application(e) || is_sqrt_tripleword_overflow_application(e) || is_sqrt_quadrupleword_application(e) || is_sqrt_quadrupleword_overflow_application(e));
        return atermpp::down_cast<application>(e)[2];
      }

      ///\brief Function for projecting out argument.
      ///        arg4 from an application.
      /// \param e A data expression.
      /// \pre arg4 is defined for e.
      /// \return The argument of e that corresponds to arg4.
      inline
      const data_expression& arg4(const data_expression& e)
      {
        assert(is_div_double_doubleword_application(e) || is_div_triple_doubleword_application(e) || is_sqrt_quadrupleword_application(e) || is_sqrt_quadrupleword_overflow_application(e));
        return atermpp::down_cast<application>(e)[3];
      }

      ///\brief Function for projecting out argument.
      ///        arg5 from an application.
      /// \param e A data expression.
      /// \pre arg5 is defined for e.
      /// \return The argument of e that corresponds to arg5.
      inline
      const data_expression& arg5(const data_expression& e)
      {
        assert(is_div_triple_doubleword_application(e));
        return atermpp::down_cast<application>(e)[4];
      }

      /// \brief Give all system defined equations for machine_word
      /// \return All system defined equations for sort machine_word
      inline
      data_equation_vector machine_word_generate_equations_code()
      {
        variable vw1("w1",machine_word());
        variable vw2("w2",machine_word());

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vw1, vw2}), equal_to(vw1, vw2), equal_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vw1, vw2}), less(vw1, vw2), less_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vw1, vw2}), less_equal(vw1, vw2), less_equal_word(vw1, vw2)));
        return result;
      }

    } // namespace sort_machine_word

  } // namespace data

} // namespace mcrl2

#include "mcrl2/data/detail/machine_word.h" // This file contains the manual implementations of rewrite functions.
#endif // MCRL2_DATA_MACHINE_WORD_H
