// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/pos.h
/// \brief The standard sort pos.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/pos.spec.

#ifndef MCRL2_DATA_POS64_H
#define MCRL2_DATA_POS64_H

#include "mcrl2/utilities/exception.h"
#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/machine_word.h"

#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort pos.
    namespace sort_pos {

      inline
      const core::identifier_string& pos_name()
      {
        static core::identifier_string pos_name = core::identifier_string("Pos");
        return pos_name;
      }

      /// \brief Constructor for sort expression Pos.
      /// \return Sort expression Pos.
      inline
      const basic_sort& pos()
      {
        static basic_sort pos = basic_sort(pos_name());
        return pos;
      }

      /// \brief Recogniser for sort expression Pos
      /// \param e A sort expression
      /// \return true iff e == pos()
      inline
      bool is_pos(const sort_expression& e)
      {
        if (is_basic_sort(e))
        {
          return basic_sort(e) == pos();
        }
        return false;
      }


      /// \brief Generate identifier \@c1.
      /// \return Identifier \@c1.
      inline
      const core::identifier_string& c1_name()
      {
        static core::identifier_string c1_name = core::identifier_string("@c1");
        return c1_name;
      }

      /// \brief Constructor for function symbol \@c1.
      
      /// \return Function symbol c1.
      inline
      const function_symbol& c1()
      {
        static function_symbol c1(c1_name(), pos());
        return c1;
      }

      /// \brief Recogniser for function \@c1.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@c1.
      inline
      bool is_c1_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == c1();
        }
        return false;
      }

      /// \brief Generate identifier \@succ_pos.
      /// \return Identifier \@succ_pos.
      inline
      const core::identifier_string& succpos_name()
      {
        static core::identifier_string succpos_name = core::identifier_string("@succ_pos");
        return succpos_name;
      }

      /// \brief Constructor for function symbol \@succ_pos.
      
      /// \return Function symbol succpos.
      inline
      const function_symbol& succpos()
      {
        static function_symbol succpos(succpos_name(), make_function_sort_(pos(), pos()));
        return succpos;
      }

      /// \brief Recogniser for function \@succ_pos.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@succ_pos.
      inline
      bool is_succpos_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == succpos();
        }
        return false;
      }

      /// \brief Application of function symbol \@succ_pos.
      
      /// \param arg0 A data expression.
      /// \return Application of \@succ_pos to a number of arguments.
      inline
      application succpos(const data_expression& arg0)
      {
        return sort_pos::succpos()(arg0);
      }

      /// \brief Make an application of function symbol \@succ_pos.
      /// \param result The data expression where the \@succ_pos expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_succpos(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_pos::succpos(),arg0);
      }

      /// \brief Recogniser for application of \@succ_pos.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol succpos to a
      ///     number of arguments.
      inline
      bool is_succpos_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_succpos_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for pos.
      /// \return All system defined constructors for pos.
      inline
      function_symbol_vector pos_generate_constructors_code()
      {
        function_symbol_vector result;
        result.push_back(sort_pos::c1());
        result.push_back(sort_pos::succpos());

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for pos.
      /// \return All system defined constructors that can be used in an mCRL2 specification for pos.
      inline
      function_symbol_vector pos_mCRL2_usable_constructors()
      {
        function_symbol_vector result;
        result.push_back(sort_pos::c1());
        result.push_back(sort_pos::succpos());

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for pos.
      /// \return All system defined constructors that are to be implemented in C++ for pos.
      inline
      implementation_map pos_cpp_implementable_constructors()
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier \@most_significant_digit.
      /// \return Identifier \@most_significant_digit.
      inline
      const core::identifier_string& most_significant_digit_name()
      {
        static core::identifier_string most_significant_digit_name = core::identifier_string("@most_significant_digit");
        return most_significant_digit_name;
      }

      /// \brief Constructor for function symbol \@most_significant_digit.
      
      /// \return Function symbol most_significant_digit.
      inline
      const function_symbol& most_significant_digit()
      {
        static function_symbol most_significant_digit(most_significant_digit_name(), make_function_sort_(sort_machine_word::machine_word(), pos()));
        return most_significant_digit;
      }

      /// \brief Recogniser for function \@most_significant_digit.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@most_significant_digit.
      inline
      bool is_most_significant_digit_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == most_significant_digit();
        }
        return false;
      }

      /// \brief Application of function symbol \@most_significant_digit.
      
      /// \param arg0 A data expression.
      /// \return Application of \@most_significant_digit to a number of arguments.
      inline
      application most_significant_digit(const data_expression& arg0)
      {
        return sort_pos::most_significant_digit()(arg0);
      }

      /// \brief Make an application of function symbol \@most_significant_digit.
      /// \param result The data expression where the \@most_significant_digit expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_most_significant_digit(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_pos::most_significant_digit(),arg0);
      }

      /// \brief Recogniser for application of \@most_significant_digit.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol most_significant_digit to a
      ///     number of arguments.
      inline
      bool is_most_significant_digit_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_most_significant_digit_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@concat_digit.
      /// \return Identifier \@concat_digit.
      inline
      const core::identifier_string& concat_digit_name()
      {
        static core::identifier_string concat_digit_name = core::identifier_string("@concat_digit");
        return concat_digit_name;
      }

      /// \brief Constructor for function symbol \@concat_digit.
      
      /// \return Function symbol concat_digit.
      inline
      const function_symbol& concat_digit()
      {
        static function_symbol concat_digit(concat_digit_name(), make_function_sort_(pos(), sort_machine_word::machine_word(), pos()));
        return concat_digit;
      }

      /// \brief Recogniser for function \@concat_digit.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@concat_digit.
      inline
      bool is_concat_digit_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == concat_digit();
        }
        return false;
      }

      /// \brief Application of function symbol \@concat_digit.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@concat_digit to a number of arguments.
      inline
      application concat_digit(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::concat_digit()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@concat_digit.
      /// \param result The data expression where the \@concat_digit expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_concat_digit(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::concat_digit(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@concat_digit.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol concat_digit to a
      ///     number of arguments.
      inline
      bool is_concat_digit_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_concat_digit_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@equals_one.
      /// \return Identifier \@equals_one.
      inline
      const core::identifier_string& equals_one_name()
      {
        static core::identifier_string equals_one_name = core::identifier_string("@equals_one");
        return equals_one_name;
      }

      /// \brief Constructor for function symbol \@equals_one.
      
      /// \return Function symbol equals_one.
      inline
      const function_symbol& equals_one()
      {
        static function_symbol equals_one(equals_one_name(), make_function_sort_(pos(), sort_bool::bool_()));
        return equals_one;
      }

      /// \brief Recogniser for function \@equals_one.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@equals_one.
      inline
      bool is_equals_one_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == equals_one();
        }
        return false;
      }

      /// \brief Application of function symbol \@equals_one.
      
      /// \param arg0 A data expression.
      /// \return Application of \@equals_one to a number of arguments.
      inline
      application equals_one(const data_expression& arg0)
      {
        return sort_pos::equals_one()(arg0);
      }

      /// \brief Make an application of function symbol \@equals_one.
      /// \param result The data expression where the \@equals_one expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_equals_one(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_pos::equals_one(),arg0);
      }

      /// \brief Recogniser for application of \@equals_one.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol equals_one to a
      ///     number of arguments.
      inline
      bool is_equals_one_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_equals_one_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier max.
      /// \return Identifier max.
      inline
      const core::identifier_string& maximum_name()
      {
        static core::identifier_string maximum_name = core::identifier_string("max");
        return maximum_name;
      }

      /// \brief Constructor for function symbol max.
      
      /// \return Function symbol maximum.
      inline
      const function_symbol& maximum()
      {
        static function_symbol maximum(maximum_name(), make_function_sort_(pos(), pos(), pos()));
        return maximum;
      }

      /// \brief Recogniser for function max.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching max.
      inline
      bool is_maximum_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == maximum();
        }
        return false;
      }

      /// \brief Application of function symbol max.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of max to a number of arguments.
      inline
      application maximum(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::maximum()(arg0, arg1);
      }

      /// \brief Make an application of function symbol max.
      /// \param result The data expression where the max expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_maximum(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::maximum(),arg0, arg1);
      }

      /// \brief Recogniser for application of max.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol maximum to a
      ///     number of arguments.
      inline
      bool is_maximum_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_maximum_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier min.
      /// \return Identifier min.
      inline
      const core::identifier_string& minimum_name()
      {
        static core::identifier_string minimum_name = core::identifier_string("min");
        return minimum_name;
      }

      /// \brief Constructor for function symbol min.
      
      /// \return Function symbol minimum.
      inline
      const function_symbol& minimum()
      {
        static function_symbol minimum(minimum_name(), make_function_sort_(pos(), pos(), pos()));
        return minimum;
      }

      /// \brief Recogniser for function min.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching min.
      inline
      bool is_minimum_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == minimum();
        }
        return false;
      }

      /// \brief Application of function symbol min.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of min to a number of arguments.
      inline
      application minimum(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::minimum()(arg0, arg1);
      }

      /// \brief Make an application of function symbol min.
      /// \param result The data expression where the min expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_minimum(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::minimum(),arg0, arg1);
      }

      /// \brief Recogniser for application of min.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol minimum to a
      ///     number of arguments.
      inline
      bool is_minimum_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_minimum_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier succ.
      /// \return Identifier succ.
      inline
      const core::identifier_string& succ_name()
      {
        static core::identifier_string succ_name = core::identifier_string("succ");
        return succ_name;
      }

      /// \brief Constructor for function symbol succ.
      
      /// \return Function symbol succ.
      inline
      const function_symbol& succ()
      {
        static function_symbol succ(succ_name(), make_function_sort_(pos(), pos()));
        return succ;
      }

      /// \brief Recogniser for function succ.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching succ.
      inline
      bool is_succ_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == succ();
        }
        return false;
      }

      /// \brief Application of function symbol succ.
      
      /// \param arg0 A data expression.
      /// \return Application of succ to a number of arguments.
      inline
      application succ(const data_expression& arg0)
      {
        return sort_pos::succ()(arg0);
      }

      /// \brief Make an application of function symbol succ.
      /// \param result The data expression where the succ expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_succ(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_pos::succ(),arg0);
      }

      /// \brief Recogniser for application of succ.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol succ to a
      ///     number of arguments.
      inline
      bool is_succ_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_succ_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@pospred.
      /// \return Identifier \@pospred.
      inline
      const core::identifier_string& pos_predecessor_name()
      {
        static core::identifier_string pos_predecessor_name = core::identifier_string("@pospred");
        return pos_predecessor_name;
      }

      /// \brief Constructor for function symbol \@pospred.
      
      /// \return Function symbol pos_predecessor.
      inline
      const function_symbol& pos_predecessor()
      {
        static function_symbol pos_predecessor(pos_predecessor_name(), make_function_sort_(pos(), pos()));
        return pos_predecessor;
      }

      /// \brief Recogniser for function \@pospred.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@pospred.
      inline
      bool is_pos_predecessor_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == pos_predecessor();
        }
        return false;
      }

      /// \brief Application of function symbol \@pospred.
      
      /// \param arg0 A data expression.
      /// \return Application of \@pospred to a number of arguments.
      inline
      application pos_predecessor(const data_expression& arg0)
      {
        return sort_pos::pos_predecessor()(arg0);
      }

      /// \brief Make an application of function symbol \@pospred.
      /// \param result The data expression where the \@pospred expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pos_predecessor(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_pos::pos_predecessor(),arg0);
      }

      /// \brief Recogniser for application of \@pospred.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol pos_predecessor to a
      ///     number of arguments.
      inline
      bool is_pos_predecessor_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_pos_predecessor_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier +.
      /// \return Identifier +.
      inline
      const core::identifier_string& plus_name()
      {
        static core::identifier_string plus_name = core::identifier_string("+");
        return plus_name;
      }

      /// \brief Constructor for function symbol +.
      
      /// \return Function symbol plus.
      inline
      const function_symbol& plus()
      {
        static function_symbol plus(plus_name(), make_function_sort_(pos(), pos(), pos()));
        return plus;
      }

      /// \brief Recogniser for function +.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching +.
      inline
      bool is_plus_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == plus();
        }
        return false;
      }

      /// \brief Application of function symbol +.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of + to a number of arguments.
      inline
      application plus(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::plus()(arg0, arg1);
      }

      /// \brief Make an application of function symbol +.
      /// \param result The data expression where the + expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_plus(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::plus(),arg0, arg1);
      }

      /// \brief Recogniser for application of +.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol plus to a
      ///     number of arguments.
      inline
      bool is_plus_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_plus_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@add_with_carry.
      /// \return Identifier \@add_with_carry.
      inline
      const core::identifier_string& add_with_carry_name()
      {
        static core::identifier_string add_with_carry_name = core::identifier_string("@add_with_carry");
        return add_with_carry_name;
      }

      /// \brief Constructor for function symbol \@add_with_carry.
      
      /// \return Function symbol add_with_carry.
      inline
      const function_symbol& add_with_carry()
      {
        static function_symbol add_with_carry(add_with_carry_name(), make_function_sort_(pos(), pos(), pos()));
        return add_with_carry;
      }

      /// \brief Recogniser for function \@add_with_carry.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@add_with_carry.
      inline
      bool is_add_with_carry_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == add_with_carry();
        }
        return false;
      }

      /// \brief Application of function symbol \@add_with_carry.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@add_with_carry to a number of arguments.
      inline
      application add_with_carry(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::add_with_carry()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@add_with_carry.
      /// \param result The data expression where the \@add_with_carry expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_add_with_carry(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::add_with_carry(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@add_with_carry.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol add_with_carry to a
      ///     number of arguments.
      inline
      bool is_add_with_carry_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_add_with_carry_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@plus_pos.
      /// \return Identifier \@plus_pos.
      inline
      const core::identifier_string& auxiliary_plus_pos_name()
      {
        static core::identifier_string auxiliary_plus_pos_name = core::identifier_string("@plus_pos");
        return auxiliary_plus_pos_name;
      }

      /// \brief Constructor for function symbol \@plus_pos.
      
      /// \return Function symbol auxiliary_plus_pos.
      inline
      const function_symbol& auxiliary_plus_pos()
      {
        static function_symbol auxiliary_plus_pos(auxiliary_plus_pos_name(), make_function_sort_(pos(), pos(), pos()));
        return auxiliary_plus_pos;
      }

      /// \brief Recogniser for function \@plus_pos.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@plus_pos.
      inline
      bool is_auxiliary_plus_pos_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == auxiliary_plus_pos();
        }
        return false;
      }

      /// \brief Application of function symbol \@plus_pos.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@plus_pos to a number of arguments.
      inline
      application auxiliary_plus_pos(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::auxiliary_plus_pos()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@plus_pos.
      /// \param result The data expression where the \@plus_pos expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_auxiliary_plus_pos(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::auxiliary_plus_pos(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@plus_pos.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol auxiliary_plus_pos to a
      ///     number of arguments.
      inline
      bool is_auxiliary_plus_pos_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_auxiliary_plus_pos_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier *.
      /// \return Identifier *.
      inline
      const core::identifier_string& times_name()
      {
        static core::identifier_string times_name = core::identifier_string("*");
        return times_name;
      }

      /// \brief Constructor for function symbol *.
      
      /// \return Function symbol times.
      inline
      const function_symbol& times()
      {
        static function_symbol times(times_name(), make_function_sort_(pos(), pos(), pos()));
        return times;
      }

      /// \brief Recogniser for function *.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching *.
      inline
      bool is_times_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == times();
        }
        return false;
      }

      /// \brief Application of function symbol *.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of * to a number of arguments.
      inline
      application times(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::times()(arg0, arg1);
      }

      /// \brief Make an application of function symbol *.
      /// \param result The data expression where the * expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::times(),arg0, arg1);
      }

      /// \brief Recogniser for application of *.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol times to a
      ///     number of arguments.
      inline
      bool is_times_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_times_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@times_overflow.
      /// \return Identifier \@times_overflow.
      inline
      const core::identifier_string& times_overflow_name()
      {
        static core::identifier_string times_overflow_name = core::identifier_string("@times_overflow");
        return times_overflow_name;
      }

      /// \brief Constructor for function symbol \@times_overflow.
      
      /// \return Function symbol times_overflow.
      inline
      const function_symbol& times_overflow()
      {
        static function_symbol times_overflow(times_overflow_name(), make_function_sort_(pos(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), pos()));
        return times_overflow;
      }

      /// \brief Recogniser for function \@times_overflow.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@times_overflow.
      inline
      bool is_times_overflow_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == times_overflow();
        }
        return false;
      }

      /// \brief Application of function symbol \@times_overflow.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@times_overflow to a number of arguments.
      inline
      application times_overflow(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_pos::times_overflow()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@times_overflow.
      /// \param result The data expression where the \@times_overflow expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_times_overflow(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_pos::times_overflow(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@times_overflow.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol times_overflow to a
      ///     number of arguments.
      inline
      bool is_times_overflow_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_times_overflow_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@times_ordered.
      /// \return Identifier \@times_ordered.
      inline
      const core::identifier_string& times_ordered_name()
      {
        static core::identifier_string times_ordered_name = core::identifier_string("@times_ordered");
        return times_ordered_name;
      }

      /// \brief Constructor for function symbol \@times_ordered.
      
      /// \return Function symbol times_ordered.
      inline
      const function_symbol& times_ordered()
      {
        static function_symbol times_ordered(times_ordered_name(), make_function_sort_(pos(), pos(), pos()));
        return times_ordered;
      }

      /// \brief Recogniser for function \@times_ordered.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@times_ordered.
      inline
      bool is_times_ordered_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == times_ordered();
        }
        return false;
      }

      /// \brief Application of function symbol \@times_ordered.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@times_ordered to a number of arguments.
      inline
      application times_ordered(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::times_ordered()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@times_ordered.
      /// \param result The data expression where the \@times_ordered expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times_ordered(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::times_ordered(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@times_ordered.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol times_ordered to a
      ///     number of arguments.
      inline
      bool is_times_ordered_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_times_ordered_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@times_whr_mult_overflow.
      /// \return Identifier \@times_whr_mult_overflow.
      inline
      const core::identifier_string& times_whr_mult_overflow_name()
      {
        static core::identifier_string times_whr_mult_overflow_name = core::identifier_string("@times_whr_mult_overflow");
        return times_whr_mult_overflow_name;
      }

      /// \brief Constructor for function symbol \@times_whr_mult_overflow.
      
      /// \return Function symbol times_whr_mult_overflow.
      inline
      const function_symbol& times_whr_mult_overflow()
      {
        static function_symbol times_whr_mult_overflow(times_whr_mult_overflow_name(), make_function_sort_(sort_machine_word::machine_word(), sort_machine_word::machine_word(), pos()));
        return times_whr_mult_overflow;
      }

      /// \brief Recogniser for function \@times_whr_mult_overflow.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@times_whr_mult_overflow.
      inline
      bool is_times_whr_mult_overflow_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == times_whr_mult_overflow();
        }
        return false;
      }

      /// \brief Application of function symbol \@times_whr_mult_overflow.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@times_whr_mult_overflow to a number of arguments.
      inline
      application times_whr_mult_overflow(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::times_whr_mult_overflow()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@times_whr_mult_overflow.
      /// \param result The data expression where the \@times_whr_mult_overflow expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times_whr_mult_overflow(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::times_whr_mult_overflow(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@times_whr_mult_overflow.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol times_whr_mult_overflow to a
      ///     number of arguments.
      inline
      bool is_times_whr_mult_overflow_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_times_whr_mult_overflow_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for pos
      /// \return All system defined mappings for pos
      inline
      function_symbol_vector pos_generate_functions_code()
      {
        function_symbol_vector result;
        result.push_back(sort_pos::most_significant_digit());
        result.push_back(sort_pos::concat_digit());
        result.push_back(sort_pos::equals_one());
        result.push_back(sort_pos::maximum());
        result.push_back(sort_pos::minimum());
        result.push_back(sort_pos::succ());
        result.push_back(sort_pos::pos_predecessor());
        result.push_back(sort_pos::plus());
        result.push_back(sort_pos::add_with_carry());
        result.push_back(sort_pos::auxiliary_plus_pos());
        result.push_back(sort_pos::times());
        result.push_back(sort_pos::times_overflow());
        result.push_back(sort_pos::times_ordered());
        result.push_back(sort_pos::times_whr_mult_overflow());
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for pos
      /// \return All system defined mappings for pos
      inline
      function_symbol_vector pos_generate_constructors_and_functions_code()
      {
        function_symbol_vector result=pos_generate_functions_code();
        for(const function_symbol& f: pos_generate_constructors_code())
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for pos
      /// \return All system defined mappings for that can be used in mCRL2 specificationis pos
      inline
      function_symbol_vector pos_mCRL2_usable_mappings()
      {
        function_symbol_vector result;
        result.push_back(sort_pos::most_significant_digit());
        result.push_back(sort_pos::concat_digit());
        result.push_back(sort_pos::equals_one());
        result.push_back(sort_pos::maximum());
        result.push_back(sort_pos::minimum());
        result.push_back(sort_pos::succ());
        result.push_back(sort_pos::pos_predecessor());
        result.push_back(sort_pos::plus());
        result.push_back(sort_pos::add_with_carry());
        result.push_back(sort_pos::auxiliary_plus_pos());
        result.push_back(sort_pos::times());
        result.push_back(sort_pos::times_overflow());
        result.push_back(sort_pos::times_ordered());
        result.push_back(sort_pos::times_whr_mult_overflow());
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for pos
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for pos
      inline
      implementation_map pos_cpp_implementable_mappings()
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
        assert(is_succpos_application(e) || is_most_significant_digit_application(e) || is_equals_one_application(e) || is_succ_application(e) || is_pos_predecessor_application(e));
        return atermpp::down_cast<application>(e)[0];
      }

      ///\brief Function for projecting out argument.
      ///        arg1 from an application.
      /// \param e A data expression.
      /// \pre arg1 is defined for e.
      /// \return The argument of e that corresponds to arg1.
      inline
      const data_expression& arg1(const data_expression& e)
      {
        assert(is_concat_digit_application(e) || is_times_overflow_application(e) || is_times_whr_mult_overflow_application(e));
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
        assert(is_concat_digit_application(e) || is_times_overflow_application(e) || is_times_whr_mult_overflow_application(e));
        return atermpp::down_cast<application>(e)[1];
      }

      ///\brief Function for projecting out argument.
      ///        left from an application.
      /// \param e A data expression.
      /// \pre left is defined for e.
      /// \return The argument of e that corresponds to left.
      inline
      const data_expression& left(const data_expression& e)
      {
        assert(is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_add_with_carry_application(e) || is_auxiliary_plus_pos_application(e) || is_times_application(e) || is_times_ordered_application(e));
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
        assert(is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_add_with_carry_application(e) || is_auxiliary_plus_pos_application(e) || is_times_application(e) || is_times_ordered_application(e));
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
        assert(is_times_overflow_application(e));
        return atermpp::down_cast<application>(e)[2];
      }

      /// \brief Give all system defined equations for pos
      /// \return All system defined equations for sort pos
      inline
      data_equation_vector pos_generate_equations_code()
      {
        variable vb("b",sort_bool::bool_());
        variable vp("p",pos());
        variable vp1("p1",pos());
        variable vp2("p2",pos());
        variable vw("w",sort_machine_word::machine_word());
        variable vw1("w1",sort_machine_word::machine_word());
        variable vw2("w2",sort_machine_word::machine_word());
        variable voverflow("overflow",sort_machine_word::machine_word());

        data_equation_vector result;
        result.push_back(data_equation(variable_list(), c1(), most_significant_digit(sort_machine_word::one_word())));
        result.push_back(data_equation(variable_list({vw}), equals_one(most_significant_digit(vw)), sort_machine_word::equals_one_word(vw)));
        result.push_back(data_equation(variable_list({vp, vw}), equals_one(concat_digit(vp, vw)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp}), equals_one(succpos(vp)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp}), succ(vp), succpos(vp)));
        result.push_back(data_equation(variable_list({vw1}), succpos(most_significant_digit(vw1)), if_(sort_machine_word::equals_max_word(vw1), concat_digit(most_significant_digit(sort_machine_word::one_word()), sort_machine_word::zero_word()), most_significant_digit(sort_machine_word::succ_word(vw1)))));
        result.push_back(data_equation(variable_list({vp, vw1}), succpos(concat_digit(vp, vw1)), if_(sort_machine_word::equals_max_word(vw1), concat_digit(succpos(vp), sort_machine_word::zero_word()), concat_digit(vp, sort_machine_word::succ_word(vw1)))));
        result.push_back(data_equation(variable_list({vw1, vw2}), equal_to(most_significant_digit(vw1), most_significant_digit(vw2)), sort_machine_word::equal_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), equal_to(concat_digit(vp, vw1), most_significant_digit(vw2)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), equal_to(most_significant_digit(vw1), concat_digit(vp, vw2)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp1, vp2, vw1, vw2}), equal_to(concat_digit(vp1, vw1), concat_digit(vp2, vw2)), sort_bool::and_(sort_machine_word::equal_word(vw1, vw2), equal_to(vp1, vp2))));
        result.push_back(data_equation(variable_list({vp1, vp2}), equal_to(succpos(vp1), vp2), sort_bool::and_(sort_bool::not_(equals_one(vp2)), equal_to(vp1, pos_predecessor(vp2)))));
        result.push_back(data_equation(variable_list({vp1, vp2}), equal_to(vp1, succpos(vp2)), sort_bool::and_(sort_bool::not_(equals_one(vp1)), equal_to(pos_predecessor(vp1), vp2))));
        result.push_back(data_equation(variable_list({vw1, vw2}), less(most_significant_digit(vw1), most_significant_digit(vw2)), sort_machine_word::less_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), less(concat_digit(vp, vw1), most_significant_digit(vw2)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), less(most_significant_digit(vw1), concat_digit(vp, vw2)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vp1, vp2, vw1, vw2}), less(concat_digit(vp1, vw1), concat_digit(vp2, vw2)), if_(sort_machine_word::less_word(vw1, vw2), less_equal(vp1, vp2), less(vp1, vp2))));
        result.push_back(data_equation(variable_list({vp1, vp2}), less(succpos(vp1), vp2), sort_bool::and_(less(most_significant_digit(sort_machine_word::two_word()), vp2), less(vp1, pos_predecessor(vp2)))));
        result.push_back(data_equation(variable_list({vp1, vp2}), less(vp1, succpos(vp2)), less_equal(vp1, vp2)));
        result.push_back(data_equation(variable_list({vp, vw1}), sort_machine_word::equals_one_word(vw1), less(vp, most_significant_digit(vw1)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vw1, vw2}), less_equal(most_significant_digit(vw1), most_significant_digit(vw2)), sort_machine_word::less_equal_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), less_equal(concat_digit(vp, vw1), most_significant_digit(vw2)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), less_equal(most_significant_digit(vw1), concat_digit(vp, vw2)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vp1, vp2, vw1, vw2}), less_equal(concat_digit(vp1, vw1), concat_digit(vp2, vw2)), if_(sort_machine_word::less_equal_word(vw1, vw2), less_equal(vp1, vp2), less(vp1, vp2))));
        result.push_back(data_equation(variable_list({vp1, vp2}), less_equal(succpos(vp1), vp2), less(vp1, vp2)));
        result.push_back(data_equation(variable_list({vp1, vp2}), less_equal(vp1, succpos(vp2)), sort_bool::or_(equals_one(vp1), less_equal(pos_predecessor(vp1), vp2))));
        result.push_back(data_equation(variable_list({vp, vw1}), sort_machine_word::equals_one_word(vw1), less_equal(most_significant_digit(vw1), vp), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vp1, vp2}), maximum(vp1, vp2), if_(less_equal(vp1, vp2), vp2, vp1)));
        result.push_back(data_equation(variable_list({vp1, vp2}), minimum(vp1, vp2), if_(less_equal(vp1, vp2), vp1, vp2)));
        result.push_back(data_equation(variable_list({vw1}), pos_predecessor(most_significant_digit(vw1)), if_(sort_machine_word::equals_one_word(vw1), most_significant_digit(sort_machine_word::one_word()), most_significant_digit(sort_machine_word::pred_word(vw1)))));
        result.push_back(data_equation(variable_list({vp, vw1}), pos_predecessor(concat_digit(vp, vw1)), if_(sort_machine_word::equals_zero_word(vw1), if_(equals_one(vp), most_significant_digit(sort_machine_word::max_word()), concat_digit(pos_predecessor(vp), sort_machine_word::max_word())), concat_digit(vp, sort_machine_word::pred_word(vw1)))));
        result.push_back(data_equation(variable_list({vp}), pos_predecessor(succpos(vp)), vp));
        result.push_back(data_equation(variable_list({vw1, vw2}), plus(most_significant_digit(vw1), most_significant_digit(vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(most_significant_digit(sort_machine_word::one_word()), sort_machine_word::add_word(vw1, vw2)), most_significant_digit(sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vw1, vw2}), add_with_carry(most_significant_digit(vw1), most_significant_digit(vw2)), if_(sort_machine_word::add_with_carry_overflow_word(vw1, vw2), concat_digit(most_significant_digit(sort_machine_word::one_word()), sort_machine_word::add_with_carry_word(vw1, vw2)), most_significant_digit(sort_machine_word::add_with_carry_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp1, vw1, vw2}), plus(concat_digit(vp1, vw1), most_significant_digit(vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(succpos(vp1), sort_machine_word::add_word(vw1, vw2)), concat_digit(vp1, sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp1, vw1, vw2}), add_with_carry(concat_digit(vp1, vw1), most_significant_digit(vw2)), if_(sort_machine_word::add_with_carry_overflow_word(vw1, vw2), concat_digit(succpos(vp1), sort_machine_word::add_with_carry_word(vw1, vw2)), concat_digit(vp1, sort_machine_word::add_with_carry_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp2, vw1, vw2}), plus(most_significant_digit(vw1), concat_digit(vp2, vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(succpos(vp2), sort_machine_word::add_word(vw1, vw2)), concat_digit(vp2, sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp2, vw1, vw2}), add_with_carry(most_significant_digit(vw1), concat_digit(vp2, vw2)), if_(sort_machine_word::add_with_carry_overflow_word(vw1, vw2), concat_digit(succpos(vp2), sort_machine_word::add_with_carry_word(vw1, vw2)), concat_digit(vp2, sort_machine_word::add_with_carry_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp1, vp2, vw1, vw2}), plus(concat_digit(vp1, vw1), concat_digit(vp2, vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(add_with_carry(vp1, vp2), sort_machine_word::add_word(vw1, vw2)), concat_digit(plus(vp1, vp2), sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp1, vp2, vw1, vw2}), add_with_carry(concat_digit(vp1, vw1), concat_digit(vp2, vw2)), if_(sort_machine_word::add_with_carry_overflow_word(vw1, vw2), concat_digit(add_with_carry(vp1, vp2), sort_machine_word::add_with_carry_word(vw1, vw2)), concat_digit(plus(vp1, vp2), sort_machine_word::add_with_carry_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp1, vp2}), plus(succpos(vp1), vp2), succpos(plus(vp1, vp2))));
        result.push_back(data_equation(variable_list({vp1, vp2}), plus(vp1, succpos(vp2)), succpos(plus(vp1, vp2))));
        result.push_back(data_equation(variable_list({vp1, vp2}), auxiliary_plus_pos(vp1, vp2), plus(vp1, vp2)));
        result.push_back(data_equation(variable_list({vw1, vw2}), times(most_significant_digit(vw1), most_significant_digit(vw2)), times_whr_mult_overflow(sort_machine_word::times_word(vw1, vw2), sort_machine_word::times_overflow_word(vw1, vw2))));
        result.push_back(data_equation(variable_list({vp2, vw1, vw2}), times(most_significant_digit(vw1), concat_digit(vp2, vw2)), concat_digit(times_overflow(vp2, vw1, sort_machine_word::times_overflow_word(vw1, vw2)), sort_machine_word::times_word(vw1, vw2))));
        result.push_back(data_equation(variable_list({vp1, vw1, vw2}), times(concat_digit(vp1, vw1), most_significant_digit(vw2)), concat_digit(times_overflow(vp1, vw2, sort_machine_word::times_overflow_word(vw1, vw2)), sort_machine_word::times_word(vw1, vw2))));
        result.push_back(data_equation(variable_list({vp1, vp2, vw1, vw2}), times(concat_digit(vp1, vw1), concat_digit(vp2, vw2)), if_(less(vp1, vp2), times_ordered(concat_digit(vp1, vw1), concat_digit(vp2, vw2)), times_ordered(concat_digit(vp2, vw2), concat_digit(vp1, vw1)))));
        result.push_back(data_equation(variable_list({vp2, vw1, vw2}), times_ordered(most_significant_digit(vw1), concat_digit(vp2, vw2)), concat_digit(times_overflow(vp2, vw1, sort_machine_word::times_overflow_word(vw1, vw2)), sort_machine_word::times_word(vw1, vw2))));
        result.push_back(data_equation(variable_list({vp1, vp2, vw1}), times_ordered(concat_digit(vp1, vw1), vp2), plus(concat_digit(times_ordered(vp1, vp2), sort_machine_word::zero_word()), times_overflow(vp2, vw1, sort_machine_word::zero_word()))));
        result.push_back(data_equation(variable_list({vw1, vw2}), times_whr_mult_overflow(vw1, vw2), if_(sort_machine_word::equals_zero_word(vw2), most_significant_digit(vw1), concat_digit(most_significant_digit(vw2), vw1))));
        result.push_back(data_equation(variable_list({voverflow, vw1, vw2}), times_overflow(most_significant_digit(vw1), vw2, voverflow), times_whr_mult_overflow(sort_machine_word::times_with_carry_word(vw1, vw2, voverflow), sort_machine_word::times_with_carry_overflow_word(vw1, vw2, voverflow))));
        result.push_back(data_equation(variable_list({voverflow, vp1, vw1, vw2}), times_overflow(concat_digit(vp1, vw1), vw2, voverflow), concat_digit(times_overflow(vp1, vw2, sort_machine_word::times_with_carry_overflow_word(vw1, vw2, voverflow)), sort_machine_word::times_with_carry_word(vw1, vw2, voverflow))));
        return result;
      }

    } // namespace sort_pos

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_POS64_H
