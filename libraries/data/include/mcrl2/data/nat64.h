// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/nat.h
/// \brief The standard sort nat.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/nat.spec.

#ifndef MCRL2_DATA_NAT64_H
#define MCRL2_DATA_NAT64_H

#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/machine_word.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/pos64.h"

#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort nat.
    namespace sort_nat {

      inline
      const core::identifier_string& nat_name()
      {
        static core::identifier_string nat_name = core::identifier_string("Nat");
        return nat_name;
      }

      /// \brief Constructor for sort expression Nat.
      /// \return Sort expression Nat.
      inline
      const basic_sort& nat()
      {
        static basic_sort nat = basic_sort(nat_name());
        return nat;
      }

      /// \brief Recogniser for sort expression Nat
      /// \param e A sort expression
      /// \return true iff e == nat()
      inline
      bool is_nat(const sort_expression& e)
      {
        if (is_basic_sort(e))
        {
          return basic_sort(e) == nat();
        }
        return false;
      }

      inline
      const core::identifier_string& natnatpair_name()
      {
        static core::identifier_string natnatpair_name = core::identifier_string("@NatNatPair");
        return natnatpair_name;
      }

      /// \brief Constructor for sort expression \@NatNatPair.
      /// \return Sort expression \@NatNatPair.
      inline
      const basic_sort& natnatpair()
      {
        static basic_sort natnatpair = basic_sort(natnatpair_name());
        return natnatpair;
      }

      /// \brief Recogniser for sort expression \@NatNatPair
      /// \param e A sort expression
      /// \return true iff e == natnatpair()
      inline
      bool is_natnatpair(const sort_expression& e)
      {
        if (is_basic_sort(e))
        {
          return basic_sort(e) == natnatpair();
        }
        return false;
      }


      /// \brief Generate identifier \@c0.
      /// \return Identifier \@c0.
      inline
      const core::identifier_string& c0_name()
      {
        static core::identifier_string c0_name = core::identifier_string("@c0");
        return c0_name;
      }

      /// \brief Constructor for function symbol \@c0.
      
      /// \return Function symbol c0.
      inline
      const function_symbol& c0()
      {
        static function_symbol c0(c0_name(), nat());
        return c0;
      }

      /// \brief Recogniser for function \@c0.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@c0.
      inline
      bool is_c0_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == c0();
        }
        return false;
      }

      /// \brief Generate identifier \@succ_nat.
      /// \return Identifier \@succ_nat.
      inline
      const core::identifier_string& succ_nat_name()
      {
        static core::identifier_string succ_nat_name = core::identifier_string("@succ_nat");
        return succ_nat_name;
      }

      /// \brief Constructor for function symbol \@succ_nat.
      
      /// \return Function symbol succ_nat.
      inline
      const function_symbol& succ_nat()
      {
        static function_symbol succ_nat(succ_nat_name(), make_function_sort_(nat(), nat()));
        return succ_nat;
      }

      /// \brief Recogniser for function \@succ_nat.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@succ_nat.
      inline
      bool is_succ_nat_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == succ_nat();
        }
        return false;
      }

      /// \brief Application of function symbol \@succ_nat.
      
      /// \param arg0 A data expression.
      /// \return Application of \@succ_nat to a number of arguments.
      inline
      application succ_nat(const data_expression& arg0)
      {
        return sort_nat::succ_nat()(arg0);
      }

      /// \brief Make an application of function symbol \@succ_nat.
      /// \param result The data expression where the \@succ_nat expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_succ_nat(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::succ_nat(),arg0);
      }

      /// \brief Recogniser for application of \@succ_nat.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol succ_nat to a
      ///     number of arguments.
      inline
      bool is_succ_nat_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_succ_nat_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@nnPair.
      /// \return Identifier \@nnPair.
      inline
      const core::identifier_string& nnpair_name()
      {
        static core::identifier_string nnpair_name = core::identifier_string("@nnPair");
        return nnpair_name;
      }

      /// \brief Constructor for function symbol \@nnPair.
      
      /// \return Function symbol nnpair.
      inline
      const function_symbol& nnpair()
      {
        static function_symbol nnpair(nnpair_name(), make_function_sort_(nat(), nat(), natnatpair()));
        return nnpair;
      }

      /// \brief Recogniser for function \@nnPair.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@nnPair.
      inline
      bool is_nnpair_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == nnpair();
        }
        return false;
      }

      /// \brief Application of function symbol \@nnPair.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@nnPair to a number of arguments.
      inline
      application nnpair(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::nnpair()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@nnPair.
      /// \param result The data expression where the \@nnPair expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_nnpair(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::nnpair(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@nnPair.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol nnpair to a
      ///     number of arguments.
      inline
      bool is_nnpair_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_nnpair_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for nat.
      /// \return All system defined constructors for nat.
      inline
      function_symbol_vector nat_generate_constructors_code()
      {
        function_symbol_vector result;
        result.push_back(sort_nat::c0());
        result.push_back(sort_nat::succ_nat());
        result.push_back(sort_nat::nnpair());

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for nat.
      /// \return All system defined constructors that can be used in an mCRL2 specification for nat.
      inline
      function_symbol_vector nat_mCRL2_usable_constructors()
      {
        function_symbol_vector result;
        result.push_back(sort_nat::c0());
        result.push_back(sort_nat::succ_nat());
        result.push_back(sort_nat::nnpair());

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for nat.
      /// \return All system defined constructors that are to be implemented in C++ for nat.
      inline
      implementation_map nat_cpp_implementable_constructors()
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier \@most_significant_digitNat.
      /// \return Identifier \@most_significant_digitNat.
      inline
      const core::identifier_string& most_significant_digit_nat_name()
      {
        static core::identifier_string most_significant_digit_nat_name = core::identifier_string("@most_significant_digitNat");
        return most_significant_digit_nat_name;
      }

      /// \brief Constructor for function symbol \@most_significant_digitNat.
      
      /// \return Function symbol most_significant_digit_nat.
      inline
      const function_symbol& most_significant_digit_nat()
      {
        static function_symbol most_significant_digit_nat(most_significant_digit_nat_name(), make_function_sort_(sort_machine_word::machine_word(), nat()));
        return most_significant_digit_nat;
      }

      /// \brief Recogniser for function \@most_significant_digitNat.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@most_significant_digitNat.
      inline
      bool is_most_significant_digit_nat_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == most_significant_digit_nat();
        }
        return false;
      }

      /// \brief Application of function symbol \@most_significant_digitNat.
      
      /// \param arg0 A data expression.
      /// \return Application of \@most_significant_digitNat to a number of arguments.
      inline
      application most_significant_digit_nat(const data_expression& arg0)
      {
        return sort_nat::most_significant_digit_nat()(arg0);
      }

      /// \brief Make an application of function symbol \@most_significant_digitNat.
      /// \param result The data expression where the \@most_significant_digitNat expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_most_significant_digit_nat(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::most_significant_digit_nat(),arg0);
      }

      /// \brief Recogniser for application of \@most_significant_digitNat.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol most_significant_digit_nat to a
      ///     number of arguments.
      inline
      bool is_most_significant_digit_nat_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_most_significant_digit_nat_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@concat_digit.
      /// \return Identifier \@concat_digit.
      inline
      const core::identifier_string& concat_digit_name()
      {
        static core::identifier_string concat_digit_name = core::identifier_string("@concat_digit");
        return concat_digit_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol concat_digit(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == nat() && s1 == sort_machine_word::machine_word())
        {
          target_sort = nat();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_machine_word::machine_word())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for concat_digit with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol concat_digit(concat_digit_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == concat_digit_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == concat_digit(nat(), sort_machine_word::machine_word()) || f == concat_digit(sort_pos::pos(), sort_machine_word::machine_word()));
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
        return sort_nat::concat_digit(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@concat_digit.
      /// \param result The data expression where the \@concat_digit expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_concat_digit(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::concat_digit(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier \@equals_zero.
      /// \return Identifier \@equals_zero.
      inline
      const core::identifier_string& equals_zero_name()
      {
        static core::identifier_string equals_zero_name = core::identifier_string("@equals_zero");
        return equals_zero_name;
      }

      /// \brief Constructor for function symbol \@equals_zero.
      
      /// \return Function symbol equals_zero.
      inline
      const function_symbol& equals_zero()
      {
        static function_symbol equals_zero(equals_zero_name(), make_function_sort_(nat(), sort_bool::bool_()));
        return equals_zero;
      }

      /// \brief Recogniser for function \@equals_zero.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@equals_zero.
      inline
      bool is_equals_zero_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == equals_zero();
        }
        return false;
      }

      /// \brief Application of function symbol \@equals_zero.
      
      /// \param arg0 A data expression.
      /// \return Application of \@equals_zero to a number of arguments.
      inline
      application equals_zero(const data_expression& arg0)
      {
        return sort_nat::equals_zero()(arg0);
      }

      /// \brief Make an application of function symbol \@equals_zero.
      /// \param result The data expression where the \@equals_zero expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_equals_zero(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::equals_zero(),arg0);
      }

      /// \brief Recogniser for application of \@equals_zero.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol equals_zero to a
      ///     number of arguments.
      inline
      bool is_equals_zero_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_equals_zero_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@not_equals_zero.
      /// \return Identifier \@not_equals_zero.
      inline
      const core::identifier_string& not_equals_zero_name()
      {
        static core::identifier_string not_equals_zero_name = core::identifier_string("@not_equals_zero");
        return not_equals_zero_name;
      }

      /// \brief Constructor for function symbol \@not_equals_zero.
      
      /// \return Function symbol not_equals_zero.
      inline
      const function_symbol& not_equals_zero()
      {
        static function_symbol not_equals_zero(not_equals_zero_name(), make_function_sort_(nat(), sort_bool::bool_()));
        return not_equals_zero;
      }

      /// \brief Recogniser for function \@not_equals_zero.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@not_equals_zero.
      inline
      bool is_not_equals_zero_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == not_equals_zero();
        }
        return false;
      }

      /// \brief Application of function symbol \@not_equals_zero.
      
      /// \param arg0 A data expression.
      /// \return Application of \@not_equals_zero to a number of arguments.
      inline
      application not_equals_zero(const data_expression& arg0)
      {
        return sort_nat::not_equals_zero()(arg0);
      }

      /// \brief Make an application of function symbol \@not_equals_zero.
      /// \param result The data expression where the \@not_equals_zero expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_not_equals_zero(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::not_equals_zero(),arg0);
      }

      /// \brief Recogniser for application of \@not_equals_zero.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol not_equals_zero to a
      ///     number of arguments.
      inline
      bool is_not_equals_zero_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_not_equals_zero_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@equals_one.
      /// \return Identifier \@equals_one.
      inline
      const core::identifier_string& equals_one_name()
      {
        static core::identifier_string equals_one_name = core::identifier_string("@equals_one");
        return equals_one_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol equals_one(const sort_expression& s0)
      {
        sort_expression target_sort(sort_bool::bool_());
        function_symbol equals_one(equals_one_name(), make_function_sort_(s0, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == equals_one_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == equals_one(nat()) || f == equals_one(sort_pos::pos()));
        }
        return false;
      }

      /// \brief Application of function symbol \@equals_one.
      
      /// \param arg0 A data expression.
      /// \return Application of \@equals_one to a number of arguments.
      inline
      application equals_one(const data_expression& arg0)
      {
        return sort_nat::equals_one(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol \@equals_one.
      /// \param result The data expression where the \@equals_one expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_equals_one(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::equals_one(arg0.sort()),arg0);
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

      /// \brief Generate identifier Pos2Nat.
      /// \return Identifier Pos2Nat.
      inline
      const core::identifier_string& pos2nat_name()
      {
        static core::identifier_string pos2nat_name = core::identifier_string("Pos2Nat");
        return pos2nat_name;
      }

      /// \brief Constructor for function symbol Pos2Nat.
      
      /// \return Function symbol pos2nat.
      inline
      const function_symbol& pos2nat()
      {
        static function_symbol pos2nat(pos2nat_name(), make_function_sort_(sort_pos::pos(), nat()));
        return pos2nat;
      }

      /// \brief Recogniser for function Pos2Nat.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Pos2Nat.
      inline
      bool is_pos2nat_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == pos2nat();
        }
        return false;
      }

      /// \brief Application of function symbol Pos2Nat.
      
      /// \param arg0 A data expression.
      /// \return Application of Pos2Nat to a number of arguments.
      inline
      application pos2nat(const data_expression& arg0)
      {
        return sort_nat::pos2nat()(arg0);
      }

      /// \brief Make an application of function symbol Pos2Nat.
      /// \param result The data expression where the Pos2Nat expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pos2nat(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::pos2nat(),arg0);
      }

      /// \brief Recogniser for application of Pos2Nat.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol pos2nat to a
      ///     number of arguments.
      inline
      bool is_pos2nat_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_pos2nat_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Nat2Pos.
      /// \return Identifier Nat2Pos.
      inline
      const core::identifier_string& nat2pos_name()
      {
        static core::identifier_string nat2pos_name = core::identifier_string("Nat2Pos");
        return nat2pos_name;
      }

      /// \brief Constructor for function symbol Nat2Pos.
      
      /// \return Function symbol nat2pos.
      inline
      const function_symbol& nat2pos()
      {
        static function_symbol nat2pos(nat2pos_name(), make_function_sort_(nat(), sort_pos::pos()));
        return nat2pos;
      }

      /// \brief Recogniser for function Nat2Pos.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Nat2Pos.
      inline
      bool is_nat2pos_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == nat2pos();
        }
        return false;
      }

      /// \brief Application of function symbol Nat2Pos.
      
      /// \param arg0 A data expression.
      /// \return Application of Nat2Pos to a number of arguments.
      inline
      application nat2pos(const data_expression& arg0)
      {
        return sort_nat::nat2pos()(arg0);
      }

      /// \brief Make an application of function symbol Nat2Pos.
      /// \param result The data expression where the Nat2Pos expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_nat2pos(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::nat2pos(),arg0);
      }

      /// \brief Recogniser for application of Nat2Pos.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol nat2pos to a
      ///     number of arguments.
      inline
      bool is_nat2pos_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_nat2pos_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier succ.
      /// \return Identifier succ.
      inline
      const core::identifier_string& succ_name()
      {
        static core::identifier_string succ_name = core::identifier_string("succ");
        return succ_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol succ(const sort_expression& s0)
      {
        sort_expression target_sort(sort_pos::pos());
        function_symbol succ(succ_name(), make_function_sort_(s0, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == succ_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == succ(nat()) || f == succ(sort_pos::pos()));
        }
        return false;
      }

      /// \brief Application of function symbol succ.
      
      /// \param arg0 A data expression.
      /// \return Application of succ to a number of arguments.
      inline
      application succ(const data_expression& arg0)
      {
        return sort_nat::succ(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol succ.
      /// \param result The data expression where the succ expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_succ(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::succ(arg0.sort()),arg0);
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

      /// \brief Generate identifier max.
      /// \return Identifier max.
      inline
      const core::identifier_string& maximum_name()
      {
        static core::identifier_string maximum_name = core::identifier_string("max");
        return maximum_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol maximum(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == sort_pos::pos() && s1 == nat())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == nat() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == nat() && s1 == nat())
        {
          target_sort = nat();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for maximum with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol maximum(maximum_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == maximum_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == maximum(sort_pos::pos(), nat()) || f == maximum(nat(), sort_pos::pos()) || f == maximum(nat(), nat()) || f == maximum(sort_pos::pos(), sort_pos::pos()));
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
        return sort_nat::maximum(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol max.
      /// \param result The data expression where the max expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_maximum(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::maximum(arg0.sort(), arg1.sort()),arg0, arg1);
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

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol minimum(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == nat() && s1 == nat())
        {
          target_sort = nat();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for minimum with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol minimum(minimum_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == minimum_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == minimum(nat(), nat()) || f == minimum(sort_pos::pos(), sort_pos::pos()));
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
        return sort_nat::minimum(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol min.
      /// \param result The data expression where the min expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_minimum(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::minimum(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier pred.
      /// \return Identifier pred.
      inline
      const core::identifier_string& pred_name()
      {
        static core::identifier_string pred_name = core::identifier_string("pred");
        return pred_name;
      }

      /// \brief Constructor for function symbol pred.
      
      /// \return Function symbol pred.
      inline
      const function_symbol& pred()
      {
        static function_symbol pred(pred_name(), make_function_sort_(sort_pos::pos(), nat()));
        return pred;
      }

      /// \brief Recogniser for function pred.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching pred.
      inline
      bool is_pred_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == pred();
        }
        return false;
      }

      /// \brief Application of function symbol pred.
      
      /// \param arg0 A data expression.
      /// \return Application of pred to a number of arguments.
      inline
      application pred(const data_expression& arg0)
      {
        return sort_nat::pred()(arg0);
      }

      /// \brief Make an application of function symbol pred.
      /// \param result The data expression where the pred expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pred(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::pred(),arg0);
      }

      /// \brief Recogniser for application of pred.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol pred to a
      ///     number of arguments.
      inline
      bool is_pred_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_pred_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@pred_whr.
      /// \return Identifier \@pred_whr.
      inline
      const core::identifier_string& pred_whr_name()
      {
        static core::identifier_string pred_whr_name = core::identifier_string("@pred_whr");
        return pred_whr_name;
      }

      /// \brief Constructor for function symbol \@pred_whr.
      
      /// \return Function symbol pred_whr.
      inline
      const function_symbol& pred_whr()
      {
        static function_symbol pred_whr(pred_whr_name(), make_function_sort_(nat(), nat()));
        return pred_whr;
      }

      /// \brief Recogniser for function \@pred_whr.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@pred_whr.
      inline
      bool is_pred_whr_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == pred_whr();
        }
        return false;
      }

      /// \brief Application of function symbol \@pred_whr.
      
      /// \param arg0 A data expression.
      /// \return Application of \@pred_whr to a number of arguments.
      inline
      application pred_whr(const data_expression& arg0)
      {
        return sort_nat::pred_whr()(arg0);
      }

      /// \brief Make an application of function symbol \@pred_whr.
      /// \param result The data expression where the \@pred_whr expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pred_whr(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::pred_whr(),arg0);
      }

      /// \brief Recogniser for application of \@pred_whr.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol pred_whr to a
      ///     number of arguments.
      inline
      bool is_pred_whr_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_pred_whr_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier +.
      /// \return Identifier +.
      inline
      const core::identifier_string& plus_name()
      {
        static core::identifier_string plus_name = core::identifier_string("+");
        return plus_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol plus(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == sort_pos::pos() && s1 == nat())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == nat() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == nat() && s1 == nat())
        {
          target_sort = nat();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for plus with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol plus(plus_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == plus_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == plus(sort_pos::pos(), nat()) || f == plus(nat(), sort_pos::pos()) || f == plus(nat(), nat()) || f == plus(sort_pos::pos(), sort_pos::pos()));
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
        return sort_nat::plus(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol +.
      /// \param result The data expression where the + expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_plus(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::plus(arg0.sort(), arg1.sort()),arg0, arg1);
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

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol add_with_carry(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == nat() && s1 == nat())
        {
          target_sort = nat();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for add_with_carry with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol add_with_carry(add_with_carry_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == add_with_carry_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == add_with_carry(nat(), nat()) || f == add_with_carry(sort_pos::pos(), sort_pos::pos()));
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
        return sort_nat::add_with_carry(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@add_with_carry.
      /// \param result The data expression where the \@add_with_carry expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_add_with_carry(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::add_with_carry(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier \@plus_nat.
      /// \return Identifier \@plus_nat.
      inline
      const core::identifier_string& auxiliary_plus_nat_name()
      {
        static core::identifier_string auxiliary_plus_nat_name = core::identifier_string("@plus_nat");
        return auxiliary_plus_nat_name;
      }

      /// \brief Constructor for function symbol \@plus_nat.
      
      /// \return Function symbol auxiliary_plus_nat.
      inline
      const function_symbol& auxiliary_plus_nat()
      {
        static function_symbol auxiliary_plus_nat(auxiliary_plus_nat_name(), make_function_sort_(nat(), nat(), nat()));
        return auxiliary_plus_nat;
      }

      /// \brief Recogniser for function \@plus_nat.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@plus_nat.
      inline
      bool is_auxiliary_plus_nat_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == auxiliary_plus_nat();
        }
        return false;
      }

      /// \brief Application of function symbol \@plus_nat.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@plus_nat to a number of arguments.
      inline
      application auxiliary_plus_nat(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::auxiliary_plus_nat()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@plus_nat.
      /// \param result The data expression where the \@plus_nat expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_auxiliary_plus_nat(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::auxiliary_plus_nat(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@plus_nat.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol auxiliary_plus_nat to a
      ///     number of arguments.
      inline
      bool is_auxiliary_plus_nat_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_auxiliary_plus_nat_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier *.
      /// \return Identifier *.
      inline
      const core::identifier_string& times_name()
      {
        static core::identifier_string times_name = core::identifier_string("*");
        return times_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol times(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == nat() && s1 == nat())
        {
          target_sort = nat();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for times with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol times(times_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == times_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == times(nat(), nat()) || f == times(sort_pos::pos(), sort_pos::pos()));
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
        return sort_nat::times(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol *.
      /// \param result The data expression where the * expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::times(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier \@times_ordered.
      /// \return Identifier \@times_ordered.
      inline
      const core::identifier_string& times_ordered_name()
      {
        static core::identifier_string times_ordered_name = core::identifier_string("@times_ordered");
        return times_ordered_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol times_ordered(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == nat() && s1 == nat())
        {
          target_sort = nat();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for times_ordered with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol times_ordered(times_ordered_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == times_ordered_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == times_ordered(nat(), nat()) || f == times_ordered(sort_pos::pos(), sort_pos::pos()));
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
        return sort_nat::times_ordered(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@times_ordered.
      /// \param result The data expression where the \@times_ordered expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times_ordered(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::times_ordered(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier \@times_overflow.
      /// \return Identifier \@times_overflow.
      inline
      const core::identifier_string& times_overflow_name()
      {
        static core::identifier_string times_overflow_name = core::identifier_string("@times_overflow");
        return times_overflow_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol times_overflow(const sort_expression& s0, const sort_expression& s1, const sort_expression& s2)
      {
        sort_expression target_sort;
        if (s0 == nat() && s1 == sort_machine_word::machine_word() && s2 == sort_machine_word::machine_word())
        {
          target_sort = nat();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_machine_word::machine_word() && s2 == sort_machine_word::machine_word())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for times_overflow with domain sorts " + pp(s0) + ", " + pp(s1) + ", " + pp(s2) + ". ");
        }

        function_symbol times_overflow(times_overflow_name(), make_function_sort_(s0, s1, s2, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == times_overflow_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 3 && (f == times_overflow(nat(), sort_machine_word::machine_word(), sort_machine_word::machine_word()) || f == times_overflow(sort_pos::pos(), sort_machine_word::machine_word(), sort_machine_word::machine_word()));
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
        return sort_nat::times_overflow(arg0.sort(), arg1.sort(), arg2.sort())(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@times_overflow.
      /// \param result The data expression where the \@times_overflow expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_times_overflow(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::times_overflow(arg0.sort(), arg1.sort(), arg2.sort()),arg0, arg1, arg2);
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

      /// \brief Generate identifier div.
      /// \return Identifier div.
      inline
      const core::identifier_string& div_name()
      {
        static core::identifier_string div_name = core::identifier_string("div");
        return div_name;
      }

      /// \brief Constructor for function symbol div.
      
      /// \return Function symbol div.
      inline
      const function_symbol& div()
      {
        static function_symbol div(div_name(), make_function_sort_(nat(), sort_pos::pos(), nat()));
        return div;
      }

      /// \brief Recogniser for function div.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching div.
      inline
      bool is_div_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div();
        }
        return false;
      }

      /// \brief Application of function symbol div.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of div to a number of arguments.
      inline
      application div(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::div()(arg0, arg1);
      }

      /// \brief Make an application of function symbol div.
      /// \param result The data expression where the div expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_div(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::div(),arg0, arg1);
      }

      /// \brief Recogniser for application of div.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div to a
      ///     number of arguments.
      inline
      bool is_div_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier mod.
      /// \return Identifier mod.
      inline
      const core::identifier_string& mod_name()
      {
        static core::identifier_string mod_name = core::identifier_string("mod");
        return mod_name;
      }

      /// \brief Constructor for function symbol mod.
      
      /// \return Function symbol mod.
      inline
      const function_symbol& mod()
      {
        static function_symbol mod(mod_name(), make_function_sort_(nat(), sort_pos::pos(), nat()));
        return mod;
      }

      /// \brief Recogniser for function mod.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching mod.
      inline
      bool is_mod_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == mod();
        }
        return false;
      }

      /// \brief Application of function symbol mod.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of mod to a number of arguments.
      inline
      application mod(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::mod()(arg0, arg1);
      }

      /// \brief Make an application of function symbol mod.
      /// \param result The data expression where the mod expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_mod(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::mod(),arg0, arg1);
      }

      /// \brief Recogniser for application of mod.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol mod to a
      ///     number of arguments.
      inline
      bool is_mod_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_mod_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier exp.
      /// \return Identifier exp.
      inline
      const core::identifier_string& exp_name()
      {
        static core::identifier_string exp_name = core::identifier_string("exp");
        return exp_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol exp(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == sort_pos::pos() && s1 == nat())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == nat() && s1 == nat())
        {
          target_sort = nat();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for exp with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol exp(exp_name(), make_function_sort_(s0, s1, target_sort));
        return exp;
      }

      /// \brief Recogniser for function exp.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching exp.
      inline
      bool is_exp_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == exp_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == exp(sort_pos::pos(), nat()) || f == exp(nat(), nat()));
        }
        return false;
      }

      /// \brief Application of function symbol exp.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of exp to a number of arguments.
      inline
      application exp(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::exp(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol exp.
      /// \param result The data expression where the exp expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_exp(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::exp(arg0.sort(), arg1.sort()),arg0, arg1);
      }

      /// \brief Recogniser for application of exp.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp to a
      ///     number of arguments.
      inline
      bool is_exp_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier sqrt.
      /// \return Identifier sqrt.
      inline
      const core::identifier_string& sqrt_name()
      {
        static core::identifier_string sqrt_name = core::identifier_string("sqrt");
        return sqrt_name;
      }

      /// \brief Constructor for function symbol sqrt.
      
      /// \return Function symbol sqrt.
      inline
      const function_symbol& sqrt()
      {
        static function_symbol sqrt(sqrt_name(), make_function_sort_(nat(), nat()));
        return sqrt;
      }

      /// \brief Recogniser for function sqrt.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching sqrt.
      inline
      bool is_sqrt_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt();
        }
        return false;
      }

      /// \brief Application of function symbol sqrt.
      
      /// \param arg0 A data expression.
      /// \return Application of sqrt to a number of arguments.
      inline
      application sqrt(const data_expression& arg0)
      {
        return sort_nat::sqrt()(arg0);
      }

      /// \brief Make an application of function symbol sqrt.
      /// \param result The data expression where the sqrt expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_sqrt(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::sqrt(),arg0);
      }

      /// \brief Recogniser for application of sqrt.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt to a
      ///     number of arguments.
      inline
      bool is_sqrt_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@natpred.
      /// \return Identifier \@natpred.
      inline
      const core::identifier_string& natpred_name()
      {
        static core::identifier_string natpred_name = core::identifier_string("@natpred");
        return natpred_name;
      }

      /// \brief Constructor for function symbol \@natpred.
      
      /// \return Function symbol natpred.
      inline
      const function_symbol& natpred()
      {
        static function_symbol natpred(natpred_name(), make_function_sort_(nat(), nat()));
        return natpred;
      }

      /// \brief Recogniser for function \@natpred.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@natpred.
      inline
      bool is_natpred_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == natpred();
        }
        return false;
      }

      /// \brief Application of function symbol \@natpred.
      
      /// \param arg0 A data expression.
      /// \return Application of \@natpred to a number of arguments.
      inline
      application natpred(const data_expression& arg0)
      {
        return sort_nat::natpred()(arg0);
      }

      /// \brief Make an application of function symbol \@natpred.
      /// \param result The data expression where the \@natpred expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_natpred(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::natpred(),arg0);
      }

      /// \brief Recogniser for application of \@natpred.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol natpred to a
      ///     number of arguments.
      inline
      bool is_natpred_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_natpred_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@is_odd.
      /// \return Identifier \@is_odd.
      inline
      const core::identifier_string& is_odd_name()
      {
        static core::identifier_string is_odd_name = core::identifier_string("@is_odd");
        return is_odd_name;
      }

      /// \brief Constructor for function symbol \@is_odd.
      
      /// \return Function symbol is_odd.
      inline
      const function_symbol& is_odd()
      {
        static function_symbol is_odd(is_odd_name(), make_function_sort_(nat(), sort_bool::bool_()));
        return is_odd;
      }

      /// \brief Recogniser for function \@is_odd.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@is_odd.
      inline
      bool is_is_odd_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == is_odd();
        }
        return false;
      }

      /// \brief Application of function symbol \@is_odd.
      
      /// \param arg0 A data expression.
      /// \return Application of \@is_odd to a number of arguments.
      inline
      application is_odd(const data_expression& arg0)
      {
        return sort_nat::is_odd()(arg0);
      }

      /// \brief Make an application of function symbol \@is_odd.
      /// \param result The data expression where the \@is_odd expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_is_odd(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::is_odd(),arg0);
      }

      /// \brief Recogniser for application of \@is_odd.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol is_odd to a
      ///     number of arguments.
      inline
      bool is_is_odd_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_is_odd_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@div2.
      /// \return Identifier \@div2.
      inline
      const core::identifier_string& div2_name()
      {
        static core::identifier_string div2_name = core::identifier_string("@div2");
        return div2_name;
      }

      /// \brief Constructor for function symbol \@div2.
      
      /// \return Function symbol div2.
      inline
      const function_symbol& div2()
      {
        static function_symbol div2(div2_name(), make_function_sort_(nat(), nat()));
        return div2;
      }

      /// \brief Recogniser for function \@div2.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div2.
      inline
      bool is_div2_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div2();
        }
        return false;
      }

      /// \brief Application of function symbol \@div2.
      
      /// \param arg0 A data expression.
      /// \return Application of \@div2 to a number of arguments.
      inline
      application div2(const data_expression& arg0)
      {
        return sort_nat::div2()(arg0);
      }

      /// \brief Make an application of function symbol \@div2.
      /// \param result The data expression where the \@div2 expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_div2(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::div2(),arg0);
      }

      /// \brief Recogniser for application of \@div2.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div2 to a
      ///     number of arguments.
      inline
      bool is_div2_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div2_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@monus.
      /// \return Identifier \@monus.
      inline
      const core::identifier_string& monus_name()
      {
        static core::identifier_string monus_name = core::identifier_string("@monus");
        return monus_name;
      }

      /// \brief Constructor for function symbol \@monus.
      
      /// \return Function symbol monus.
      inline
      const function_symbol& monus()
      {
        static function_symbol monus(monus_name(), make_function_sort_(nat(), nat(), nat()));
        return monus;
      }

      /// \brief Recogniser for function \@monus.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@monus.
      inline
      bool is_monus_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == monus();
        }
        return false;
      }

      /// \brief Application of function symbol \@monus.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@monus to a number of arguments.
      inline
      application monus(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::monus()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@monus.
      /// \param result The data expression where the \@monus expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_monus(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::monus(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@monus.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol monus to a
      ///     number of arguments.
      inline
      bool is_monus_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_monus_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@monus_whr.
      /// \return Identifier \@monus_whr.
      inline
      const core::identifier_string& monus_whr_name()
      {
        static core::identifier_string monus_whr_name = core::identifier_string("@monus_whr");
        return monus_whr_name;
      }

      /// \brief Constructor for function symbol \@monus_whr.
      
      /// \return Function symbol monus_whr.
      inline
      const function_symbol& monus_whr()
      {
        static function_symbol monus_whr(monus_whr_name(), make_function_sort_(nat(), sort_machine_word::machine_word(), nat(), sort_machine_word::machine_word(), nat(), nat()));
        return monus_whr;
      }

      /// \brief Recogniser for function \@monus_whr.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@monus_whr.
      inline
      bool is_monus_whr_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == monus_whr();
        }
        return false;
      }

      /// \brief Application of function symbol \@monus_whr.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@monus_whr to a number of arguments.
      inline
      application monus_whr(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_nat::monus_whr()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@monus_whr.
      /// \param result The data expression where the \@monus_whr expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_monus_whr(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_nat::monus_whr(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@monus_whr.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol monus_whr to a
      ///     number of arguments.
      inline
      bool is_monus_whr_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_monus_whr_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@exp_aux3p.
      /// \return Identifier \@exp_aux3p.
      inline
      const core::identifier_string& exp_aux3p_name()
      {
        static core::identifier_string exp_aux3p_name = core::identifier_string("@exp_aux3p");
        return exp_aux3p_name;
      }

      /// \brief Constructor for function symbol \@exp_aux3p.
      
      /// \return Function symbol exp_aux3p.
      inline
      const function_symbol& exp_aux3p()
      {
        static function_symbol exp_aux3p(exp_aux3p_name(), make_function_sort_(sort_bool::bool_(), sort_pos::pos(), sort_machine_word::machine_word(), sort_pos::pos()));
        return exp_aux3p;
      }

      /// \brief Recogniser for function \@exp_aux3p.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@exp_aux3p.
      inline
      bool is_exp_aux3p_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == exp_aux3p();
        }
        return false;
      }

      /// \brief Application of function symbol \@exp_aux3p.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@exp_aux3p to a number of arguments.
      inline
      application exp_aux3p(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_nat::exp_aux3p()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@exp_aux3p.
      /// \param result The data expression where the \@exp_aux3p expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_exp_aux3p(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::exp_aux3p(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@exp_aux3p.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp_aux3p to a
      ///     number of arguments.
      inline
      bool is_exp_aux3p_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_aux3p_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@exp_aux4p.
      /// \return Identifier \@exp_aux4p.
      inline
      const core::identifier_string& exp_aux4p_name()
      {
        static core::identifier_string exp_aux4p_name = core::identifier_string("@exp_aux4p");
        return exp_aux4p_name;
      }

      /// \brief Constructor for function symbol \@exp_aux4p.
      
      /// \return Function symbol exp_aux4p.
      inline
      const function_symbol& exp_aux4p()
      {
        static function_symbol exp_aux4p(exp_aux4p_name(), make_function_sort_(sort_bool::bool_(), sort_pos::pos(), nat(), sort_machine_word::machine_word(), sort_pos::pos()));
        return exp_aux4p;
      }

      /// \brief Recogniser for function \@exp_aux4p.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@exp_aux4p.
      inline
      bool is_exp_aux4p_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == exp_aux4p();
        }
        return false;
      }

      /// \brief Application of function symbol \@exp_aux4p.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@exp_aux4p to a number of arguments.
      inline
      application exp_aux4p(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::exp_aux4p()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@exp_aux4p.
      /// \param result The data expression where the \@exp_aux4p expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_exp_aux4p(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::exp_aux4p(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@exp_aux4p.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp_aux4p to a
      ///     number of arguments.
      inline
      bool is_exp_aux4p_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_aux4p_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@exp_aux3n.
      /// \return Identifier \@exp_aux3n.
      inline
      const core::identifier_string& exp_aux3n_name()
      {
        static core::identifier_string exp_aux3n_name = core::identifier_string("@exp_aux3n");
        return exp_aux3n_name;
      }

      /// \brief Constructor for function symbol \@exp_aux3n.
      
      /// \return Function symbol exp_aux3n.
      inline
      const function_symbol& exp_aux3n()
      {
        static function_symbol exp_aux3n(exp_aux3n_name(), make_function_sort_(sort_bool::bool_(), nat(), sort_machine_word::machine_word(), nat()));
        return exp_aux3n;
      }

      /// \brief Recogniser for function \@exp_aux3n.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@exp_aux3n.
      inline
      bool is_exp_aux3n_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == exp_aux3n();
        }
        return false;
      }

      /// \brief Application of function symbol \@exp_aux3n.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@exp_aux3n to a number of arguments.
      inline
      application exp_aux3n(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_nat::exp_aux3n()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@exp_aux3n.
      /// \param result The data expression where the \@exp_aux3n expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_exp_aux3n(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::exp_aux3n(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@exp_aux3n.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp_aux3n to a
      ///     number of arguments.
      inline
      bool is_exp_aux3n_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_aux3n_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@exp_aux4n.
      /// \return Identifier \@exp_aux4n.
      inline
      const core::identifier_string& exp_aux4n_name()
      {
        static core::identifier_string exp_aux4n_name = core::identifier_string("@exp_aux4n");
        return exp_aux4n_name;
      }

      /// \brief Constructor for function symbol \@exp_aux4n.
      
      /// \return Function symbol exp_aux4n.
      inline
      const function_symbol& exp_aux4n()
      {
        static function_symbol exp_aux4n(exp_aux4n_name(), make_function_sort_(sort_bool::bool_(), nat(), nat(), sort_machine_word::machine_word(), nat()));
        return exp_aux4n;
      }

      /// \brief Recogniser for function \@exp_aux4n.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@exp_aux4n.
      inline
      bool is_exp_aux4n_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == exp_aux4n();
        }
        return false;
      }

      /// \brief Application of function symbol \@exp_aux4n.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@exp_aux4n to a number of arguments.
      inline
      application exp_aux4n(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::exp_aux4n()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@exp_aux4n.
      /// \param result The data expression where the \@exp_aux4n expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_exp_aux4n(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::exp_aux4n(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@exp_aux4n.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp_aux4n to a
      ///     number of arguments.
      inline
      bool is_exp_aux4n_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_aux4n_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@exp_auxtruep.
      /// \return Identifier \@exp_auxtruep.
      inline
      const core::identifier_string& exp_auxtruep_name()
      {
        static core::identifier_string exp_auxtruep_name = core::identifier_string("@exp_auxtruep");
        return exp_auxtruep_name;
      }

      /// \brief Constructor for function symbol \@exp_auxtruep.
      
      /// \return Function symbol exp_auxtruep.
      inline
      const function_symbol& exp_auxtruep()
      {
        static function_symbol exp_auxtruep(exp_auxtruep_name(), make_function_sort_(sort_pos::pos(), nat(), sort_machine_word::machine_word(), nat()));
        return exp_auxtruep;
      }

      /// \brief Recogniser for function \@exp_auxtruep.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@exp_auxtruep.
      inline
      bool is_exp_auxtruep_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == exp_auxtruep();
        }
        return false;
      }

      /// \brief Application of function symbol \@exp_auxtruep.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@exp_auxtruep to a number of arguments.
      inline
      application exp_auxtruep(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_nat::exp_auxtruep()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@exp_auxtruep.
      /// \param result The data expression where the \@exp_auxtruep expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_exp_auxtruep(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::exp_auxtruep(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@exp_auxtruep.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp_auxtruep to a
      ///     number of arguments.
      inline
      bool is_exp_auxtruep_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_auxtruep_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@exp_auxtruen.
      /// \return Identifier \@exp_auxtruen.
      inline
      const core::identifier_string& exp_auxtruen_name()
      {
        static core::identifier_string exp_auxtruen_name = core::identifier_string("@exp_auxtruen");
        return exp_auxtruen_name;
      }

      /// \brief Constructor for function symbol \@exp_auxtruen.
      
      /// \return Function symbol exp_auxtruen.
      inline
      const function_symbol& exp_auxtruen()
      {
        static function_symbol exp_auxtruen(exp_auxtruen_name(), make_function_sort_(nat(), nat(), sort_machine_word::machine_word(), nat()));
        return exp_auxtruen;
      }

      /// \brief Recogniser for function \@exp_auxtruen.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@exp_auxtruen.
      inline
      bool is_exp_auxtruen_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == exp_auxtruen();
        }
        return false;
      }

      /// \brief Application of function symbol \@exp_auxtruen.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@exp_auxtruen to a number of arguments.
      inline
      application exp_auxtruen(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_nat::exp_auxtruen()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@exp_auxtruen.
      /// \param result The data expression where the \@exp_auxtruen expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_exp_auxtruen(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::exp_auxtruen(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@exp_auxtruen.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp_auxtruen to a
      ///     number of arguments.
      inline
      bool is_exp_auxtruen_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_auxtruen_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@exp_auxfalsep.
      /// \return Identifier \@exp_auxfalsep.
      inline
      const core::identifier_string& exp_auxfalsep_name()
      {
        static core::identifier_string exp_auxfalsep_name = core::identifier_string("@exp_auxfalsep");
        return exp_auxfalsep_name;
      }

      /// \brief Constructor for function symbol \@exp_auxfalsep.
      
      /// \return Function symbol exp_auxfalsep.
      inline
      const function_symbol& exp_auxfalsep()
      {
        static function_symbol exp_auxfalsep(exp_auxfalsep_name(), make_function_sort_(sort_pos::pos(), nat(), sort_machine_word::machine_word(), nat()));
        return exp_auxfalsep;
      }

      /// \brief Recogniser for function \@exp_auxfalsep.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@exp_auxfalsep.
      inline
      bool is_exp_auxfalsep_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == exp_auxfalsep();
        }
        return false;
      }

      /// \brief Application of function symbol \@exp_auxfalsep.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@exp_auxfalsep to a number of arguments.
      inline
      application exp_auxfalsep(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_nat::exp_auxfalsep()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@exp_auxfalsep.
      /// \param result The data expression where the \@exp_auxfalsep expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_exp_auxfalsep(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::exp_auxfalsep(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@exp_auxfalsep.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp_auxfalsep to a
      ///     number of arguments.
      inline
      bool is_exp_auxfalsep_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_auxfalsep_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@exp_auxfalsen.
      /// \return Identifier \@exp_auxfalsen.
      inline
      const core::identifier_string& exp_auxfalsen_name()
      {
        static core::identifier_string exp_auxfalsen_name = core::identifier_string("@exp_auxfalsen");
        return exp_auxfalsen_name;
      }

      /// \brief Constructor for function symbol \@exp_auxfalsen.
      
      /// \return Function symbol exp_auxfalsen.
      inline
      const function_symbol& exp_auxfalsen()
      {
        static function_symbol exp_auxfalsen(exp_auxfalsen_name(), make_function_sort_(nat(), nat(), sort_machine_word::machine_word(), nat()));
        return exp_auxfalsen;
      }

      /// \brief Recogniser for function \@exp_auxfalsen.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@exp_auxfalsen.
      inline
      bool is_exp_auxfalsen_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == exp_auxfalsen();
        }
        return false;
      }

      /// \brief Application of function symbol \@exp_auxfalsen.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@exp_auxfalsen to a number of arguments.
      inline
      application exp_auxfalsen(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_nat::exp_auxfalsen()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@exp_auxfalsen.
      /// \param result The data expression where the \@exp_auxfalsen expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_exp_auxfalsen(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::exp_auxfalsen(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@exp_auxfalsen.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol exp_auxfalsen to a
      ///     number of arguments.
      inline
      bool is_exp_auxfalsen_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_exp_auxfalsen_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@div_bold.
      /// \return Identifier \@div_bold.
      inline
      const core::identifier_string& div_bold_name()
      {
        static core::identifier_string div_bold_name = core::identifier_string("@div_bold");
        return div_bold_name;
      }

      /// \brief Constructor for function symbol \@div_bold.
      
      /// \return Function symbol div_bold.
      inline
      const function_symbol& div_bold()
      {
        static function_symbol div_bold(div_bold_name(), make_function_sort_(nat(), sort_pos::pos(), sort_machine_word::machine_word()));
        return div_bold;
      }

      /// \brief Recogniser for function \@div_bold.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div_bold.
      inline
      bool is_div_bold_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div_bold();
        }
        return false;
      }

      /// \brief Application of function symbol \@div_bold.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@div_bold to a number of arguments.
      inline
      application div_bold(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::div_bold()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@div_bold.
      /// \param result The data expression where the \@div_bold expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_div_bold(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::div_bold(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@div_bold.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div_bold to a
      ///     number of arguments.
      inline
      bool is_div_bold_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_bold_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@div_bold_whr.
      /// \return Identifier \@div_bold_whr.
      inline
      const core::identifier_string& div_bold_whr_name()
      {
        static core::identifier_string div_bold_whr_name = core::identifier_string("@div_bold_whr");
        return div_bold_whr_name;
      }

      /// \brief Constructor for function symbol \@div_bold_whr.
      
      /// \return Function symbol div_bold_whr.
      inline
      const function_symbol& div_bold_whr()
      {
        static function_symbol div_bold_whr(div_bold_whr_name(), make_function_sort_(nat(), sort_machine_word::machine_word(), sort_pos::pos(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word()));
        return div_bold_whr;
      }

      /// \brief Recogniser for function \@div_bold_whr.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div_bold_whr.
      inline
      bool is_div_bold_whr_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div_bold_whr();
        }
        return false;
      }

      /// \brief Application of function symbol \@div_bold_whr.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \param arg5 A data expression.
      /// \return Application of \@div_bold_whr to a number of arguments.
      inline
      application div_bold_whr(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4, const data_expression& arg5)
      {
        return sort_nat::div_bold_whr()(arg0, arg1, arg2, arg3, arg4, arg5);
      }

      /// \brief Make an application of function symbol \@div_bold_whr.
      /// \param result The data expression where the \@div_bold_whr expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \param arg5 A data expression.
      inline
      void make_div_bold_whr(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4, const data_expression& arg5)
      {
        make_application(result, sort_nat::div_bold_whr(),arg0, arg1, arg2, arg3, arg4, arg5);
      }

      /// \brief Recogniser for application of \@div_bold_whr.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div_bold_whr to a
      ///     number of arguments.
      inline
      bool is_div_bold_whr_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_bold_whr_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@div_whr1.
      /// \return Identifier \@div_whr1.
      inline
      const core::identifier_string& div_whr1_name()
      {
        static core::identifier_string div_whr1_name = core::identifier_string("@div_whr1");
        return div_whr1_name;
      }

      /// \brief Constructor for function symbol \@div_whr1.
      
      /// \return Function symbol div_whr1.
      inline
      const function_symbol& div_whr1()
      {
        static function_symbol div_whr1(div_whr1_name(), make_function_sort_(nat(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), natnatpair(), nat()));
        return div_whr1;
      }

      /// \brief Recogniser for function \@div_whr1.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div_whr1.
      inline
      bool is_div_whr1_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div_whr1();
        }
        return false;
      }

      /// \brief Application of function symbol \@div_whr1.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@div_whr1 to a number of arguments.
      inline
      application div_whr1(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::div_whr1()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@div_whr1.
      /// \param result The data expression where the \@div_whr1 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_div_whr1(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::div_whr1(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@div_whr1.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div_whr1 to a
      ///     number of arguments.
      inline
      bool is_div_whr1_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_whr1_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@div_whr2.
      /// \return Identifier \@div_whr2.
      inline
      const core::identifier_string& div_whr2_name()
      {
        static core::identifier_string div_whr2_name = core::identifier_string("@div_whr2");
        return div_whr2_name;
      }

      /// \brief Constructor for function symbol \@div_whr2.
      
      /// \return Function symbol div_whr2.
      inline
      const function_symbol& div_whr2()
      {
        static function_symbol div_whr2(div_whr2_name(), make_function_sort_(nat(), sort_machine_word::machine_word(), sort_pos::pos(), sort_machine_word::machine_word(), natnatpair(), nat()));
        return div_whr2;
      }

      /// \brief Recogniser for function \@div_whr2.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@div_whr2.
      inline
      bool is_div_whr2_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == div_whr2();
        }
        return false;
      }

      /// \brief Application of function symbol \@div_whr2.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@div_whr2 to a number of arguments.
      inline
      application div_whr2(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_nat::div_whr2()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@div_whr2.
      /// \param result The data expression where the \@div_whr2 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_div_whr2(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_nat::div_whr2(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@div_whr2.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol div_whr2 to a
      ///     number of arguments.
      inline
      bool is_div_whr2_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_div_whr2_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@mod_whr1.
      /// \return Identifier \@mod_whr1.
      inline
      const core::identifier_string& mod_whr1_name()
      {
        static core::identifier_string mod_whr1_name = core::identifier_string("@mod_whr1");
        return mod_whr1_name;
      }

      /// \brief Constructor for function symbol \@mod_whr1.
      
      /// \return Function symbol mod_whr1.
      inline
      const function_symbol& mod_whr1()
      {
        static function_symbol mod_whr1(mod_whr1_name(), make_function_sort_(sort_machine_word::machine_word(), sort_pos::pos(), sort_machine_word::machine_word(), nat(), nat()));
        return mod_whr1;
      }

      /// \brief Recogniser for function \@mod_whr1.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@mod_whr1.
      inline
      bool is_mod_whr1_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == mod_whr1();
        }
        return false;
      }

      /// \brief Application of function symbol \@mod_whr1.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@mod_whr1 to a number of arguments.
      inline
      application mod_whr1(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::mod_whr1()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@mod_whr1.
      /// \param result The data expression where the \@mod_whr1 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_mod_whr1(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::mod_whr1(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@mod_whr1.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol mod_whr1 to a
      ///     number of arguments.
      inline
      bool is_mod_whr1_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_mod_whr1_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@divmod_aux.
      /// \return Identifier \@divmod_aux.
      inline
      const core::identifier_string& divmod_aux_name()
      {
        static core::identifier_string divmod_aux_name = core::identifier_string("@divmod_aux");
        return divmod_aux_name;
      }

      /// \brief Constructor for function symbol \@divmod_aux.
      
      /// \return Function symbol divmod_aux.
      inline
      const function_symbol& divmod_aux()
      {
        static function_symbol divmod_aux(divmod_aux_name(), make_function_sort_(nat(), sort_pos::pos(), natnatpair()));
        return divmod_aux;
      }

      /// \brief Recogniser for function \@divmod_aux.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@divmod_aux.
      inline
      bool is_divmod_aux_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == divmod_aux();
        }
        return false;
      }

      /// \brief Application of function symbol \@divmod_aux.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@divmod_aux to a number of arguments.
      inline
      application divmod_aux(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::divmod_aux()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@divmod_aux.
      /// \param result The data expression where the \@divmod_aux expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_divmod_aux(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::divmod_aux(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@divmod_aux.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol divmod_aux to a
      ///     number of arguments.
      inline
      bool is_divmod_aux_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_divmod_aux_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@divmod_aux_whr1.
      /// \return Identifier \@divmod_aux_whr1.
      inline
      const core::identifier_string& divmod_aux_whr1_name()
      {
        static core::identifier_string divmod_aux_whr1_name = core::identifier_string("@divmod_aux_whr1");
        return divmod_aux_whr1_name;
      }

      /// \brief Constructor for function symbol \@divmod_aux_whr1.
      
      /// \return Function symbol divmod_aux_whr1.
      inline
      const function_symbol& divmod_aux_whr1()
      {
        static function_symbol divmod_aux_whr1(divmod_aux_whr1_name(), make_function_sort_(nat(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), natnatpair(), natnatpair()));
        return divmod_aux_whr1;
      }

      /// \brief Recogniser for function \@divmod_aux_whr1.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@divmod_aux_whr1.
      inline
      bool is_divmod_aux_whr1_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == divmod_aux_whr1();
        }
        return false;
      }

      /// \brief Application of function symbol \@divmod_aux_whr1.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@divmod_aux_whr1 to a number of arguments.
      inline
      application divmod_aux_whr1(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::divmod_aux_whr1()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@divmod_aux_whr1.
      /// \param result The data expression where the \@divmod_aux_whr1 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_divmod_aux_whr1(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::divmod_aux_whr1(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@divmod_aux_whr1.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol divmod_aux_whr1 to a
      ///     number of arguments.
      inline
      bool is_divmod_aux_whr1_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_divmod_aux_whr1_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@divmod_aux_whr2.
      /// \return Identifier \@divmod_aux_whr2.
      inline
      const core::identifier_string& divmod_aux_whr2_name()
      {
        static core::identifier_string divmod_aux_whr2_name = core::identifier_string("@divmod_aux_whr2");
        return divmod_aux_whr2_name;
      }

      /// \brief Constructor for function symbol \@divmod_aux_whr2.
      
      /// \return Function symbol divmod_aux_whr2.
      inline
      const function_symbol& divmod_aux_whr2()
      {
        static function_symbol divmod_aux_whr2(divmod_aux_whr2_name(), make_function_sort_(nat(), sort_machine_word::machine_word(), sort_pos::pos(), sort_machine_word::machine_word(), nat(), natnatpair()));
        return divmod_aux_whr2;
      }

      /// \brief Recogniser for function \@divmod_aux_whr2.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@divmod_aux_whr2.
      inline
      bool is_divmod_aux_whr2_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == divmod_aux_whr2();
        }
        return false;
      }

      /// \brief Application of function symbol \@divmod_aux_whr2.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@divmod_aux_whr2 to a number of arguments.
      inline
      application divmod_aux_whr2(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_nat::divmod_aux_whr2()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@divmod_aux_whr2.
      /// \param result The data expression where the \@divmod_aux_whr2 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_divmod_aux_whr2(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_nat::divmod_aux_whr2(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@divmod_aux_whr2.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol divmod_aux_whr2 to a
      ///     number of arguments.
      inline
      bool is_divmod_aux_whr2_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_divmod_aux_whr2_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@divmod_aux_whr3.
      /// \return Identifier \@divmod_aux_whr3.
      inline
      const core::identifier_string& divmod_aux_whr3_name()
      {
        static core::identifier_string divmod_aux_whr3_name = core::identifier_string("@divmod_aux_whr3");
        return divmod_aux_whr3_name;
      }

      /// \brief Constructor for function symbol \@divmod_aux_whr3.
      
      /// \return Function symbol divmod_aux_whr3.
      inline
      const function_symbol& divmod_aux_whr3()
      {
        static function_symbol divmod_aux_whr3(divmod_aux_whr3_name(), make_function_sort_(nat(), sort_machine_word::machine_word(), sort_pos::pos(), sort_machine_word::machine_word(), nat(), natnatpair()));
        return divmod_aux_whr3;
      }

      /// \brief Recogniser for function \@divmod_aux_whr3.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@divmod_aux_whr3.
      inline
      bool is_divmod_aux_whr3_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == divmod_aux_whr3();
        }
        return false;
      }

      /// \brief Application of function symbol \@divmod_aux_whr3.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@divmod_aux_whr3 to a number of arguments.
      inline
      application divmod_aux_whr3(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_nat::divmod_aux_whr3()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@divmod_aux_whr3.
      /// \param result The data expression where the \@divmod_aux_whr3 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_divmod_aux_whr3(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_nat::divmod_aux_whr3(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@divmod_aux_whr3.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol divmod_aux_whr3 to a
      ///     number of arguments.
      inline
      bool is_divmod_aux_whr3_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_divmod_aux_whr3_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@divmod_aux_whr4.
      /// \return Identifier \@divmod_aux_whr4.
      inline
      const core::identifier_string& divmod_aux_whr4_name()
      {
        static core::identifier_string divmod_aux_whr4_name = core::identifier_string("@divmod_aux_whr4");
        return divmod_aux_whr4_name;
      }

      /// \brief Constructor for function symbol \@divmod_aux_whr4.
      
      /// \return Function symbol divmod_aux_whr4.
      inline
      const function_symbol& divmod_aux_whr4()
      {
        static function_symbol divmod_aux_whr4(divmod_aux_whr4_name(), make_function_sort_(sort_machine_word::machine_word(), sort_pos::pos(), sort_machine_word::machine_word(), natnatpair(), natnatpair()));
        return divmod_aux_whr4;
      }

      /// \brief Recogniser for function \@divmod_aux_whr4.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@divmod_aux_whr4.
      inline
      bool is_divmod_aux_whr4_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == divmod_aux_whr4();
        }
        return false;
      }

      /// \brief Application of function symbol \@divmod_aux_whr4.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@divmod_aux_whr4 to a number of arguments.
      inline
      application divmod_aux_whr4(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::divmod_aux_whr4()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@divmod_aux_whr4.
      /// \param result The data expression where the \@divmod_aux_whr4 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_divmod_aux_whr4(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::divmod_aux_whr4(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@divmod_aux_whr4.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol divmod_aux_whr4 to a
      ///     number of arguments.
      inline
      bool is_divmod_aux_whr4_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_divmod_aux_whr4_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@divmod_aux_whr5.
      /// \return Identifier \@divmod_aux_whr5.
      inline
      const core::identifier_string& divmod_aux_whr5_name()
      {
        static core::identifier_string divmod_aux_whr5_name = core::identifier_string("@divmod_aux_whr5");
        return divmod_aux_whr5_name;
      }

      /// \brief Constructor for function symbol \@divmod_aux_whr5.
      
      /// \return Function symbol divmod_aux_whr5.
      inline
      const function_symbol& divmod_aux_whr5()
      {
        static function_symbol divmod_aux_whr5(divmod_aux_whr5_name(), make_function_sort_(sort_pos::pos(), sort_machine_word::machine_word(), natnatpair(), nat(), natnatpair()));
        return divmod_aux_whr5;
      }

      /// \brief Recogniser for function \@divmod_aux_whr5.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@divmod_aux_whr5.
      inline
      bool is_divmod_aux_whr5_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == divmod_aux_whr5();
        }
        return false;
      }

      /// \brief Application of function symbol \@divmod_aux_whr5.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@divmod_aux_whr5 to a number of arguments.
      inline
      application divmod_aux_whr5(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::divmod_aux_whr5()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@divmod_aux_whr5.
      /// \param result The data expression where the \@divmod_aux_whr5 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_divmod_aux_whr5(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::divmod_aux_whr5(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@divmod_aux_whr5.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol divmod_aux_whr5 to a
      ///     number of arguments.
      inline
      bool is_divmod_aux_whr5_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_divmod_aux_whr5_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@divmod_aux_whr6.
      /// \return Identifier \@divmod_aux_whr6.
      inline
      const core::identifier_string& divmod_aux_whr6_name()
      {
        static core::identifier_string divmod_aux_whr6_name = core::identifier_string("@divmod_aux_whr6");
        return divmod_aux_whr6_name;
      }

      /// \brief Constructor for function symbol \@divmod_aux_whr6.
      
      /// \return Function symbol divmod_aux_whr6.
      inline
      const function_symbol& divmod_aux_whr6()
      {
        static function_symbol divmod_aux_whr6(divmod_aux_whr6_name(), make_function_sort_(sort_pos::pos(), sort_machine_word::machine_word(), natnatpair(), nat(), nat(), natnatpair()));
        return divmod_aux_whr6;
      }

      /// \brief Recogniser for function \@divmod_aux_whr6.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@divmod_aux_whr6.
      inline
      bool is_divmod_aux_whr6_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == divmod_aux_whr6();
        }
        return false;
      }

      /// \brief Application of function symbol \@divmod_aux_whr6.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@divmod_aux_whr6 to a number of arguments.
      inline
      application divmod_aux_whr6(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_nat::divmod_aux_whr6()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@divmod_aux_whr6.
      /// \param result The data expression where the \@divmod_aux_whr6 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_divmod_aux_whr6(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_nat::divmod_aux_whr6(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@divmod_aux_whr6.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol divmod_aux_whr6 to a
      ///     number of arguments.
      inline
      bool is_divmod_aux_whr6_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_divmod_aux_whr6_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@msd.
      /// \return Identifier \@msd.
      inline
      const core::identifier_string& msd_name()
      {
        static core::identifier_string msd_name = core::identifier_string("@msd");
        return msd_name;
      }

      /// \brief Constructor for function symbol \@msd.
      
      /// \return Function symbol msd.
      inline
      const function_symbol& msd()
      {
        static function_symbol msd(msd_name(), make_function_sort_(nat(), sort_machine_word::machine_word()));
        return msd;
      }

      /// \brief Recogniser for function \@msd.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@msd.
      inline
      bool is_msd_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == msd();
        }
        return false;
      }

      /// \brief Application of function symbol \@msd.
      
      /// \param arg0 A data expression.
      /// \return Application of \@msd to a number of arguments.
      inline
      application msd(const data_expression& arg0)
      {
        return sort_nat::msd()(arg0);
      }

      /// \brief Make an application of function symbol \@msd.
      /// \param result The data expression where the \@msd expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_msd(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::msd(),arg0);
      }

      /// \brief Recogniser for application of \@msd.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol msd to a
      ///     number of arguments.
      inline
      bool is_msd_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_msd_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@swap_zero.
      /// \return Identifier \@swap_zero.
      inline
      const core::identifier_string& swap_zero_name()
      {
        static core::identifier_string swap_zero_name = core::identifier_string("@swap_zero");
        return swap_zero_name;
      }

      /// \brief Constructor for function symbol \@swap_zero.
      
      /// \return Function symbol swap_zero.
      inline
      const function_symbol& swap_zero()
      {
        static function_symbol swap_zero(swap_zero_name(), make_function_sort_(nat(), nat(), nat()));
        return swap_zero;
      }

      /// \brief Recogniser for function \@swap_zero.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@swap_zero.
      inline
      bool is_swap_zero_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == swap_zero();
        }
        return false;
      }

      /// \brief Application of function symbol \@swap_zero.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@swap_zero to a number of arguments.
      inline
      application swap_zero(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_nat::swap_zero()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@swap_zero.
      /// \param result The data expression where the \@swap_zero expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_swap_zero(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_nat::swap_zero(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@swap_zero.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol swap_zero to a
      ///     number of arguments.
      inline
      bool is_swap_zero_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_swap_zero_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@swap_zero_add.
      /// \return Identifier \@swap_zero_add.
      inline
      const core::identifier_string& swap_zero_add_name()
      {
        static core::identifier_string swap_zero_add_name = core::identifier_string("@swap_zero_add");
        return swap_zero_add_name;
      }

      /// \brief Constructor for function symbol \@swap_zero_add.
      
      /// \return Function symbol swap_zero_add.
      inline
      const function_symbol& swap_zero_add()
      {
        static function_symbol swap_zero_add(swap_zero_add_name(), make_function_sort_(nat(), nat(), nat(), nat(), nat()));
        return swap_zero_add;
      }

      /// \brief Recogniser for function \@swap_zero_add.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@swap_zero_add.
      inline
      bool is_swap_zero_add_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == swap_zero_add();
        }
        return false;
      }

      /// \brief Application of function symbol \@swap_zero_add.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@swap_zero_add to a number of arguments.
      inline
      application swap_zero_add(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::swap_zero_add()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@swap_zero_add.
      /// \param result The data expression where the \@swap_zero_add expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_swap_zero_add(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::swap_zero_add(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@swap_zero_add.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol swap_zero_add to a
      ///     number of arguments.
      inline
      bool is_swap_zero_add_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_swap_zero_add_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@swap_zero_min.
      /// \return Identifier \@swap_zero_min.
      inline
      const core::identifier_string& swap_zero_min_name()
      {
        static core::identifier_string swap_zero_min_name = core::identifier_string("@swap_zero_min");
        return swap_zero_min_name;
      }

      /// \brief Constructor for function symbol \@swap_zero_min.
      
      /// \return Function symbol swap_zero_min.
      inline
      const function_symbol& swap_zero_min()
      {
        static function_symbol swap_zero_min(swap_zero_min_name(), make_function_sort_(nat(), nat(), nat(), nat(), nat()));
        return swap_zero_min;
      }

      /// \brief Recogniser for function \@swap_zero_min.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@swap_zero_min.
      inline
      bool is_swap_zero_min_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == swap_zero_min();
        }
        return false;
      }

      /// \brief Application of function symbol \@swap_zero_min.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@swap_zero_min to a number of arguments.
      inline
      application swap_zero_min(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::swap_zero_min()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@swap_zero_min.
      /// \param result The data expression where the \@swap_zero_min expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_swap_zero_min(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::swap_zero_min(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@swap_zero_min.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol swap_zero_min to a
      ///     number of arguments.
      inline
      bool is_swap_zero_min_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_swap_zero_min_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@swap_zero_monus.
      /// \return Identifier \@swap_zero_monus.
      inline
      const core::identifier_string& swap_zero_monus_name()
      {
        static core::identifier_string swap_zero_monus_name = core::identifier_string("@swap_zero_monus");
        return swap_zero_monus_name;
      }

      /// \brief Constructor for function symbol \@swap_zero_monus.
      
      /// \return Function symbol swap_zero_monus.
      inline
      const function_symbol& swap_zero_monus()
      {
        static function_symbol swap_zero_monus(swap_zero_monus_name(), make_function_sort_(nat(), nat(), nat(), nat(), nat()));
        return swap_zero_monus;
      }

      /// \brief Recogniser for function \@swap_zero_monus.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@swap_zero_monus.
      inline
      bool is_swap_zero_monus_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == swap_zero_monus();
        }
        return false;
      }

      /// \brief Application of function symbol \@swap_zero_monus.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@swap_zero_monus to a number of arguments.
      inline
      application swap_zero_monus(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::swap_zero_monus()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@swap_zero_monus.
      /// \param result The data expression where the \@swap_zero_monus expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_swap_zero_monus(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::swap_zero_monus(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@swap_zero_monus.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol swap_zero_monus to a
      ///     number of arguments.
      inline
      bool is_swap_zero_monus_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_swap_zero_monus_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_whr1.
      /// \return Identifier \@sqrt_whr1.
      inline
      const core::identifier_string& sqrt_whr1_name()
      {
        static core::identifier_string sqrt_whr1_name = core::identifier_string("@sqrt_whr1");
        return sqrt_whr1_name;
      }

      /// \brief Constructor for function symbol \@sqrt_whr1.
      
      /// \return Function symbol sqrt_whr1.
      inline
      const function_symbol& sqrt_whr1()
      {
        static function_symbol sqrt_whr1(sqrt_whr1_name(), make_function_sort_(sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), nat()));
        return sqrt_whr1;
      }

      /// \brief Recogniser for function \@sqrt_whr1.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_whr1.
      inline
      bool is_sqrt_whr1_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_whr1();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_whr1.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@sqrt_whr1 to a number of arguments.
      inline
      application sqrt_whr1(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::sqrt_whr1()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@sqrt_whr1.
      /// \param result The data expression where the \@sqrt_whr1 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_sqrt_whr1(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::sqrt_whr1(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@sqrt_whr1.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_whr1 to a
      ///     number of arguments.
      inline
      bool is_sqrt_whr1_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_whr1_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_whr2.
      /// \return Identifier \@sqrt_whr2.
      inline
      const core::identifier_string& sqrt_whr2_name()
      {
        static core::identifier_string sqrt_whr2_name = core::identifier_string("@sqrt_whr2");
        return sqrt_whr2_name;
      }

      /// \brief Constructor for function symbol \@sqrt_whr2.
      
      /// \return Function symbol sqrt_whr2.
      inline
      const function_symbol& sqrt_whr2()
      {
        static function_symbol sqrt_whr2(sqrt_whr2_name(), make_function_sort_(sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), nat()));
        return sqrt_whr2;
      }

      /// \brief Recogniser for function \@sqrt_whr2.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_whr2.
      inline
      bool is_sqrt_whr2_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_whr2();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_whr2.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@sqrt_whr2 to a number of arguments.
      inline
      application sqrt_whr2(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_nat::sqrt_whr2()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@sqrt_whr2.
      /// \param result The data expression where the \@sqrt_whr2 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_sqrt_whr2(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_nat::sqrt_whr2(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@sqrt_whr2.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_whr2 to a
      ///     number of arguments.
      inline
      bool is_sqrt_whr2_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_whr2_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_pair.
      /// \return Identifier \@sqrt_pair.
      inline
      const core::identifier_string& sqrt_pair_name()
      {
        static core::identifier_string sqrt_pair_name = core::identifier_string("@sqrt_pair");
        return sqrt_pair_name;
      }

      /// \brief Constructor for function symbol \@sqrt_pair.
      
      /// \return Function symbol sqrt_pair.
      inline
      const function_symbol& sqrt_pair()
      {
        static function_symbol sqrt_pair(sqrt_pair_name(), make_function_sort_(nat(), natnatpair()));
        return sqrt_pair;
      }

      /// \brief Recogniser for function \@sqrt_pair.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_pair.
      inline
      bool is_sqrt_pair_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_pair();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_pair.
      
      /// \param arg0 A data expression.
      /// \return Application of \@sqrt_pair to a number of arguments.
      inline
      application sqrt_pair(const data_expression& arg0)
      {
        return sort_nat::sqrt_pair()(arg0);
      }

      /// \brief Make an application of function symbol \@sqrt_pair.
      /// \param result The data expression where the \@sqrt_pair expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_sqrt_pair(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::sqrt_pair(),arg0);
      }

      /// \brief Recogniser for application of \@sqrt_pair.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_pair to a
      ///     number of arguments.
      inline
      bool is_sqrt_pair_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_pair_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_pair_whr1.
      /// \return Identifier \@sqrt_pair_whr1.
      inline
      const core::identifier_string& sqrt_pair_whr1_name()
      {
        static core::identifier_string sqrt_pair_whr1_name = core::identifier_string("@sqrt_pair_whr1");
        return sqrt_pair_whr1_name;
      }

      /// \brief Constructor for function symbol \@sqrt_pair_whr1.
      
      /// \return Function symbol sqrt_pair_whr1.
      inline
      const function_symbol& sqrt_pair_whr1()
      {
        static function_symbol sqrt_pair_whr1(sqrt_pair_whr1_name(), make_function_sort_(sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), nat(), natnatpair()));
        return sqrt_pair_whr1;
      }

      /// \brief Recogniser for function \@sqrt_pair_whr1.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_pair_whr1.
      inline
      bool is_sqrt_pair_whr1_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_pair_whr1();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_pair_whr1.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@sqrt_pair_whr1 to a number of arguments.
      inline
      application sqrt_pair_whr1(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_nat::sqrt_pair_whr1()(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@sqrt_pair_whr1.
      /// \param result The data expression where the \@sqrt_pair_whr1 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_sqrt_pair_whr1(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_nat::sqrt_pair_whr1(),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@sqrt_pair_whr1.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_pair_whr1 to a
      ///     number of arguments.
      inline
      bool is_sqrt_pair_whr1_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_pair_whr1_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_pair_whr2.
      /// \return Identifier \@sqrt_pair_whr2.
      inline
      const core::identifier_string& sqrt_pair_whr2_name()
      {
        static core::identifier_string sqrt_pair_whr2_name = core::identifier_string("@sqrt_pair_whr2");
        return sqrt_pair_whr2_name;
      }

      /// \brief Constructor for function symbol \@sqrt_pair_whr2.
      
      /// \return Function symbol sqrt_pair_whr2.
      inline
      const function_symbol& sqrt_pair_whr2()
      {
        static function_symbol sqrt_pair_whr2(sqrt_pair_whr2_name(), make_function_sort_(sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), sort_machine_word::machine_word(), nat(), natnatpair()));
        return sqrt_pair_whr2;
      }

      /// \brief Recogniser for function \@sqrt_pair_whr2.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_pair_whr2.
      inline
      bool is_sqrt_pair_whr2_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_pair_whr2();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_pair_whr2.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@sqrt_pair_whr2 to a number of arguments.
      inline
      application sqrt_pair_whr2(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_nat::sqrt_pair_whr2()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@sqrt_pair_whr2.
      /// \param result The data expression where the \@sqrt_pair_whr2 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_sqrt_pair_whr2(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_nat::sqrt_pair_whr2(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@sqrt_pair_whr2.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_pair_whr2 to a
      ///     number of arguments.
      inline
      bool is_sqrt_pair_whr2_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_pair_whr2_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_pair_whr3.
      /// \return Identifier \@sqrt_pair_whr3.
      inline
      const core::identifier_string& sqrt_pair_whr3_name()
      {
        static core::identifier_string sqrt_pair_whr3_name = core::identifier_string("@sqrt_pair_whr3");
        return sqrt_pair_whr3_name;
      }

      /// \brief Constructor for function symbol \@sqrt_pair_whr3.
      
      /// \return Function symbol sqrt_pair_whr3.
      inline
      const function_symbol& sqrt_pair_whr3()
      {
        static function_symbol sqrt_pair_whr3(sqrt_pair_whr3_name(), make_function_sort_(sort_machine_word::machine_word(), sort_machine_word::machine_word(), natnatpair(), natnatpair()));
        return sqrt_pair_whr3;
      }

      /// \brief Recogniser for function \@sqrt_pair_whr3.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_pair_whr3.
      inline
      bool is_sqrt_pair_whr3_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_pair_whr3();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_pair_whr3.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@sqrt_pair_whr3 to a number of arguments.
      inline
      application sqrt_pair_whr3(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_nat::sqrt_pair_whr3()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@sqrt_pair_whr3.
      /// \param result The data expression where the \@sqrt_pair_whr3 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_sqrt_pair_whr3(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::sqrt_pair_whr3(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@sqrt_pair_whr3.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_pair_whr3 to a
      ///     number of arguments.
      inline
      bool is_sqrt_pair_whr3_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_pair_whr3_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_pair_whr4.
      /// \return Identifier \@sqrt_pair_whr4.
      inline
      const core::identifier_string& sqrt_pair_whr4_name()
      {
        static core::identifier_string sqrt_pair_whr4_name = core::identifier_string("@sqrt_pair_whr4");
        return sqrt_pair_whr4_name;
      }

      /// \brief Constructor for function symbol \@sqrt_pair_whr4.
      
      /// \return Function symbol sqrt_pair_whr4.
      inline
      const function_symbol& sqrt_pair_whr4()
      {
        static function_symbol sqrt_pair_whr4(sqrt_pair_whr4_name(), make_function_sort_(nat(), sort_machine_word::machine_word(), natnatpair(), nat(), nat(), nat(), natnatpair()));
        return sqrt_pair_whr4;
      }

      /// \brief Recogniser for function \@sqrt_pair_whr4.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_pair_whr4.
      inline
      bool is_sqrt_pair_whr4_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_pair_whr4();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_pair_whr4.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \param arg5 A data expression.
      /// \return Application of \@sqrt_pair_whr4 to a number of arguments.
      inline
      application sqrt_pair_whr4(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4, const data_expression& arg5)
      {
        return sort_nat::sqrt_pair_whr4()(arg0, arg1, arg2, arg3, arg4, arg5);
      }

      /// \brief Make an application of function symbol \@sqrt_pair_whr4.
      /// \param result The data expression where the \@sqrt_pair_whr4 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \param arg5 A data expression.
      inline
      void make_sqrt_pair_whr4(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4, const data_expression& arg5)
      {
        make_application(result, sort_nat::sqrt_pair_whr4(),arg0, arg1, arg2, arg3, arg4, arg5);
      }

      /// \brief Recogniser for application of \@sqrt_pair_whr4.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_pair_whr4 to a
      ///     number of arguments.
      inline
      bool is_sqrt_pair_whr4_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_pair_whr4_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_pair_whr5.
      /// \return Identifier \@sqrt_pair_whr5.
      inline
      const core::identifier_string& sqrt_pair_whr5_name()
      {
        static core::identifier_string sqrt_pair_whr5_name = core::identifier_string("@sqrt_pair_whr5");
        return sqrt_pair_whr5_name;
      }

      /// \brief Constructor for function symbol \@sqrt_pair_whr5.
      
      /// \return Function symbol sqrt_pair_whr5.
      inline
      const function_symbol& sqrt_pair_whr5()
      {
        static function_symbol sqrt_pair_whr5(sqrt_pair_whr5_name(), make_function_sort_(natnatpair(), nat(), nat(), nat(), nat(), natnatpair()));
        return sqrt_pair_whr5;
      }

      /// \brief Recogniser for function \@sqrt_pair_whr5.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_pair_whr5.
      inline
      bool is_sqrt_pair_whr5_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_pair_whr5();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_pair_whr5.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      /// \return Application of \@sqrt_pair_whr5 to a number of arguments.
      inline
      application sqrt_pair_whr5(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        return sort_nat::sqrt_pair_whr5()(arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Make an application of function symbol \@sqrt_pair_whr5.
      /// \param result The data expression where the \@sqrt_pair_whr5 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \param arg4 A data expression.
      inline
      void make_sqrt_pair_whr5(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3, const data_expression& arg4)
      {
        make_application(result, sort_nat::sqrt_pair_whr5(),arg0, arg1, arg2, arg3, arg4);
      }

      /// \brief Recogniser for application of \@sqrt_pair_whr5.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_pair_whr5 to a
      ///     number of arguments.
      inline
      bool is_sqrt_pair_whr5_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_pair_whr5_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@sqrt_pair_whr6.
      /// \return Identifier \@sqrt_pair_whr6.
      inline
      const core::identifier_string& sqrt_pair_whr6_name()
      {
        static core::identifier_string sqrt_pair_whr6_name = core::identifier_string("@sqrt_pair_whr6");
        return sqrt_pair_whr6_name;
      }

      /// \brief Constructor for function symbol \@sqrt_pair_whr6.
      
      /// \return Function symbol sqrt_pair_whr6.
      inline
      const function_symbol& sqrt_pair_whr6()
      {
        static function_symbol sqrt_pair_whr6(sqrt_pair_whr6_name(), make_function_sort_(nat(), nat(), nat(), natnatpair()));
        return sqrt_pair_whr6;
      }

      /// \brief Recogniser for function \@sqrt_pair_whr6.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@sqrt_pair_whr6.
      inline
      bool is_sqrt_pair_whr6_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == sqrt_pair_whr6();
        }
        return false;
      }

      /// \brief Application of function symbol \@sqrt_pair_whr6.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@sqrt_pair_whr6 to a number of arguments.
      inline
      application sqrt_pair_whr6(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_nat::sqrt_pair_whr6()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@sqrt_pair_whr6.
      /// \param result The data expression where the \@sqrt_pair_whr6 expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_sqrt_pair_whr6(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_nat::sqrt_pair_whr6(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@sqrt_pair_whr6.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol sqrt_pair_whr6 to a
      ///     number of arguments.
      inline
      bool is_sqrt_pair_whr6_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_sqrt_pair_whr6_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@first.
      /// \return Identifier \@first.
      inline
      const core::identifier_string& first_name()
      {
        static core::identifier_string first_name = core::identifier_string("@first");
        return first_name;
      }

      /// \brief Constructor for function symbol \@first.
      
      /// \return Function symbol first.
      inline
      const function_symbol& first()
      {
        static function_symbol first(first_name(), make_function_sort_(natnatpair(), nat()));
        return first;
      }

      /// \brief Recogniser for function \@first.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@first.
      inline
      bool is_first_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == first();
        }
        return false;
      }

      /// \brief Application of function symbol \@first.
      
      /// \param arg0 A data expression.
      /// \return Application of \@first to a number of arguments.
      inline
      application first(const data_expression& arg0)
      {
        return sort_nat::first()(arg0);
      }

      /// \brief Make an application of function symbol \@first.
      /// \param result The data expression where the \@first expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_first(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::first(),arg0);
      }

      /// \brief Recogniser for application of \@first.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol first to a
      ///     number of arguments.
      inline
      bool is_first_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_first_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@last.
      /// \return Identifier \@last.
      inline
      const core::identifier_string& last_name()
      {
        static core::identifier_string last_name = core::identifier_string("@last");
        return last_name;
      }

      /// \brief Constructor for function symbol \@last.
      
      /// \return Function symbol last.
      inline
      const function_symbol& last()
      {
        static function_symbol last(last_name(), make_function_sort_(natnatpair(), nat()));
        return last;
      }

      /// \brief Recogniser for function \@last.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@last.
      inline
      bool is_last_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == last();
        }
        return false;
      }

      /// \brief Application of function symbol \@last.
      
      /// \param arg0 A data expression.
      /// \return Application of \@last to a number of arguments.
      inline
      application last(const data_expression& arg0)
      {
        return sort_nat::last()(arg0);
      }

      /// \brief Make an application of function symbol \@last.
      /// \param result The data expression where the \@last expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_last(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_nat::last(),arg0);
      }

      /// \brief Recogniser for application of \@last.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol last to a
      ///     number of arguments.
      inline
      bool is_last_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_last_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for nat
      /// \return All system defined mappings for nat
      inline
      function_symbol_vector nat_generate_functions_code()
      {
        function_symbol_vector result;
        result.push_back(sort_nat::most_significant_digit_nat());
        result.push_back(sort_nat::concat_digit(nat(), sort_machine_word::machine_word()));
        result.push_back(sort_nat::equals_zero());
        result.push_back(sort_nat::not_equals_zero());
        result.push_back(sort_nat::equals_one(nat()));
        result.push_back(sort_nat::pos2nat());
        result.push_back(sort_nat::nat2pos());
        result.push_back(sort_nat::succ(nat()));
        result.push_back(sort_nat::maximum(sort_pos::pos(), nat()));
        result.push_back(sort_nat::maximum(nat(), sort_pos::pos()));
        result.push_back(sort_nat::maximum(nat(), nat()));
        result.push_back(sort_nat::minimum(nat(), nat()));
        result.push_back(sort_nat::pred());
        result.push_back(sort_nat::pred_whr());
        result.push_back(sort_nat::plus(sort_pos::pos(), nat()));
        result.push_back(sort_nat::plus(nat(), sort_pos::pos()));
        result.push_back(sort_nat::plus(nat(), nat()));
        result.push_back(sort_nat::add_with_carry(nat(), nat()));
        result.push_back(sort_nat::auxiliary_plus_nat());
        result.push_back(sort_nat::times(nat(), nat()));
        result.push_back(sort_nat::times_ordered(nat(), nat()));
        result.push_back(sort_nat::times_overflow(nat(), sort_machine_word::machine_word(), sort_machine_word::machine_word()));
        result.push_back(sort_nat::div());
        result.push_back(sort_nat::mod());
        result.push_back(sort_nat::exp(sort_pos::pos(), nat()));
        result.push_back(sort_nat::exp(nat(), nat()));
        result.push_back(sort_nat::sqrt());
        result.push_back(sort_nat::natpred());
        result.push_back(sort_nat::is_odd());
        result.push_back(sort_nat::div2());
        result.push_back(sort_nat::monus());
        result.push_back(sort_nat::monus_whr());
        result.push_back(sort_nat::exp_aux3p());
        result.push_back(sort_nat::exp_aux4p());
        result.push_back(sort_nat::exp_aux3n());
        result.push_back(sort_nat::exp_aux4n());
        result.push_back(sort_nat::exp_auxtruep());
        result.push_back(sort_nat::exp_auxtruen());
        result.push_back(sort_nat::exp_auxfalsep());
        result.push_back(sort_nat::exp_auxfalsen());
        result.push_back(sort_nat::div_bold());
        result.push_back(sort_nat::div_bold_whr());
        result.push_back(sort_nat::div_whr1());
        result.push_back(sort_nat::div_whr2());
        result.push_back(sort_nat::mod_whr1());
        result.push_back(sort_nat::divmod_aux());
        result.push_back(sort_nat::divmod_aux_whr1());
        result.push_back(sort_nat::divmod_aux_whr2());
        result.push_back(sort_nat::divmod_aux_whr3());
        result.push_back(sort_nat::divmod_aux_whr4());
        result.push_back(sort_nat::divmod_aux_whr5());
        result.push_back(sort_nat::divmod_aux_whr6());
        result.push_back(sort_nat::msd());
        result.push_back(sort_nat::swap_zero());
        result.push_back(sort_nat::swap_zero_add());
        result.push_back(sort_nat::swap_zero_min());
        result.push_back(sort_nat::swap_zero_monus());
        result.push_back(sort_nat::sqrt_whr1());
        result.push_back(sort_nat::sqrt_whr2());
        result.push_back(sort_nat::sqrt_pair());
        result.push_back(sort_nat::sqrt_pair_whr1());
        result.push_back(sort_nat::sqrt_pair_whr2());
        result.push_back(sort_nat::sqrt_pair_whr3());
        result.push_back(sort_nat::sqrt_pair_whr4());
        result.push_back(sort_nat::sqrt_pair_whr5());
        result.push_back(sort_nat::sqrt_pair_whr6());
        result.push_back(sort_nat::first());
        result.push_back(sort_nat::last());
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for nat
      /// \return All system defined mappings for nat
      inline
      function_symbol_vector nat_generate_constructors_and_functions_code()
      {
        function_symbol_vector result=nat_generate_functions_code();
        for(const function_symbol& f: nat_generate_constructors_code())
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for nat
      /// \return All system defined mappings for that can be used in mCRL2 specificationis nat
      inline
      function_symbol_vector nat_mCRL2_usable_mappings()
      {
        function_symbol_vector result;
        result.push_back(sort_nat::most_significant_digit_nat());
        result.push_back(sort_nat::concat_digit(nat(), sort_machine_word::machine_word()));
        result.push_back(sort_nat::equals_zero());
        result.push_back(sort_nat::not_equals_zero());
        result.push_back(sort_nat::equals_one(nat()));
        result.push_back(sort_nat::pos2nat());
        result.push_back(sort_nat::nat2pos());
        result.push_back(sort_nat::succ(nat()));
        result.push_back(sort_nat::maximum(sort_pos::pos(), nat()));
        result.push_back(sort_nat::maximum(nat(), sort_pos::pos()));
        result.push_back(sort_nat::maximum(nat(), nat()));
        result.push_back(sort_nat::minimum(nat(), nat()));
        result.push_back(sort_nat::pred());
        result.push_back(sort_nat::pred_whr());
        result.push_back(sort_nat::plus(sort_pos::pos(), nat()));
        result.push_back(sort_nat::plus(nat(), sort_pos::pos()));
        result.push_back(sort_nat::plus(nat(), nat()));
        result.push_back(sort_nat::add_with_carry(nat(), nat()));
        result.push_back(sort_nat::auxiliary_plus_nat());
        result.push_back(sort_nat::times(nat(), nat()));
        result.push_back(sort_nat::times_ordered(nat(), nat()));
        result.push_back(sort_nat::times_overflow(nat(), sort_machine_word::machine_word(), sort_machine_word::machine_word()));
        result.push_back(sort_nat::div());
        result.push_back(sort_nat::mod());
        result.push_back(sort_nat::exp(sort_pos::pos(), nat()));
        result.push_back(sort_nat::exp(nat(), nat()));
        result.push_back(sort_nat::sqrt());
        result.push_back(sort_nat::natpred());
        result.push_back(sort_nat::is_odd());
        result.push_back(sort_nat::div2());
        result.push_back(sort_nat::monus());
        result.push_back(sort_nat::monus_whr());
        result.push_back(sort_nat::exp_aux3p());
        result.push_back(sort_nat::exp_aux4p());
        result.push_back(sort_nat::exp_aux3n());
        result.push_back(sort_nat::exp_aux4n());
        result.push_back(sort_nat::exp_auxtruep());
        result.push_back(sort_nat::exp_auxtruen());
        result.push_back(sort_nat::exp_auxfalsep());
        result.push_back(sort_nat::exp_auxfalsen());
        result.push_back(sort_nat::div_bold());
        result.push_back(sort_nat::div_bold_whr());
        result.push_back(sort_nat::div_whr1());
        result.push_back(sort_nat::div_whr2());
        result.push_back(sort_nat::mod_whr1());
        result.push_back(sort_nat::divmod_aux());
        result.push_back(sort_nat::divmod_aux_whr1());
        result.push_back(sort_nat::divmod_aux_whr2());
        result.push_back(sort_nat::divmod_aux_whr3());
        result.push_back(sort_nat::divmod_aux_whr4());
        result.push_back(sort_nat::divmod_aux_whr5());
        result.push_back(sort_nat::divmod_aux_whr6());
        result.push_back(sort_nat::msd());
        result.push_back(sort_nat::swap_zero());
        result.push_back(sort_nat::swap_zero_add());
        result.push_back(sort_nat::swap_zero_min());
        result.push_back(sort_nat::swap_zero_monus());
        result.push_back(sort_nat::sqrt_whr1());
        result.push_back(sort_nat::sqrt_whr2());
        result.push_back(sort_nat::sqrt_pair());
        result.push_back(sort_nat::sqrt_pair_whr1());
        result.push_back(sort_nat::sqrt_pair_whr2());
        result.push_back(sort_nat::sqrt_pair_whr3());
        result.push_back(sort_nat::sqrt_pair_whr4());
        result.push_back(sort_nat::sqrt_pair_whr5());
        result.push_back(sort_nat::sqrt_pair_whr6());
        result.push_back(sort_nat::first());
        result.push_back(sort_nat::last());
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for nat
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for nat
      inline
      implementation_map nat_cpp_implementable_mappings()
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
        assert(is_succ_nat_application(e) || is_most_significant_digit_nat_application(e) || is_equals_zero_application(e) || is_not_equals_zero_application(e) || is_equals_one_application(e) || is_pos2nat_application(e) || is_nat2pos_application(e) || is_succ_application(e) || is_pred_application(e) || is_pred_whr_application(e) || is_sqrt_application(e) || is_natpred_application(e) || is_is_odd_application(e) || is_div2_application(e) || is_msd_application(e) || is_sqrt_pair_application(e) || is_first_application(e) || is_last_application(e));
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
        assert(is_nnpair_application(e) || is_concat_digit_application(e) || is_times_overflow_application(e) || is_monus_whr_application(e) || is_exp_aux3p_application(e) || is_exp_aux4p_application(e) || is_exp_aux3n_application(e) || is_exp_aux4n_application(e) || is_exp_auxtruep_application(e) || is_exp_auxtruen_application(e) || is_exp_auxfalsep_application(e) || is_exp_auxfalsen_application(e) || is_div_bold_application(e) || is_div_bold_whr_application(e) || is_div_whr1_application(e) || is_div_whr2_application(e) || is_mod_whr1_application(e) || is_divmod_aux_application(e) || is_divmod_aux_whr1_application(e) || is_divmod_aux_whr2_application(e) || is_divmod_aux_whr3_application(e) || is_divmod_aux_whr4_application(e) || is_divmod_aux_whr5_application(e) || is_divmod_aux_whr6_application(e) || is_swap_zero_add_application(e) || is_swap_zero_min_application(e) || is_swap_zero_monus_application(e) || is_sqrt_whr1_application(e) || is_sqrt_whr2_application(e) || is_sqrt_pair_whr1_application(e) || is_sqrt_pair_whr2_application(e) || is_sqrt_pair_whr3_application(e) || is_sqrt_pair_whr4_application(e) || is_sqrt_pair_whr5_application(e) || is_sqrt_pair_whr6_application(e));
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
        assert(is_nnpair_application(e) || is_concat_digit_application(e) || is_times_overflow_application(e) || is_monus_whr_application(e) || is_exp_aux3p_application(e) || is_exp_aux4p_application(e) || is_exp_aux3n_application(e) || is_exp_aux4n_application(e) || is_exp_auxtruep_application(e) || is_exp_auxtruen_application(e) || is_exp_auxfalsep_application(e) || is_exp_auxfalsen_application(e) || is_div_bold_application(e) || is_div_bold_whr_application(e) || is_div_whr1_application(e) || is_div_whr2_application(e) || is_mod_whr1_application(e) || is_divmod_aux_application(e) || is_divmod_aux_whr1_application(e) || is_divmod_aux_whr2_application(e) || is_divmod_aux_whr3_application(e) || is_divmod_aux_whr4_application(e) || is_divmod_aux_whr5_application(e) || is_divmod_aux_whr6_application(e) || is_swap_zero_add_application(e) || is_swap_zero_min_application(e) || is_swap_zero_monus_application(e) || is_sqrt_whr1_application(e) || is_sqrt_whr2_application(e) || is_sqrt_pair_whr1_application(e) || is_sqrt_pair_whr2_application(e) || is_sqrt_pair_whr3_application(e) || is_sqrt_pair_whr4_application(e) || is_sqrt_pair_whr5_application(e) || is_sqrt_pair_whr6_application(e));
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
        assert(is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_add_with_carry_application(e) || is_auxiliary_plus_nat_application(e) || is_times_application(e) || is_times_ordered_application(e) || is_div_application(e) || is_mod_application(e) || is_exp_application(e) || is_monus_application(e) || is_swap_zero_application(e));
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
        assert(is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_add_with_carry_application(e) || is_auxiliary_plus_nat_application(e) || is_times_application(e) || is_times_ordered_application(e) || is_div_application(e) || is_mod_application(e) || is_exp_application(e) || is_monus_application(e) || is_swap_zero_application(e));
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
        assert(is_times_overflow_application(e) || is_monus_whr_application(e) || is_exp_aux3p_application(e) || is_exp_aux4p_application(e) || is_exp_aux3n_application(e) || is_exp_aux4n_application(e) || is_exp_auxtruep_application(e) || is_exp_auxtruen_application(e) || is_exp_auxfalsep_application(e) || is_exp_auxfalsen_application(e) || is_div_bold_whr_application(e) || is_div_whr1_application(e) || is_div_whr2_application(e) || is_mod_whr1_application(e) || is_divmod_aux_whr1_application(e) || is_divmod_aux_whr2_application(e) || is_divmod_aux_whr3_application(e) || is_divmod_aux_whr4_application(e) || is_divmod_aux_whr5_application(e) || is_divmod_aux_whr6_application(e) || is_swap_zero_add_application(e) || is_swap_zero_min_application(e) || is_swap_zero_monus_application(e) || is_sqrt_whr1_application(e) || is_sqrt_whr2_application(e) || is_sqrt_pair_whr1_application(e) || is_sqrt_pair_whr2_application(e) || is_sqrt_pair_whr3_application(e) || is_sqrt_pair_whr4_application(e) || is_sqrt_pair_whr5_application(e) || is_sqrt_pair_whr6_application(e));
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
        assert(is_monus_whr_application(e) || is_exp_aux4p_application(e) || is_exp_aux4n_application(e) || is_div_bold_whr_application(e) || is_div_whr1_application(e) || is_div_whr2_application(e) || is_mod_whr1_application(e) || is_divmod_aux_whr1_application(e) || is_divmod_aux_whr2_application(e) || is_divmod_aux_whr3_application(e) || is_divmod_aux_whr4_application(e) || is_divmod_aux_whr5_application(e) || is_divmod_aux_whr6_application(e) || is_swap_zero_add_application(e) || is_swap_zero_min_application(e) || is_swap_zero_monus_application(e) || is_sqrt_whr1_application(e) || is_sqrt_whr2_application(e) || is_sqrt_pair_whr1_application(e) || is_sqrt_pair_whr2_application(e) || is_sqrt_pair_whr4_application(e) || is_sqrt_pair_whr5_application(e));
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
        assert(is_monus_whr_application(e) || is_div_bold_whr_application(e) || is_div_whr2_application(e) || is_divmod_aux_whr2_application(e) || is_divmod_aux_whr3_application(e) || is_divmod_aux_whr6_application(e) || is_sqrt_whr2_application(e) || is_sqrt_pair_whr2_application(e) || is_sqrt_pair_whr4_application(e) || is_sqrt_pair_whr5_application(e));
        return atermpp::down_cast<application>(e)[4];
      }

      ///\brief Function for projecting out argument.
      ///        arg6 from an application.
      /// \param e A data expression.
      /// \pre arg6 is defined for e.
      /// \return The argument of e that corresponds to arg6.
      inline
      const data_expression& arg6(const data_expression& e)
      {
        assert(is_div_bold_whr_application(e) || is_sqrt_pair_whr4_application(e));
        return atermpp::down_cast<application>(e)[5];
      }

      /// \brief Give all system defined equations for nat
      /// \return All system defined equations for sort nat
      inline
      data_equation_vector nat_generate_equations_code()
      {
        variable vb("b",sort_bool::bool_());
        variable vp("p",sort_pos::pos());
        variable vp1("p1",sort_pos::pos());
        variable vp2("p2",sort_pos::pos());
        variable vn("n",nat());
        variable vn1("n1",nat());
        variable vn2("n2",nat());
        variable vm("m",nat());
        variable vm1("m1",nat());
        variable vm2("m2",nat());
        variable vm3("m3",nat());
        variable vm4("m4",nat());
        variable vm5("m5",nat());
        variable vpredp("predp",nat());
        variable vdiff("diff",nat());
        variable vshift_n1("shift_n1",nat());
        variable vsolution("solution",nat());
        variable vpq("pq",nat());
        variable vy("y",nat());
        variable vy_guess("y_guess",nat());
        variable vpair_("pair_",natnatpair());
        variable vlp("lp",nat());
        variable vw("w",sort_machine_word::machine_word());
        variable vw1("w1",sort_machine_word::machine_word());
        variable vw2("w2",sort_machine_word::machine_word());
        variable vw3("w3",sort_machine_word::machine_word());
        variable vw4("w4",sort_machine_word::machine_word());
        variable vshift_w("shift_w",sort_machine_word::machine_word());
        variable voverflow("overflow",sort_machine_word::machine_word());

        data_equation_vector result;
        result.push_back(data_equation(variable_list(), c0(), most_significant_digit_nat(sort_machine_word::zero_word())));
        result.push_back(data_equation(variable_list({vw}), equals_zero(most_significant_digit_nat(vw)), sort_machine_word::equals_zero_word(vw)));
        result.push_back(data_equation(variable_list({vn, vw}), equals_zero(concat_digit(vn, vw)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vn}), equals_zero(succ_nat(vn)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vw}), not_equals_zero(most_significant_digit_nat(vw)), sort_machine_word::not_equals_zero_word(vw)));
        result.push_back(data_equation(variable_list({vn, vw}), not_equals_zero(concat_digit(vn, vw)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vn}), not_equals_zero(succ_nat(vn)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vw}), equals_one(most_significant_digit_nat(vw)), sort_machine_word::equals_one_word(vw)));
        result.push_back(data_equation(variable_list({vn, vw}), equals_one(concat_digit(vn, vw)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vw}), succ_nat(most_significant_digit_nat(vw)), if_(sort_machine_word::equals_max_word(vw), concat_digit(most_significant_digit_nat(sort_machine_word::one_word()), sort_machine_word::zero_word()), most_significant_digit_nat(sort_machine_word::succ_word(vw)))));
        result.push_back(data_equation(variable_list({vn, vw}), succ_nat(concat_digit(vn, vw)), if_(sort_machine_word::equals_max_word(vw), concat_digit(succ_nat(vn), sort_machine_word::zero_word()), concat_digit(vn, sort_machine_word::succ_word(vw)))));
        result.push_back(data_equation(variable_list({vw}), succ(most_significant_digit_nat(vw)), if_(sort_machine_word::equals_max_word(vw), concat_digit(sort_pos::most_significant_digit(sort_machine_word::one_word()), sort_machine_word::zero_word()), sort_pos::most_significant_digit(sort_machine_word::succ_word(vw)))));
        result.push_back(data_equation(variable_list({vn, vw}), succ(concat_digit(vn, vw)), if_(sort_machine_word::equals_max_word(vw), concat_digit(succ(vn), sort_machine_word::zero_word()), concat_digit(nat2pos(vn), sort_machine_word::succ_word(vw)))));
        result.push_back(data_equation(variable_list({vw1, vw2}), equal_to(most_significant_digit_nat(vw1), most_significant_digit_nat(vw2)), sort_machine_word::equal_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), equal_to(concat_digit(vn, vw1), most_significant_digit_nat(vw2)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), equal_to(most_significant_digit_nat(vw1), concat_digit(vn, vw2)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vn1, vn2, vw1, vw2}), equal_to(concat_digit(vn1, vw1), concat_digit(vn2, vw2)), sort_bool::and_(sort_machine_word::equal_word(vw1, vw2), equal_to(vn1, vn2))));
        result.push_back(data_equation(variable_list({vn1, vn2}), equal_to(succ_nat(vn1), vn2), sort_bool::and_(not_equals_zero(vn2), equal_to(vn1, natpred(vn2)))));
        result.push_back(data_equation(variable_list({vn1, vn2}), equal_to(vn1, succ_nat(vn2)), sort_bool::and_(not_equals_zero(vn1), equal_to(natpred(vn1), vn2))));
        result.push_back(data_equation(variable_list({vw1, vw2}), less(most_significant_digit_nat(vw1), most_significant_digit_nat(vw2)), sort_machine_word::less_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), less(concat_digit(vn, vw1), most_significant_digit_nat(vw2)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), less(most_significant_digit_nat(vw1), concat_digit(vn, vw2)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vn1, vn2, vw1, vw2}), less(concat_digit(vn1, vw1), concat_digit(vn2, vw2)), if_(sort_machine_word::less_word(vw1, vw2), less_equal(vn1, vn2), less(vn1, vn2))));
        result.push_back(data_equation(variable_list({vn1, vn2}), less(succ_nat(vn1), vn2), sort_bool::and_(less(most_significant_digit_nat(sort_machine_word::one_word()), vn2), less(vn1, natpred(vn2)))));
        result.push_back(data_equation(variable_list({vn1, vn2}), less(vn1, succ_nat(vn2)), less_equal(vn1, vn2)));
        result.push_back(data_equation(variable_list({vn, vw1}), sort_machine_word::equals_zero_word(vw1), less(vn, most_significant_digit_nat(vw1)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vw1, vw2}), less_equal(most_significant_digit_nat(vw1), most_significant_digit_nat(vw2)), sort_machine_word::less_equal_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), less_equal(concat_digit(vn, vw1), most_significant_digit_nat(vw2)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), less_equal(most_significant_digit_nat(vw1), concat_digit(vn, vw2)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vn1, vn2, vw1, vw2}), less_equal(concat_digit(vn1, vw1), concat_digit(vn2, vw2)), if_(sort_machine_word::less_equal_word(vw1, vw2), less_equal(vn1, vn2), less(vn1, vn2))));
        result.push_back(data_equation(variable_list({vn1, vn2}), less_equal(succ_nat(vn1), vn2), less(vn1, vn2)));
        result.push_back(data_equation(variable_list({vn1, vn2}), less_equal(vn1, succ_nat(vn2)), less_equal(natpred(vn1), vn2)));
        result.push_back(data_equation(variable_list({vn, vw1}), sort_machine_word::equals_zero_word(vw1), less_equal(most_significant_digit_nat(vw1), vn), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vw}), pos2nat(sort_pos::most_significant_digit(vw)), most_significant_digit_nat(vw)));
        result.push_back(data_equation(variable_list({vp, vw}), pos2nat(concat_digit(vp, vw)), concat_digit(pos2nat(vp), vw)));
        result.push_back(data_equation(variable_list({vp}), pos2nat(sort_pos::succpos(vp)), succ_nat(pos2nat(vp))));
        result.push_back(data_equation(variable_list({vw}), sort_machine_word::not_equals_zero_word(vw), nat2pos(most_significant_digit_nat(vw)), sort_pos::most_significant_digit(vw)));
        result.push_back(data_equation(variable_list({vn, vw}), nat2pos(concat_digit(vn, vw)), concat_digit(nat2pos(vn), vw)));
        result.push_back(data_equation(variable_list({vn, vp}), maximum(vp, vn), if_(less_equal(vn, pos2nat(vp)), vp, nat2pos(vn))));
        result.push_back(data_equation(variable_list({vn, vp}), maximum(vn, vp), if_(less_equal(vn, pos2nat(vp)), vp, nat2pos(vn))));
        result.push_back(data_equation(variable_list({vm, vn}), maximum(vm, vn), if_(less_equal(vm, vn), vn, vm)));
        result.push_back(data_equation(variable_list({vm, vn}), minimum(vm, vn), if_(less_equal(vm, vn), vm, vn)));
        result.push_back(data_equation(variable_list({vw}), pred(sort_pos::most_significant_digit(vw)), most_significant_digit_nat(sort_machine_word::pred_word(vw))));
        result.push_back(data_equation(variable_list({vp, vw}), pred(concat_digit(vp, vw)), if_(sort_machine_word::equals_zero_word(vw), pred_whr(pred(vp)), concat_digit(pos2nat(vp), sort_machine_word::pred_word(vw)))));
        result.push_back(data_equation(variable_list({vpredp}), pred_whr(vpredp), if_(equals_zero(vpredp), most_significant_digit_nat(sort_machine_word::max_word()), concat_digit(vpredp, sort_machine_word::max_word()))));
        result.push_back(data_equation(variable_list({vn, vp}), plus(vp, vn), plus(vn, vp)));
        result.push_back(data_equation(variable_list({vw1, vw2}), plus(most_significant_digit_nat(vw1), sort_pos::most_significant_digit(vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(sort_pos::most_significant_digit(sort_machine_word::one_word()), sort_machine_word::add_word(vw1, vw2)), sort_pos::most_significant_digit(sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vw1, vw2}), plus(concat_digit(vn1, vw1), sort_pos::most_significant_digit(vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(succ(vn1), sort_machine_word::add_word(vw1, vw2)), concat_digit(nat2pos(vn1), sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), plus(most_significant_digit_nat(vw1), concat_digit(vp, vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(sort_pos::succpos(vp), sort_machine_word::add_word(vw1, vw2)), concat_digit(vp, sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vp, vw1, vw2}), plus(concat_digit(vn1, vw1), concat_digit(vp, vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(plus(succ(vn1), vp), sort_machine_word::add_word(vw1, vw2)), concat_digit(plus(vn1, vp), sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vw1, vw2}), plus(most_significant_digit_nat(vw1), most_significant_digit_nat(vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(most_significant_digit_nat(sort_machine_word::one_word()), sort_machine_word::add_word(vw1, vw2)), most_significant_digit_nat(sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vw1, vw2}), add_with_carry(most_significant_digit_nat(vw1), most_significant_digit_nat(vw2)), if_(sort_machine_word::add_with_carry_overflow_word(vw1, vw2), concat_digit(most_significant_digit_nat(sort_machine_word::one_word()), sort_machine_word::add_with_carry_word(vw1, vw2)), most_significant_digit_nat(sort_machine_word::add_with_carry_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vw1, vw2}), plus(concat_digit(vn1, vw1), most_significant_digit_nat(vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(succ_nat(vn1), sort_machine_word::add_word(vw1, vw2)), concat_digit(vn1, sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vw1, vw2}), add_with_carry(concat_digit(vn1, vw1), most_significant_digit_nat(vw2)), if_(sort_machine_word::add_with_carry_overflow_word(vw1, vw2), concat_digit(succ_nat(vn1), sort_machine_word::add_with_carry_word(vw1, vw2)), concat_digit(vn1, sort_machine_word::add_with_carry_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn2, vw1, vw2}), plus(most_significant_digit_nat(vw1), concat_digit(vn2, vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(succ_nat(vn2), sort_machine_word::add_word(vw1, vw2)), concat_digit(vn2, sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn2, vw1, vw2}), add_with_carry(most_significant_digit_nat(vw1), concat_digit(vn2, vw2)), if_(sort_machine_word::add_with_carry_overflow_word(vw1, vw2), concat_digit(succ_nat(vn2), sort_machine_word::add_with_carry_word(vw1, vw2)), concat_digit(vn2, sort_machine_word::add_with_carry_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vn2, vw1, vw2}), plus(concat_digit(vn1, vw1), concat_digit(vn2, vw2)), if_(sort_machine_word::add_overflow_word(vw1, vw2), concat_digit(add_with_carry(vn1, vn2), sort_machine_word::add_word(vw1, vw2)), concat_digit(plus(vn1, vn2), sort_machine_word::add_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vn2, vw1, vw2}), add_with_carry(concat_digit(vn1, vw1), concat_digit(vn2, vw2)), if_(sort_machine_word::add_with_carry_overflow_word(vw1, vw2), concat_digit(add_with_carry(vn1, vn2), sort_machine_word::add_with_carry_word(vw1, vw2)), concat_digit(plus(vn1, vn2), sort_machine_word::add_with_carry_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vn2}), plus(succ_nat(vn1), vn2), succ_nat(plus(vn1, vn2))));
        result.push_back(data_equation(variable_list({vn1, vn2}), plus(vn1, succ_nat(vn2)), succ_nat(plus(vn1, vn2))));
        result.push_back(data_equation(variable_list({vn1, vp2}), plus(succ_nat(vn1), vp2), sort_pos::succpos(plus(vn1, vp2))));
        result.push_back(data_equation(variable_list({vn1, vp2}), plus(vn1, sort_pos::succpos(vp2)), sort_pos::succpos(plus(vn1, vp2))));
        result.push_back(data_equation(variable_list({vn2, vp1}), plus(sort_pos::succpos(vp1), vn2), sort_pos::succpos(plus(vp1, vn2))));
        result.push_back(data_equation(variable_list({vn2, vp1}), plus(vp1, succ_nat(vn2)), sort_pos::succpos(plus(vp1, vn2))));
        result.push_back(data_equation(variable_list({vn1, vn2}), auxiliary_plus_nat(vn1, vn2), plus(vn1, vn2)));
        result.push_back(data_equation(variable_list({vw}), natpred(most_significant_digit_nat(vw)), if_(sort_machine_word::equals_zero_word(vw), most_significant_digit_nat(sort_machine_word::zero_word()), most_significant_digit_nat(sort_machine_word::pred_word(vw)))));
        result.push_back(data_equation(variable_list({vn, vw}), natpred(concat_digit(vn, vw)), if_(sort_machine_word::equals_zero_word(vw), if_(equals_one(vn), most_significant_digit_nat(sort_machine_word::max_word()), concat_digit(natpred(vn), sort_machine_word::max_word())), concat_digit(vn, sort_machine_word::pred_word(vw)))));
        result.push_back(data_equation(variable_list({vn}), natpred(succ_nat(vn)), vn));
        result.push_back(data_equation(variable_list({vw1, vw2}), monus(most_significant_digit_nat(vw1), most_significant_digit_nat(vw2)), most_significant_digit_nat(sort_machine_word::monus_word(vw1, vw2))));
        result.push_back(data_equation(variable_list({vn1, vw1, vw2}), monus(concat_digit(vn1, vw1), most_significant_digit_nat(vw2)), if_(sort_machine_word::less_word(vw1, vw2), if_(equals_one(vn1), most_significant_digit_nat(sort_machine_word::minus_word(vw1, vw2)), concat_digit(natpred(vn1), sort_machine_word::minus_word(vw1, vw2))), concat_digit(vn1, sort_machine_word::minus_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn2, vw1, vw2}), monus(most_significant_digit_nat(vw1), concat_digit(vn2, vw2)), most_significant_digit_nat(sort_machine_word::zero_word())));
        result.push_back(data_equation(variable_list({vn1, vn2, vw1, vw2}), monus(concat_digit(vn1, vw1), concat_digit(vn2, vw2)), monus_whr(vn1, vw1, vn2, vw2, monus(vn1, vn2))));
        result.push_back(data_equation(variable_list({vdiff, vn1, vn2, vw1, vw2}), monus_whr(vn1, vw1, vn2, vw2, vdiff), if_(sort_machine_word::less_word(vw1, vw2), if_(equals_zero(vdiff), most_significant_digit_nat(sort_machine_word::zero_word()), if_(equals_one(vdiff), most_significant_digit_nat(sort_machine_word::minus_word(vw1, vw2)), concat_digit(natpred(vdiff), sort_machine_word::minus_word(vw1, vw2)))), if_(equals_zero(vdiff), most_significant_digit_nat(sort_machine_word::minus_word(vw1, vw2)), concat_digit(vdiff, sort_machine_word::minus_word(vw1, vw2))))));
        result.push_back(data_equation(variable_list({vw1, vw2}), times(most_significant_digit_nat(vw1), most_significant_digit_nat(vw2)), if_(sort_machine_word::equals_zero_word(sort_machine_word::times_overflow_word(vw1, vw2)), most_significant_digit_nat(sort_machine_word::times_word(vw1, vw2)), concat_digit(most_significant_digit_nat(sort_machine_word::times_overflow_word(vw1, vw2)), sort_machine_word::times_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn2, vw1, vw2}), times(most_significant_digit_nat(vw1), concat_digit(vn2, vw2)), if_(sort_machine_word::equals_zero_word(vw1), most_significant_digit_nat(sort_machine_word::zero_word()), concat_digit(times_overflow(vn2, vw1, sort_machine_word::times_overflow_word(vw1, vw2)), sort_machine_word::times_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vw1, vw2}), times(concat_digit(vn1, vw1), most_significant_digit_nat(vw2)), if_(sort_machine_word::equals_zero_word(vw2), most_significant_digit_nat(sort_machine_word::zero_word()), concat_digit(times_overflow(vn1, vw2, sort_machine_word::times_overflow_word(vw1, vw2)), sort_machine_word::times_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn1, vn2, vw1, vw2}), times(concat_digit(vn1, vw1), concat_digit(vn2, vw2)), if_(less(vn1, vn2), times_ordered(concat_digit(vn1, vw1), concat_digit(vn2, vw2)), times_ordered(concat_digit(vn2, vw2), concat_digit(vn1, vw1)))));
        result.push_back(data_equation(variable_list({vn2, vw1, vw2}), times_ordered(most_significant_digit_nat(vw1), concat_digit(vn2, vw2)), concat_digit(times_overflow(vn2, vw1, sort_machine_word::times_overflow_word(vw1, vw2)), sort_machine_word::times_word(vw1, vw2))));
        result.push_back(data_equation(variable_list({vn1, vn2, vw1}), times_ordered(concat_digit(vn1, vw1), vn2), plus(concat_digit(times_ordered(vn1, vn2), sort_machine_word::zero_word()), times_overflow(vn2, vw1, sort_machine_word::zero_word()))));
        result.push_back(data_equation(variable_list({voverflow, vw1, vw2}), times_overflow(most_significant_digit_nat(vw1), vw2, voverflow), if_(sort_machine_word::equals_zero_word(sort_machine_word::times_with_carry_overflow_word(vw1, vw2, voverflow)), most_significant_digit_nat(sort_machine_word::times_with_carry_word(vw1, vw2, voverflow)), concat_digit(most_significant_digit_nat(sort_machine_word::times_with_carry_overflow_word(vw1, vw2, voverflow)), sort_machine_word::times_with_carry_word(vw1, vw2, voverflow)))));
        result.push_back(data_equation(variable_list({vn1, voverflow, vw1, vw2}), times_overflow(concat_digit(vn1, vw1), vw2, voverflow), if_(sort_machine_word::equals_zero_word(vw2), most_significant_digit_nat(voverflow), concat_digit(times_overflow(vn1, vw2, sort_machine_word::times_with_carry_overflow_word(vw1, vw2, voverflow)), sort_machine_word::times_with_carry_word(vw1, vw2, voverflow)))));
        result.push_back(data_equation(variable_list({vw}), is_odd(most_significant_digit_nat(vw)), sort_machine_word::rightmost_bit(vw)));
        result.push_back(data_equation(variable_list({vn, vw}), is_odd(concat_digit(vn, vw)), sort_machine_word::rightmost_bit(vw)));
        result.push_back(data_equation(variable_list({vw}), div2(most_significant_digit_nat(vw)), most_significant_digit_nat(sort_machine_word::shift_right(sort_bool::false_(), vw))));
        result.push_back(data_equation(variable_list({vn, vw}), div2(concat_digit(vn, vw)), if_(equals_zero(vn), most_significant_digit_nat(sort_machine_word::shift_right(is_odd(vn), vw)), concat_digit(div2(vn), sort_machine_word::shift_right(is_odd(vn), vw)))));
        result.push_back(data_equation(variable_list({vw}), msd(most_significant_digit_nat(vw)), vw));
        result.push_back(data_equation(variable_list({vn, vw}), msd(concat_digit(vn, vw)), msd(vn)));
        result.push_back(data_equation(variable_list({vn, vw}), exp(vn, most_significant_digit_nat(vw)), exp_aux3n(sort_machine_word::rightmost_bit(vw), vn, vw)));
        result.push_back(data_equation(variable_list({vn, vn1, vw1}), exp(vn, concat_digit(vn1, vw1)), exp_aux4n(sort_machine_word::rightmost_bit(vw1), vn, vn1, vw1)));
        result.push_back(data_equation(variable_list({vn, vw}), exp_aux3n(sort_bool::true_(), vn, vw), if_(sort_machine_word::equals_one_word(vw), vn, times(vn, exp_aux3n(sort_machine_word::rightmost_bit(sort_machine_word::shift_right(sort_bool::false_(), vw)), times(vn, vn), sort_machine_word::shift_right(sort_bool::false_(), vw))))));
        result.push_back(data_equation(variable_list({vn, vw}), exp_aux3n(sort_bool::false_(), vn, vw), if_(sort_machine_word::equals_zero_word(vw), most_significant_digit_nat(sort_machine_word::one_word()), exp_aux3n(sort_machine_word::rightmost_bit(sort_machine_word::shift_right(sort_bool::false_(), vw)), times(vn, vn), sort_machine_word::shift_right(sort_bool::false_(), vw)))));
        result.push_back(data_equation(variable_list({vn, vn1, vw}), exp_aux4n(sort_bool::true_(), vn, vn1, vw), exp_auxtruen(vn, div2(vn1), sort_machine_word::shift_right(is_odd(vn1), vw))));
        result.push_back(data_equation(variable_list({vn, vshift_n1, vshift_w}), exp_auxtruen(vn, vshift_n1, vshift_w), if_(equals_zero(vshift_n1), times(vn, exp_aux3n(sort_machine_word::rightmost_bit(vshift_w), times(vn, vn), vshift_w)), times(vn, exp_aux4n(sort_machine_word::rightmost_bit(vshift_w), times(vn, vn), vshift_n1, vshift_w)))));
        result.push_back(data_equation(variable_list({vn, vn1, vw}), exp_aux4n(sort_bool::false_(), vn, vn1, vw), exp_auxfalsen(vn, div2(vn1), sort_machine_word::shift_right(is_odd(vn1), vw))));
        result.push_back(data_equation(variable_list({vn, vshift_n1, vshift_w}), exp_auxfalsen(vn, vshift_n1, vshift_w), if_(equals_zero(vshift_n1), exp_aux3n(sort_machine_word::rightmost_bit(vshift_w), times(vn, vn), vshift_w), exp_aux4n(sort_machine_word::rightmost_bit(vshift_w), times(vn, vn), vshift_n1, vshift_w))));
        result.push_back(data_equation(variable_list({vp, vw}), exp(vp, most_significant_digit_nat(vw)), exp_aux3p(sort_machine_word::rightmost_bit(vw), vp, vw)));
        result.push_back(data_equation(variable_list({vn1, vp, vw1}), exp(vp, concat_digit(vn1, vw1)), exp_aux4p(sort_machine_word::rightmost_bit(vw1), vp, vn1, vw1)));
        result.push_back(data_equation(variable_list({vp, vw}), exp_aux3p(sort_bool::true_(), vp, vw), if_(sort_machine_word::equals_one_word(vw), vp, times(vp, exp_aux3p(sort_machine_word::rightmost_bit(sort_machine_word::shift_right(sort_bool::false_(), vw)), times(vp, vp), sort_machine_word::shift_right(sort_bool::false_(), vw))))));
        result.push_back(data_equation(variable_list({vp, vw}), exp_aux3p(sort_bool::false_(), vp, vw), if_(sort_machine_word::equals_zero_word(vw), sort_pos::most_significant_digit(sort_machine_word::one_word()), exp_aux3p(sort_machine_word::rightmost_bit(sort_machine_word::shift_right(sort_bool::false_(), vw)), times(vp, vp), sort_machine_word::shift_right(sort_bool::false_(), vw)))));
        result.push_back(data_equation(variable_list({vn1, vp, vw}), exp_aux4p(sort_bool::true_(), vp, vn1, vw), exp_auxtruep(vp, div2(vn1), sort_machine_word::shift_right(is_odd(vn1), vw))));
        result.push_back(data_equation(variable_list({vp, vshift_n1, vshift_w}), exp_auxtruep(vp, vshift_n1, vshift_w), if_(equals_zero(vshift_n1), times(vp, exp_aux3p(sort_machine_word::rightmost_bit(vshift_w), times(vp, vp), vshift_w)), times(vp, exp_aux4p(sort_machine_word::rightmost_bit(vshift_w), times(vp, vp), vshift_n1, vshift_w)))));
        result.push_back(data_equation(variable_list({vn1, vp, vw}), exp_aux4p(sort_bool::false_(), vp, vn1, vw), exp_auxfalsep(vp, div2(vn1), sort_machine_word::shift_right(is_odd(vn1), vw))));
        result.push_back(data_equation(variable_list({vp, vshift_n1, vshift_w}), exp_auxfalsep(vp, vshift_n1, vshift_w), if_(equals_zero(vshift_n1), exp_aux3p(sort_machine_word::rightmost_bit(vshift_w), times(vp, vp), vshift_w), exp_aux4p(sort_machine_word::rightmost_bit(vshift_w), times(vp, vp), vshift_n1, vshift_w))));
        result.push_back(data_equation(variable_list({vw1, vw2}), div(most_significant_digit_nat(vw1), sort_pos::most_significant_digit(vw2)), most_significant_digit_nat(sort_machine_word::div_word(vw1, vw2))));
        result.push_back(data_equation(variable_list({vw1, vw2}), mod(most_significant_digit_nat(vw1), sort_pos::most_significant_digit(vw2)), most_significant_digit_nat(sort_machine_word::mod_word(vw1, vw2))));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), div(most_significant_digit_nat(vw1), concat_digit(vp, vw2)), most_significant_digit_nat(sort_machine_word::zero_word())));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), mod(most_significant_digit_nat(vw1), concat_digit(vp, vw2)), most_significant_digit_nat(vw1)));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), div(concat_digit(vn, vw1), sort_pos::most_significant_digit(vw2)), div_whr1(vn, vw1, vw2, divmod_aux(vn, sort_pos::most_significant_digit(vw2)))));
        result.push_back(data_equation(variable_list({vn, vpair_, vw1, vw2}), div_whr1(vn, vw1, vw2, vpair_), if_(less(vn, most_significant_digit_nat(vw2)), most_significant_digit_nat(div_bold(concat_digit(vn, vw1), sort_pos::most_significant_digit(vw2))), if_(equals_zero(first(vpair_)), most_significant_digit_nat(div_bold(if_(equals_zero(last(vpair_)), most_significant_digit_nat(vw1), concat_digit(last(vpair_), vw1)), sort_pos::most_significant_digit(vw2))), concat_digit(first(vpair_), div_bold(if_(equals_zero(last(vpair_)), most_significant_digit_nat(vw1), concat_digit(last(vpair_), vw1)), sort_pos::most_significant_digit(vw2)))))));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), mod(concat_digit(vn, vw1), sort_pos::most_significant_digit(vw2)), most_significant_digit_nat(sort_machine_word::mod_doubleword(msd(mod(vn, sort_pos::most_significant_digit(vw2))), vw1, vw2))));
        result.push_back(data_equation(variable_list({vn, vp, vw1, vw2}), div(concat_digit(vn, vw1), concat_digit(vp, vw2)), if_(less(vn, pos2nat(concat_digit(vp, vw2))), most_significant_digit_nat(div_bold(concat_digit(vn, vw1), concat_digit(vp, vw2))), div_whr2(vn, vw1, vp, vw2, divmod_aux(vn, concat_digit(vp, vw2))))));
        result.push_back(data_equation(variable_list({vn, vp, vpair_, vw1, vw2}), div_whr2(vn, vw1, vp, vw2, vpair_), plus(if_(equals_zero(first(vpair_)), most_significant_digit_nat(sort_machine_word::zero_word()), concat_digit(first(vpair_), sort_machine_word::zero_word())), most_significant_digit_nat(div_bold(if_(equals_zero(last(vpair_)), most_significant_digit_nat(vw1), concat_digit(last(vpair_), vw1)), concat_digit(vp, vw2))))));
        result.push_back(data_equation(variable_list({vn, vp, vw1, vw2}), mod(concat_digit(vn, vw1), concat_digit(vp, vw2)), mod_whr1(vw1, vp, vw2, mod(vn, concat_digit(vp, vw2)))));
        result.push_back(data_equation(variable_list({vm1, vp, vw1, vw2}), mod_whr1(vw1, vp, vw2, vm1), monus(if_(less(most_significant_digit_nat(sort_machine_word::zero_word()), vm1), concat_digit(vm1, vw1), most_significant_digit_nat(vw1)), times(concat_digit(pos2nat(vp), vw2), most_significant_digit_nat(div_bold(if_(less(most_significant_digit_nat(sort_machine_word::zero_word()), vm1), concat_digit(vm1, vw1), most_significant_digit_nat(vw1)), concat_digit(vp, vw2)))))));
        result.push_back(data_equation(variable_list({vw1, vw2}), divmod_aux(most_significant_digit_nat(vw1), sort_pos::most_significant_digit(vw2)), nnpair(most_significant_digit_nat(sort_machine_word::div_word(vw1, vw2)), most_significant_digit_nat(sort_machine_word::mod_word(vw1, vw2)))));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), divmod_aux(most_significant_digit_nat(vw1), concat_digit(vp, vw2)), nnpair(most_significant_digit_nat(sort_machine_word::zero_word()), most_significant_digit_nat(vw1))));
        result.push_back(data_equation(variable_list({vn, vw1, vw2}), divmod_aux(concat_digit(vn, vw1), sort_pos::most_significant_digit(vw2)), divmod_aux_whr1(vn, vw1, vw2, divmod_aux(vn, sort_pos::most_significant_digit(vw2)))));
        result.push_back(data_equation(variable_list({vn, vpair_, vw1, vw2}), divmod_aux_whr1(vn, vw1, vw2, vpair_), nnpair(if_(less(vn, most_significant_digit_nat(vw2)), most_significant_digit_nat(div_bold(concat_digit(vn, vw1), sort_pos::most_significant_digit(vw2))), if_(equals_zero(first(vpair_)), most_significant_digit_nat(div_bold(if_(equals_zero(last(vpair_)), most_significant_digit_nat(vw1), concat_digit(last(vpair_), vw1)), sort_pos::most_significant_digit(vw2))), concat_digit(first(vpair_), div_bold(if_(equals_zero(last(vpair_)), most_significant_digit_nat(vw1), concat_digit(last(vpair_), vw1)), sort_pos::most_significant_digit(vw2))))), most_significant_digit_nat(sort_machine_word::mod_doubleword(msd(last(vpair_)), vw1, vw2)))));
        result.push_back(data_equation(variable_list({vn, vp, vw1, vw2}), divmod_aux(concat_digit(vn, vw1), concat_digit(vp, vw2)), if_(less(vn, pos2nat(concat_digit(vp, vw2))), divmod_aux_whr2(vn, vw1, vp, vw2, mod(vn, concat_digit(vp, vw2))), divmod_aux_whr4(vw1, vp, vw2, divmod_aux(vn, concat_digit(vp, vw2))))));
        result.push_back(data_equation(variable_list({vlp, vn, vp, vw1, vw2}), divmod_aux_whr2(vn, vw1, vp, vw2, vlp), divmod_aux_whr3(vn, vw1, vp, vw2, if_(equals_zero(vlp), most_significant_digit_nat(vw1), concat_digit(vlp, vw1)))));
        result.push_back(data_equation(variable_list({vm, vn, vp, vw1, vw2}), divmod_aux_whr3(vn, vw1, vp, vw2, vm), nnpair(most_significant_digit_nat(div_bold(concat_digit(vn, vw1), concat_digit(vp, vw2))), monus(vm, times(concat_digit(pos2nat(vp), vw2), most_significant_digit_nat(div_bold(vm, concat_digit(vp, vw2))))))));
        result.push_back(data_equation(variable_list({vp, vpair_, vw1, vw2}), divmod_aux_whr4(vw1, vp, vw2, vpair_), divmod_aux_whr5(vp, vw2, vpair_, if_(equals_zero(last(vpair_)), most_significant_digit_nat(vw1), concat_digit(last(vpair_), vw1)))));
        result.push_back(data_equation(variable_list({vm, vp, vpair_, vw2}), divmod_aux_whr5(vp, vw2, vpair_, vm), divmod_aux_whr6(vp, vw2, vpair_, vm, most_significant_digit_nat(div_bold(vm, concat_digit(vp, vw2))))));
        result.push_back(data_equation(variable_list({vm, vm1, vp, vpair_, vw2}), divmod_aux_whr6(vp, vw2, vpair_, vm, vm1), nnpair(plus(concat_digit(first(vpair_), sort_machine_word::zero_word()), vm1), monus(vm, times(concat_digit(pos2nat(vp), vw2), vm1)))));
        result.push_back(data_equation(variable_list({vw1, vw2}), div_bold(most_significant_digit_nat(vw1), sort_pos::most_significant_digit(vw2)), sort_machine_word::div_word(vw1, vw2)));
        result.push_back(data_equation(variable_list({vp, vw1, vw2}), div_bold(most_significant_digit_nat(vw1), concat_digit(vp, vw2)), sort_machine_word::zero_word()));
        result.push_back(data_equation(variable_list({vw1, vw2, vw3}), div_bold(concat_digit(most_significant_digit_nat(vw1), vw2), sort_pos::most_significant_digit(vw3)), sort_machine_word::div_doubleword(vw1, vw2, vw3)));
        result.push_back(data_equation(variable_list({vw, vw1, vw2, vw3}), div_bold(concat_digit(most_significant_digit_nat(vw), vw1), concat_digit(sort_pos::most_significant_digit(vw2), vw3)), sort_machine_word::div_double_doubleword(vw, vw1, vw2, vw3)));
        result.push_back(data_equation(variable_list({vw, vw1, vw2, vw3, vw4}), div_bold(concat_digit(concat_digit(most_significant_digit_nat(vw), vw1), vw2), concat_digit(sort_pos::most_significant_digit(vw3), vw4)), sort_machine_word::div_triple_doubleword(vw, vw1, vw2, vw3, vw4)));
        result.push_back(data_equation(variable_list({vn, vp, vw1, vw2, vw3}), div_bold(concat_digit(vn, vw1), concat_digit(concat_digit(vp, vw2), vw3)), if_(less_equal(concat_digit(pos2nat(concat_digit(vp, vw2)), sort_machine_word::zero_word()), vn), sort_machine_word::max_word(), div_bold_whr(vn, vw1, vp, vw2, vw3, div_bold(vn, concat_digit(vp, vw2))))));
        result.push_back(data_equation(variable_list({vn, vp, vw, vw1, vw2, vw3}), div_bold_whr(vn, vw1, vp, vw2, vw3, vw), if_(less(concat_digit(vn, vw1), times(most_significant_digit_nat(vw), pos2nat(concat_digit(concat_digit(vp, vw2), vw3)))), sort_machine_word::pred_word(vw), vw)));
        result.push_back(data_equation(variable_list({vm1, vm2, vn1, vn2}), equal_to(nnpair(vn1, vn2), nnpair(vm1, vm2)), sort_bool::and_(equal_to(vn1, vm1), equal_to(vn2, vm2))));
        result.push_back(data_equation(variable_list({vm1, vm2, vn1, vn2}), less(nnpair(vn1, vn2), nnpair(vm1, vm2)), sort_bool::or_(less(vn1, vm1), sort_bool::and_(equal_to(vn1, vm1), less(vn2, vm2)))));
        result.push_back(data_equation(variable_list({vm1, vm2, vn1, vn2}), less_equal(nnpair(vn1, vn2), nnpair(vm1, vm2)), sort_bool::or_(less(vn1, vm1), sort_bool::and_(equal_to(vn1, vm1), less_equal(vn2, vm2)))));
        result.push_back(data_equation(variable_list({vm, vn}), first(nnpair(vm, vn)), vm));
        result.push_back(data_equation(variable_list({vm, vn}), last(nnpair(vm, vn)), vn));
        result.push_back(data_equation(variable_list({vm, vn}), swap_zero(vm, vn), if_(equals_zero(vn), vm, if_(equals_zero(vm), vn, if_(equal_to(vn, vm), most_significant_digit_nat(sort_machine_word::zero_word()), vn)))));
        result.push_back(data_equation(variable_list({vm1, vm2, vn1, vn2}), swap_zero_add(vn1, vn2, vm1, vm2), if_(equals_zero(vn1), if_(equals_zero(vn2), plus(vm1, vm2), if_(equals_zero(vm1), vm2, swap_zero(vn2, plus(vm1, swap_zero(vn2, vm2))))), if_(equals_zero(vn2), if_(equals_zero(vm2), vn1, swap_zero(vn1, plus(swap_zero(vn1, vm1), vm2))), swap_zero(plus(vn1, vn2), plus(swap_zero(vn1, vm1), swap_zero(vn2, vm2)))))));
        result.push_back(data_equation(variable_list({vm1, vm2, vn1, vn2}), swap_zero_min(vn1, vn2, vm1, vm2), if_(equals_zero(vn1), if_(equals_zero(vn2), minimum(vm1, vm2), if_(equals_zero(vm1), most_significant_digit_nat(sort_machine_word::zero_word()), minimum(vm1, swap_zero(vn2, vm2)))), if_(equals_zero(vn2), if_(equals_zero(vm2), most_significant_digit_nat(sort_machine_word::zero_word()), minimum(swap_zero(vn1, vm1), vm2)), swap_zero(minimum(vn1, vn2), minimum(swap_zero(vn1, vm1), swap_zero(vn2, vm2)))))));
        result.push_back(data_equation(variable_list({vm1, vm2, vn1, vn2}), swap_zero_monus(vn1, vn2, vm1, vm2), if_(equals_zero(vn1), if_(equals_zero(vn2), monus(vm1, vm2), if_(equals_zero(vm1), most_significant_digit_nat(sort_machine_word::zero_word()), monus(vm1, swap_zero(vn2, vm2)))), if_(equals_zero(vn2), if_(equals_zero(vm2), vm1, swap_zero(vn1, monus(swap_zero(vn1, vm1), vm2))), swap_zero(monus(vn1, vn2), monus(swap_zero(vn1, vm1), swap_zero(vn2, vm2)))))));
        result.push_back(data_equation(variable_list({vw}), sqrt(most_significant_digit_nat(vw)), most_significant_digit_nat(sort_machine_word::sqrt_word(vw))));
        result.push_back(data_equation(variable_list({vw1, vw2}), sqrt(concat_digit(most_significant_digit_nat(vw1), vw2)), most_significant_digit_nat(sort_machine_word::sqrt_doubleword(vw1, vw2))));
        result.push_back(data_equation(variable_list({vw1, vw2, vw3}), sqrt(concat_digit(concat_digit(most_significant_digit_nat(vw1), vw2), vw3)), sqrt_whr1(vw1, vw2, vw3, sort_machine_word::sqrt_tripleword_overflow(vw1, vw2, vw3))));
        result.push_back(data_equation(variable_list({voverflow, vw1, vw2, vw3}), sqrt_whr1(vw1, vw2, vw3, voverflow), if_(sort_machine_word::equals_zero_word(voverflow), most_significant_digit_nat(sort_machine_word::sqrt_tripleword(vw1, vw2, vw3)), concat_digit(most_significant_digit_nat(voverflow), sort_machine_word::sqrt_tripleword(vw1, vw2, vw3)))));
        result.push_back(data_equation(variable_list({vw1, vw2, vw3, vw4}), sqrt(concat_digit(concat_digit(concat_digit(most_significant_digit_nat(vw1), vw2), vw3), vw4)), sqrt_whr2(vw1, vw2, vw3, vw4, sort_machine_word::sqrt_quadrupleword_overflow(vw1, vw2, vw3, vw4))));
        result.push_back(data_equation(variable_list({voverflow, vw1, vw2, vw3, vw4}), sqrt_whr2(vw1, vw2, vw3, vw4, voverflow), if_(sort_machine_word::equals_zero_word(voverflow), most_significant_digit_nat(sort_machine_word::sqrt_quadrupleword(vw1, vw2, vw3, vw4)), concat_digit(most_significant_digit_nat(voverflow), sort_machine_word::sqrt_quadrupleword(vw1, vw2, vw3, vw4)))));
        result.push_back(data_equation(variable_list({vn, vw1, vw2, vw3, vw4}), sqrt(concat_digit(concat_digit(concat_digit(concat_digit(vn, vw1), vw2), vw3), vw4)), first(sqrt_pair(concat_digit(concat_digit(concat_digit(concat_digit(vn, vw1), vw2), vw3), vw4)))));
        result.push_back(data_equation(variable_list({vw}), sqrt_pair(most_significant_digit_nat(vw)), nnpair(most_significant_digit_nat(sort_machine_word::sqrt_word(vw)), most_significant_digit_nat(sort_machine_word::minus_word(vw, sort_machine_word::times_word(sort_machine_word::sqrt_word(vw), sort_machine_word::sqrt_word(vw)))))));
        result.push_back(data_equation(variable_list({vw1, vw2}), sqrt_pair(concat_digit(most_significant_digit_nat(vw1), vw2)), nnpair(most_significant_digit_nat(sort_machine_word::sqrt_doubleword(vw1, vw2)), monus(concat_digit(most_significant_digit_nat(vw1), vw2), exp(most_significant_digit_nat(sort_machine_word::sqrt_doubleword(vw1, vw2)), most_significant_digit_nat(sort_machine_word::two_word()))))));
        result.push_back(data_equation(variable_list({vw1, vw2, vw3}), sqrt_pair(concat_digit(concat_digit(most_significant_digit_nat(vw1), vw2), vw3)), sqrt_pair_whr1(vw1, vw2, vw3, sqrt(concat_digit(concat_digit(most_significant_digit_nat(vw1), vw2), vw3)))));
        result.push_back(data_equation(variable_list({vsolution, vw1, vw2, vw3}), sqrt_pair_whr1(vw1, vw2, vw3, vsolution), nnpair(vsolution, monus(concat_digit(concat_digit(most_significant_digit_nat(vw1), vw2), vw3), times(vsolution, vsolution)))));
        result.push_back(data_equation(variable_list({vw1, vw2, vw3, vw4}), sqrt_pair(concat_digit(concat_digit(concat_digit(most_significant_digit_nat(vw1), vw2), vw3), vw4)), sqrt_pair_whr2(vw1, vw2, vw3, vw4, sqrt(concat_digit(concat_digit(concat_digit(most_significant_digit_nat(vw1), vw2), vw3), vw4)))));
        result.push_back(data_equation(variable_list({vsolution, vw1, vw2, vw3, vw4}), sqrt_pair_whr2(vw1, vw2, vw3, vw4, vsolution), nnpair(vsolution, monus(concat_digit(concat_digit(concat_digit(most_significant_digit_nat(vw1), vw2), vw3), vw4), times(vsolution, vsolution)))));
        result.push_back(data_equation(variable_list({vn, vw1, vw2, vw3, vw4}), sqrt_pair(concat_digit(concat_digit(concat_digit(concat_digit(vn, vw1), vw2), vw3), vw4)), sqrt_pair_whr3(vw3, vw4, sqrt_pair(concat_digit(concat_digit(vn, vw1), vw2)))));
        result.push_back(data_equation(variable_list({vpq, vw3, vw4}), sqrt_pair_whr3(vw3, vw4, vpq), sqrt_pair_whr4(vw3, vw4, vpq, if_(equals_zero(first(vpq)), most_significant_digit_nat(sort_machine_word::zero_word()), concat_digit(first(vpq), sort_machine_word::zero_word())), if_(less(most_significant_digit_nat(sort_machine_word::zero_word()), last(vpq)), concat_digit(concat_digit(last(vpq), vw3), vw4), if_(sort_machine_word::not_equals_zero_word(vw3), concat_digit(most_significant_digit_nat(vw3), vw4), most_significant_digit_nat(vw4))), plus(times(first(vpq), first(vpq)), last(vpq)))));
        result.push_back(data_equation(variable_list({vm2, vm3, vm5, vpq, vw3, vw4}), sqrt_pair_whr4(vw3, vw4, vpq, vm3, vm2, vm5), sqrt_pair_whr5(vpq, vm3, vm2, if_(less(most_significant_digit_nat(sort_machine_word::zero_word()), vm5), concat_digit(concat_digit(vm5, vw3), vw4), if_(sort_machine_word::not_equals_zero_word(vw3), concat_digit(most_significant_digit_nat(vw3), vw4), most_significant_digit_nat(vw4))), div(vm2, nat2pos(times(most_significant_digit_nat(sort_machine_word::two_word()), vm3))))));
        result.push_back(data_equation(variable_list({vm2, vm3, vm4, vpq, vy_guess}), sqrt_pair_whr5(vpq, vm3, vm2, vm4, vy_guess), sqrt_pair_whr6(vm2, if_(less(times(most_significant_digit_nat(sort_machine_word::four_word()), first(vpq)), concat_digit(most_significant_digit_nat(sort_machine_word::three_word()), sort_machine_word::zero_word())), monus(sqrt(vm4), vm3), if_(less(vm2, times(plus(times(most_significant_digit_nat(sort_machine_word::two_word()), vm3), vy_guess), vy_guess)), natpred(vy_guess), vy_guess)), if_(equals_zero(first(vpq)), most_significant_digit_nat(sort_machine_word::zero_word()), concat_digit(first(vpq), sort_machine_word::zero_word())))));
        result.push_back(data_equation(variable_list({vm1, vm2, vy}), sqrt_pair_whr6(vm2, vy, vm1), nnpair(plus(vm1, vy), monus(vm2, times(plus(times(vm1, most_significant_digit_nat(sort_machine_word::two_word())), vy), vy)))));
        return result;
      }

    } // namespace sort_nat

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_NAT64_H
