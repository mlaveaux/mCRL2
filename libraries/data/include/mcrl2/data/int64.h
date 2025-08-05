// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/int.h
/// \brief The standard sort int_.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/int.spec.

#ifndef MCRL2_DATA_INT64_H
#define MCRL2_DATA_INT64_H

#include "mcrl2/utilities/exception.h"
#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/pos64.h"
#include "mcrl2/data/nat64.h"

#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort int_.
    namespace sort_int {

      inline
      const core::identifier_string& int_name()
      {
        static core::identifier_string int_name = core::identifier_string("Int");
        return int_name;
      }

      /// \brief Constructor for sort expression Int.
      /// \return Sort expression Int.
      inline
      const basic_sort& int_()
      {
        static basic_sort int_ = basic_sort(int_name());
        return int_;
      }

      /// \brief Recogniser for sort expression Int
      /// \param e A sort expression
      /// \return true iff e == int_()
      inline
      bool is_int(const sort_expression& e)
      {
        if (is_basic_sort(e))
        {
          return basic_sort(e) == int_();
        }
        return false;
      }


      /// \brief Generate identifier \@cInt.
      /// \return Identifier \@cInt.
      inline
      const core::identifier_string& cint_name()
      {
        static core::identifier_string cint_name = core::identifier_string("@cInt");
        return cint_name;
      }

      /// \brief Constructor for function symbol \@cInt.
      
      /// \return Function symbol cint.
      inline
      const function_symbol& cint()
      {
        static function_symbol cint(cint_name(), make_function_sort_(sort_nat::nat(), int_()));
        return cint;
      }

      /// \brief Recogniser for function \@cInt.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@cInt.
      inline
      bool is_cint_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == cint();
        }
        return false;
      }

      /// \brief Application of function symbol \@cInt.
      
      /// \param arg0 A data expression.
      /// \return Application of \@cInt to a number of arguments.
      inline
      application cint(const data_expression& arg0)
      {
        return sort_int::cint()(arg0);
      }

      /// \brief Make an application of function symbol \@cInt.
      /// \param result The data expression where the \@cInt expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_cint(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::cint(),arg0);
      }

      /// \brief Recogniser for application of \@cInt.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol cint to a
      ///     number of arguments.
      inline
      bool is_cint_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_cint_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@cNeg.
      /// \return Identifier \@cNeg.
      inline
      const core::identifier_string& cneg_name()
      {
        static core::identifier_string cneg_name = core::identifier_string("@cNeg");
        return cneg_name;
      }

      /// \brief Constructor for function symbol \@cNeg.
      
      /// \return Function symbol cneg.
      inline
      const function_symbol& cneg()
      {
        static function_symbol cneg(cneg_name(), make_function_sort_(sort_pos::pos(), int_()));
        return cneg;
      }

      /// \brief Recogniser for function \@cNeg.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@cNeg.
      inline
      bool is_cneg_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == cneg();
        }
        return false;
      }

      /// \brief Application of function symbol \@cNeg.
      
      /// \param arg0 A data expression.
      /// \return Application of \@cNeg to a number of arguments.
      inline
      application cneg(const data_expression& arg0)
      {
        return sort_int::cneg()(arg0);
      }

      /// \brief Make an application of function symbol \@cNeg.
      /// \param result The data expression where the \@cNeg expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_cneg(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::cneg(),arg0);
      }

      /// \brief Recogniser for application of \@cNeg.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol cneg to a
      ///     number of arguments.
      inline
      bool is_cneg_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_cneg_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for int_.
      /// \return All system defined constructors for int_.
      inline
      function_symbol_vector int_generate_constructors_code()
      {
        function_symbol_vector result;
        result.push_back(sort_int::cint());
        result.push_back(sort_int::cneg());

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for int_.
      /// \return All system defined constructors that can be used in an mCRL2 specification for int_.
      inline
      function_symbol_vector int_mCRL2_usable_constructors()
      {
        function_symbol_vector result;
        result.push_back(sort_int::cint());
        result.push_back(sort_int::cneg());

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for int_.
      /// \return All system defined constructors that are to be implemented in C++ for int_.
      inline
      implementation_map int_cpp_implementable_constructors()
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier Nat2Int.
      /// \return Identifier Nat2Int.
      inline
      const core::identifier_string& nat2int_name()
      {
        static core::identifier_string nat2int_name = core::identifier_string("Nat2Int");
        return nat2int_name;
      }

      /// \brief Constructor for function symbol Nat2Int.
      
      /// \return Function symbol nat2int.
      inline
      const function_symbol& nat2int()
      {
        static function_symbol nat2int(nat2int_name(), make_function_sort_(sort_nat::nat(), int_()));
        return nat2int;
      }

      /// \brief Recogniser for function Nat2Int.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Nat2Int.
      inline
      bool is_nat2int_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == nat2int();
        }
        return false;
      }

      /// \brief Application of function symbol Nat2Int.
      
      /// \param arg0 A data expression.
      /// \return Application of Nat2Int to a number of arguments.
      inline
      application nat2int(const data_expression& arg0)
      {
        return sort_int::nat2int()(arg0);
      }

      /// \brief Make an application of function symbol Nat2Int.
      /// \param result The data expression where the Nat2Int expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_nat2int(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::nat2int(),arg0);
      }

      /// \brief Recogniser for application of Nat2Int.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol nat2int to a
      ///     number of arguments.
      inline
      bool is_nat2int_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_nat2int_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Int2Nat.
      /// \return Identifier Int2Nat.
      inline
      const core::identifier_string& int2nat_name()
      {
        static core::identifier_string int2nat_name = core::identifier_string("Int2Nat");
        return int2nat_name;
      }

      /// \brief Constructor for function symbol Int2Nat.
      
      /// \return Function symbol int2nat.
      inline
      const function_symbol& int2nat()
      {
        static function_symbol int2nat(int2nat_name(), make_function_sort_(int_(), sort_nat::nat()));
        return int2nat;
      }

      /// \brief Recogniser for function Int2Nat.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Int2Nat.
      inline
      bool is_int2nat_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == int2nat();
        }
        return false;
      }

      /// \brief Application of function symbol Int2Nat.
      
      /// \param arg0 A data expression.
      /// \return Application of Int2Nat to a number of arguments.
      inline
      application int2nat(const data_expression& arg0)
      {
        return sort_int::int2nat()(arg0);
      }

      /// \brief Make an application of function symbol Int2Nat.
      /// \param result The data expression where the Int2Nat expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_int2nat(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::int2nat(),arg0);
      }

      /// \brief Recogniser for application of Int2Nat.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol int2nat to a
      ///     number of arguments.
      inline
      bool is_int2nat_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_int2nat_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Pos2Int.
      /// \return Identifier Pos2Int.
      inline
      const core::identifier_string& pos2int_name()
      {
        static core::identifier_string pos2int_name = core::identifier_string("Pos2Int");
        return pos2int_name;
      }

      /// \brief Constructor for function symbol Pos2Int.
      
      /// \return Function symbol pos2int.
      inline
      const function_symbol& pos2int()
      {
        static function_symbol pos2int(pos2int_name(), make_function_sort_(sort_pos::pos(), int_()));
        return pos2int;
      }

      /// \brief Recogniser for function Pos2Int.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Pos2Int.
      inline
      bool is_pos2int_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == pos2int();
        }
        return false;
      }

      /// \brief Application of function symbol Pos2Int.
      
      /// \param arg0 A data expression.
      /// \return Application of Pos2Int to a number of arguments.
      inline
      application pos2int(const data_expression& arg0)
      {
        return sort_int::pos2int()(arg0);
      }

      /// \brief Make an application of function symbol Pos2Int.
      /// \param result The data expression where the Pos2Int expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pos2int(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::pos2int(),arg0);
      }

      /// \brief Recogniser for application of Pos2Int.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol pos2int to a
      ///     number of arguments.
      inline
      bool is_pos2int_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_pos2int_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Int2Pos.
      /// \return Identifier Int2Pos.
      inline
      const core::identifier_string& int2pos_name()
      {
        static core::identifier_string int2pos_name = core::identifier_string("Int2Pos");
        return int2pos_name;
      }

      /// \brief Constructor for function symbol Int2Pos.
      
      /// \return Function symbol int2pos.
      inline
      const function_symbol& int2pos()
      {
        static function_symbol int2pos(int2pos_name(), make_function_sort_(int_(), sort_pos::pos()));
        return int2pos;
      }

      /// \brief Recogniser for function Int2Pos.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Int2Pos.
      inline
      bool is_int2pos_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == int2pos();
        }
        return false;
      }

      /// \brief Application of function symbol Int2Pos.
      
      /// \param arg0 A data expression.
      /// \return Application of Int2Pos to a number of arguments.
      inline
      application int2pos(const data_expression& arg0)
      {
        return sort_int::int2pos()(arg0);
      }

      /// \brief Make an application of function symbol Int2Pos.
      /// \param result The data expression where the Int2Pos expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_int2pos(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::int2pos(),arg0);
      }

      /// \brief Recogniser for application of Int2Pos.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol int2pos to a
      ///     number of arguments.
      inline
      bool is_int2pos_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_int2pos_function_symbol(atermpp::down_cast<application>(e).head());
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
        if (s0 == sort_pos::pos() && s1 == int_())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == int_() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_nat::nat() && s1 == int_())
        {
          target_sort = sort_nat::nat();
        }
        else if (s0 == int_() && s1 == sort_nat::nat())
        {
          target_sort = sort_nat::nat();
        }
        else if (s0 == int_() && s1 == int_())
        {
          target_sort = int_();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_nat::nat())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_nat::nat())
        {
          target_sort = sort_nat::nat();
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
          return f.name() == maximum_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == maximum(sort_pos::pos(), int_()) || f == maximum(int_(), sort_pos::pos()) || f == maximum(sort_nat::nat(), int_()) || f == maximum(int_(), sort_nat::nat()) || f == maximum(int_(), int_()) || f == maximum(sort_pos::pos(), sort_nat::nat()) || f == maximum(sort_nat::nat(), sort_pos::pos()) || f == maximum(sort_nat::nat(), sort_nat::nat()) || f == maximum(sort_pos::pos(), sort_pos::pos()));
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
        return sort_int::maximum(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol max.
      /// \param result The data expression where the max expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_maximum(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_int::maximum(arg0.sort(), arg1.sort()),arg0, arg1);
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
        if (s0 == int_() && s1 == int_())
        {
          target_sort = int_();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_nat::nat())
        {
          target_sort = sort_nat::nat();
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
          return f.name() == minimum_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == minimum(int_(), int_()) || f == minimum(sort_nat::nat(), sort_nat::nat()) || f == minimum(sort_pos::pos(), sort_pos::pos()));
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
        return sort_int::minimum(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol min.
      /// \param result The data expression where the min expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_minimum(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_int::minimum(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier abs.
      /// \return Identifier abs.
      inline
      const core::identifier_string& abs_name()
      {
        static core::identifier_string abs_name = core::identifier_string("abs");
        return abs_name;
      }

      /// \brief Constructor for function symbol abs.
      
      /// \return Function symbol abs.
      inline
      const function_symbol& abs()
      {
        static function_symbol abs(abs_name(), make_function_sort_(int_(), sort_nat::nat()));
        return abs;
      }

      /// \brief Recogniser for function abs.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching abs.
      inline
      bool is_abs_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == abs();
        }
        return false;
      }

      /// \brief Application of function symbol abs.
      
      /// \param arg0 A data expression.
      /// \return Application of abs to a number of arguments.
      inline
      application abs(const data_expression& arg0)
      {
        return sort_int::abs()(arg0);
      }

      /// \brief Make an application of function symbol abs.
      /// \param result The data expression where the abs expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_abs(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::abs(),arg0);
      }

      /// \brief Recogniser for application of abs.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol abs to a
      ///     number of arguments.
      inline
      bool is_abs_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_abs_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier -.
      /// \return Identifier -.
      inline
      const core::identifier_string& negate_name()
      {
        static core::identifier_string negate_name = core::identifier_string("-");
        return negate_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol negate(const sort_expression& s0)
      {
        sort_expression target_sort(int_());
        function_symbol negate(negate_name(), make_function_sort_(s0, target_sort));
        return negate;
      }

      /// \brief Recogniser for function -.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching -.
      inline
      bool is_negate_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == negate_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == negate(sort_pos::pos()) || f == negate(sort_nat::nat()) || f == negate(int_()));
        }
        return false;
      }

      /// \brief Application of function symbol -.
      
      /// \param arg0 A data expression.
      /// \return Application of - to a number of arguments.
      inline
      application negate(const data_expression& arg0)
      {
        return sort_int::negate(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol -.
      /// \param result The data expression where the - expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_negate(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::negate(arg0.sort()),arg0);
      }

      /// \brief Recogniser for application of -.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol negate to a
      ///     number of arguments.
      inline
      bool is_negate_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_negate_function_symbol(atermpp::down_cast<application>(e).head());
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
        sort_expression target_sort;
        if (s0 == int_())
        {
          target_sort = int_();
        }
        else if (s0 == sort_nat::nat())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for succ with domain sorts " + pp(s0) + ". ");
        }

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
          return f.name() == succ_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == succ(int_()) || f == succ(sort_nat::nat()) || f == succ(sort_pos::pos()));
        }
        return false;
      }

      /// \brief Application of function symbol succ.
      
      /// \param arg0 A data expression.
      /// \return Application of succ to a number of arguments.
      inline
      application succ(const data_expression& arg0)
      {
        return sort_int::succ(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol succ.
      /// \param result The data expression where the succ expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_succ(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::succ(arg0.sort()),arg0);
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

      /// \brief Generate identifier pred.
      /// \return Identifier pred.
      inline
      const core::identifier_string& pred_name()
      {
        static core::identifier_string pred_name = core::identifier_string("pred");
        return pred_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol pred(const sort_expression& s0)
      {
        sort_expression target_sort;
        if (s0 == sort_nat::nat())
        {
          target_sort = int_();
        }
        else if (s0 == int_())
        {
          target_sort = int_();
        }
        else if (s0 == sort_pos::pos())
        {
          target_sort = sort_nat::nat();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for pred with domain sorts " + pp(s0) + ". ");
        }

        function_symbol pred(pred_name(), make_function_sort_(s0, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == pred_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == pred(sort_nat::nat()) || f == pred(int_()) || f == pred(sort_pos::pos()));
        }
        return false;
      }

      /// \brief Application of function symbol pred.
      
      /// \param arg0 A data expression.
      /// \return Application of pred to a number of arguments.
      inline
      application pred(const data_expression& arg0)
      {
        return sort_int::pred(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol pred.
      /// \param result The data expression where the pred expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pred(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_int::pred(arg0.sort()),arg0);
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
        if (s0 == int_() && s1 == int_())
        {
          target_sort = int_();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_nat::nat())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_nat::nat())
        {
          target_sort = sort_nat::nat();
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
          return f.name() == plus_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == plus(int_(), int_()) || f == plus(sort_pos::pos(), sort_nat::nat()) || f == plus(sort_nat::nat(), sort_pos::pos()) || f == plus(sort_nat::nat(), sort_nat::nat()) || f == plus(sort_pos::pos(), sort_pos::pos()));
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
        return sort_int::plus(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol +.
      /// \param result The data expression where the + expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_plus(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_int::plus(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier -.
      /// \return Identifier -.
      inline
      const core::identifier_string& minus_name()
      {
        static core::identifier_string minus_name = core::identifier_string("-");
        return minus_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol minus(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort(int_());
        function_symbol minus(minus_name(), make_function_sort_(s0, s1, target_sort));
        return minus;
      }

      /// \brief Recogniser for function -.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching -.
      inline
      bool is_minus_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == minus_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == minus(sort_pos::pos(), sort_pos::pos()) || f == minus(sort_nat::nat(), sort_nat::nat()) || f == minus(int_(), int_()));
        }
        return false;
      }

      /// \brief Application of function symbol -.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of - to a number of arguments.
      inline
      application minus(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_int::minus(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol -.
      /// \param result The data expression where the - expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_minus(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_int::minus(arg0.sort(), arg1.sort()),arg0, arg1);
      }

      /// \brief Recogniser for application of -.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol minus to a
      ///     number of arguments.
      inline
      bool is_minus_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_minus_function_symbol(atermpp::down_cast<application>(e).head());
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
        if (s0 == int_() && s1 == int_())
        {
          target_sort = int_();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_nat::nat())
        {
          target_sort = sort_nat::nat();
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
          return f.name() == times_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == times(int_(), int_()) || f == times(sort_nat::nat(), sort_nat::nat()) || f == times(sort_pos::pos(), sort_pos::pos()));
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
        return sort_int::times(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol *.
      /// \param result The data expression where the * expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_int::times(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier div.
      /// \return Identifier div.
      inline
      const core::identifier_string& div_name()
      {
        static core::identifier_string div_name = core::identifier_string("div");
        return div_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol div(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == int_() && s1 == sort_pos::pos())
        {
          target_sort = int_();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_pos::pos())
        {
          target_sort = sort_nat::nat();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for div with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol div(div_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == div_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == div(int_(), sort_pos::pos()) || f == div(sort_nat::nat(), sort_pos::pos()));
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
        return sort_int::div(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol div.
      /// \param result The data expression where the div expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_div(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_int::div(arg0.sort(), arg1.sort()),arg0, arg1);
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

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol mod(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort(sort_nat::nat());
        function_symbol mod(mod_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == mod_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == mod(int_(), sort_pos::pos()) || f == mod(sort_nat::nat(), sort_pos::pos()));
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
        return sort_int::mod(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol mod.
      /// \param result The data expression where the mod expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_mod(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_int::mod(arg0.sort(), arg1.sort()),arg0, arg1);
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
        if (s0 == int_() && s1 == sort_nat::nat())
        {
          target_sort = int_();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_nat::nat())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_nat::nat())
        {
          target_sort = sort_nat::nat();
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
          return f.name() == exp_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == exp(int_(), sort_nat::nat()) || f == exp(sort_pos::pos(), sort_nat::nat()) || f == exp(sort_nat::nat(), sort_nat::nat()));
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
        return sort_int::exp(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol exp.
      /// \param result The data expression where the exp expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_exp(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_int::exp(arg0.sort(), arg1.sort()),arg0, arg1);
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
      /// \brief Give all system defined mappings for int_
      /// \return All system defined mappings for int_
      inline
      function_symbol_vector int_generate_functions_code()
      {
        function_symbol_vector result;
        result.push_back(sort_int::nat2int());
        result.push_back(sort_int::int2nat());
        result.push_back(sort_int::pos2int());
        result.push_back(sort_int::int2pos());
        result.push_back(sort_int::maximum(sort_pos::pos(), int_()));
        result.push_back(sort_int::maximum(int_(), sort_pos::pos()));
        result.push_back(sort_int::maximum(sort_nat::nat(), int_()));
        result.push_back(sort_int::maximum(int_(), sort_nat::nat()));
        result.push_back(sort_int::maximum(int_(), int_()));
        result.push_back(sort_int::minimum(int_(), int_()));
        result.push_back(sort_int::abs());
        result.push_back(sort_int::negate(sort_pos::pos()));
        result.push_back(sort_int::negate(sort_nat::nat()));
        result.push_back(sort_int::negate(int_()));
        result.push_back(sort_int::succ(int_()));
        result.push_back(sort_int::pred(sort_nat::nat()));
        result.push_back(sort_int::pred(int_()));
        result.push_back(sort_int::plus(int_(), int_()));
        result.push_back(sort_int::minus(sort_pos::pos(), sort_pos::pos()));
        result.push_back(sort_int::minus(sort_nat::nat(), sort_nat::nat()));
        result.push_back(sort_int::minus(int_(), int_()));
        result.push_back(sort_int::times(int_(), int_()));
        result.push_back(sort_int::div(int_(), sort_pos::pos()));
        result.push_back(sort_int::mod(int_(), sort_pos::pos()));
        result.push_back(sort_int::exp(int_(), sort_nat::nat()));
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for int_
      /// \return All system defined mappings for int_
      inline
      function_symbol_vector int_generate_constructors_and_functions_code()
      {
        function_symbol_vector result=int_generate_functions_code();
        for(const function_symbol& f: int_generate_constructors_code())
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for int_
      /// \return All system defined mappings for that can be used in mCRL2 specificationis int_
      inline
      function_symbol_vector int_mCRL2_usable_mappings()
      {
        function_symbol_vector result;
        result.push_back(sort_int::nat2int());
        result.push_back(sort_int::int2nat());
        result.push_back(sort_int::pos2int());
        result.push_back(sort_int::int2pos());
        result.push_back(sort_int::maximum(sort_pos::pos(), int_()));
        result.push_back(sort_int::maximum(int_(), sort_pos::pos()));
        result.push_back(sort_int::maximum(sort_nat::nat(), int_()));
        result.push_back(sort_int::maximum(int_(), sort_nat::nat()));
        result.push_back(sort_int::maximum(int_(), int_()));
        result.push_back(sort_int::minimum(int_(), int_()));
        result.push_back(sort_int::abs());
        result.push_back(sort_int::negate(sort_pos::pos()));
        result.push_back(sort_int::negate(sort_nat::nat()));
        result.push_back(sort_int::negate(int_()));
        result.push_back(sort_int::succ(int_()));
        result.push_back(sort_int::pred(sort_nat::nat()));
        result.push_back(sort_int::pred(int_()));
        result.push_back(sort_int::plus(int_(), int_()));
        result.push_back(sort_int::minus(sort_pos::pos(), sort_pos::pos()));
        result.push_back(sort_int::minus(sort_nat::nat(), sort_nat::nat()));
        result.push_back(sort_int::minus(int_(), int_()));
        result.push_back(sort_int::times(int_(), int_()));
        result.push_back(sort_int::div(int_(), sort_pos::pos()));
        result.push_back(sort_int::mod(int_(), sort_pos::pos()));
        result.push_back(sort_int::exp(int_(), sort_nat::nat()));
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for int_
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for int_
      inline
      implementation_map int_cpp_implementable_mappings()
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
        assert(is_cint_application(e) || is_cneg_application(e) || is_nat2int_application(e) || is_int2nat_application(e) || is_pos2int_application(e) || is_int2pos_application(e) || is_abs_application(e) || is_negate_application(e) || is_succ_application(e) || is_pred_application(e));
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
        assert(is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_minus_application(e) || is_times_application(e) || is_div_application(e) || is_mod_application(e) || is_exp_application(e));
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
        assert(is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_minus_application(e) || is_times_application(e) || is_div_application(e) || is_mod_application(e) || is_exp_application(e));
        return atermpp::down_cast<application>(e)[1];
      }

      /// \brief Give all system defined equations for int_
      /// \return All system defined equations for sort int_
      inline
      data_equation_vector int_generate_equations_code()
      {
        variable vb("b",sort_bool::bool_());
        variable vn("n",sort_nat::nat());
        variable vm("m",sort_nat::nat());
        variable vp("p",sort_pos::pos());
        variable vq("q",sort_pos::pos());
        variable vx("x",int_());
        variable vy("y",int_());

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vm, vn}), equal_to(cint(vm), cint(vn)), equal_to(vm, vn)));
        result.push_back(data_equation(variable_list({vn, vp}), equal_to(cint(vn), cneg(vp)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vn, vp}), equal_to(cneg(vp), cint(vn)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp, vq}), equal_to(cneg(vp), cneg(vq)), equal_to(vp, vq)));
        result.push_back(data_equation(variable_list({vm, vn}), less(cint(vm), cint(vn)), less(vm, vn)));
        result.push_back(data_equation(variable_list({vn, vp}), less(cint(vn), cneg(vp)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vn, vp}), less(cneg(vp), cint(vn)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vp, vq}), less(cneg(vp), cneg(vq)), less(vq, vp)));
        result.push_back(data_equation(variable_list({vm, vn}), less_equal(cint(vm), cint(vn)), less_equal(vm, vn)));
        result.push_back(data_equation(variable_list({vn, vp}), less_equal(cint(vn), cneg(vp)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vn, vp}), less_equal(cneg(vp), cint(vn)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vp, vq}), less_equal(cneg(vp), cneg(vq)), less_equal(vq, vp)));
        result.push_back(data_equation(variable_list({vn}), nat2int(vn), cint(vn)));
        result.push_back(data_equation(variable_list({vn}), int2nat(cint(vn)), vn));
        result.push_back(data_equation(variable_list({vp}), pos2int(vp), cint(sort_nat::pos2nat(vp))));
        result.push_back(data_equation(variable_list({vn}), int2pos(cint(vn)), sort_nat::nat2pos(vn)));
        result.push_back(data_equation(variable_list({vn, vp}), maximum(vp, cint(vn)), maximum(vp, vn)));
        result.push_back(data_equation(variable_list({vp, vq}), maximum(vp, cneg(vq)), vp));
        result.push_back(data_equation(variable_list({vn, vp}), maximum(cint(vn), vp), maximum(vn, vp)));
        result.push_back(data_equation(variable_list({vp, vq}), maximum(cneg(vq), vp), vp));
        result.push_back(data_equation(variable_list({vm, vn}), maximum(vm, cint(vn)), if_(less_equal(vm, vn), vn, vm)));
        result.push_back(data_equation(variable_list({vn, vp}), maximum(vn, cneg(vp)), vn));
        result.push_back(data_equation(variable_list({vm, vn}), maximum(cint(vm), vn), if_(less_equal(vm, vn), vn, vm)));
        result.push_back(data_equation(variable_list({vn, vp}), maximum(cneg(vp), vn), vn));
        result.push_back(data_equation(variable_list({vx, vy}), maximum(vx, vy), if_(less_equal(vx, vy), vy, vx)));
        result.push_back(data_equation(variable_list({vx, vy}), minimum(vx, vy), if_(less_equal(vx, vy), vx, vy)));
        result.push_back(data_equation(variable_list({vn}), abs(cint(vn)), vn));
        result.push_back(data_equation(variable_list({vp}), abs(cneg(vp)), sort_nat::pos2nat(vp)));
        result.push_back(data_equation(variable_list({vp}), negate(vp), cneg(vp)));
        result.push_back(data_equation(variable_list(), negate(sort_nat::c0()), cint(sort_nat::c0())));
        result.push_back(data_equation(variable_list({vn}), negate(vn), if_(sort_nat::equals_zero(vn), cint(sort_nat::most_significant_digit_nat(sort_machine_word::zero_word())), cneg(sort_nat::nat2pos(vn)))));
        result.push_back(data_equation(variable_list({vn}), negate(cint(vn)), negate(vn)));
        result.push_back(data_equation(variable_list({vp}), negate(cneg(vp)), cint(sort_nat::pos2nat(vp))));
        result.push_back(data_equation(variable_list({vn}), succ(cint(vn)), cint(sort_nat::succ_nat(vn))));
        result.push_back(data_equation(variable_list({vp}), succ(cneg(vp)), negate(pred(vp))));
        result.push_back(data_equation(variable_list({vn}), pred(vn), if_(sort_nat::equals_zero(vn), cneg(sort_pos::most_significant_digit(sort_machine_word::one_word())), cint(sort_nat::natpred(vn)))));
        result.push_back(data_equation(variable_list({vn}), pred(cint(vn)), pred(vn)));
        result.push_back(data_equation(variable_list({vp}), pred(cneg(vp)), cneg(succ(vp))));
        result.push_back(data_equation(variable_list({vm, vn}), plus(cint(vm), cint(vn)), cint(plus(vm, vn))));
        result.push_back(data_equation(variable_list({vn, vp}), plus(cint(vn), cneg(vp)), minus(vn, sort_nat::pos2nat(vp))));
        result.push_back(data_equation(variable_list({vn, vp}), plus(cneg(vp), cint(vn)), minus(vn, sort_nat::pos2nat(vp))));
        result.push_back(data_equation(variable_list({vp, vq}), plus(cneg(vp), cneg(vq)), cneg(plus(vp, vq))));
        result.push_back(data_equation(variable_list({vp, vq}), less_equal(vq, vp), minus(vp, vq), cint(sort_nat::monus(sort_nat::pos2nat(vp), sort_nat::pos2nat(vq)))));
        result.push_back(data_equation(variable_list({vp, vq}), less(vp, vq), minus(vp, vq), negate(sort_nat::monus(sort_nat::pos2nat(vq), sort_nat::pos2nat(vp)))));
        result.push_back(data_equation(variable_list({vm, vn}), less_equal(vn, vm), minus(vm, vn), cint(sort_nat::monus(vm, vn))));
        result.push_back(data_equation(variable_list({vm, vn}), less(vm, vn), minus(vm, vn), negate(sort_nat::monus(vn, vm))));
        result.push_back(data_equation(variable_list({vx, vy}), minus(vx, vy), plus(vx, negate(vy))));
        result.push_back(data_equation(variable_list({vm, vn}), times(cint(vm), cint(vn)), cint(times(vm, vn))));
        result.push_back(data_equation(variable_list({vn, vp}), times(cint(vn), cneg(vp)), negate(times(sort_nat::pos2nat(vp), vn))));
        result.push_back(data_equation(variable_list({vn, vp}), times(cneg(vp), cint(vn)), negate(times(sort_nat::pos2nat(vp), vn))));
        result.push_back(data_equation(variable_list({vp, vq}), times(cneg(vp), cneg(vq)), cint(sort_nat::pos2nat(times(vp, vq)))));
        result.push_back(data_equation(variable_list({vn, vp}), div(cint(vn), vp), cint(div(vn, vp))));
        result.push_back(data_equation(variable_list({vp, vq}), div(cneg(vp), vq), cneg(succ(div(pred(vp), vq)))));
        result.push_back(data_equation(variable_list({vn, vp}), mod(cint(vn), vp), mod(vn, vp)));
        result.push_back(data_equation(variable_list({vp, vq}), mod(cneg(vp), vq), int2nat(minus(vq, succ(mod(pred(vp), vq))))));
        result.push_back(data_equation(variable_list({vm, vn}), exp(cint(vm), vn), cint(exp(vm, vn))));
        result.push_back(data_equation(variable_list({vn, vp}), exp(cneg(vp), vn), if_(sort_nat::is_odd(vn), cneg(exp(vp, vn)), cint(sort_nat::pos2nat(exp(vp, vn))))));
        return result;
      }

    } // namespace sort_int_

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_INT64_H
