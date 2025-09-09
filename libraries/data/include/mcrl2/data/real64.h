// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/real.h
/// \brief The standard sort real_.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/real.spec.

#ifndef MCRL2_DATA_REAL64_H
#define MCRL2_DATA_REAL64_H

#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/pos64.h"
#include "mcrl2/data/nat64.h"
#include "mcrl2/data/int64.h"

#include <functional>
#include <map>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort real_.
    namespace sort_real {

      inline
      const core::identifier_string& real_name()
      {
        static core::identifier_string real_name = core::identifier_string("Real");
        return real_name;
      }

      /// \brief Constructor for sort expression Real.
      /// \return Sort expression Real.
      inline
      const basic_sort& real_()
      {
        static basic_sort real_ = basic_sort(real_name());
        return real_;
      }

      /// \brief Recogniser for sort expression Real
      /// \param e A sort expression
      /// \return true iff e == real_()
      inline
      bool is_real(const sort_expression& e)
      {
        if (is_basic_sort(e))
        {
          return basic_sort(e) == real_();
        }
        return false;
      }

      /// \brief Give all system defined constructors for real_.
      /// \return All system defined constructors for real_.
      inline
      function_symbol_vector real_generate_constructors_code()
      {
        function_symbol_vector result;
        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for real_.
      /// \return All system defined constructors that can be used in an mCRL2 specification for real_.
      inline
      function_symbol_vector real_mCRL2_usable_constructors()
      {
        function_symbol_vector result;
        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for real_.
      /// \return All system defined constructors that are to be implemented in C++ for real_.
      inline
      implementation_map real_cpp_implementable_constructors()
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier \@cReal.
      /// \return Identifier \@cReal.
      inline
      const core::identifier_string& creal_name()
      {
        static core::identifier_string creal_name = core::identifier_string("@cReal");
        return creal_name;
      }

      /// \brief Constructor for function symbol \@cReal.
      
      /// \return Function symbol creal.
      inline
      const function_symbol& creal()
      {
        static function_symbol creal(creal_name(), make_function_sort_(sort_int::int_(), sort_pos::pos(), real_()));
        return creal;
      }

      /// \brief Recogniser for function \@cReal.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@cReal.
      inline
      bool is_creal_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == creal();
        }
        return false;
      }

      /// \brief Application of function symbol \@cReal.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@cReal to a number of arguments.
      inline
      application creal(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_real::creal()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@cReal.
      /// \param result The data expression where the \@cReal expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_creal(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::creal(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@cReal.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol creal to a
      ///     number of arguments.
      inline
      bool is_creal_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_creal_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Pos2Real.
      /// \return Identifier Pos2Real.
      inline
      const core::identifier_string& pos2real_name()
      {
        static core::identifier_string pos2real_name = core::identifier_string("Pos2Real");
        return pos2real_name;
      }

      /// \brief Constructor for function symbol Pos2Real.
      
      /// \return Function symbol pos2real.
      inline
      const function_symbol& pos2real()
      {
        static function_symbol pos2real(pos2real_name(), make_function_sort_(sort_pos::pos(), real_()));
        return pos2real;
      }

      /// \brief Recogniser for function Pos2Real.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Pos2Real.
      inline
      bool is_pos2real_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == pos2real();
        }
        return false;
      }

      /// \brief Application of function symbol Pos2Real.
      
      /// \param arg0 A data expression.
      /// \return Application of Pos2Real to a number of arguments.
      inline
      application pos2real(const data_expression& arg0)
      {
        return sort_real::pos2real()(arg0);
      }

      /// \brief Make an application of function symbol Pos2Real.
      /// \param result The data expression where the Pos2Real expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pos2real(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::pos2real(),arg0);
      }

      /// \brief Recogniser for application of Pos2Real.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol pos2real to a
      ///     number of arguments.
      inline
      bool is_pos2real_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_pos2real_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Nat2Real.
      /// \return Identifier Nat2Real.
      inline
      const core::identifier_string& nat2real_name()
      {
        static core::identifier_string nat2real_name = core::identifier_string("Nat2Real");
        return nat2real_name;
      }

      /// \brief Constructor for function symbol Nat2Real.
      
      /// \return Function symbol nat2real.
      inline
      const function_symbol& nat2real()
      {
        static function_symbol nat2real(nat2real_name(), make_function_sort_(sort_nat::nat(), real_()));
        return nat2real;
      }

      /// \brief Recogniser for function Nat2Real.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Nat2Real.
      inline
      bool is_nat2real_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == nat2real();
        }
        return false;
      }

      /// \brief Application of function symbol Nat2Real.
      
      /// \param arg0 A data expression.
      /// \return Application of Nat2Real to a number of arguments.
      inline
      application nat2real(const data_expression& arg0)
      {
        return sort_real::nat2real()(arg0);
      }

      /// \brief Make an application of function symbol Nat2Real.
      /// \param result The data expression where the Nat2Real expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_nat2real(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::nat2real(),arg0);
      }

      /// \brief Recogniser for application of Nat2Real.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol nat2real to a
      ///     number of arguments.
      inline
      bool is_nat2real_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_nat2real_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Int2Real.
      /// \return Identifier Int2Real.
      inline
      const core::identifier_string& int2real_name()
      {
        static core::identifier_string int2real_name = core::identifier_string("Int2Real");
        return int2real_name;
      }

      /// \brief Constructor for function symbol Int2Real.
      
      /// \return Function symbol int2real.
      inline
      const function_symbol& int2real()
      {
        static function_symbol int2real(int2real_name(), make_function_sort_(sort_int::int_(), real_()));
        return int2real;
      }

      /// \brief Recogniser for function Int2Real.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Int2Real.
      inline
      bool is_int2real_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == int2real();
        }
        return false;
      }

      /// \brief Application of function symbol Int2Real.
      
      /// \param arg0 A data expression.
      /// \return Application of Int2Real to a number of arguments.
      inline
      application int2real(const data_expression& arg0)
      {
        return sort_real::int2real()(arg0);
      }

      /// \brief Make an application of function symbol Int2Real.
      /// \param result The data expression where the Int2Real expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_int2real(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::int2real(),arg0);
      }

      /// \brief Recogniser for application of Int2Real.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol int2real to a
      ///     number of arguments.
      inline
      bool is_int2real_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_int2real_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Real2Pos.
      /// \return Identifier Real2Pos.
      inline
      const core::identifier_string& real2pos_name()
      {
        static core::identifier_string real2pos_name = core::identifier_string("Real2Pos");
        return real2pos_name;
      }

      /// \brief Constructor for function symbol Real2Pos.
      
      /// \return Function symbol real2pos.
      inline
      const function_symbol& real2pos()
      {
        static function_symbol real2pos(real2pos_name(), make_function_sort_(real_(), sort_pos::pos()));
        return real2pos;
      }

      /// \brief Recogniser for function Real2Pos.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Real2Pos.
      inline
      bool is_real2pos_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == real2pos();
        }
        return false;
      }

      /// \brief Application of function symbol Real2Pos.
      
      /// \param arg0 A data expression.
      /// \return Application of Real2Pos to a number of arguments.
      inline
      application real2pos(const data_expression& arg0)
      {
        return sort_real::real2pos()(arg0);
      }

      /// \brief Make an application of function symbol Real2Pos.
      /// \param result The data expression where the Real2Pos expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_real2pos(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::real2pos(),arg0);
      }

      /// \brief Recogniser for application of Real2Pos.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol real2pos to a
      ///     number of arguments.
      inline
      bool is_real2pos_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_real2pos_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Real2Nat.
      /// \return Identifier Real2Nat.
      inline
      const core::identifier_string& real2nat_name()
      {
        static core::identifier_string real2nat_name = core::identifier_string("Real2Nat");
        return real2nat_name;
      }

      /// \brief Constructor for function symbol Real2Nat.
      
      /// \return Function symbol real2nat.
      inline
      const function_symbol& real2nat()
      {
        static function_symbol real2nat(real2nat_name(), make_function_sort_(real_(), sort_nat::nat()));
        return real2nat;
      }

      /// \brief Recogniser for function Real2Nat.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Real2Nat.
      inline
      bool is_real2nat_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == real2nat();
        }
        return false;
      }

      /// \brief Application of function symbol Real2Nat.
      
      /// \param arg0 A data expression.
      /// \return Application of Real2Nat to a number of arguments.
      inline
      application real2nat(const data_expression& arg0)
      {
        return sort_real::real2nat()(arg0);
      }

      /// \brief Make an application of function symbol Real2Nat.
      /// \param result The data expression where the Real2Nat expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_real2nat(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::real2nat(),arg0);
      }

      /// \brief Recogniser for application of Real2Nat.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol real2nat to a
      ///     number of arguments.
      inline
      bool is_real2nat_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_real2nat_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Real2Int.
      /// \return Identifier Real2Int.
      inline
      const core::identifier_string& real2int_name()
      {
        static core::identifier_string real2int_name = core::identifier_string("Real2Int");
        return real2int_name;
      }

      /// \brief Constructor for function symbol Real2Int.
      
      /// \return Function symbol real2int.
      inline
      const function_symbol& real2int()
      {
        static function_symbol real2int(real2int_name(), make_function_sort_(real_(), sort_int::int_()));
        return real2int;
      }

      /// \brief Recogniser for function Real2Int.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Real2Int.
      inline
      bool is_real2int_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == real2int();
        }
        return false;
      }

      /// \brief Application of function symbol Real2Int.
      
      /// \param arg0 A data expression.
      /// \return Application of Real2Int to a number of arguments.
      inline
      application real2int(const data_expression& arg0)
      {
        return sort_real::real2int()(arg0);
      }

      /// \brief Make an application of function symbol Real2Int.
      /// \param result The data expression where the Real2Int expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_real2int(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::real2int(),arg0);
      }

      /// \brief Recogniser for application of Real2Int.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol real2int to a
      ///     number of arguments.
      inline
      bool is_real2int_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_real2int_function_symbol(atermpp::down_cast<application>(e).head());
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
        if (s0 == real_() && s1 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_int::int_())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_int::int_() && s1 == sort_pos::pos())
        {
          target_sort = sort_pos::pos();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_int::int_())
        {
          target_sort = sort_nat::nat();
        }
        else if (s0 == sort_int::int_() && s1 == sort_nat::nat())
        {
          target_sort = sort_nat::nat();
        }
        else if (s0 == sort_int::int_() && s1 == sort_int::int_())
        {
          target_sort = sort_int::int_();
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
          return f.name() == maximum_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == maximum(real_(), real_()) || f == maximum(sort_pos::pos(), sort_int::int_()) || f == maximum(sort_int::int_(), sort_pos::pos()) || f == maximum(sort_nat::nat(), sort_int::int_()) || f == maximum(sort_int::int_(), sort_nat::nat()) || f == maximum(sort_int::int_(), sort_int::int_()) || f == maximum(sort_pos::pos(), sort_nat::nat()) || f == maximum(sort_nat::nat(), sort_pos::pos()) || f == maximum(sort_nat::nat(), sort_nat::nat()) || f == maximum(sort_pos::pos(), sort_pos::pos()));
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
        return sort_real::maximum(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol max.
      /// \param result The data expression where the max expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_maximum(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::maximum(arg0.sort(), arg1.sort()),arg0, arg1);
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
        if (s0 == real_() && s1 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_int::int_() && s1 == sort_int::int_())
        {
          target_sort = sort_int::int_();
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
          return f.name() == minimum_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == minimum(real_(), real_()) || f == minimum(sort_int::int_(), sort_int::int_()) || f == minimum(sort_nat::nat(), sort_nat::nat()) || f == minimum(sort_pos::pos(), sort_pos::pos()));
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
        return sort_real::minimum(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol min.
      /// \param result The data expression where the min expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_minimum(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::minimum(arg0.sort(), arg1.sort()),arg0, arg1);
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

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol abs(const sort_expression& s0)
      {
        sort_expression target_sort;
        if (s0 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_int::int_())
        {
          target_sort = sort_nat::nat();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for abs with domain sorts " + pp(s0) + ". ");
        }

        function_symbol abs(abs_name(), make_function_sort_(s0, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == abs_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == abs(real_()) || f == abs(sort_int::int_()));
        }
        return false;
      }

      /// \brief Application of function symbol abs.
      
      /// \param arg0 A data expression.
      /// \return Application of abs to a number of arguments.
      inline
      application abs(const data_expression& arg0)
      {
        return sort_real::abs(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol abs.
      /// \param result The data expression where the abs expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_abs(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::abs(arg0.sort()),arg0);
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
        sort_expression target_sort;
        if (s0 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_pos::pos())
        {
          target_sort = sort_int::int_();
        }
        else if (s0 == sort_nat::nat())
        {
          target_sort = sort_int::int_();
        }
        else if (s0 == sort_int::int_())
        {
          target_sort = sort_int::int_();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for negate with domain sorts " + pp(s0) + ". ");
        }

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
          return f.name() == negate_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == negate(real_()) || f == negate(sort_pos::pos()) || f == negate(sort_nat::nat()) || f == negate(sort_int::int_()));
        }
        return false;
      }

      /// \brief Application of function symbol -.
      
      /// \param arg0 A data expression.
      /// \return Application of - to a number of arguments.
      inline
      application negate(const data_expression& arg0)
      {
        return sort_real::negate(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol -.
      /// \param result The data expression where the - expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_negate(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::negate(arg0.sort()),arg0);
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
        if (s0 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_int::int_())
        {
          target_sort = sort_int::int_();
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
          return f.name() == succ_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == succ(real_()) || f == succ(sort_int::int_()) || f == succ(sort_nat::nat()) || f == succ(sort_pos::pos()));
        }
        return false;
      }

      /// \brief Application of function symbol succ.
      
      /// \param arg0 A data expression.
      /// \return Application of succ to a number of arguments.
      inline
      application succ(const data_expression& arg0)
      {
        return sort_real::succ(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol succ.
      /// \param result The data expression where the succ expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_succ(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::succ(arg0.sort()),arg0);
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
        if (s0 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_nat::nat())
        {
          target_sort = sort_int::int_();
        }
        else if (s0 == sort_int::int_())
        {
          target_sort = sort_int::int_();
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
          return f.name() == pred_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 1 && (f == pred(real_()) || f == pred(sort_nat::nat()) || f == pred(sort_int::int_()) || f == pred(sort_pos::pos()));
        }
        return false;
      }

      /// \brief Application of function symbol pred.
      
      /// \param arg0 A data expression.
      /// \return Application of pred to a number of arguments.
      inline
      application pred(const data_expression& arg0)
      {
        return sort_real::pred(arg0.sort())(arg0);
      }

      /// \brief Make an application of function symbol pred.
      /// \param result The data expression where the pred expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_pred(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::pred(arg0.sort()),arg0);
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
        if (s0 == real_() && s1 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_int::int_() && s1 == sort_int::int_())
        {
          target_sort = sort_int::int_();
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
          return f.name() == plus_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == plus(real_(), real_()) || f == plus(sort_int::int_(), sort_int::int_()) || f == plus(sort_pos::pos(), sort_nat::nat()) || f == plus(sort_nat::nat(), sort_pos::pos()) || f == plus(sort_nat::nat(), sort_nat::nat()) || f == plus(sort_pos::pos(), sort_pos::pos()));
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
        return sort_real::plus(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol +.
      /// \param result The data expression where the + expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_plus(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::plus(arg0.sort(), arg1.sort()),arg0, arg1);
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
        sort_expression target_sort;
        if (s0 == real_() && s1 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_pos::pos() && s1 == sort_pos::pos())
        {
          target_sort = sort_int::int_();
        }
        else if (s0 == sort_nat::nat() && s1 == sort_nat::nat())
        {
          target_sort = sort_int::int_();
        }
        else if (s0 == sort_int::int_() && s1 == sort_int::int_())
        {
          target_sort = sort_int::int_();
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for minus with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

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
          return f.name() == minus_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == minus(real_(), real_()) || f == minus(sort_pos::pos(), sort_pos::pos()) || f == minus(sort_nat::nat(), sort_nat::nat()) || f == minus(sort_int::int_(), sort_int::int_()));
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
        return sort_real::minus(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol -.
      /// \param result The data expression where the - expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_minus(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::minus(arg0.sort(), arg1.sort()),arg0, arg1);
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
        if (s0 == real_() && s1 == real_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_int::int_() && s1 == sort_int::int_())
        {
          target_sort = sort_int::int_();
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
          return f.name() == times_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == times(real_(), real_()) || f == times(sort_int::int_(), sort_int::int_()) || f == times(sort_nat::nat(), sort_nat::nat()) || f == times(sort_pos::pos(), sort_pos::pos()));
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
        return sort_real::times(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol *.
      /// \param result The data expression where the * expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_times(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::times(arg0.sort(), arg1.sort()),arg0, arg1);
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
        if (s0 == real_() && s1 == sort_int::int_())
        {
          target_sort = real_();
        }
        else if (s0 == sort_int::int_() && s1 == sort_nat::nat())
        {
          target_sort = sort_int::int_();
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
          return f.name() == exp_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == exp(real_(), sort_int::int_()) || f == exp(sort_int::int_(), sort_nat::nat()) || f == exp(sort_pos::pos(), sort_nat::nat()) || f == exp(sort_nat::nat(), sort_nat::nat()));
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
        return sort_real::exp(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol exp.
      /// \param result The data expression where the exp expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_exp(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::exp(arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier /.
      /// \return Identifier /.
      inline
      const core::identifier_string& divides_name()
      {
        static core::identifier_string divides_name = core::identifier_string("/");
        return divides_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol divides(const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort(real_());
        function_symbol divides(divides_name(), make_function_sort_(s0, s1, target_sort));
        return divides;
      }

      /// \brief Recogniser for function /.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching /.
      inline
      bool is_divides_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == divides_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2 && (f == divides(sort_pos::pos(), sort_pos::pos()) || f == divides(sort_nat::nat(), sort_nat::nat()) || f == divides(sort_int::int_(), sort_int::int_()) || f == divides(real_(), real_()));
        }
        return false;
      }

      /// \brief Application of function symbol /.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of / to a number of arguments.
      inline
      application divides(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_real::divides(arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol /.
      /// \param result The data expression where the / expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_divides(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::divides(arg0.sort(), arg1.sort()),arg0, arg1);
      }

      /// \brief Recogniser for application of /.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol divides to a
      ///     number of arguments.
      inline
      bool is_divides_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_divides_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier floor.
      /// \return Identifier floor.
      inline
      const core::identifier_string& floor_name()
      {
        static core::identifier_string floor_name = core::identifier_string("floor");
        return floor_name;
      }

      /// \brief Constructor for function symbol floor.
      
      /// \return Function symbol floor.
      inline
      const function_symbol& floor()
      {
        static function_symbol floor(floor_name(), make_function_sort_(real_(), sort_int::int_()));
        return floor;
      }

      /// \brief Recogniser for function floor.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching floor.
      inline
      bool is_floor_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == floor();
        }
        return false;
      }

      /// \brief Application of function symbol floor.
      
      /// \param arg0 A data expression.
      /// \return Application of floor to a number of arguments.
      inline
      application floor(const data_expression& arg0)
      {
        return sort_real::floor()(arg0);
      }

      /// \brief Make an application of function symbol floor.
      /// \param result The data expression where the floor expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_floor(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::floor(),arg0);
      }

      /// \brief Recogniser for application of floor.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol floor to a
      ///     number of arguments.
      inline
      bool is_floor_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_floor_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier ceil.
      /// \return Identifier ceil.
      inline
      const core::identifier_string& ceil_name()
      {
        static core::identifier_string ceil_name = core::identifier_string("ceil");
        return ceil_name;
      }

      /// \brief Constructor for function symbol ceil.
      
      /// \return Function symbol ceil.
      inline
      const function_symbol& ceil()
      {
        static function_symbol ceil(ceil_name(), make_function_sort_(real_(), sort_int::int_()));
        return ceil;
      }

      /// \brief Recogniser for function ceil.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching ceil.
      inline
      bool is_ceil_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == ceil();
        }
        return false;
      }

      /// \brief Application of function symbol ceil.
      
      /// \param arg0 A data expression.
      /// \return Application of ceil to a number of arguments.
      inline
      application ceil(const data_expression& arg0)
      {
        return sort_real::ceil()(arg0);
      }

      /// \brief Make an application of function symbol ceil.
      /// \param result The data expression where the ceil expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_ceil(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::ceil(),arg0);
      }

      /// \brief Recogniser for application of ceil.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol ceil to a
      ///     number of arguments.
      inline
      bool is_ceil_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_ceil_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier round.
      /// \return Identifier round.
      inline
      const core::identifier_string& round_name()
      {
        static core::identifier_string round_name = core::identifier_string("round");
        return round_name;
      }

      /// \brief Constructor for function symbol round.
      
      /// \return Function symbol round.
      inline
      const function_symbol& round()
      {
        static function_symbol round(round_name(), make_function_sort_(real_(), sort_int::int_()));
        return round;
      }

      /// \brief Recogniser for function round.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching round.
      inline
      bool is_round_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == round();
        }
        return false;
      }

      /// \brief Application of function symbol round.
      
      /// \param arg0 A data expression.
      /// \return Application of round to a number of arguments.
      inline
      application round(const data_expression& arg0)
      {
        return sort_real::round()(arg0);
      }

      /// \brief Make an application of function symbol round.
      /// \param result The data expression where the round expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_round(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_real::round(),arg0);
      }

      /// \brief Recogniser for application of round.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol round to a
      ///     number of arguments.
      inline
      bool is_round_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_round_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@redfrac.
      /// \return Identifier \@redfrac.
      inline
      const core::identifier_string& reduce_fraction_name()
      {
        static core::identifier_string reduce_fraction_name = core::identifier_string("@redfrac");
        return reduce_fraction_name;
      }

      /// \brief Constructor for function symbol \@redfrac.
      
      /// \return Function symbol reduce_fraction.
      inline
      const function_symbol& reduce_fraction()
      {
        static function_symbol reduce_fraction(reduce_fraction_name(), make_function_sort_(sort_int::int_(), sort_int::int_(), real_()));
        return reduce_fraction;
      }

      /// \brief Recogniser for function \@redfrac.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@redfrac.
      inline
      bool is_reduce_fraction_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == reduce_fraction();
        }
        return false;
      }

      /// \brief Application of function symbol \@redfrac.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@redfrac to a number of arguments.
      inline
      application reduce_fraction(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_real::reduce_fraction()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@redfrac.
      /// \param result The data expression where the \@redfrac expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_reduce_fraction(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::reduce_fraction(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@redfrac.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol reduce_fraction to a
      ///     number of arguments.
      inline
      bool is_reduce_fraction_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_reduce_fraction_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@redfracwhr.
      /// \return Identifier \@redfracwhr.
      inline
      const core::identifier_string& reduce_fraction_where_name()
      {
        static core::identifier_string reduce_fraction_where_name = core::identifier_string("@redfracwhr");
        return reduce_fraction_where_name;
      }

      /// \brief Constructor for function symbol \@redfracwhr.
      
      /// \return Function symbol reduce_fraction_where.
      inline
      const function_symbol& reduce_fraction_where()
      {
        static function_symbol reduce_fraction_where(reduce_fraction_where_name(), make_function_sort_(sort_nat::nat(), sort_int::int_(), sort_nat::nat(), real_()));
        return reduce_fraction_where;
      }

      /// \brief Recogniser for function \@redfracwhr.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@redfracwhr.
      inline
      bool is_reduce_fraction_where_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == reduce_fraction_where();
        }
        return false;
      }

      /// \brief Application of function symbol \@redfracwhr.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@redfracwhr to a number of arguments.
      inline
      application reduce_fraction_where(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_real::reduce_fraction_where()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@redfracwhr.
      /// \param result The data expression where the \@redfracwhr expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_reduce_fraction_where(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_real::reduce_fraction_where(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@redfracwhr.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol reduce_fraction_where to a
      ///     number of arguments.
      inline
      bool is_reduce_fraction_where_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_reduce_fraction_where_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@redfrachlp.
      /// \return Identifier \@redfrachlp.
      inline
      const core::identifier_string& reduce_fraction_helper_name()
      {
        static core::identifier_string reduce_fraction_helper_name = core::identifier_string("@redfrachlp");
        return reduce_fraction_helper_name;
      }

      /// \brief Constructor for function symbol \@redfrachlp.
      
      /// \return Function symbol reduce_fraction_helper.
      inline
      const function_symbol& reduce_fraction_helper()
      {
        static function_symbol reduce_fraction_helper(reduce_fraction_helper_name(), make_function_sort_(real_(), sort_int::int_(), real_()));
        return reduce_fraction_helper;
      }

      /// \brief Recogniser for function \@redfrachlp.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@redfrachlp.
      inline
      bool is_reduce_fraction_helper_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == reduce_fraction_helper();
        }
        return false;
      }

      /// \brief Application of function symbol \@redfrachlp.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@redfrachlp to a number of arguments.
      inline
      application reduce_fraction_helper(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_real::reduce_fraction_helper()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@redfrachlp.
      /// \param result The data expression where the \@redfrachlp expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_reduce_fraction_helper(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_real::reduce_fraction_helper(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@redfrachlp.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol reduce_fraction_helper to a
      ///     number of arguments.
      inline
      bool is_reduce_fraction_helper_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_reduce_fraction_helper_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for real_
      /// \return All system defined mappings for real_
      inline
      function_symbol_vector real_generate_functions_code()
      {
        function_symbol_vector result;
        result.push_back(sort_real::creal());
        result.push_back(sort_real::pos2real());
        result.push_back(sort_real::nat2real());
        result.push_back(sort_real::int2real());
        result.push_back(sort_real::real2pos());
        result.push_back(sort_real::real2nat());
        result.push_back(sort_real::real2int());
        result.push_back(sort_real::maximum(real_(), real_()));
        result.push_back(sort_real::minimum(real_(), real_()));
        result.push_back(sort_real::abs(real_()));
        result.push_back(sort_real::negate(real_()));
        result.push_back(sort_real::succ(real_()));
        result.push_back(sort_real::pred(real_()));
        result.push_back(sort_real::plus(real_(), real_()));
        result.push_back(sort_real::minus(real_(), real_()));
        result.push_back(sort_real::times(real_(), real_()));
        result.push_back(sort_real::exp(real_(), sort_int::int_()));
        result.push_back(sort_real::divides(sort_pos::pos(), sort_pos::pos()));
        result.push_back(sort_real::divides(sort_nat::nat(), sort_nat::nat()));
        result.push_back(sort_real::divides(sort_int::int_(), sort_int::int_()));
        result.push_back(sort_real::divides(real_(), real_()));
        result.push_back(sort_real::floor());
        result.push_back(sort_real::ceil());
        result.push_back(sort_real::round());
        result.push_back(sort_real::reduce_fraction());
        result.push_back(sort_real::reduce_fraction_where());
        result.push_back(sort_real::reduce_fraction_helper());
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for real_
      /// \return All system defined mappings for real_
      inline
      function_symbol_vector real_generate_constructors_and_functions_code()
      {
        function_symbol_vector result=real_generate_functions_code();
        for(const function_symbol& f: real_generate_constructors_code())
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for real_
      /// \return All system defined mappings for that can be used in mCRL2 specificationis real_
      inline
      function_symbol_vector real_mCRL2_usable_mappings()
      {
        function_symbol_vector result;
        result.push_back(sort_real::creal());
        result.push_back(sort_real::pos2real());
        result.push_back(sort_real::nat2real());
        result.push_back(sort_real::int2real());
        result.push_back(sort_real::real2pos());
        result.push_back(sort_real::real2nat());
        result.push_back(sort_real::real2int());
        result.push_back(sort_real::maximum(real_(), real_()));
        result.push_back(sort_real::minimum(real_(), real_()));
        result.push_back(sort_real::abs(real_()));
        result.push_back(sort_real::negate(real_()));
        result.push_back(sort_real::succ(real_()));
        result.push_back(sort_real::pred(real_()));
        result.push_back(sort_real::plus(real_(), real_()));
        result.push_back(sort_real::minus(real_(), real_()));
        result.push_back(sort_real::times(real_(), real_()));
        result.push_back(sort_real::exp(real_(), sort_int::int_()));
        result.push_back(sort_real::divides(sort_pos::pos(), sort_pos::pos()));
        result.push_back(sort_real::divides(sort_nat::nat(), sort_nat::nat()));
        result.push_back(sort_real::divides(sort_int::int_(), sort_int::int_()));
        result.push_back(sort_real::divides(real_(), real_()));
        result.push_back(sort_real::floor());
        result.push_back(sort_real::ceil());
        result.push_back(sort_real::round());
        result.push_back(sort_real::reduce_fraction());
        result.push_back(sort_real::reduce_fraction_where());
        result.push_back(sort_real::reduce_fraction_helper());
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for real_
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for real_
      inline
      implementation_map real_cpp_implementable_mappings()
      {
        implementation_map result;
        return result;
      }
      ///\brief Function for projecting out argument.
      ///        left from an application.
      /// \param e A data expression.
      /// \pre left is defined for e.
      /// \return The argument of e that corresponds to left.
      inline
      const data_expression& left(const data_expression& e)
      {
        assert(is_creal_application(e) || is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_minus_application(e) || is_times_application(e) || is_exp_application(e) || is_divides_application(e) || is_reduce_fraction_application(e) || is_reduce_fraction_helper_application(e));
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
        assert(is_creal_application(e) || is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_minus_application(e) || is_times_application(e) || is_exp_application(e) || is_divides_application(e) || is_reduce_fraction_application(e) || is_reduce_fraction_helper_application(e));
        return atermpp::down_cast<application>(e)[1];
      }

      ///\brief Function for projecting out argument.
      ///        arg from an application.
      /// \param e A data expression.
      /// \pre arg is defined for e.
      /// \return The argument of e that corresponds to arg.
      inline
      const data_expression& arg(const data_expression& e)
      {
        assert(is_pos2real_application(e) || is_nat2real_application(e) || is_int2real_application(e) || is_real2pos_application(e) || is_real2nat_application(e) || is_real2int_application(e) || is_abs_application(e) || is_negate_application(e) || is_succ_application(e) || is_pred_application(e) || is_floor_application(e) || is_ceil_application(e) || is_round_application(e));
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
        assert(is_reduce_fraction_where_application(e));
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
        assert(is_reduce_fraction_where_application(e));
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
        assert(is_reduce_fraction_where_application(e));
        return atermpp::down_cast<application>(e)[2];
      }

      /// \brief Give all system defined equations for real_
      /// \return All system defined equations for sort real_
      inline
      data_equation_vector real_generate_equations_code()
      {
        variable vm("m",sort_nat::nat());
        variable vn("n",sort_nat::nat());
        variable vp("p",sort_pos::pos());
        variable vq("q",sort_pos::pos());
        variable vx("x",sort_int::int_());
        variable vy("y",sort_int::int_());
        variable vr("r",real_());
        variable vs("s",real_());

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vp, vq, vx, vy}), equal_to(creal(vx, vp), creal(vy, vq)), equal_to(times(vx, sort_int::cint(sort_nat::pos2nat(vq))), times(vy, sort_int::cint(sort_nat::pos2nat(vp))))));
        result.push_back(data_equation(variable_list({vp, vq, vx, vy}), less(creal(vx, vp), creal(vy, vq)), less(times(vx, sort_int::cint(sort_nat::pos2nat(vq))), times(vy, sort_int::cint(sort_nat::pos2nat(vp))))));
        result.push_back(data_equation(variable_list({vp, vq, vx, vy}), less_equal(creal(vx, vp), creal(vy, vq)), less_equal(times(vx, sort_int::cint(sort_nat::pos2nat(vq))), times(vy, sort_int::cint(sort_nat::pos2nat(vp))))));
        result.push_back(data_equation(variable_list({vx}), int2real(vx), creal(vx, sort_pos::c1())));
        result.push_back(data_equation(variable_list({vn}), nat2real(vn), creal(sort_int::cint(vn), sort_pos::c1())));
        result.push_back(data_equation(variable_list({vp}), pos2real(vp), creal(sort_int::cint(sort_nat::pos2nat(vp)), sort_pos::c1())));
        result.push_back(data_equation(variable_list({vp, vx}), equal_to(vp, sort_pos::c1()), real2int(creal(vx, vp)), vx));
        result.push_back(data_equation(variable_list({vp, vx}), equal_to(vp, sort_pos::c1()), real2nat(creal(vx, vp)), sort_int::int2nat(vx)));
        result.push_back(data_equation(variable_list({vp, vx}), equal_to(vp, sort_pos::c1()), real2pos(creal(vx, vp)), sort_int::int2pos(vx)));
        result.push_back(data_equation(variable_list({vr, vs}), minimum(vr, vs), if_(less(vr, vs), vr, vs)));
        result.push_back(data_equation(variable_list({vr, vs}), maximum(vr, vs), if_(less(vr, vs), vs, vr)));
        result.push_back(data_equation(variable_list({vr}), abs(vr), if_(less(vr, creal(sort_int::cint(sort_nat::c0()), sort_pos::c1())), negate(vr), vr)));
        result.push_back(data_equation(variable_list({vp, vx}), negate(creal(vx, vp)), creal(negate(vx), vp)));
        result.push_back(data_equation(variable_list({vp, vx}), succ(creal(vx, vp)), creal(plus(vx, sort_int::cint(sort_nat::pos2nat(vp))), vp)));
        result.push_back(data_equation(variable_list({vp, vx}), pred(creal(vx, vp)), creal(minus(vx, sort_int::cint(sort_nat::pos2nat(vp))), vp)));
        result.push_back(data_equation(variable_list({vp, vq, vx, vy}), plus(creal(vx, vp), creal(vy, vq)), reduce_fraction(plus(times(vx, sort_int::cint(sort_nat::pos2nat(vq))), times(vy, sort_int::cint(sort_nat::pos2nat(vp)))), sort_int::cint(sort_nat::pos2nat(times(vp, vq))))));
        result.push_back(data_equation(variable_list({vp, vq, vx, vy}), minus(creal(vx, vp), creal(vy, vq)), reduce_fraction(minus(times(vx, sort_int::cint(sort_nat::pos2nat(vq))), times(vy, sort_int::cint(sort_nat::pos2nat(vp)))), sort_int::cint(sort_nat::pos2nat(times(vp, vq))))));
        result.push_back(data_equation(variable_list({vp, vq, vx, vy}), times(creal(vx, vp), creal(vy, vq)), reduce_fraction(times(vx, vy), sort_int::cint(sort_nat::pos2nat(times(vp, vq))))));
        result.push_back(data_equation(variable_list({vm, vp, vr}), sort_nat::equals_zero(vm), times(vr, creal(sort_int::cint(vm), vp)), creal(sort_int::cint(sort_nat::c0()), sort_pos::c1())));
        result.push_back(data_equation(variable_list({vm, vp, vr}), sort_nat::equals_zero(vm), times(creal(sort_int::cint(vm), vp), vr), creal(sort_int::cint(sort_nat::c0()), sort_pos::c1())));
        result.push_back(data_equation(variable_list({vp, vq, vx, vy}), not_equal_to(vy, sort_int::cint(sort_nat::c0())), divides(creal(vx, vp), creal(vy, vq)), reduce_fraction(times(vx, sort_int::cint(sort_nat::pos2nat(vq))), times(vy, sort_int::cint(sort_nat::pos2nat(vp))))));
        result.push_back(data_equation(variable_list({vp, vq}), divides(vp, vq), reduce_fraction(sort_int::cint(sort_nat::pos2nat(vp)), sort_int::cint(sort_nat::pos2nat(vq)))));
        result.push_back(data_equation(variable_list({vm, vn}), not_equal_to(vn, sort_nat::c0()), divides(vm, vn), reduce_fraction(sort_int::cint(vm), sort_int::cint(vn))));
        result.push_back(data_equation(variable_list({vx, vy}), not_equal_to(vy, sort_int::cint(sort_nat::c0())), divides(vx, vy), reduce_fraction(vx, vy)));
        result.push_back(data_equation(variable_list({vn, vp, vx}), exp(creal(vx, vp), sort_int::cint(vn)), reduce_fraction(exp(vx, vn), sort_int::cint(sort_nat::pos2nat(exp(vp, vn))))));
        result.push_back(data_equation(variable_list({vp, vq, vx}), not_equal_to(vx, sort_int::cint(sort_nat::c0())), exp(creal(vx, vp), sort_int::cneg(vq)), reduce_fraction(sort_int::cint(sort_nat::pos2nat(exp(vp, sort_nat::pos2nat(vq)))), exp(vx, sort_nat::pos2nat(vq)))));
        result.push_back(data_equation(variable_list({vp, vx}), floor(creal(vx, vp)), sort_int::div(vx, vp)));
        result.push_back(data_equation(variable_list({vr}), ceil(vr), negate(floor(negate(vr)))));
        result.push_back(data_equation(variable_list({vr}), round(vr), floor(plus(vr, creal(sort_int::cint(sort_nat::pos2nat(sort_pos::c1())), plus(sort_pos::c1(), sort_pos::c1()))))));
        result.push_back(data_equation(variable_list({vp, vx}), reduce_fraction(vx, sort_int::cneg(vp)), reduce_fraction(negate(vx), sort_int::cint(sort_nat::pos2nat(vp)))));
        result.push_back(data_equation(variable_list({vn, vx}), reduce_fraction(vx, sort_int::cint(vn)), reduce_fraction_where(vn, sort_int::div(vx, sort_nat::nat2pos(vn)), sort_int::mod(vx, sort_nat::nat2pos(vn)))));
        result.push_back(data_equation(variable_list({vm, vn, vx}), sort_nat::equals_zero(vm), reduce_fraction_where(vn, vx, vm), creal(vx, sort_pos::c1())));
        result.push_back(data_equation(variable_list({vm, vn, vx}), sort_nat::not_equals_zero(vm), reduce_fraction_where(vn, vx, vm), reduce_fraction_helper(reduce_fraction(sort_int::cint(vn), sort_int::cint(vm)), vx)));
        result.push_back(data_equation(variable_list({vp, vx, vy}), reduce_fraction_helper(creal(vx, vp), vy), creal(plus(sort_int::cint(sort_nat::pos2nat(vp)), times(vy, vx)), sort_int::int2pos(vx))));
        return result;
      }

    } // namespace sort_real_

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_REAL64_H
