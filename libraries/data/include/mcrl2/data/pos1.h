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

#ifndef MCRL2_DATA_POS1_H
#define MCRL2_DATA_POS1_H

#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/bool.h"

#include <functional>
#include <map>

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

      /// \brief Generate identifier \@cDub.
      /// \return Identifier \@cDub.
      inline
      const core::identifier_string& cdub_name()
      {
        static core::identifier_string cdub_name = core::identifier_string("@cDub");
        return cdub_name;
      }

      /// \brief Constructor for function symbol \@cDub.
      
      /// \return Function symbol cdub.
      inline
      const function_symbol& cdub()
      {
        static function_symbol cdub(cdub_name(), make_function_sort_(sort_bool::bool_(), pos(), pos()));
        return cdub;
      }

      /// \brief Recogniser for function \@cDub.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@cDub.
      inline
      bool is_cdub_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == cdub();
        }
        return false;
      }

      /// \brief Application of function symbol \@cDub.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@cDub to a number of arguments.
      inline
      application cdub(const data_expression& arg0, const data_expression& arg1)
      {
        return sort_pos::cdub()(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@cDub.
      /// \param result The data expression where the \@cDub expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_cdub(data_expression& result, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_pos::cdub(),arg0, arg1);
      }

      /// \brief Recogniser for application of \@cDub.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol cdub to a
      ///     number of arguments.
      inline
      bool is_cdub_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_cdub_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for pos.
      /// \return All system defined constructors for pos.
      inline
      function_symbol_vector pos_generate_constructors_code()
      {
        function_symbol_vector result;
        result.push_back(sort_pos::c1());
        result.push_back(sort_pos::cdub());

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for pos.
      /// \return All system defined constructors that can be used in an mCRL2 specification for pos.
      inline
      function_symbol_vector pos_mCRL2_usable_constructors()
      {
        function_symbol_vector result;
        result.push_back(sort_pos::c1());
        result.push_back(sort_pos::cdub());

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

      /// \brief Generate identifier \@addc.
      /// \return Identifier \@addc.
      inline
      const core::identifier_string& add_with_carry_name()
      {
        static core::identifier_string add_with_carry_name = core::identifier_string("@addc");
        return add_with_carry_name;
      }

      /// \brief Constructor for function symbol \@addc.
      
      /// \return Function symbol add_with_carry.
      inline
      const function_symbol& add_with_carry()
      {
        static function_symbol add_with_carry(add_with_carry_name(), make_function_sort_(sort_bool::bool_(), pos(), pos(), pos()));
        return add_with_carry;
      }

      /// \brief Recogniser for function \@addc.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@addc.
      inline
      bool is_add_with_carry_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == add_with_carry();
        }
        return false;
      }

      /// \brief Application of function symbol \@addc.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@addc to a number of arguments.
      inline
      application add_with_carry(const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_pos::add_with_carry()(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@addc.
      /// \param result The data expression where the \@addc expression is put.
      
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_add_with_carry(data_expression& result, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_pos::add_with_carry(),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@addc.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol add_with_carry to a
      ///     number of arguments.
      inline
      bool is_add_with_carry_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_add_with_carry_function_symbol(atermpp::down_cast<application>(e).head());
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

      /// \brief Generate identifier \@powerlog2.
      /// \return Identifier \@powerlog2.
      inline
      const core::identifier_string& powerlog2_pos_name()
      {
        static core::identifier_string powerlog2_pos_name = core::identifier_string("@powerlog2");
        return powerlog2_pos_name;
      }

      /// \brief Constructor for function symbol \@powerlog2.
      
      /// \return Function symbol powerlog2_pos.
      inline
      const function_symbol& powerlog2_pos()
      {
        static function_symbol powerlog2_pos(powerlog2_pos_name(), make_function_sort_(pos(), pos()));
        return powerlog2_pos;
      }

      /// \brief Recogniser for function \@powerlog2.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@powerlog2.
      inline
      bool is_powerlog2_pos_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e) == powerlog2_pos();
        }
        return false;
      }

      /// \brief Application of function symbol \@powerlog2.
      
      /// \param arg0 A data expression.
      /// \return Application of \@powerlog2 to a number of arguments.
      inline
      application powerlog2_pos(const data_expression& arg0)
      {
        return sort_pos::powerlog2_pos()(arg0);
      }

      /// \brief Make an application of function symbol \@powerlog2.
      /// \param result The data expression where the \@powerlog2 expression is put.
      
      /// \param arg0 A data expression.
      inline
      void make_powerlog2_pos(data_expression& result, const data_expression& arg0)
      {
        make_application(result, sort_pos::powerlog2_pos(),arg0);
      }

      /// \brief Recogniser for application of \@powerlog2.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol powerlog2_pos to a
      ///     number of arguments.
      inline
      bool is_powerlog2_pos_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_powerlog2_pos_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for pos
      /// \return All system defined mappings for pos
      inline
      function_symbol_vector pos_generate_functions_code()
      {
        function_symbol_vector result;
        result.push_back(sort_pos::maximum());
        result.push_back(sort_pos::minimum());
        result.push_back(sort_pos::succ());
        result.push_back(sort_pos::pos_predecessor());
        result.push_back(sort_pos::plus());
        result.push_back(sort_pos::add_with_carry());
        result.push_back(sort_pos::times());
        result.push_back(sort_pos::powerlog2_pos());
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
        result.push_back(sort_pos::maximum());
        result.push_back(sort_pos::minimum());
        result.push_back(sort_pos::succ());
        result.push_back(sort_pos::pos_predecessor());
        result.push_back(sort_pos::plus());
        result.push_back(sort_pos::add_with_carry());
        result.push_back(sort_pos::times());
        result.push_back(sort_pos::powerlog2_pos());
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
      ///        left from an application.
      /// \param e A data expression.
      /// \pre left is defined for e.
      /// \return The argument of e that corresponds to left.
      inline
      const data_expression& left(const data_expression& e)
      {
        assert(is_cdub_application(e) || is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_times_application(e));
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
        assert(is_cdub_application(e) || is_maximum_application(e) || is_minimum_application(e) || is_plus_application(e) || is_times_application(e));
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
        assert(is_succ_application(e) || is_pos_predecessor_application(e) || is_powerlog2_pos_application(e));
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
        assert(is_add_with_carry_application(e));
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
        assert(is_add_with_carry_application(e));
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
        assert(is_add_with_carry_application(e));
        return atermpp::down_cast<application>(e)[2];
      }

      /// \brief Give all system defined equations for pos
      /// \return All system defined equations for sort pos
      inline
      data_equation_vector pos_generate_equations_code()
      {
        variable vb("b",sort_bool::bool_());
        variable vc("c",sort_bool::bool_());
        variable vp("p",pos());
        variable vq("q",pos());
        variable vp1("p1",pos());
        variable vq1("q1",pos());

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vb, vp}), equal_to(c1(), cdub(vb, vp)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vb, vp}), equal_to(cdub(vb, vp), c1()), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vb, vp, vq}), equal_to(cdub(vb, vp), cdub(vb, vq)), equal_to(vp, vq)));
        result.push_back(data_equation(variable_list({vp, vq}), equal_to(cdub(sort_bool::false_(), vp), cdub(sort_bool::true_(), vq)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp, vq}), equal_to(cdub(sort_bool::true_(), vp), cdub(sort_bool::false_(), vq)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp}), equal_to(succ(vp), c1()), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vq}), equal_to(c1(), succ(vq)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vc, vp, vq}), equal_to(succ(vp), cdub(vc, vq)), equal_to(vp, pos_predecessor(cdub(vc, vq)))));
        result.push_back(data_equation(variable_list({vb, vp, vq}), equal_to(cdub(vb, vp), succ(vq)), equal_to(pos_predecessor(cdub(vb, vp)), vq)));
        result.push_back(data_equation(variable_list({vp}), less(vp, c1()), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vb, vp}), less(c1(), cdub(vb, vp)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vb, vc, vp, vq}), less(cdub(vb, vp), cdub(vc, vq)), if_(sort_bool::implies(vc, vb), less(vp, vq), less_equal(vp, vq))));
        result.push_back(data_equation(variable_list({vc, vp, vq}), less(succ(vp), cdub(vc, vq)), less(vp, pos_predecessor(cdub(vc, vq)))));
        result.push_back(data_equation(variable_list({vb, vp, vq}), less(cdub(vb, vp), succ(vq)), less_equal(cdub(vb, vp), vq)));
        result.push_back(data_equation(variable_list({vq}), less(c1(), succ(vq)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vp}), less_equal(c1(), vp), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vb, vp}), less_equal(cdub(vb, vp), c1()), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vb, vc, vp, vq}), less_equal(cdub(vb, vp), cdub(vc, vq)), if_(sort_bool::implies(vb, vc), less_equal(vp, vq), less(vp, vq))));
        result.push_back(data_equation(variable_list({vc, vp, vq}), less_equal(succ(vp), cdub(vc, vq)), less(vp, cdub(vc, vq))));
        result.push_back(data_equation(variable_list({vb, vp, vq}), less_equal(cdub(vb, vp), succ(vq)), less_equal(pos_predecessor(cdub(vb, vp)), vq)));
        result.push_back(data_equation(variable_list({vp}), less_equal(succ(vp), c1()), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vp, vq}), maximum(vp, vq), if_(less_equal(vp, vq), vq, vp)));
        result.push_back(data_equation(variable_list({vp, vq}), minimum(vp, vq), if_(less_equal(vp, vq), vp, vq)));
        result.push_back(data_equation(variable_list(), succ(c1()), cdub(sort_bool::false_(), c1())));
        result.push_back(data_equation(variable_list({vp}), succ(cdub(sort_bool::false_(), vp)), cdub(sort_bool::true_(), vp)));
        result.push_back(data_equation(variable_list({vp}), succ(cdub(sort_bool::true_(), vp)), cdub(sort_bool::false_(), succ(vp))));
        result.push_back(data_equation(variable_list(), pos_predecessor(c1()), c1()));
        result.push_back(data_equation(variable_list(), pos_predecessor(cdub(sort_bool::false_(), c1())), c1()));
        result.push_back(data_equation(variable_list({vb, vp}), pos_predecessor(cdub(sort_bool::false_(), cdub(vb, vp))), cdub(sort_bool::true_(), pos_predecessor(cdub(vb, vp)))));
        result.push_back(data_equation(variable_list({vp}), pos_predecessor(cdub(sort_bool::true_(), vp)), cdub(sort_bool::false_(), vp)));
        result.push_back(data_equation(variable_list({vp, vq}), plus(vp, vq), add_with_carry(sort_bool::false_(), vp, vq)));
        result.push_back(data_equation(variable_list({vp}), add_with_carry(sort_bool::false_(), c1(), vp), succ(vp)));
        result.push_back(data_equation(variable_list({vp}), add_with_carry(sort_bool::true_(), c1(), vp), succ(succ(vp))));
        result.push_back(data_equation(variable_list({vp}), add_with_carry(sort_bool::false_(), vp, c1()), succ(vp)));
        result.push_back(data_equation(variable_list({vp}), add_with_carry(sort_bool::true_(), vp, c1()), succ(succ(vp))));
        result.push_back(data_equation(variable_list({vb, vc, vp, vq}), add_with_carry(vb, cdub(vc, vp), cdub(vc, vq)), cdub(vb, add_with_carry(vc, vp, vq))));
        result.push_back(data_equation(variable_list({vb, vp, vq}), add_with_carry(vb, cdub(sort_bool::false_(), vp), cdub(sort_bool::true_(), vq)), cdub(sort_bool::not_(vb), add_with_carry(vb, vp, vq))));
        result.push_back(data_equation(variable_list({vb, vp, vq}), add_with_carry(vb, cdub(sort_bool::true_(), vp), cdub(sort_bool::false_(), vq)), cdub(sort_bool::not_(vb), add_with_carry(vb, vp, vq))));
        result.push_back(data_equation(variable_list({vp}), times(c1(), vp), vp));
        result.push_back(data_equation(variable_list({vp}), times(vp, c1()), vp));
        result.push_back(data_equation(variable_list({vp, vq}), times(cdub(sort_bool::false_(), vp), vq), cdub(sort_bool::false_(), times(vp, vq))));
        result.push_back(data_equation(variable_list({vp, vq}), times(vp, cdub(sort_bool::false_(), vq)), cdub(sort_bool::false_(), times(vp, vq))));
        result.push_back(data_equation(variable_list({vp, vq}), times(cdub(sort_bool::true_(), vp), cdub(sort_bool::true_(), vq)), cdub(sort_bool::true_(), add_with_carry(sort_bool::false_(), vp, add_with_carry(sort_bool::false_(), vq, cdub(sort_bool::false_(), times(vp, vq)))))));
        result.push_back(data_equation(variable_list(), powerlog2_pos(c1()), c1()));
        result.push_back(data_equation(variable_list({vb}), powerlog2_pos(cdub(vb, c1())), c1()));
        result.push_back(data_equation(variable_list({vb, vc, vp}), powerlog2_pos(cdub(vb, cdub(vc, vp))), cdub(sort_bool::false_(), powerlog2_pos(vp))));
        return result;
      }

    } // namespace sort_pos

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_POS1_H
