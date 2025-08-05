// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/fbag.h
/// \brief The standard sort fbag.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/fbag.spec.

#ifndef MCRL2_DATA_FBAG64_H
#define MCRL2_DATA_FBAG64_H

#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/container_sort.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/pos64.h"
#include "mcrl2/data/nat64.h"
#include "mcrl2/data/fset64.h"

#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort fbag.
    namespace sort_fbag {

      /// \brief Constructor for sort expression FBag(S)
      /// \param s A sort expression
      /// \return Sort expression fbag(s)
      inline
      container_sort fbag(const sort_expression& s)
      {
        container_sort fbag(fbag_container(), s);
        return fbag;
      }

      /// \brief Recogniser for sort expression FBag(s)
      /// \param e A sort expression
      /// \return true iff e is a container sort of which the name matches
      ///      fbag
      inline
      bool is_fbag(const sort_expression& e)
      {
        if (is_container_sort(e))
        {
          return container_sort(e).container_name() == fbag_container();
        }
        return false;
      }


      /// \brief Generate identifier {:}.
      /// \return Identifier {:}.
      inline
      const core::identifier_string& empty_name()
      {
        static core::identifier_string empty_name = core::identifier_string("{:}");
        return empty_name;
      }

      /// \brief Constructor for function symbol {:}.
      /// \param s A sort expression.
      /// \return Function symbol empty.
      inline
      function_symbol empty(const sort_expression& s)
      {
        function_symbol empty(empty_name(), fbag(s));
        return empty;
      }

      /// \brief Recogniser for function {:}.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching {:}.
      inline
      bool is_empty_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == empty_name();
        }
        return false;
      }

      /// \brief Generate identifier \@fbag_insert.
      /// \return Identifier \@fbag_insert.
      inline
      const core::identifier_string& insert_name()
      {
        static core::identifier_string insert_name = core::identifier_string("@fbag_insert");
        return insert_name;
      }

      /// \brief Constructor for function symbol \@fbag_insert.
      /// \param s A sort expression.
      /// \return Function symbol insert.
      inline
      function_symbol insert(const sort_expression& s)
      {
        function_symbol insert(insert_name(), make_function_sort_(s, sort_pos::pos(), fbag(s), fbag(s)));
        return insert;
      }

      /// \brief Recogniser for function \@fbag_insert.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fbag_insert.
      inline
      bool is_insert_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == insert_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fbag_insert.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@fbag_insert to a number of arguments.
      inline
      application insert(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_fbag::insert(s)(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@fbag_insert.
      /// \param result The data expression where the \@fbag_insert expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_insert(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_fbag::insert(s),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@fbag_insert.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol insert to a
      ///     number of arguments.
      inline
      bool is_insert_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_insert_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for fbag.
      /// \param s A sort expression.
      /// \return All system defined constructors for fbag.
      inline
      function_symbol_vector fbag_generate_constructors_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_fbag::empty(s));
        result.push_back(sort_fbag::insert(s));

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for fbag.
      /// \param s A sort expression.
      /// \return All system defined constructors that can be used in an mCRL2 specification for fbag.
      inline
      function_symbol_vector fbag_mCRL2_usable_constructors(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_fbag::empty(s));
        result.push_back(sort_fbag::insert(s));

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for fbag.
      /// \param s A sort expression.
      /// \return All system defined constructors that are to be implemented in C++ for fbag.
      inline
      implementation_map fbag_cpp_implementable_constructors(const sort_expression& s)
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier \@fbag_cons.
      /// \return Identifier \@fbag_cons.
      inline
      const core::identifier_string& cons_name()
      {
        static core::identifier_string cons_name = core::identifier_string("@fbag_cons");
        return cons_name;
      }

      /// \brief Constructor for function symbol \@fbag_cons.
      /// \param s A sort expression.
      /// \return Function symbol cons_.
      inline
      function_symbol cons_(const sort_expression& s)
      {
        function_symbol cons_(cons_name(), make_function_sort_(s, sort_pos::pos(), fbag(s), fbag(s)));
        return cons_;
      }

      /// \brief Recogniser for function \@fbag_cons.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fbag_cons.
      inline
      bool is_cons_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == cons_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fbag_cons.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@fbag_cons to a number of arguments.
      inline
      application cons_(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_fbag::cons_(s)(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@fbag_cons.
      /// \param result The data expression where the \@fbag_cons expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_cons_(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_fbag::cons_(s),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@fbag_cons.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol cons_ to a
      ///     number of arguments.
      inline
      bool is_cons_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_cons_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fbag_cinsert.
      /// \return Identifier \@fbag_cinsert.
      inline
      const core::identifier_string& cinsert_name()
      {
        static core::identifier_string cinsert_name = core::identifier_string("@fbag_cinsert");
        return cinsert_name;
      }

      /// \brief Constructor for function symbol \@fbag_cinsert.
      /// \param s A sort expression.
      /// \return Function symbol cinsert.
      inline
      function_symbol cinsert(const sort_expression& s)
      {
        function_symbol cinsert(cinsert_name(), make_function_sort_(s, sort_nat::nat(), fbag(s), fbag(s)));
        return cinsert;
      }

      /// \brief Recogniser for function \@fbag_cinsert.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fbag_cinsert.
      inline
      bool is_cinsert_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == cinsert_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fbag_cinsert.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@fbag_cinsert to a number of arguments.
      inline
      application cinsert(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_fbag::cinsert(s)(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@fbag_cinsert.
      /// \param result The data expression where the \@fbag_cinsert expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_cinsert(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_fbag::cinsert(s),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@fbag_cinsert.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol cinsert to a
      ///     number of arguments.
      inline
      bool is_cinsert_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_cinsert_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier count.
      /// \return Identifier count.
      inline
      const core::identifier_string& count_name()
      {
        static core::identifier_string count_name = core::identifier_string("count");
        return count_name;
      }

      /// \brief Constructor for function symbol count.
      /// \param s A sort expression.
      /// \return Function symbol count.
      inline
      function_symbol count(const sort_expression& s)
      {
        function_symbol count(count_name(), make_function_sort_(s, fbag(s), sort_nat::nat()));
        return count;
      }

      /// \brief Recogniser for function count.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching count.
      inline
      bool is_count_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == count_name();
        }
        return false;
      }

      /// \brief Application of function symbol count.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of count to a number of arguments.
      inline
      application count(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_fbag::count(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol count.
      /// \param result The data expression where the count expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_count(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fbag::count(s),arg0, arg1);
      }

      /// \brief Recogniser for application of count.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol count to a
      ///     number of arguments.
      inline
      bool is_count_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_count_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier in.
      /// \return Identifier in.
      inline
      const core::identifier_string& in_name()
      {
        static core::identifier_string in_name = core::identifier_string("in");
        return in_name;
      }

      /// \brief Constructor for function symbol in.
      /// \param s A sort expression.
      /// \return Function symbol in.
      inline
      function_symbol in(const sort_expression& s)
      {
        function_symbol in(in_name(), make_function_sort_(s, fbag(s), sort_bool::bool_()));
        return in;
      }

      /// \brief Recogniser for function in.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching in.
      inline
      bool is_in_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == in_name();
        }
        return false;
      }

      /// \brief Application of function symbol in.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of in to a number of arguments.
      inline
      application in(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_fbag::in(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol in.
      /// \param result The data expression where the in expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_in(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fbag::in(s),arg0, arg1);
      }

      /// \brief Recogniser for application of in.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol in to a
      ///     number of arguments.
      inline
      bool is_in_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_in_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier +.
      /// \return Identifier +.
      inline
      const core::identifier_string& union_name()
      {
        static core::identifier_string union_name = core::identifier_string("+");
        return union_name;
      }

      /// \brief Constructor for function symbol +.
      /// \param s A sort expression.
      /// \return Function symbol union_.
      inline
      function_symbol union_(const sort_expression& s)
      {
        function_symbol union_(union_name(), make_function_sort_(fbag(s), fbag(s), fbag(s)));
        return union_;
      }

      /// \brief Recogniser for function +.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching +.
      inline
      bool is_union_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == union_name();
        }
        return false;
      }

      /// \brief Application of function symbol +.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of + to a number of arguments.
      inline
      application union_(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_fbag::union_(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol +.
      /// \param result The data expression where the + expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_union_(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fbag::union_(s),arg0, arg1);
      }

      /// \brief Recogniser for application of +.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol union_ to a
      ///     number of arguments.
      inline
      bool is_union_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_union_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier *.
      /// \return Identifier *.
      inline
      const core::identifier_string& intersection_name()
      {
        static core::identifier_string intersection_name = core::identifier_string("*");
        return intersection_name;
      }

      /// \brief Constructor for function symbol *.
      /// \param s A sort expression.
      /// \return Function symbol intersection.
      inline
      function_symbol intersection(const sort_expression& s)
      {
        function_symbol intersection(intersection_name(), make_function_sort_(fbag(s), fbag(s), fbag(s)));
        return intersection;
      }

      /// \brief Recogniser for function *.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching *.
      inline
      bool is_intersection_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == intersection_name();
        }
        return false;
      }

      /// \brief Application of function symbol *.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of * to a number of arguments.
      inline
      application intersection(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_fbag::intersection(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol *.
      /// \param result The data expression where the * expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_intersection(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fbag::intersection(s),arg0, arg1);
      }

      /// \brief Recogniser for application of *.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol intersection to a
      ///     number of arguments.
      inline
      bool is_intersection_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_intersection_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier -.
      /// \return Identifier -.
      inline
      const core::identifier_string& difference_name()
      {
        static core::identifier_string difference_name = core::identifier_string("-");
        return difference_name;
      }

      /// \brief Constructor for function symbol -.
      /// \param s A sort expression.
      /// \return Function symbol difference.
      inline
      function_symbol difference(const sort_expression& s)
      {
        function_symbol difference(difference_name(), make_function_sort_(fbag(s), fbag(s), fbag(s)));
        return difference;
      }

      /// \brief Recogniser for function -.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching -.
      inline
      bool is_difference_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == difference_name();
        }
        return false;
      }

      /// \brief Application of function symbol -.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of - to a number of arguments.
      inline
      application difference(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_fbag::difference(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol -.
      /// \param result The data expression where the - expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_difference(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fbag::difference(s),arg0, arg1);
      }

      /// \brief Recogniser for application of -.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol difference to a
      ///     number of arguments.
      inline
      bool is_difference_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_difference_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier #.
      /// \return Identifier #.
      inline
      const core::identifier_string& count_all_name()
      {
        static core::identifier_string count_all_name = core::identifier_string("#");
        return count_all_name;
      }

      /// \brief Constructor for function symbol #.
      /// \param s A sort expression.
      /// \return Function symbol count_all.
      inline
      function_symbol count_all(const sort_expression& s)
      {
        function_symbol count_all(count_all_name(), make_function_sort_(fbag(s), sort_nat::nat()));
        return count_all;
      }

      /// \brief Recogniser for function #.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching #.
      inline
      bool is_count_all_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == count_all_name();
        }
        return false;
      }

      /// \brief Application of function symbol #.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of # to a number of arguments.
      inline
      application count_all(const sort_expression& s, const data_expression& arg0)
      {
        return sort_fbag::count_all(s)(arg0);
      }

      /// \brief Make an application of function symbol #.
      /// \param result The data expression where the # expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_count_all(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_fbag::count_all(s),arg0);
      }

      /// \brief Recogniser for application of #.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol count_all to a
      ///     number of arguments.
      inline
      bool is_count_all_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_count_all_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier pick.
      /// \return Identifier pick.
      inline
      const core::identifier_string& pick_name()
      {
        static core::identifier_string pick_name = core::identifier_string("pick");
        return pick_name;
      }

      /// \brief Constructor for function symbol pick.
      /// \param s A sort expression.
      /// \return Function symbol pick.
      inline
      function_symbol pick(const sort_expression& s)
      {
        function_symbol pick(pick_name(), make_function_sort_(fbag(s), s));
        return pick;
      }

      /// \brief Recogniser for function pick.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching pick.
      inline
      bool is_pick_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == pick_name();
        }
        return false;
      }

      /// \brief Application of function symbol pick.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of pick to a number of arguments.
      inline
      application pick(const sort_expression& s, const data_expression& arg0)
      {
        return sort_fbag::pick(s)(arg0);
      }

      /// \brief Make an application of function symbol pick.
      /// \param result The data expression where the pick expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_pick(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_fbag::pick(s),arg0);
      }

      /// \brief Recogniser for application of pick.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol pick to a
      ///     number of arguments.
      inline
      bool is_pick_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_pick_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for fbag
      /// \param s A sort expression
      /// \return All system defined mappings for fbag
      inline
      function_symbol_vector fbag_generate_functions_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_fbag::cons_(s));
        result.push_back(sort_fbag::cinsert(s));
        result.push_back(sort_fbag::count(s));
        result.push_back(sort_fbag::in(s));
        result.push_back(sort_fbag::union_(s));
        result.push_back(sort_fbag::intersection(s));
        result.push_back(sort_fbag::difference(s));
        result.push_back(sort_fbag::count_all(s));
        result.push_back(sort_fbag::pick(s));
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for fbag
      /// \param s A sort expression
      /// \return All system defined mappings for fbag
      inline
      function_symbol_vector fbag_generate_constructors_and_functions_code(const sort_expression& s)
      {
        function_symbol_vector result=fbag_generate_functions_code(s);
        for(const function_symbol& f: fbag_generate_constructors_code(s))
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for fbag
      /// \param s A sort expression
      /// \return All system defined mappings for that can be used in mCRL2 specificationis fbag
      inline
      function_symbol_vector fbag_mCRL2_usable_mappings(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_fbag::cons_(s));
        result.push_back(sort_fbag::cinsert(s));
        result.push_back(sort_fbag::count(s));
        result.push_back(sort_fbag::in(s));
        result.push_back(sort_fbag::union_(s));
        result.push_back(sort_fbag::intersection(s));
        result.push_back(sort_fbag::difference(s));
        result.push_back(sort_fbag::count_all(s));
        result.push_back(sort_fbag::pick(s));
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for fbag
      /// \param s A sort expression
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for fbag
      inline
      implementation_map fbag_cpp_implementable_mappings(const sort_expression& s)
      {
        implementation_map result;
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
        assert(is_insert_application(e) || is_cons_application(e) || is_cinsert_application(e));
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
        assert(is_insert_application(e) || is_cons_application(e) || is_cinsert_application(e));
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
        assert(is_insert_application(e) || is_cons_application(e) || is_cinsert_application(e));
        return atermpp::down_cast<application>(e)[2];
      }

      ///\brief Function for projecting out argument.
      ///        left from an application.
      /// \param e A data expression.
      /// \pre left is defined for e.
      /// \return The argument of e that corresponds to left.
      inline
      const data_expression& left(const data_expression& e)
      {
        assert(is_count_application(e) || is_in_application(e) || is_union_application(e) || is_intersection_application(e) || is_difference_application(e));
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
        assert(is_count_application(e) || is_in_application(e) || is_union_application(e) || is_intersection_application(e) || is_difference_application(e));
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
        assert(is_count_all_application(e) || is_pick_application(e));
        return atermpp::down_cast<application>(e)[0];
      }

      /// \brief Give all system defined equations for fbag
      /// \param s A sort expression
      /// \return All system defined equations for sort fbag
      inline
      data_equation_vector fbag_generate_equations_code(const sort_expression& s)
      {
        variable vd("d",s);
        variable ve("e",s);
        variable vp("p",sort_pos::pos());
        variable vq("q",sort_pos::pos());
        variable vn("n",sort_nat::nat());
        variable vb("b",fbag(s));
        variable vc("c",fbag(s));
        variable vf("f",make_function_sort_(s, sort_nat::nat()));
        variable vg("g",make_function_sort_(s, sort_nat::nat()));

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vb, vd, vp}), equal_to(cons_(s, vd, vp, vb), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vb, vd, vp}), equal_to(empty(s), cons_(s, vd, vp, vb)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), equal_to(cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), sort_bool::and_(equal_to(vp, vq), sort_bool::and_(equal_to(vd, ve), equal_to(vb, vc)))));
        result.push_back(data_equation(variable_list({vb, vd, vp}), less_equal(cons_(s, vd, vp, vb), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vb, vd, vp}), less_equal(empty(s), cons_(s, vd, vp, vb)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), less_equal(cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), if_(less(vd, ve), sort_bool::false_(), if_(equal_to(vd, ve), sort_bool::and_(less_equal(vp, vq), less_equal(vb, vc)), less_equal(cons_(s, vd, vp, vb), vc)))));
        result.push_back(data_equation(variable_list({vb, vd, vp}), less(cons_(s, vd, vp, vb), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vb, vd, vp}), less(empty(s), cons_(s, vd, vp, vb)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), less(cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), if_(less(vd, ve), sort_bool::false_(), if_(equal_to(vd, ve), sort_bool::or_(sort_bool::and_(equal_to(vp, vq), less(vb, vc)), sort_bool::and_(less(vp, vq), less_equal(vb, vc))), less_equal(cons_(s, vd, vp, vb), vc)))));
        result.push_back(data_equation(variable_list({vd, vp}), insert(s, vd, vp, empty(s)), cons_(s, vd, vp, empty(s))));
        result.push_back(data_equation(variable_list({vb, vd, vp, vq}), insert(s, vd, vp, cons_(s, vd, vq, vb)), cons_(s, vd, sort_pos::auxiliary_plus_pos(vp, vq), vb)));
        result.push_back(data_equation(variable_list({vb, vd, ve, vp, vq}), less(vd, ve), insert(s, vd, vp, cons_(s, ve, vq, vb)), cons_(s, vd, vp, cons_(s, ve, vq, vb))));
        result.push_back(data_equation(variable_list({vb, vd, ve, vp, vq}), less(ve, vd), insert(s, vd, vp, cons_(s, ve, vq, vb)), cons_(s, ve, vq, insert(s, vd, vp, vb))));
        result.push_back(data_equation(variable_list({vb, vd, vn}), cinsert(s, vd, vn, vb), if_(equal_to(vn, sort_nat::most_significant_digit_nat(sort_machine_word::zero_word())), vb, insert(s, vd, sort_nat::nat2pos(vn), vb))));
        result.push_back(data_equation(variable_list({vd}), count(s, vd, empty(s)), sort_nat::c0()));
        result.push_back(data_equation(variable_list({vb, vd, vp}), count(s, vd, cons_(s, vd, vp, vb)), sort_nat::pos2nat(vp)));
        result.push_back(data_equation(variable_list({vb, vd, ve, vp}), less(vd, ve), count(s, vd, cons_(s, ve, vp, vb)), sort_nat::c0()));
        result.push_back(data_equation(variable_list({vb, vd, ve, vp}), less(ve, vd), count(s, vd, cons_(s, ve, vp, vb)), count(s, vd, vb)));
        result.push_back(data_equation(variable_list({vb, vd}), in(s, vd, vb), greater(count(s, vd, vb), sort_nat::c0())));
        result.push_back(data_equation(variable_list({vb}), difference(s, vb, empty(s)), vb));
        result.push_back(data_equation(variable_list({vc}), difference(s, empty(s), vc), empty(s)));
        result.push_back(data_equation(variable_list({vb, vc, vd, vp}), difference(s, cons_(s, vd, vp, vb), cons_(s, vd, vp, vc)), difference(s, vb, vc)));
        result.push_back(data_equation(variable_list({vb, vc, vd, vp, vq}), less(vp, vq), difference(s, cons_(s, vd, vp, vb), cons_(s, vd, vq, vc)), difference(s, vb, vc)));
        result.push_back(data_equation(variable_list({vb, vc, vd, vp, vq}), less(vq, vp), difference(s, cons_(s, vd, vp, vb), cons_(s, vd, vq, vc)), cons_(s, vd, sort_nat::nat2pos(sort_nat::monus(sort_nat::pos2nat(vp), sort_nat::pos2nat(vq))), difference(s, vb, vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), less(vd, ve), difference(s, cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), cons_(s, vd, vp, difference(s, vb, cons_(s, ve, vq, vc)))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), less(ve, vd), difference(s, cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), difference(s, cons_(s, vd, vp, vb), vc)));
        result.push_back(data_equation(variable_list({vb}), union_(s, vb, empty(s)), vb));
        result.push_back(data_equation(variable_list({vc}), union_(s, empty(s), vc), vc));
        result.push_back(data_equation(variable_list({vb, vc, vd, vp, vq}), union_(s, cons_(s, vd, vp, vb), cons_(s, vd, vq, vc)), cons_(s, vd, sort_pos::auxiliary_plus_pos(vp, vq), union_(s, vb, vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), less(vd, ve), union_(s, cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), cons_(s, vd, vp, union_(s, vb, cons_(s, ve, vq, vc)))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), less(ve, vd), union_(s, cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), cons_(s, ve, vq, union_(s, cons_(s, vd, vp, vb), vc))));
        result.push_back(data_equation(variable_list({vb}), intersection(s, vb, empty(s)), empty(s)));
        result.push_back(data_equation(variable_list({vc}), intersection(s, empty(s), vc), empty(s)));
        result.push_back(data_equation(variable_list({vb, vc, vd, vp, vq}), intersection(s, cons_(s, vd, vp, vb), cons_(s, vd, vq, vc)), cons_(s, vd, sort_nat::minimum(vp, vq), intersection(s, vb, vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), less(vd, ve), intersection(s, cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), intersection(s, vb, cons_(s, ve, vq, vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vp, vq}), less(ve, vd), intersection(s, cons_(s, vd, vp, vb), cons_(s, ve, vq, vc)), intersection(s, cons_(s, vd, vp, vb), vc)));
        result.push_back(data_equation(variable_list(), count_all(s, empty(s)), sort_nat::c0()));
        result.push_back(data_equation(variable_list({vd, vp}), count_all(s, cons_(s, vd, vp, empty(s))), sort_nat::pos2nat(vp)));
        result.push_back(data_equation(variable_list({vb, vd, ve, vp, vq}), count_all(s, cons_(s, vd, vp, cons_(s, ve, vq, vb))), sort_nat::pos2nat(sort_pos::auxiliary_plus_pos(vp, sort_nat::nat2pos(count_all(s, cons_(s, ve, vq, vb)))))));
        result.push_back(data_equation(variable_list({vb, vd, vp}), pick(s, cons_(s, vd, vp, vb)), vd));
        return result;
      }

    } // namespace sort_fbag

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_FBAG64_H
