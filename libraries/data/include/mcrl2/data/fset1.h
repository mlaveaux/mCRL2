// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/fset.h
/// \brief The standard sort fset.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/fset.spec.

#ifndef MCRL2_DATA_FSET1_H
#define MCRL2_DATA_FSET1_H

#include "mcrl2/utilities/exception.h"
#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/container_sort.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/nat1.h"

#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort fset.
    namespace sort_fset {

      /// \brief Constructor for sort expression FSet(S)
      /// \param s A sort expression
      /// \return Sort expression fset(s)
      inline
      container_sort fset(const sort_expression& s)
      {
        container_sort fset(fset_container(), s);
        return fset;
      }

      /// \brief Recogniser for sort expression FSet(s)
      /// \param e A sort expression
      /// \return true iff e is a container sort of which the name matches
      ///      fset
      inline
      bool is_fset(const sort_expression& e)
      {
        if (is_container_sort(e))
        {
          return container_sort(e).container_name() == fset_container();
        }
        return false;
      }


      /// \brief Generate identifier {}.
      /// \return Identifier {}.
      inline
      const core::identifier_string& empty_name()
      {
        static core::identifier_string empty_name = core::identifier_string("{}");
        return empty_name;
      }

      /// \brief Constructor for function symbol {}.
      /// \param s A sort expression.
      /// \return Function symbol empty.
      inline
      function_symbol empty(const sort_expression& s)
      {
        function_symbol empty(empty_name(), fset(s));
        return empty;
      }

      /// \brief Recogniser for function {}.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching {}.
      inline
      bool is_empty_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == empty_name();
        }
        return false;
      }

      /// \brief Generate identifier \@fset_insert.
      /// \return Identifier \@fset_insert.
      inline
      const core::identifier_string& insert_name()
      {
        static core::identifier_string insert_name = core::identifier_string("@fset_insert");
        return insert_name;
      }

      /// \brief Constructor for function symbol \@fset_insert.
      /// \param s A sort expression.
      /// \return Function symbol insert.
      inline
      function_symbol insert(const sort_expression& s)
      {
        function_symbol insert(insert_name(), make_function_sort_(s, fset(s), fset(s)));
        return insert;
      }

      /// \brief Recogniser for function \@fset_insert.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fset_insert.
      inline
      bool is_insert_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == insert_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fset_insert.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@fset_insert to a number of arguments.
      inline
      application insert(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_fset::insert(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@fset_insert.
      /// \param result The data expression where the \@fset_insert expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_insert(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fset::insert(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@fset_insert.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol insert to a
      ///     number of arguments.
      inline
      bool is_insert_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_insert_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for fset.
      /// \param s A sort expression.
      /// \return All system defined constructors for fset.
      inline
      function_symbol_vector fset_generate_constructors_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_fset::empty(s));
        result.push_back(sort_fset::insert(s));

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for fset.
      /// \param s A sort expression.
      /// \return All system defined constructors that can be used in an mCRL2 specification for fset.
      inline
      function_symbol_vector fset_mCRL2_usable_constructors(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_fset::empty(s));
        result.push_back(sort_fset::insert(s));

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for fset.
      /// \param s A sort expression.
      /// \return All system defined constructors that are to be implemented in C++ for fset.
      inline
      implementation_map fset_cpp_implementable_constructors(const sort_expression& s)
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier \@fset_cons.
      /// \return Identifier \@fset_cons.
      inline
      const core::identifier_string& cons_name()
      {
        static core::identifier_string cons_name = core::identifier_string("@fset_cons");
        return cons_name;
      }

      /// \brief Constructor for function symbol \@fset_cons.
      /// \param s A sort expression.
      /// \return Function symbol cons_.
      inline
      function_symbol cons_(const sort_expression& s)
      {
        function_symbol cons_(cons_name(), make_function_sort_(s, fset(s), fset(s)));
        return cons_;
      }

      /// \brief Recogniser for function \@fset_cons.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fset_cons.
      inline
      bool is_cons_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == cons_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fset_cons.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@fset_cons to a number of arguments.
      inline
      application cons_(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_fset::cons_(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@fset_cons.
      /// \param result The data expression where the \@fset_cons expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_cons_(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fset::cons_(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@fset_cons.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol cons_ to a
      ///     number of arguments.
      inline
      bool is_cons_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_cons_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fset_cinsert.
      /// \return Identifier \@fset_cinsert.
      inline
      const core::identifier_string& cinsert_name()
      {
        static core::identifier_string cinsert_name = core::identifier_string("@fset_cinsert");
        return cinsert_name;
      }

      /// \brief Constructor for function symbol \@fset_cinsert.
      /// \param s A sort expression.
      /// \return Function symbol cinsert.
      inline
      function_symbol cinsert(const sort_expression& s)
      {
        function_symbol cinsert(cinsert_name(), make_function_sort_(s, sort_bool::bool_(), fset(s), fset(s)));
        return cinsert;
      }

      /// \brief Recogniser for function \@fset_cinsert.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fset_cinsert.
      inline
      bool is_cinsert_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == cinsert_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fset_cinsert.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \return Application of \@fset_cinsert to a number of arguments.
      inline
      application cinsert(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        return sort_fset::cinsert(s)(arg0, arg1, arg2);
      }

      /// \brief Make an application of function symbol \@fset_cinsert.
      /// \param result The data expression where the \@fset_cinsert expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      inline
      void make_cinsert(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2)
      {
        make_application(result, sort_fset::cinsert(s),arg0, arg1, arg2);
      }

      /// \brief Recogniser for application of \@fset_cinsert.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol cinsert to a
      ///     number of arguments.
      inline
      bool is_cinsert_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_cinsert_function_symbol(atermpp::down_cast<application>(e).head());
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
        function_symbol in(in_name(), make_function_sort_(s, fset(s), sort_bool::bool_()));
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
        return sort_fset::in(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol in.
      /// \param result The data expression where the in expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_in(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fset::in(s),arg0, arg1);
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
        function_symbol difference(difference_name(), make_function_sort_(fset(s), fset(s), fset(s)));
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
        return sort_fset::difference(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol -.
      /// \param result The data expression where the - expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_difference(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fset::difference(s),arg0, arg1);
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
        function_symbol union_(union_name(), make_function_sort_(fset(s), fset(s), fset(s)));
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
        return sort_fset::union_(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol +.
      /// \param result The data expression where the + expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_union_(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fset::union_(s),arg0, arg1);
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
        function_symbol intersection(intersection_name(), make_function_sort_(fset(s), fset(s), fset(s)));
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
        return sort_fset::intersection(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol *.
      /// \param result The data expression where the * expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_intersection(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_fset::intersection(s),arg0, arg1);
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

      /// \brief Generate identifier #.
      /// \return Identifier #.
      inline
      const core::identifier_string& count_name()
      {
        static core::identifier_string count_name = core::identifier_string("#");
        return count_name;
      }

      /// \brief Constructor for function symbol #.
      /// \param s A sort expression.
      /// \return Function symbol count.
      inline
      function_symbol count(const sort_expression& s)
      {
        function_symbol count(count_name(), make_function_sort_(fset(s), sort_nat::nat()));
        return count;
      }

      /// \brief Recogniser for function #.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching #.
      inline
      bool is_count_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == count_name();
        }
        return false;
      }

      /// \brief Application of function symbol #.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of # to a number of arguments.
      inline
      application count(const sort_expression& s, const data_expression& arg0)
      {
        return sort_fset::count(s)(arg0);
      }

      /// \brief Make an application of function symbol #.
      /// \param result The data expression where the # expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_count(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_fset::count(s),arg0);
      }

      /// \brief Recogniser for application of #.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol count to a
      ///     number of arguments.
      inline
      bool is_count_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_count_function_symbol(atermpp::down_cast<application>(e).head());
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
        function_symbol pick(pick_name(), make_function_sort_(fset(s), s));
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
        return sort_fset::pick(s)(arg0);
      }

      /// \brief Make an application of function symbol pick.
      /// \param result The data expression where the pick expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_pick(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_fset::pick(s),arg0);
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
      /// \brief Give all system defined mappings for fset
      /// \param s A sort expression
      /// \return All system defined mappings for fset
      inline
      function_symbol_vector fset_generate_functions_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_fset::cons_(s));
        result.push_back(sort_fset::cinsert(s));
        result.push_back(sort_fset::in(s));
        result.push_back(sort_fset::difference(s));
        result.push_back(sort_fset::union_(s));
        result.push_back(sort_fset::intersection(s));
        result.push_back(sort_fset::count(s));
        result.push_back(sort_fset::pick(s));
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for fset
      /// \param s A sort expression
      /// \return All system defined mappings for fset
      inline
      function_symbol_vector fset_generate_constructors_and_functions_code(const sort_expression& s)
      {
        function_symbol_vector result=fset_generate_functions_code(s);
        for(const function_symbol& f: fset_generate_constructors_code(s))
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for fset
      /// \param s A sort expression
      /// \return All system defined mappings for that can be used in mCRL2 specificationis fset
      inline
      function_symbol_vector fset_mCRL2_usable_mappings(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_fset::cons_(s));
        result.push_back(sort_fset::cinsert(s));
        result.push_back(sort_fset::in(s));
        result.push_back(sort_fset::difference(s));
        result.push_back(sort_fset::union_(s));
        result.push_back(sort_fset::intersection(s));
        result.push_back(sort_fset::count(s));
        result.push_back(sort_fset::pick(s));
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for fset
      /// \param s A sort expression
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for fset
      inline
      implementation_map fset_cpp_implementable_mappings(const sort_expression& s)
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
        assert(is_insert_application(e) || is_cons_application(e) || is_in_application(e) || is_difference_application(e) || is_union_application(e) || is_intersection_application(e));
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
        assert(is_insert_application(e) || is_cons_application(e) || is_in_application(e) || is_difference_application(e) || is_union_application(e) || is_intersection_application(e));
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
        assert(is_cinsert_application(e));
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
        assert(is_cinsert_application(e));
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
        assert(is_cinsert_application(e));
        return atermpp::down_cast<application>(e)[2];
      }

      ///\brief Function for projecting out argument.
      ///        arg from an application.
      /// \param e A data expression.
      /// \pre arg is defined for e.
      /// \return The argument of e that corresponds to arg.
      inline
      const data_expression& arg(const data_expression& e)
      {
        assert(is_count_application(e) || is_pick_application(e));
        return atermpp::down_cast<application>(e)[0];
      }

      /// \brief Give all system defined equations for fset
      /// \param s A sort expression
      /// \return All system defined equations for sort fset
      inline
      data_equation_vector fset_generate_equations_code(const sort_expression& s)
      {
        variable vd("d",s);
        variable ve("e",s);
        variable vs("s",fset(s));
        variable vt("t",fset(s));

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vd, vs}), equal_to(empty(s), cons_(s, vd, vs)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, vs}), equal_to(cons_(s, vd, vs), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), equal_to(cons_(s, vd, vs), cons_(s, ve, vt)), sort_bool::and_(equal_to(vd, ve), equal_to(vs, vt))));
        result.push_back(data_equation(variable_list({vd, vs}), less_equal(empty(s), cons_(s, vd, vs)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vd, vs}), less_equal(cons_(s, vd, vs), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less_equal(cons_(s, vd, vs), cons_(s, ve, vt)), if_(less(vd, ve), sort_bool::false_(), if_(equal_to(vd, ve), less_equal(vs, vt), less_equal(cons_(s, vd, vs), vt)))));
        result.push_back(data_equation(variable_list({vd, vs}), less(empty(s), cons_(s, vd, vs)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vd, vs}), less(cons_(s, vd, vs), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less(cons_(s, vd, vs), cons_(s, ve, vt)), if_(less(vd, ve), sort_bool::false_(), if_(equal_to(vd, ve), less(vs, vt), less_equal(cons_(s, vd, vs), vt)))));
        result.push_back(data_equation(variable_list({vd}), insert(s, vd, empty(s)), cons_(s, vd, empty(s))));
        result.push_back(data_equation(variable_list({vd, vs}), insert(s, vd, cons_(s, vd, vs)), cons_(s, vd, vs)));
        result.push_back(data_equation(variable_list({vd, ve, vs}), less(vd, ve), insert(s, vd, cons_(s, ve, vs)), cons_(s, vd, cons_(s, ve, vs))));
        result.push_back(data_equation(variable_list({vd, ve, vs}), less(ve, vd), insert(s, vd, cons_(s, ve, vs)), cons_(s, ve, insert(s, vd, vs))));
        result.push_back(data_equation(variable_list({vd, vs}), cinsert(s, vd, sort_bool::false_(), vs), vs));
        result.push_back(data_equation(variable_list({vd, vs}), cinsert(s, vd, sort_bool::true_(), vs), insert(s, vd, vs)));
        result.push_back(data_equation(variable_list({vd}), in(s, vd, empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, ve, vs}), in(s, vd, cons_(s, ve, vs)), sort_bool::or_(equal_to(vd, ve), in(s, vd, vs))));
        result.push_back(data_equation(variable_list({vd, ve, vs}), in(s, vd, insert(s, ve, vs)), sort_bool::or_(equal_to(vd, ve), in(s, vd, vs))));
        result.push_back(data_equation(variable_list({vs}), difference(s, vs, empty(s)), vs));
        result.push_back(data_equation(variable_list({vt}), difference(s, empty(s), vt), empty(s)));
        result.push_back(data_equation(variable_list({vd, vs, vt}), difference(s, cons_(s, vd, vs), cons_(s, vd, vt)), difference(s, vs, vt)));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less(vd, ve), difference(s, cons_(s, vd, vs), cons_(s, ve, vt)), cons_(s, vd, difference(s, vs, cons_(s, ve, vt)))));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less(ve, vd), difference(s, cons_(s, vd, vs), cons_(s, ve, vt)), difference(s, cons_(s, vd, vs), vt)));
        result.push_back(data_equation(variable_list({vs}), union_(s, vs, empty(s)), vs));
        result.push_back(data_equation(variable_list({vt}), union_(s, empty(s), vt), vt));
        result.push_back(data_equation(variable_list({vd, vs, vt}), union_(s, cons_(s, vd, vs), cons_(s, vd, vt)), cons_(s, vd, union_(s, vs, vt))));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less(vd, ve), union_(s, cons_(s, vd, vs), cons_(s, ve, vt)), cons_(s, vd, union_(s, vs, cons_(s, ve, vt)))));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less(ve, vd), union_(s, cons_(s, vd, vs), cons_(s, ve, vt)), cons_(s, ve, union_(s, cons_(s, vd, vs), vt))));
        result.push_back(data_equation(variable_list({vs}), intersection(s, vs, empty(s)), empty(s)));
        result.push_back(data_equation(variable_list({vt}), intersection(s, empty(s), vt), empty(s)));
        result.push_back(data_equation(variable_list({vd, vs, vt}), intersection(s, cons_(s, vd, vs), cons_(s, vd, vt)), cons_(s, vd, intersection(s, vs, vt))));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less(vd, ve), intersection(s, cons_(s, vd, vs), cons_(s, ve, vt)), intersection(s, vs, cons_(s, ve, vt))));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less(ve, vd), intersection(s, cons_(s, vd, vs), cons_(s, ve, vt)), intersection(s, cons_(s, vd, vs), vt)));
        result.push_back(data_equation(variable_list(), count(s, empty(s)), sort_nat::c0()));
        result.push_back(data_equation(variable_list({vd, vs}), count(s, cons_(s, vd, vs)), sort_nat::cnat(sort_nat::succ(count(s, vs)))));
        result.push_back(data_equation(variable_list({vs, vt}), not_equal_to(vs, vt), sort_bool::not_(equal_to(vs, vt))));
        result.push_back(data_equation(variable_list({vd, vs}), pick(s, cons_(s, vd, vs)), vd));
        return result;
      }

    } // namespace sort_fset

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_FSET1_H
