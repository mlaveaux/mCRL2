// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/list.h
/// \brief The standard sort list.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/list.spec.

#ifndef MCRL2_DATA_LIST64_H
#define MCRL2_DATA_LIST64_H

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

#include <functional>
#include <map>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort list.
    namespace sort_list {

      /// \brief Constructor for sort expression List(S)
      /// \param s A sort expression
      /// \return Sort expression list(s)
      inline
      container_sort list(const sort_expression& s)
      {
        container_sort list(list_container(), s);
        return list;
      }

      /// \brief Recogniser for sort expression List(s)
      /// \param e A sort expression
      /// \return true iff e is a container sort of which the name matches
      ///      list
      inline
      bool is_list(const sort_expression& e)
      {
        if (is_container_sort(e))
        {
          return container_sort(e).container_name() == list_container();
        }
        return false;
      }


      /// \brief Generate identifier [].
      /// \return Identifier [].
      inline
      const core::identifier_string& empty_name()
      {
        static core::identifier_string empty_name = core::identifier_string("[]");
        return empty_name;
      }

      /// \brief Constructor for function symbol [].
      /// \param s A sort expression.
      /// \return Function symbol empty.
      inline
      function_symbol empty(const sort_expression& s)
      {
        function_symbol empty(empty_name(), list(s));
        return empty;
      }

      /// \brief Recogniser for function [].
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching [].
      inline
      bool is_empty_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == empty_name();
        }
        return false;
      }

      /// \brief Generate identifier |>.
      /// \return Identifier |>.
      inline
      const core::identifier_string& cons_name()
      {
        static core::identifier_string cons_name = core::identifier_string("|>");
        return cons_name;
      }

      /// \brief Constructor for function symbol |>.
      /// \param s A sort expression.
      /// \return Function symbol cons_.
      inline
      function_symbol cons_(const sort_expression& s)
      {
        function_symbol cons_(cons_name(), make_function_sort_(s, list(s), list(s)));
        return cons_;
      }

      /// \brief Recogniser for function |>.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching |>.
      inline
      bool is_cons_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == cons_name();
        }
        return false;
      }

      /// \brief Application of function symbol |>.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of |> to a number of arguments.
      inline
      application cons_(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_list::cons_(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol |>.
      /// \param result The data expression where the |> expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_cons_(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_list::cons_(s),arg0, arg1);
      }

      /// \brief Recogniser for application of |>.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol cons_ to a
      ///     number of arguments.
      inline
      bool is_cons_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_cons_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for list.
      /// \param s A sort expression.
      /// \return All system defined constructors for list.
      inline
      function_symbol_vector list_generate_constructors_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_list::empty(s));
        result.push_back(sort_list::cons_(s));

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for list.
      /// \param s A sort expression.
      /// \return All system defined constructors that can be used in an mCRL2 specification for list.
      inline
      function_symbol_vector list_mCRL2_usable_constructors(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_list::empty(s));
        result.push_back(sort_list::cons_(s));

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for list.
      /// \param s A sort expression.
      /// \return All system defined constructors that are to be implemented in C++ for list.
      inline
      implementation_map list_cpp_implementable_constructors(const sort_expression& s)
      {
        implementation_map result;
        return result;
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
        function_symbol in(in_name(), make_function_sort_(s, list(s), sort_bool::bool_()));
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
        return sort_list::in(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol in.
      /// \param result The data expression where the in expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_in(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_list::in(s),arg0, arg1);
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
        function_symbol count(count_name(), make_function_sort_(list(s), sort_nat::nat()));
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
        return sort_list::count(s)(arg0);
      }

      /// \brief Make an application of function symbol #.
      /// \param result The data expression where the # expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_count(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_list::count(s),arg0);
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

      /// \brief Generate identifier <|.
      /// \return Identifier <|.
      inline
      const core::identifier_string& snoc_name()
      {
        static core::identifier_string snoc_name = core::identifier_string("<|");
        return snoc_name;
      }

      /// \brief Constructor for function symbol <|.
      /// \param s A sort expression.
      /// \return Function symbol snoc.
      inline
      function_symbol snoc(const sort_expression& s)
      {
        function_symbol snoc(snoc_name(), make_function_sort_(list(s), s, list(s)));
        return snoc;
      }

      /// \brief Recogniser for function <|.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching <|.
      inline
      bool is_snoc_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == snoc_name();
        }
        return false;
      }

      /// \brief Application of function symbol <|.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of <| to a number of arguments.
      inline
      application snoc(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_list::snoc(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol <|.
      /// \param result The data expression where the <| expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_snoc(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_list::snoc(s),arg0, arg1);
      }

      /// \brief Recogniser for application of <|.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol snoc to a
      ///     number of arguments.
      inline
      bool is_snoc_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_snoc_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier ++.
      /// \return Identifier ++.
      inline
      const core::identifier_string& concat_name()
      {
        static core::identifier_string concat_name = core::identifier_string("++");
        return concat_name;
      }

      /// \brief Constructor for function symbol ++.
      /// \param s A sort expression.
      /// \return Function symbol concat.
      inline
      function_symbol concat(const sort_expression& s)
      {
        function_symbol concat(concat_name(), make_function_sort_(list(s), list(s), list(s)));
        return concat;
      }

      /// \brief Recogniser for function ++.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching ++.
      inline
      bool is_concat_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == concat_name();
        }
        return false;
      }

      /// \brief Application of function symbol ++.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of ++ to a number of arguments.
      inline
      application concat(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_list::concat(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol ++.
      /// \param result The data expression where the ++ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_concat(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_list::concat(s),arg0, arg1);
      }

      /// \brief Recogniser for application of ++.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol concat to a
      ///     number of arguments.
      inline
      bool is_concat_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_concat_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier ..
      /// \return Identifier ..
      inline
      const core::identifier_string& element_at_name()
      {
        static core::identifier_string element_at_name = core::identifier_string(".");
        return element_at_name;
      }

      /// \brief Constructor for function symbol ..
      /// \param s A sort expression.
      /// \return Function symbol element_at.
      inline
      function_symbol element_at(const sort_expression& s)
      {
        function_symbol element_at(element_at_name(), make_function_sort_(list(s), sort_nat::nat(), s));
        return element_at;
      }

      /// \brief Recogniser for function ..
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching ..
      inline
      bool is_element_at_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == element_at_name();
        }
        return false;
      }

      /// \brief Application of function symbol ..
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of . to a number of arguments.
      inline
      application element_at(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_list::element_at(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol ..
      /// \param result The data expression where the . expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_element_at(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_list::element_at(s),arg0, arg1);
      }

      /// \brief Recogniser for application of ..
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol element_at to a
      ///     number of arguments.
      inline
      bool is_element_at_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_element_at_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier head.
      /// \return Identifier head.
      inline
      const core::identifier_string& head_name()
      {
        static core::identifier_string head_name = core::identifier_string("head");
        return head_name;
      }

      /// \brief Constructor for function symbol head.
      /// \param s A sort expression.
      /// \return Function symbol head.
      inline
      function_symbol head(const sort_expression& s)
      {
        function_symbol head(head_name(), make_function_sort_(list(s), s));
        return head;
      }

      /// \brief Recogniser for function head.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching head.
      inline
      bool is_head_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == head_name();
        }
        return false;
      }

      /// \brief Application of function symbol head.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of head to a number of arguments.
      inline
      application head(const sort_expression& s, const data_expression& arg0)
      {
        return sort_list::head(s)(arg0);
      }

      /// \brief Make an application of function symbol head.
      /// \param result The data expression where the head expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_head(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_list::head(s),arg0);
      }

      /// \brief Recogniser for application of head.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol head to a
      ///     number of arguments.
      inline
      bool is_head_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_head_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier tail.
      /// \return Identifier tail.
      inline
      const core::identifier_string& tail_name()
      {
        static core::identifier_string tail_name = core::identifier_string("tail");
        return tail_name;
      }

      /// \brief Constructor for function symbol tail.
      /// \param s A sort expression.
      /// \return Function symbol tail.
      inline
      function_symbol tail(const sort_expression& s)
      {
        function_symbol tail(tail_name(), make_function_sort_(list(s), list(s)));
        return tail;
      }

      /// \brief Recogniser for function tail.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching tail.
      inline
      bool is_tail_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == tail_name();
        }
        return false;
      }

      /// \brief Application of function symbol tail.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of tail to a number of arguments.
      inline
      application tail(const sort_expression& s, const data_expression& arg0)
      {
        return sort_list::tail(s)(arg0);
      }

      /// \brief Make an application of function symbol tail.
      /// \param result The data expression where the tail expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_tail(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_list::tail(s),arg0);
      }

      /// \brief Recogniser for application of tail.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol tail to a
      ///     number of arguments.
      inline
      bool is_tail_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_tail_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier rhead.
      /// \return Identifier rhead.
      inline
      const core::identifier_string& rhead_name()
      {
        static core::identifier_string rhead_name = core::identifier_string("rhead");
        return rhead_name;
      }

      /// \brief Constructor for function symbol rhead.
      /// \param s A sort expression.
      /// \return Function symbol rhead.
      inline
      function_symbol rhead(const sort_expression& s)
      {
        function_symbol rhead(rhead_name(), make_function_sort_(list(s), s));
        return rhead;
      }

      /// \brief Recogniser for function rhead.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching rhead.
      inline
      bool is_rhead_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == rhead_name();
        }
        return false;
      }

      /// \brief Application of function symbol rhead.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of rhead to a number of arguments.
      inline
      application rhead(const sort_expression& s, const data_expression& arg0)
      {
        return sort_list::rhead(s)(arg0);
      }

      /// \brief Make an application of function symbol rhead.
      /// \param result The data expression where the rhead expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_rhead(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_list::rhead(s),arg0);
      }

      /// \brief Recogniser for application of rhead.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol rhead to a
      ///     number of arguments.
      inline
      bool is_rhead_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_rhead_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier rtail.
      /// \return Identifier rtail.
      inline
      const core::identifier_string& rtail_name()
      {
        static core::identifier_string rtail_name = core::identifier_string("rtail");
        return rtail_name;
      }

      /// \brief Constructor for function symbol rtail.
      /// \param s A sort expression.
      /// \return Function symbol rtail.
      inline
      function_symbol rtail(const sort_expression& s)
      {
        function_symbol rtail(rtail_name(), make_function_sort_(list(s), list(s)));
        return rtail;
      }

      /// \brief Recogniser for function rtail.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching rtail.
      inline
      bool is_rtail_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == rtail_name();
        }
        return false;
      }

      /// \brief Application of function symbol rtail.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of rtail to a number of arguments.
      inline
      application rtail(const sort_expression& s, const data_expression& arg0)
      {
        return sort_list::rtail(s)(arg0);
      }

      /// \brief Make an application of function symbol rtail.
      /// \param result The data expression where the rtail expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_rtail(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_list::rtail(s),arg0);
      }

      /// \brief Recogniser for application of rtail.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol rtail to a
      ///     number of arguments.
      inline
      bool is_rtail_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_rtail_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for list
      /// \param s A sort expression
      /// \return All system defined mappings for list
      inline
      function_symbol_vector list_generate_functions_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_list::in(s));
        result.push_back(sort_list::count(s));
        result.push_back(sort_list::snoc(s));
        result.push_back(sort_list::concat(s));
        result.push_back(sort_list::element_at(s));
        result.push_back(sort_list::head(s));
        result.push_back(sort_list::tail(s));
        result.push_back(sort_list::rhead(s));
        result.push_back(sort_list::rtail(s));
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for list
      /// \param s A sort expression
      /// \return All system defined mappings for list
      inline
      function_symbol_vector list_generate_constructors_and_functions_code(const sort_expression& s)
      {
        function_symbol_vector result=list_generate_functions_code(s);
        for(const function_symbol& f: list_generate_constructors_code(s))
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for list
      /// \param s A sort expression
      /// \return All system defined mappings for that can be used in mCRL2 specificationis list
      inline
      function_symbol_vector list_mCRL2_usable_mappings(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_list::in(s));
        result.push_back(sort_list::count(s));
        result.push_back(sort_list::snoc(s));
        result.push_back(sort_list::concat(s));
        result.push_back(sort_list::element_at(s));
        result.push_back(sort_list::head(s));
        result.push_back(sort_list::tail(s));
        result.push_back(sort_list::rhead(s));
        result.push_back(sort_list::rtail(s));
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for list
      /// \param s A sort expression
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for list
      inline
      implementation_map list_cpp_implementable_mappings(const sort_expression& s)
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
        assert(is_cons_application(e) || is_in_application(e) || is_snoc_application(e) || is_concat_application(e) || is_element_at_application(e));
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
        assert(is_cons_application(e) || is_in_application(e) || is_snoc_application(e) || is_concat_application(e) || is_element_at_application(e));
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
        assert(is_count_application(e) || is_head_application(e) || is_tail_application(e) || is_rhead_application(e) || is_rtail_application(e));
        return atermpp::down_cast<application>(e)[0];
      }

      /// \brief Give all system defined equations for list
      /// \param s A sort expression
      /// \return All system defined equations for sort list
      inline
      data_equation_vector list_generate_equations_code(const sort_expression& s)
      {
        variable vd("d",s);
        variable ve("e",s);
        variable vs("s",list(s));
        variable vt("t",list(s));
        variable vp("p",sort_pos::pos());
        variable vn("n",sort_nat::nat());

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vd, vs}), equal_to(empty(s), cons_(s, vd, vs)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, vs}), equal_to(cons_(s, vd, vs), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), equal_to(cons_(s, vd, vs), cons_(s, ve, vt)), sort_bool::and_(equal_to(vd, ve), equal_to(vs, vt))));
        result.push_back(data_equation(variable_list({vd, vs}), less(empty(s), cons_(s, vd, vs)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vd, vs}), less(cons_(s, vd, vs), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less(cons_(s, vd, vs), cons_(s, ve, vt)), sort_bool::or_(sort_bool::and_(equal_to(vd, ve), less(vs, vt)), less(vd, ve))));
        result.push_back(data_equation(variable_list({vd, vs}), less_equal(empty(s), cons_(s, vd, vs)), sort_bool::true_()));
        result.push_back(data_equation(variable_list({vd, vs}), less_equal(cons_(s, vd, vs), empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, ve, vs, vt}), less_equal(cons_(s, vd, vs), cons_(s, ve, vt)), sort_bool::or_(sort_bool::and_(equal_to(vd, ve), less_equal(vs, vt)), less(vd, ve))));
        result.push_back(data_equation(variable_list({vd}), in(s, vd, empty(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({vd, ve, vs}), in(s, vd, cons_(s, ve, vs)), sort_bool::or_(equal_to(vd, ve), in(s, vd, vs))));
        result.push_back(data_equation(variable_list(), count(s, empty(s)), sort_nat::c0()));
        result.push_back(data_equation(variable_list({vd, vs}), count(s, cons_(s, vd, vs)), sort_nat::succ_nat(count(s, vs))));
        result.push_back(data_equation(variable_list({vd}), snoc(s, empty(s), vd), cons_(s, vd, empty(s))));
        result.push_back(data_equation(variable_list({vd, ve, vs}), snoc(s, cons_(s, vd, vs), ve), cons_(s, vd, snoc(s, vs, ve))));
        result.push_back(data_equation(variable_list({vs}), concat(s, empty(s), vs), vs));
        result.push_back(data_equation(variable_list({vd, vs, vt}), concat(s, cons_(s, vd, vs), vt), cons_(s, vd, concat(s, vs, vt))));
        result.push_back(data_equation(variable_list({vs}), concat(s, vs, empty(s)), vs));
        result.push_back(data_equation(variable_list({vd, vn, vs}), element_at(s, cons_(s, vd, vs), vn), if_(equal_to(vn, sort_nat::most_significant_digit_nat(sort_machine_word::zero_word())), vd, element_at(s, vs, sort_nat::natpred(vn)))));
        result.push_back(data_equation(variable_list({vd, vs}), head(s, cons_(s, vd, vs)), vd));
        result.push_back(data_equation(variable_list({vd, vs}), tail(s, cons_(s, vd, vs)), vs));
        result.push_back(data_equation(variable_list({vd}), rhead(s, cons_(s, vd, empty(s))), vd));
        result.push_back(data_equation(variable_list({vd, ve, vs}), rhead(s, cons_(s, vd, cons_(s, ve, vs))), rhead(s, cons_(s, ve, vs))));
        result.push_back(data_equation(variable_list({vd}), rtail(s, cons_(s, vd, empty(s))), empty(s)));
        result.push_back(data_equation(variable_list({vd, ve, vs}), rtail(s, cons_(s, vd, cons_(s, ve, vs))), cons_(s, vd, rtail(s, cons_(s, ve, vs)))));
        return result;
      }

    } // namespace sort_list

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_LIST64_H
