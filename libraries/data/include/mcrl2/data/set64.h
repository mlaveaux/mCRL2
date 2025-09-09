// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/set.h
/// \brief The standard sort set_.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/set.spec.

#ifndef MCRL2_DATA_SET64_H
#define MCRL2_DATA_SET64_H

#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/forall.h"
#include "mcrl2/data/container_sort.h"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/fset64.h"

#include <functional>
#include <map>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort set_.
    namespace sort_set {

      /// \brief Constructor for sort expression Set(S)
      /// \param s A sort expression
      /// \return Sort expression set_(s)
      inline
      container_sort set_(const sort_expression& s)
      {
        container_sort set_(set_container(), s);
        return set_;
      }

      /// \brief Recogniser for sort expression Set(s)
      /// \param e A sort expression
      /// \return true iff e is a container sort of which the name matches
      ///      set_
      inline
      bool is_set(const sort_expression& e)
      {
        if (is_container_sort(e))
        {
          return container_sort(e).container_name() == set_container();
        }
        return false;
      }


      /// \brief Generate identifier \@set.
      /// \return Identifier \@set.
      inline
      const core::identifier_string& constructor_name()
      {
        static core::identifier_string constructor_name = core::identifier_string("@set");
        return constructor_name;
      }

      /// \brief Constructor for function symbol \@set.
      /// \param s A sort expression.
      /// \return Function symbol constructor.
      inline
      function_symbol constructor(const sort_expression& s)
      {
        function_symbol constructor(constructor_name(), make_function_sort_(make_function_sort_(s, sort_bool::bool_()), sort_fset::fset(s), set_(s)));
        return constructor;
      }

      /// \brief Recogniser for function \@set.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@set.
      inline
      bool is_constructor_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == constructor_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@set.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@set to a number of arguments.
      inline
      application constructor(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_set::constructor(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@set.
      /// \param result The data expression where the \@set expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_constructor(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_set::constructor(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@set.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol constructor to a
      ///     number of arguments.
      inline
      bool is_constructor_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_constructor_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for set_.
      /// \param s A sort expression.
      /// \return All system defined constructors for set_.
      inline
      function_symbol_vector set_generate_constructors_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_set::constructor(s));

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for set_.
      /// \param s A sort expression.
      /// \return All system defined constructors that can be used in an mCRL2 specification for set_.
      inline
      function_symbol_vector set_mCRL2_usable_constructors(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_set::constructor(s));

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for set_.
      /// \param s A sort expression.
      /// \return All system defined constructors that are to be implemented in C++ for set_.
      inline
      implementation_map set_cpp_implementable_constructors(const sort_expression& s)
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier \@setfset.
      /// \return Identifier \@setfset.
      inline
      const core::identifier_string& set_fset_name()
      {
        static core::identifier_string set_fset_name = core::identifier_string("@setfset");
        return set_fset_name;
      }

      /// \brief Constructor for function symbol \@setfset.
      /// \param s A sort expression.
      /// \return Function symbol set_fset.
      inline
      function_symbol set_fset(const sort_expression& s)
      {
        function_symbol set_fset(set_fset_name(), make_function_sort_(sort_fset::fset(s), set_(s)));
        return set_fset;
      }

      /// \brief Recogniser for function \@setfset.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@setfset.
      inline
      bool is_set_fset_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == set_fset_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@setfset.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@setfset to a number of arguments.
      inline
      application set_fset(const sort_expression& s, const data_expression& arg0)
      {
        return sort_set::set_fset(s)(arg0);
      }

      /// \brief Make an application of function symbol \@setfset.
      /// \param result The data expression where the \@setfset expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_set_fset(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_set::set_fset(s),arg0);
      }

      /// \brief Recogniser for application of \@setfset.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol set_fset to a
      ///     number of arguments.
      inline
      bool is_set_fset_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_set_fset_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@setcomp.
      /// \return Identifier \@setcomp.
      inline
      const core::identifier_string& set_comprehension_name()
      {
        static core::identifier_string set_comprehension_name = core::identifier_string("@setcomp");
        return set_comprehension_name;
      }

      /// \brief Constructor for function symbol \@setcomp.
      /// \param s A sort expression.
      /// \return Function symbol set_comprehension.
      inline
      function_symbol set_comprehension(const sort_expression& s)
      {
        function_symbol set_comprehension(set_comprehension_name(), make_function_sort_(make_function_sort_(s, sort_bool::bool_()), set_(s)));
        return set_comprehension;
      }

      /// \brief Recogniser for function \@setcomp.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@setcomp.
      inline
      bool is_set_comprehension_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == set_comprehension_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@setcomp.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@setcomp to a number of arguments.
      inline
      application set_comprehension(const sort_expression& s, const data_expression& arg0)
      {
        return sort_set::set_comprehension(s)(arg0);
      }

      /// \brief Make an application of function symbol \@setcomp.
      /// \param result The data expression where the \@setcomp expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_set_comprehension(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_set::set_comprehension(s),arg0);
      }

      /// \brief Recogniser for application of \@setcomp.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol set_comprehension to a
      ///     number of arguments.
      inline
      bool is_set_comprehension_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_set_comprehension_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier in.
      /// \return Identifier in.
      inline
      const core::identifier_string& in_name()
      {
        static core::identifier_string in_name = core::identifier_string("in");
        return in_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol in(const sort_expression& , const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort(sort_bool::bool_());
        function_symbol in(in_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == in_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2;
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
        return sort_set::in(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol in.
      /// \param result The data expression where the in expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_in(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_set::in(s, arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier !.
      /// \return Identifier !.
      inline
      const core::identifier_string& complement_name()
      {
        static core::identifier_string complement_name = core::identifier_string("!");
        return complement_name;
      }

      /// \brief Constructor for function symbol !.
      /// \param s A sort expression.
      /// \return Function symbol complement.
      inline
      function_symbol complement(const sort_expression& s)
      {
        function_symbol complement(complement_name(), make_function_sort_(set_(s), set_(s)));
        return complement;
      }

      /// \brief Recogniser for function !.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching !.
      inline
      bool is_complement_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == complement_name();
        }
        return false;
      }

      /// \brief Application of function symbol !.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of ! to a number of arguments.
      inline
      application complement(const sort_expression& s, const data_expression& arg0)
      {
        return sort_set::complement(s)(arg0);
      }

      /// \brief Make an application of function symbol !.
      /// \param result The data expression where the ! expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_complement(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_set::complement(s),arg0);
      }

      /// \brief Recogniser for application of !.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol complement to a
      ///     number of arguments.
      inline
      bool is_complement_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_complement_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier +.
      /// \return Identifier +.
      inline
      const core::identifier_string& union_name()
      {
        static core::identifier_string union_name = core::identifier_string("+");
        return union_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol union_(const sort_expression& s, const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == set_(s) && s1 == set_(s))
        {
          target_sort = set_(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == sort_fset::fset(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for union_ with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol union_(union_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == union_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2;
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
        return sort_set::union_(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol +.
      /// \param result The data expression where the + expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_union_(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_set::union_(s, arg0.sort(), arg1.sort()),arg0, arg1);
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

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol intersection(const sort_expression& s, const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == set_(s) && s1 == set_(s))
        {
          target_sort = set_(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == set_(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == set_(s) && s1 == sort_fset::fset(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == sort_fset::fset(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for intersection with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol intersection(intersection_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == intersection_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2;
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
        return sort_set::intersection(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol *.
      /// \param result The data expression where the * expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_intersection(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_set::intersection(s, arg0.sort(), arg1.sort()),arg0, arg1);
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

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol difference(const sort_expression& s, const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == set_(s) && s1 == set_(s))
        {
          target_sort = set_(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == set_(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == sort_fset::fset(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else
        {
          throw mcrl2::runtime_error("Cannot compute target sort for difference with domain sorts " + pp(s0) + ", " + pp(s1) + ". ");
        }

        function_symbol difference(difference_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == difference_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2;
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
        return sort_set::difference(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol -.
      /// \param result The data expression where the - expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_difference(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_set::difference(s, arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier \@false_.
      /// \return Identifier \@false_.
      inline
      const core::identifier_string& false_function_name()
      {
        static core::identifier_string false_function_name = core::identifier_string("@false_");
        return false_function_name;
      }

      /// \brief Constructor for function symbol \@false_.
      /// \param s A sort expression.
      /// \return Function symbol false_function.
      inline
      function_symbol false_function(const sort_expression& s)
      {
        function_symbol false_function(false_function_name(), make_function_sort_(s, sort_bool::bool_()));
        return false_function;
      }

      /// \brief Recogniser for function \@false_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@false_.
      inline
      bool is_false_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == false_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@false_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@false_ to a number of arguments.
      inline
      application false_function(const sort_expression& s, const data_expression& arg0)
      {
        return sort_set::false_function(s)(arg0);
      }

      /// \brief Make an application of function symbol \@false_.
      /// \param result The data expression where the \@false_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_false_function(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_set::false_function(s),arg0);
      }

      /// \brief Recogniser for application of \@false_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol false_function to a
      ///     number of arguments.
      inline
      bool is_false_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_false_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@true_.
      /// \return Identifier \@true_.
      inline
      const core::identifier_string& true_function_name()
      {
        static core::identifier_string true_function_name = core::identifier_string("@true_");
        return true_function_name;
      }

      /// \brief Constructor for function symbol \@true_.
      /// \param s A sort expression.
      /// \return Function symbol true_function.
      inline
      function_symbol true_function(const sort_expression& s)
      {
        function_symbol true_function(true_function_name(), make_function_sort_(s, sort_bool::bool_()));
        return true_function;
      }

      /// \brief Recogniser for function \@true_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@true_.
      inline
      bool is_true_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == true_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@true_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@true_ to a number of arguments.
      inline
      application true_function(const sort_expression& s, const data_expression& arg0)
      {
        return sort_set::true_function(s)(arg0);
      }

      /// \brief Make an application of function symbol \@true_.
      /// \param result The data expression where the \@true_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_true_function(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_set::true_function(s),arg0);
      }

      /// \brief Recogniser for application of \@true_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol true_function to a
      ///     number of arguments.
      inline
      bool is_true_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_true_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@not_.
      /// \return Identifier \@not_.
      inline
      const core::identifier_string& not_function_name()
      {
        static core::identifier_string not_function_name = core::identifier_string("@not_");
        return not_function_name;
      }

      /// \brief Constructor for function symbol \@not_.
      /// \param s A sort expression.
      /// \return Function symbol not_function.
      inline
      function_symbol not_function(const sort_expression& s)
      {
        function_symbol not_function(not_function_name(), make_function_sort_(make_function_sort_(s, sort_bool::bool_()), make_function_sort_(s, sort_bool::bool_())));
        return not_function;
      }

      /// \brief Recogniser for function \@not_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@not_.
      inline
      bool is_not_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == not_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@not_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@not_ to a number of arguments.
      inline
      application not_function(const sort_expression& s, const data_expression& arg0)
      {
        return sort_set::not_function(s)(arg0);
      }

      /// \brief Make an application of function symbol \@not_.
      /// \param result The data expression where the \@not_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_not_function(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_set::not_function(s),arg0);
      }

      /// \brief Recogniser for application of \@not_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol not_function to a
      ///     number of arguments.
      inline
      bool is_not_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_not_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@and_.
      /// \return Identifier \@and_.
      inline
      const core::identifier_string& and_function_name()
      {
        static core::identifier_string and_function_name = core::identifier_string("@and_");
        return and_function_name;
      }

      /// \brief Constructor for function symbol \@and_.
      /// \param s A sort expression.
      /// \return Function symbol and_function.
      inline
      function_symbol and_function(const sort_expression& s)
      {
        function_symbol and_function(and_function_name(), make_function_sort_(make_function_sort_(s, sort_bool::bool_()), make_function_sort_(s, sort_bool::bool_()), make_function_sort_(s, sort_bool::bool_())));
        return and_function;
      }

      /// \brief Recogniser for function \@and_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@and_.
      inline
      bool is_and_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == and_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@and_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@and_ to a number of arguments.
      inline
      application and_function(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_set::and_function(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@and_.
      /// \param result The data expression where the \@and_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_and_function(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_set::and_function(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@and_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol and_function to a
      ///     number of arguments.
      inline
      bool is_and_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_and_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@or_.
      /// \return Identifier \@or_.
      inline
      const core::identifier_string& or_function_name()
      {
        static core::identifier_string or_function_name = core::identifier_string("@or_");
        return or_function_name;
      }

      /// \brief Constructor for function symbol \@or_.
      /// \param s A sort expression.
      /// \return Function symbol or_function.
      inline
      function_symbol or_function(const sort_expression& s)
      {
        function_symbol or_function(or_function_name(), make_function_sort_(make_function_sort_(s, sort_bool::bool_()), make_function_sort_(s, sort_bool::bool_()), make_function_sort_(s, sort_bool::bool_())));
        return or_function;
      }

      /// \brief Recogniser for function \@or_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@or_.
      inline
      bool is_or_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == or_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@or_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@or_ to a number of arguments.
      inline
      application or_function(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_set::or_function(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@or_.
      /// \param result The data expression where the \@or_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_or_function(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_set::or_function(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@or_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol or_function to a
      ///     number of arguments.
      inline
      bool is_or_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_or_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fset_union.
      /// \return Identifier \@fset_union.
      inline
      const core::identifier_string& fset_union_name()
      {
        static core::identifier_string fset_union_name = core::identifier_string("@fset_union");
        return fset_union_name;
      }

      /// \brief Constructor for function symbol \@fset_union.
      /// \param s A sort expression.
      /// \return Function symbol fset_union.
      inline
      function_symbol fset_union(const sort_expression& s)
      {
        function_symbol fset_union(fset_union_name(), make_function_sort_(make_function_sort_(s, sort_bool::bool_()), make_function_sort_(s, sort_bool::bool_()), sort_fset::fset(s), sort_fset::fset(s), sort_fset::fset(s)));
        return fset_union;
      }

      /// \brief Recogniser for function \@fset_union.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fset_union.
      inline
      bool is_fset_union_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == fset_union_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fset_union.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@fset_union to a number of arguments.
      inline
      application fset_union(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_set::fset_union(s)(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@fset_union.
      /// \param result The data expression where the \@fset_union expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_fset_union(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_set::fset_union(s),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@fset_union.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol fset_union to a
      ///     number of arguments.
      inline
      bool is_fset_union_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_fset_union_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fset_inter.
      /// \return Identifier \@fset_inter.
      inline
      const core::identifier_string& fset_intersection_name()
      {
        static core::identifier_string fset_intersection_name = core::identifier_string("@fset_inter");
        return fset_intersection_name;
      }

      /// \brief Constructor for function symbol \@fset_inter.
      /// \param s A sort expression.
      /// \return Function symbol fset_intersection.
      inline
      function_symbol fset_intersection(const sort_expression& s)
      {
        function_symbol fset_intersection(fset_intersection_name(), make_function_sort_(make_function_sort_(s, sort_bool::bool_()), make_function_sort_(s, sort_bool::bool_()), sort_fset::fset(s), sort_fset::fset(s), sort_fset::fset(s)));
        return fset_intersection;
      }

      /// \brief Recogniser for function \@fset_inter.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fset_inter.
      inline
      bool is_fset_intersection_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == fset_intersection_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fset_inter.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@fset_inter to a number of arguments.
      inline
      application fset_intersection(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_set::fset_intersection(s)(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@fset_inter.
      /// \param result The data expression where the \@fset_inter expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_fset_intersection(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_set::fset_intersection(s),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@fset_inter.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol fset_intersection to a
      ///     number of arguments.
      inline
      bool is_fset_intersection_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_fset_intersection_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for set_
      /// \param s A sort expression
      /// \return All system defined mappings for set_
      inline
      function_symbol_vector set_generate_functions_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_set::set_fset(s));
        result.push_back(sort_set::set_comprehension(s));
        result.push_back(sort_set::in(s, s, set_(s)));
        result.push_back(sort_set::complement(s));
        result.push_back(sort_set::union_(s, set_(s), set_(s)));
        result.push_back(sort_set::intersection(s, set_(s), set_(s)));
        result.push_back(sort_set::intersection(s, sort_fset::fset(s), set_(s)));
        result.push_back(sort_set::intersection(s, set_(s), sort_fset::fset(s)));
        result.push_back(sort_set::difference(s, set_(s), set_(s)));
        result.push_back(sort_set::difference(s, sort_fset::fset(s), set_(s)));
        result.push_back(sort_set::false_function(s));
        result.push_back(sort_set::true_function(s));
        result.push_back(sort_set::not_function(s));
        result.push_back(sort_set::and_function(s));
        result.push_back(sort_set::or_function(s));
        result.push_back(sort_set::fset_union(s));
        result.push_back(sort_set::fset_intersection(s));
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for set_
      /// \param s A sort expression
      /// \return All system defined mappings for set_
      inline
      function_symbol_vector set_generate_constructors_and_functions_code(const sort_expression& s)
      {
        function_symbol_vector result=set_generate_functions_code(s);
        for(const function_symbol& f: set_generate_constructors_code(s))
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for set_
      /// \param s A sort expression
      /// \return All system defined mappings for that can be used in mCRL2 specificationis set_
      inline
      function_symbol_vector set_mCRL2_usable_mappings(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_set::set_fset(s));
        result.push_back(sort_set::set_comprehension(s));
        result.push_back(sort_set::in(s, s, set_(s)));
        result.push_back(sort_set::complement(s));
        result.push_back(sort_set::union_(s, set_(s), set_(s)));
        result.push_back(sort_set::intersection(s, set_(s), set_(s)));
        result.push_back(sort_set::intersection(s, sort_fset::fset(s), set_(s)));
        result.push_back(sort_set::intersection(s, set_(s), sort_fset::fset(s)));
        result.push_back(sort_set::difference(s, set_(s), set_(s)));
        result.push_back(sort_set::difference(s, sort_fset::fset(s), set_(s)));
        result.push_back(sort_set::false_function(s));
        result.push_back(sort_set::true_function(s));
        result.push_back(sort_set::not_function(s));
        result.push_back(sort_set::and_function(s));
        result.push_back(sort_set::or_function(s));
        result.push_back(sort_set::fset_union(s));
        result.push_back(sort_set::fset_intersection(s));
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for set_
      /// \param s A sort expression
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for set_
      inline
      implementation_map set_cpp_implementable_mappings(const sort_expression& s)
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
        assert(is_constructor_application(e) || is_in_application(e) || is_union_application(e) || is_intersection_application(e) || is_difference_application(e) || is_and_function_application(e) || is_or_function_application(e));
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
        assert(is_constructor_application(e) || is_in_application(e) || is_union_application(e) || is_intersection_application(e) || is_difference_application(e) || is_and_function_application(e) || is_or_function_application(e));
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
        assert(is_set_fset_application(e) || is_set_comprehension_application(e) || is_complement_application(e) || is_false_function_application(e) || is_true_function_application(e) || is_not_function_application(e));
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
        assert(is_fset_union_application(e) || is_fset_intersection_application(e));
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
        assert(is_fset_union_application(e) || is_fset_intersection_application(e));
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
        assert(is_fset_union_application(e) || is_fset_intersection_application(e));
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
        assert(is_fset_union_application(e) || is_fset_intersection_application(e));
        return atermpp::down_cast<application>(e)[3];
      }

      /// \brief Give all system defined equations for set_
      /// \param s A sort expression
      /// \return All system defined equations for sort set_
      inline
      data_equation_vector set_generate_equations_code(const sort_expression& s)
      {
        variable vd("d",s);
        variable ve("e",s);
        variable vs("s",sort_fset::fset(s));
        variable vt("t",sort_fset::fset(s));
        variable vf("f",make_function_sort_(s, sort_bool::bool_()));
        variable vg("g",make_function_sort_(s, sort_bool::bool_()));
        variable vx("x",set_(s));
        variable vy("y",set_(s));
        variable vc("c",s);

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vs}), set_fset(s, vs), constructor(s, false_function(s), vs)));
        result.push_back(data_equation(variable_list({vf}), sort_set::set_comprehension(s, vf), constructor(s, vf, sort_fset::empty(s))));
        result.push_back(data_equation(variable_list({ve, vf, vs}), in(s, ve, constructor(s, vf, vs)), not_equal_to(vf(ve), in(s, ve, vs))));
        result.push_back(data_equation(variable_list({vf, vg, vs, vt}), equal_to(constructor(s, vf, vs), constructor(s, vg, vt)), forall(variable_list({vc}), equal_to(equal_to(vf(vc), vg(vc)), equal_to(in(s, vc, vs), in(s, vc, vt))))));
        result.push_back(data_equation(variable_list({vx, vy}), less(vx, vy), sort_bool::and_(less_equal(vx, vy), not_equal_to(vx, vy))));
        result.push_back(data_equation(variable_list({vx, vy}), less_equal(vx, vy), equal_to(intersection(s, vx, vy), vx)));
        result.push_back(data_equation(variable_list({vf, vs}), complement(s, constructor(s, vf, vs)), constructor(s, not_function(s, vf), vs)));
        result.push_back(data_equation(variable_list({vx}), union_(s, vx, vx), vx));
        result.push_back(data_equation(variable_list({vx, vy}), union_(s, vx, union_(s, vx, vy)), union_(s, vx, vy)));
        result.push_back(data_equation(variable_list({vx, vy}), union_(s, vx, union_(s, vy, vx)), union_(s, vy, vx)));
        result.push_back(data_equation(variable_list({vx, vy}), union_(s, union_(s, vx, vy), vx), union_(s, vx, vy)));
        result.push_back(data_equation(variable_list({vx, vy}), union_(s, union_(s, vy, vx), vx), union_(s, vy, vx)));
        result.push_back(data_equation(variable_list({vf, vg, vs, vt}), union_(s, constructor(s, vf, vs), constructor(s, vg, vt)), constructor(s, or_function(s, vf, vg), fset_union(s, vf, vg, vs, vt))));
        result.push_back(data_equation(variable_list({vx}), intersection(s, vx, vx), vx));
        result.push_back(data_equation(variable_list({vx, vy}), intersection(s, vx, intersection(s, vx, vy)), intersection(s, vx, vy)));
        result.push_back(data_equation(variable_list({vx, vy}), intersection(s, vx, intersection(s, vy, vx)), intersection(s, vy, vx)));
        result.push_back(data_equation(variable_list({vx, vy}), intersection(s, intersection(s, vx, vy), vx), intersection(s, vx, vy)));
        result.push_back(data_equation(variable_list({vx, vy}), intersection(s, intersection(s, vy, vx), vx), intersection(s, vy, vx)));
        result.push_back(data_equation(variable_list({vf, vg, vs, vt}), intersection(s, constructor(s, vf, vs), constructor(s, vg, vt)), constructor(s, and_function(s, vf, vg), fset_intersection(s, vf, vg, vs, vt))));
        result.push_back(data_equation(variable_list({vx}), intersection(s, sort_fset::empty(s), vx), sort_fset::empty(s)));
        result.push_back(data_equation(variable_list({vd, vs, vx}), intersection(s, sort_fset::cons_(s, vd, vs), vx), if_(in(s, vd, vx), sort_fset::cons_(s, vd, intersection(s, vs, vx)), intersection(s, vs, vx))));
        result.push_back(data_equation(variable_list({vs, vx}), intersection(s, vx, vs), intersection(s, vs, vx)));
        result.push_back(data_equation(variable_list({vx, vy}), difference(s, vx, vy), intersection(s, vx, complement(s, vy))));
        result.push_back(data_equation(variable_list({vs, vx}), difference(s, vs, vx), intersection(s, vs, complement(s, vx))));
        result.push_back(data_equation(variable_list({ve}), false_function(s, ve), sort_bool::false_()));
        result.push_back(data_equation(variable_list({ve}), true_function(s, ve), sort_bool::true_()));
        result.push_back(data_equation(variable_list(), equal_to(false_function(s), true_function(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list(), equal_to(true_function(s), false_function(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({ve, vf}), not_function(s, vf)(ve), sort_bool::not_(vf(ve))));
        result.push_back(data_equation(variable_list(), not_function(s, false_function(s)), true_function(s)));
        result.push_back(data_equation(variable_list(), not_function(s, true_function(s)), false_function(s)));
        result.push_back(data_equation(variable_list({ve, vf, vg}), and_function(s, vf, vg)(ve), sort_bool::and_(vf(ve), vg(ve))));
        result.push_back(data_equation(variable_list({vf}), and_function(s, vf, vf), vf));
        result.push_back(data_equation(variable_list({vf}), and_function(s, vf, false_function(s)), false_function(s)));
        result.push_back(data_equation(variable_list({vf}), and_function(s, false_function(s), vf), false_function(s)));
        result.push_back(data_equation(variable_list({vf}), and_function(s, vf, true_function(s)), vf));
        result.push_back(data_equation(variable_list({vf}), and_function(s, true_function(s), vf), vf));
        result.push_back(data_equation(variable_list({vf}), or_function(s, vf, vf), vf));
        result.push_back(data_equation(variable_list({vf}), or_function(s, vf, false_function(s)), vf));
        result.push_back(data_equation(variable_list({vf}), or_function(s, false_function(s), vf), vf));
        result.push_back(data_equation(variable_list({vf}), or_function(s, vf, true_function(s)), true_function(s)));
        result.push_back(data_equation(variable_list({vf}), or_function(s, true_function(s), vf), true_function(s)));
        result.push_back(data_equation(variable_list({ve, vf, vg}), or_function(s, vf, vg)(ve), sort_bool::or_(vf(ve), vg(ve))));
        result.push_back(data_equation(variable_list({vs, vt}), fset_union(s, false_function(s), false_function(s), vs, vt), union_(s, vs, vt)));
        result.push_back(data_equation(variable_list({vf, vg}), fset_union(s, vf, vg, sort_fset::empty(s), sort_fset::empty(s)), sort_fset::empty(s)));
        result.push_back(data_equation(variable_list({vd, vf, vg, vs}), fset_union(s, vf, vg, sort_fset::cons_(s, vd, vs), sort_fset::empty(s)), sort_fset::cinsert(s, vd, sort_bool::not_(vg(vd)), fset_union(s, vf, vg, vs, sort_fset::empty(s)))));
        result.push_back(data_equation(variable_list({ve, vf, vg, vt}), fset_union(s, vf, vg, sort_fset::empty(s), sort_fset::cons_(s, ve, vt)), sort_fset::cinsert(s, ve, sort_bool::not_(vf(ve)), fset_union(s, vf, vg, sort_fset::empty(s), vt))));
        result.push_back(data_equation(variable_list({vd, vf, vg, vs, vt}), fset_union(s, vf, vg, sort_fset::cons_(s, vd, vs), sort_fset::cons_(s, vd, vt)), sort_fset::cinsert(s, vd, equal_to(vf(vd), vg(vd)), fset_union(s, vf, vg, vs, vt))));
        result.push_back(data_equation(variable_list({vd, ve, vf, vg, vs, vt}), less(vd, ve), fset_union(s, vf, vg, sort_fset::cons_(s, vd, vs), sort_fset::cons_(s, ve, vt)), sort_fset::cinsert(s, vd, sort_bool::not_(vg(vd)), fset_union(s, vf, vg, vs, sort_fset::cons_(s, ve, vt)))));
        result.push_back(data_equation(variable_list({vd, ve, vf, vg, vs, vt}), less(ve, vd), fset_union(s, vf, vg, sort_fset::cons_(s, vd, vs), sort_fset::cons_(s, ve, vt)), sort_fset::cinsert(s, ve, sort_bool::not_(vf(ve)), fset_union(s, vf, vg, sort_fset::cons_(s, vd, vs), vt))));
        result.push_back(data_equation(variable_list({vf, vg}), fset_intersection(s, vf, vg, sort_fset::empty(s), sort_fset::empty(s)), sort_fset::empty(s)));
        result.push_back(data_equation(variable_list({vd, vf, vg, vs}), fset_intersection(s, vf, vg, sort_fset::cons_(s, vd, vs), sort_fset::empty(s)), sort_fset::cinsert(s, vd, vg(vd), fset_intersection(s, vf, vg, vs, sort_fset::empty(s)))));
        result.push_back(data_equation(variable_list({ve, vf, vg, vt}), fset_intersection(s, vf, vg, sort_fset::empty(s), sort_fset::cons_(s, ve, vt)), sort_fset::cinsert(s, ve, vf(ve), fset_intersection(s, vf, vg, sort_fset::empty(s), vt))));
        result.push_back(data_equation(variable_list({vd, vf, vg, vs, vt}), fset_intersection(s, vf, vg, sort_fset::cons_(s, vd, vs), sort_fset::cons_(s, vd, vt)), sort_fset::cinsert(s, vd, equal_to(vf(vd), vg(vd)), fset_intersection(s, vf, vg, vs, vt))));
        result.push_back(data_equation(variable_list({vd, ve, vf, vg, vs, vt}), less(vd, ve), fset_intersection(s, vf, vg, sort_fset::cons_(s, vd, vs), sort_fset::cons_(s, ve, vt)), sort_fset::cinsert(s, vd, vg(vd), fset_intersection(s, vf, vg, vs, sort_fset::cons_(s, ve, vt)))));
        result.push_back(data_equation(variable_list({vd, ve, vf, vg, vs, vt}), less(ve, vd), fset_intersection(s, vf, vg, sort_fset::cons_(s, vd, vs), sort_fset::cons_(s, ve, vt)), sort_fset::cinsert(s, ve, vf(ve), fset_intersection(s, vf, vg, sort_fset::cons_(s, vd, vs), vt))));
        return result;
      }

    } // namespace sort_set_

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_SET64_H
