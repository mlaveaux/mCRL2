// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/bag.h
/// \brief The standard sort bag.
///
/// This file was generated from the data sort specification
/// mcrl2/data/build/bag.spec.

#ifndef MCRL2_DATA_BAG64_H
#define MCRL2_DATA_BAG64_H

#include "mcrl2/data/basic_sort.h"
#include "mcrl2/data/function_sort.h"
#include "mcrl2/data/function_symbol.h"
#include "mcrl2/data/application.h"
#include "mcrl2/data/data_equation.h"
#include "mcrl2/data/standard.h"
#include "mcrl2/data/forall.h"
#include "mcrl2/data/container_sort.h"
#include "mcrl2/data/nat64.h"
#include "mcrl2/data/fbag64.h"
#include "mcrl2/data/fset64.h"
#include "mcrl2/data/set64.h"

#include <functional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
 import utilities;
#endif

namespace mcrl2 {

  namespace data {

    /// \brief Namespace for system defined sort bag.
    namespace sort_bag {

      /// \brief Constructor for sort expression Bag(S)
      /// \param s A sort expression
      /// \return Sort expression bag(s)
      inline
      container_sort bag(const sort_expression& s)
      {
        container_sort bag(bag_container(), s);
        return bag;
      }

      /// \brief Recogniser for sort expression Bag(s)
      /// \param e A sort expression
      /// \return true iff e is a container sort of which the name matches
      ///      bag
      inline
      bool is_bag(const sort_expression& e)
      {
        if (is_container_sort(e))
        {
          return container_sort(e).container_name() == bag_container();
        }
        return false;
      }


      /// \brief Generate identifier \@bag.
      /// \return Identifier \@bag.
      inline
      const core::identifier_string& constructor_name()
      {
        static core::identifier_string constructor_name = core::identifier_string("@bag");
        return constructor_name;
      }

      /// \brief Constructor for function symbol \@bag.
      /// \param s A sort expression.
      /// \return Function symbol constructor.
      inline
      function_symbol constructor(const sort_expression& s)
      {
        function_symbol constructor(constructor_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), sort_fbag::fbag(s), bag(s)));
        return constructor;
      }

      /// \brief Recogniser for function \@bag.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@bag.
      inline
      bool is_constructor_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == constructor_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@bag.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@bag to a number of arguments.
      inline
      application constructor(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_bag::constructor(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@bag.
      /// \param result The data expression where the \@bag expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_constructor(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::constructor(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@bag.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol constructor to a
      ///     number of arguments.
      inline
      bool is_constructor_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_constructor_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined constructors for bag.
      /// \param s A sort expression.
      /// \return All system defined constructors for bag.
      inline
      function_symbol_vector bag_generate_constructors_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_bag::constructor(s));

        return result;
      }
      /// \brief Give all defined constructors which can be used in mCRL2 specs for bag.
      /// \param s A sort expression.
      /// \return All system defined constructors that can be used in an mCRL2 specification for bag.
      inline
      function_symbol_vector bag_mCRL2_usable_constructors(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_bag::constructor(s));

        return result;
      }
      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined constructors which have an implementation in C++ and not in rewrite rules for bag.
      /// \param s A sort expression.
      /// \return All system defined constructors that are to be implemented in C++ for bag.
      inline
      implementation_map bag_cpp_implementable_constructors(const sort_expression& s)
      {
        implementation_map result;
        return result;
      }

      /// \brief Generate identifier \@bagfbag.
      /// \return Identifier \@bagfbag.
      inline
      const core::identifier_string& bag_fbag_name()
      {
        static core::identifier_string bag_fbag_name = core::identifier_string("@bagfbag");
        return bag_fbag_name;
      }

      /// \brief Constructor for function symbol \@bagfbag.
      /// \param s A sort expression.
      /// \return Function symbol bag_fbag.
      inline
      function_symbol bag_fbag(const sort_expression& s)
      {
        function_symbol bag_fbag(bag_fbag_name(), make_function_sort_(sort_fbag::fbag(s), bag(s)));
        return bag_fbag;
      }

      /// \brief Recogniser for function \@bagfbag.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@bagfbag.
      inline
      bool is_bag_fbag_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == bag_fbag_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@bagfbag.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@bagfbag to a number of arguments.
      inline
      application bag_fbag(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::bag_fbag(s)(arg0);
      }

      /// \brief Make an application of function symbol \@bagfbag.
      /// \param result The data expression where the \@bagfbag expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_bag_fbag(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::bag_fbag(s),arg0);
      }

      /// \brief Recogniser for application of \@bagfbag.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol bag_fbag to a
      ///     number of arguments.
      inline
      bool is_bag_fbag_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_bag_fbag_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@bagcomp.
      /// \return Identifier \@bagcomp.
      inline
      const core::identifier_string& bag_comprehension_name()
      {
        static core::identifier_string bag_comprehension_name = core::identifier_string("@bagcomp");
        return bag_comprehension_name;
      }

      /// \brief Constructor for function symbol \@bagcomp.
      /// \param s A sort expression.
      /// \return Function symbol bag_comprehension.
      inline
      function_symbol bag_comprehension(const sort_expression& s)
      {
        function_symbol bag_comprehension(bag_comprehension_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), bag(s)));
        return bag_comprehension;
      }

      /// \brief Recogniser for function \@bagcomp.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@bagcomp.
      inline
      bool is_bag_comprehension_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == bag_comprehension_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@bagcomp.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@bagcomp to a number of arguments.
      inline
      application bag_comprehension(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::bag_comprehension(s)(arg0);
      }

      /// \brief Make an application of function symbol \@bagcomp.
      /// \param result The data expression where the \@bagcomp expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_bag_comprehension(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::bag_comprehension(s),arg0);
      }

      /// \brief Recogniser for application of \@bagcomp.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol bag_comprehension to a
      ///     number of arguments.
      inline
      bool is_bag_comprehension_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_bag_comprehension_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier count.
      /// \return Identifier count.
      inline
      const core::identifier_string& count_name()
      {
        static core::identifier_string count_name = core::identifier_string("count");
        return count_name;
      }

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol count(const sort_expression& , const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort(sort_nat::nat());
        function_symbol count(count_name(), make_function_sort_(s0, s1, target_sort));
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
          const function_symbol& f = atermpp::down_cast<function_symbol>(e);
          return f.name() == count_name() && atermpp::down_cast<function_sort>(f.sort()).domain().size() == 2;
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
        return sort_bag::count(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol count.
      /// \param result The data expression where the count expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_count(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::count(s, arg0.sort(), arg1.sort()),arg0, arg1);
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
        return sort_bag::in(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol in.
      /// \param result The data expression where the in expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_in(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::in(s, arg0.sort(), arg1.sort()),arg0, arg1);
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

      // This function is not intended for public use and therefore not documented in Doxygen.
      inline
      function_symbol union_(const sort_expression& s, const sort_expression& s0, const sort_expression& s1)
      {
        sort_expression target_sort;
        if (s0 == bag(s) && s1 == bag(s))
        {
          target_sort = bag(s);
        }
        else if (s0 == sort_set::set_(s) && s1 == sort_set::set_(s))
        {
          target_sort = sort_set::set_(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == sort_fset::fset(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == sort_fbag::fbag(s) && s1 == sort_fbag::fbag(s))
        {
          target_sort = sort_fbag::fbag(s);
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
        return sort_bag::union_(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol +.
      /// \param result The data expression where the + expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_union_(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::union_(s, arg0.sort(), arg1.sort()),arg0, arg1);
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
        if (s0 == bag(s) && s1 == bag(s))
        {
          target_sort = bag(s);
        }
        else if (s0 == sort_fbag::fbag(s) && s1 == bag(s))
        {
          target_sort = sort_fbag::fbag(s);
        }
        else if (s0 == bag(s) && s1 == sort_fbag::fbag(s))
        {
          target_sort = sort_fbag::fbag(s);
        }
        else if (s0 == sort_set::set_(s) && s1 == sort_set::set_(s))
        {
          target_sort = sort_set::set_(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == sort_set::set_(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == sort_set::set_(s) && s1 == sort_fset::fset(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == sort_fset::fset(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == sort_fbag::fbag(s) && s1 == sort_fbag::fbag(s))
        {
          target_sort = sort_fbag::fbag(s);
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
        return sort_bag::intersection(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol *.
      /// \param result The data expression where the * expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_intersection(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::intersection(s, arg0.sort(), arg1.sort()),arg0, arg1);
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
        if (s0 == bag(s) && s1 == bag(s))
        {
          target_sort = bag(s);
        }
        else if (s0 == sort_fbag::fbag(s) && s1 == bag(s))
        {
          target_sort = sort_fbag::fbag(s);
        }
        else if (s0 == sort_set::set_(s) && s1 == sort_set::set_(s))
        {
          target_sort = sort_set::set_(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == sort_set::set_(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == sort_fset::fset(s) && s1 == sort_fset::fset(s))
        {
          target_sort = sort_fset::fset(s);
        }
        else if (s0 == sort_fbag::fbag(s) && s1 == sort_fbag::fbag(s))
        {
          target_sort = sort_fbag::fbag(s);
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
        return sort_bag::difference(s, arg0.sort(), arg1.sort())(arg0, arg1);
      }

      /// \brief Make an application of function symbol -.
      /// \param result The data expression where the - expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_difference(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::difference(s, arg0.sort(), arg1.sort()),arg0, arg1);
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

      /// \brief Generate identifier Bag2Set.
      /// \return Identifier Bag2Set.
      inline
      const core::identifier_string& bag2set_name()
      {
        static core::identifier_string bag2set_name = core::identifier_string("Bag2Set");
        return bag2set_name;
      }

      /// \brief Constructor for function symbol Bag2Set.
      /// \param s A sort expression.
      /// \return Function symbol bag2set.
      inline
      function_symbol bag2set(const sort_expression& s)
      {
        function_symbol bag2set(bag2set_name(), make_function_sort_(bag(s), sort_set::set_(s)));
        return bag2set;
      }

      /// \brief Recogniser for function Bag2Set.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Bag2Set.
      inline
      bool is_bag2set_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == bag2set_name();
        }
        return false;
      }

      /// \brief Application of function symbol Bag2Set.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of Bag2Set to a number of arguments.
      inline
      application bag2set(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::bag2set(s)(arg0);
      }

      /// \brief Make an application of function symbol Bag2Set.
      /// \param result The data expression where the Bag2Set expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_bag2set(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::bag2set(s),arg0);
      }

      /// \brief Recogniser for application of Bag2Set.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol bag2set to a
      ///     number of arguments.
      inline
      bool is_bag2set_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_bag2set_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier Set2Bag.
      /// \return Identifier Set2Bag.
      inline
      const core::identifier_string& set2bag_name()
      {
        static core::identifier_string set2bag_name = core::identifier_string("Set2Bag");
        return set2bag_name;
      }

      /// \brief Constructor for function symbol Set2Bag.
      /// \param s A sort expression.
      /// \return Function symbol set2bag.
      inline
      function_symbol set2bag(const sort_expression& s)
      {
        function_symbol set2bag(set2bag_name(), make_function_sort_(sort_set::set_(s), bag(s)));
        return set2bag;
      }

      /// \brief Recogniser for function Set2Bag.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching Set2Bag.
      inline
      bool is_set2bag_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == set2bag_name();
        }
        return false;
      }

      /// \brief Application of function symbol Set2Bag.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of Set2Bag to a number of arguments.
      inline
      application set2bag(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::set2bag(s)(arg0);
      }

      /// \brief Make an application of function symbol Set2Bag.
      /// \param result The data expression where the Set2Bag expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_set2bag(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::set2bag(s),arg0);
      }

      /// \brief Recogniser for application of Set2Bag.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol set2bag to a
      ///     number of arguments.
      inline
      bool is_set2bag_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_set2bag_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@zero_.
      /// \return Identifier \@zero_.
      inline
      const core::identifier_string& zero_function_name()
      {
        static core::identifier_string zero_function_name = core::identifier_string("@zero_");
        return zero_function_name;
      }

      /// \brief Constructor for function symbol \@zero_.
      /// \param s A sort expression.
      /// \return Function symbol zero_function.
      inline
      function_symbol zero_function(const sort_expression& s)
      {
        function_symbol zero_function(zero_function_name(), make_function_sort_(s, sort_nat::nat()));
        return zero_function;
      }

      /// \brief Recogniser for function \@zero_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@zero_.
      inline
      bool is_zero_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == zero_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@zero_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@zero_ to a number of arguments.
      inline
      application zero_function(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::zero_function(s)(arg0);
      }

      /// \brief Make an application of function symbol \@zero_.
      /// \param result The data expression where the \@zero_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_zero_function(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::zero_function(s),arg0);
      }

      /// \brief Recogniser for application of \@zero_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol zero_function to a
      ///     number of arguments.
      inline
      bool is_zero_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_zero_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@one_.
      /// \return Identifier \@one_.
      inline
      const core::identifier_string& one_function_name()
      {
        static core::identifier_string one_function_name = core::identifier_string("@one_");
        return one_function_name;
      }

      /// \brief Constructor for function symbol \@one_.
      /// \param s A sort expression.
      /// \return Function symbol one_function.
      inline
      function_symbol one_function(const sort_expression& s)
      {
        function_symbol one_function(one_function_name(), make_function_sort_(s, sort_nat::nat()));
        return one_function;
      }

      /// \brief Recogniser for function \@one_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@one_.
      inline
      bool is_one_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == one_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@one_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@one_ to a number of arguments.
      inline
      application one_function(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::one_function(s)(arg0);
      }

      /// \brief Make an application of function symbol \@one_.
      /// \param result The data expression where the \@one_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_one_function(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::one_function(s),arg0);
      }

      /// \brief Recogniser for application of \@one_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol one_function to a
      ///     number of arguments.
      inline
      bool is_one_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_one_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@add_.
      /// \return Identifier \@add_.
      inline
      const core::identifier_string& add_function_name()
      {
        static core::identifier_string add_function_name = core::identifier_string("@add_");
        return add_function_name;
      }

      /// \brief Constructor for function symbol \@add_.
      /// \param s A sort expression.
      /// \return Function symbol add_function.
      inline
      function_symbol add_function(const sort_expression& s)
      {
        function_symbol add_function(add_function_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat())));
        return add_function;
      }

      /// \brief Recogniser for function \@add_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@add_.
      inline
      bool is_add_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == add_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@add_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@add_ to a number of arguments.
      inline
      application add_function(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_bag::add_function(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@add_.
      /// \param result The data expression where the \@add_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_add_function(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::add_function(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@add_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol add_function to a
      ///     number of arguments.
      inline
      bool is_add_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_add_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@min_.
      /// \return Identifier \@min_.
      inline
      const core::identifier_string& min_function_name()
      {
        static core::identifier_string min_function_name = core::identifier_string("@min_");
        return min_function_name;
      }

      /// \brief Constructor for function symbol \@min_.
      /// \param s A sort expression.
      /// \return Function symbol min_function.
      inline
      function_symbol min_function(const sort_expression& s)
      {
        function_symbol min_function(min_function_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat())));
        return min_function;
      }

      /// \brief Recogniser for function \@min_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@min_.
      inline
      bool is_min_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == min_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@min_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@min_ to a number of arguments.
      inline
      application min_function(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_bag::min_function(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@min_.
      /// \param result The data expression where the \@min_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_min_function(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::min_function(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@min_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol min_function to a
      ///     number of arguments.
      inline
      bool is_min_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_min_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@monus_.
      /// \return Identifier \@monus_.
      inline
      const core::identifier_string& monus_function_name()
      {
        static core::identifier_string monus_function_name = core::identifier_string("@monus_");
        return monus_function_name;
      }

      /// \brief Constructor for function symbol \@monus_.
      /// \param s A sort expression.
      /// \return Function symbol monus_function.
      inline
      function_symbol monus_function(const sort_expression& s)
      {
        function_symbol monus_function(monus_function_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat())));
        return monus_function;
      }

      /// \brief Recogniser for function \@monus_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@monus_.
      inline
      bool is_monus_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == monus_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@monus_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@monus_ to a number of arguments.
      inline
      application monus_function(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_bag::monus_function(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@monus_.
      /// \param result The data expression where the \@monus_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_monus_function(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::monus_function(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@monus_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol monus_function to a
      ///     number of arguments.
      inline
      bool is_monus_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_monus_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@Nat2Bool_.
      /// \return Identifier \@Nat2Bool_.
      inline
      const core::identifier_string& nat2bool_function_name()
      {
        static core::identifier_string nat2bool_function_name = core::identifier_string("@Nat2Bool_");
        return nat2bool_function_name;
      }

      /// \brief Constructor for function symbol \@Nat2Bool_.
      /// \param s A sort expression.
      /// \return Function symbol nat2bool_function.
      inline
      function_symbol nat2bool_function(const sort_expression& s)
      {
        function_symbol nat2bool_function(nat2bool_function_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_bool::bool_())));
        return nat2bool_function;
      }

      /// \brief Recogniser for function \@Nat2Bool_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@Nat2Bool_.
      inline
      bool is_nat2bool_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == nat2bool_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@Nat2Bool_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@Nat2Bool_ to a number of arguments.
      inline
      application nat2bool_function(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::nat2bool_function(s)(arg0);
      }

      /// \brief Make an application of function symbol \@Nat2Bool_.
      /// \param result The data expression where the \@Nat2Bool_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_nat2bool_function(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::nat2bool_function(s),arg0);
      }

      /// \brief Recogniser for application of \@Nat2Bool_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol nat2bool_function to a
      ///     number of arguments.
      inline
      bool is_nat2bool_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_nat2bool_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@Bool2Nat_.
      /// \return Identifier \@Bool2Nat_.
      inline
      const core::identifier_string& bool2nat_function_name()
      {
        static core::identifier_string bool2nat_function_name = core::identifier_string("@Bool2Nat_");
        return bool2nat_function_name;
      }

      /// \brief Constructor for function symbol \@Bool2Nat_.
      /// \param s A sort expression.
      /// \return Function symbol bool2nat_function.
      inline
      function_symbol bool2nat_function(const sort_expression& s)
      {
        function_symbol bool2nat_function(bool2nat_function_name(), make_function_sort_(make_function_sort_(s, sort_bool::bool_()), make_function_sort_(s, sort_nat::nat())));
        return bool2nat_function;
      }

      /// \brief Recogniser for function \@Bool2Nat_.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@Bool2Nat_.
      inline
      bool is_bool2nat_function_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == bool2nat_function_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@Bool2Nat_.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@Bool2Nat_ to a number of arguments.
      inline
      application bool2nat_function(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::bool2nat_function(s)(arg0);
      }

      /// \brief Make an application of function symbol \@Bool2Nat_.
      /// \param result The data expression where the \@Bool2Nat_ expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_bool2nat_function(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::bool2nat_function(s),arg0);
      }

      /// \brief Recogniser for application of \@Bool2Nat_.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol bool2nat_function to a
      ///     number of arguments.
      inline
      bool is_bool2nat_function_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_bool2nat_function_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fbag_join.
      /// \return Identifier \@fbag_join.
      inline
      const core::identifier_string& join_name()
      {
        static core::identifier_string join_name = core::identifier_string("@fbag_join");
        return join_name;
      }

      /// \brief Constructor for function symbol \@fbag_join.
      /// \param s A sort expression.
      /// \return Function symbol join.
      inline
      function_symbol join(const sort_expression& s)
      {
        function_symbol join(join_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat()), sort_fbag::fbag(s), sort_fbag::fbag(s), sort_fbag::fbag(s)));
        return join;
      }

      /// \brief Recogniser for function \@fbag_join.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fbag_join.
      inline
      bool is_join_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == join_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fbag_join.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@fbag_join to a number of arguments.
      inline
      application join(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_bag::join(s)(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@fbag_join.
      /// \param result The data expression where the \@fbag_join expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_join(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_bag::join(s),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@fbag_join.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol join to a
      ///     number of arguments.
      inline
      bool is_join_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_join_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fbag_inter.
      /// \return Identifier \@fbag_inter.
      inline
      const core::identifier_string& fbag_intersect_name()
      {
        static core::identifier_string fbag_intersect_name = core::identifier_string("@fbag_inter");
        return fbag_intersect_name;
      }

      /// \brief Constructor for function symbol \@fbag_inter.
      /// \param s A sort expression.
      /// \return Function symbol fbag_intersect.
      inline
      function_symbol fbag_intersect(const sort_expression& s)
      {
        function_symbol fbag_intersect(fbag_intersect_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat()), sort_fbag::fbag(s), sort_fbag::fbag(s), sort_fbag::fbag(s)));
        return fbag_intersect;
      }

      /// \brief Recogniser for function \@fbag_inter.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fbag_inter.
      inline
      bool is_fbag_intersect_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == fbag_intersect_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fbag_inter.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@fbag_inter to a number of arguments.
      inline
      application fbag_intersect(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_bag::fbag_intersect(s)(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@fbag_inter.
      /// \param result The data expression where the \@fbag_inter expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_fbag_intersect(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_bag::fbag_intersect(s),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@fbag_inter.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol fbag_intersect to a
      ///     number of arguments.
      inline
      bool is_fbag_intersect_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_fbag_intersect_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fbag_diff.
      /// \return Identifier \@fbag_diff.
      inline
      const core::identifier_string& fbag_difference_name()
      {
        static core::identifier_string fbag_difference_name = core::identifier_string("@fbag_diff");
        return fbag_difference_name;
      }

      /// \brief Constructor for function symbol \@fbag_diff.
      /// \param s A sort expression.
      /// \return Function symbol fbag_difference.
      inline
      function_symbol fbag_difference(const sort_expression& s)
      {
        function_symbol fbag_difference(fbag_difference_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), make_function_sort_(s, sort_nat::nat()), sort_fbag::fbag(s), sort_fbag::fbag(s), sort_fbag::fbag(s)));
        return fbag_difference;
      }

      /// \brief Recogniser for function \@fbag_diff.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fbag_diff.
      inline
      bool is_fbag_difference_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == fbag_difference_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fbag_diff.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      /// \return Application of \@fbag_diff to a number of arguments.
      inline
      application fbag_difference(const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        return sort_bag::fbag_difference(s)(arg0, arg1, arg2, arg3);
      }

      /// \brief Make an application of function symbol \@fbag_diff.
      /// \param result The data expression where the \@fbag_diff expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \param arg2 A data expression.
      /// \param arg3 A data expression.
      inline
      void make_fbag_difference(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1, const data_expression& arg2, const data_expression& arg3)
      {
        make_application(result, sort_bag::fbag_difference(s),arg0, arg1, arg2, arg3);
      }

      /// \brief Recogniser for application of \@fbag_diff.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol fbag_difference to a
      ///     number of arguments.
      inline
      bool is_fbag_difference_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_fbag_difference_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fbag2fset.
      /// \return Identifier \@fbag2fset.
      inline
      const core::identifier_string& fbag2fset_name()
      {
        static core::identifier_string fbag2fset_name = core::identifier_string("@fbag2fset");
        return fbag2fset_name;
      }

      /// \brief Constructor for function symbol \@fbag2fset.
      /// \param s A sort expression.
      /// \return Function symbol fbag2fset.
      inline
      function_symbol fbag2fset(const sort_expression& s)
      {
        function_symbol fbag2fset(fbag2fset_name(), make_function_sort_(make_function_sort_(s, sort_nat::nat()), sort_fbag::fbag(s), sort_fset::fset(s)));
        return fbag2fset;
      }

      /// \brief Recogniser for function \@fbag2fset.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fbag2fset.
      inline
      bool is_fbag2fset_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == fbag2fset_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fbag2fset.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      /// \return Application of \@fbag2fset to a number of arguments.
      inline
      application fbag2fset(const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        return sort_bag::fbag2fset(s)(arg0, arg1);
      }

      /// \brief Make an application of function symbol \@fbag2fset.
      /// \param result The data expression where the \@fbag2fset expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \param arg1 A data expression.
      inline
      void make_fbag2fset(data_expression& result, const sort_expression& s, const data_expression& arg0, const data_expression& arg1)
      {
        make_application(result, sort_bag::fbag2fset(s),arg0, arg1);
      }

      /// \brief Recogniser for application of \@fbag2fset.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol fbag2fset to a
      ///     number of arguments.
      inline
      bool is_fbag2fset_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_fbag2fset_function_symbol(atermpp::down_cast<application>(e).head());
      }

      /// \brief Generate identifier \@fset2fbag.
      /// \return Identifier \@fset2fbag.
      inline
      const core::identifier_string& fset2fbag_name()
      {
        static core::identifier_string fset2fbag_name = core::identifier_string("@fset2fbag");
        return fset2fbag_name;
      }

      /// \brief Constructor for function symbol \@fset2fbag.
      /// \param s A sort expression.
      /// \return Function symbol fset2fbag.
      inline
      function_symbol fset2fbag(const sort_expression& s)
      {
        function_symbol fset2fbag(fset2fbag_name(), make_function_sort_(sort_fset::fset(s), sort_fbag::fbag(s)));
        return fset2fbag;
      }

      /// \brief Recogniser for function \@fset2fbag.
      /// \param e A data expression.
      /// \return true iff e is the function symbol matching \@fset2fbag.
      inline
      bool is_fset2fbag_function_symbol(const atermpp::aterm& e)
      {
        if (is_function_symbol(e))
        {
          return atermpp::down_cast<function_symbol>(e).name() == fset2fbag_name();
        }
        return false;
      }

      /// \brief Application of function symbol \@fset2fbag.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      /// \return Application of \@fset2fbag to a number of arguments.
      inline
      application fset2fbag(const sort_expression& s, const data_expression& arg0)
      {
        return sort_bag::fset2fbag(s)(arg0);
      }

      /// \brief Make an application of function symbol \@fset2fbag.
      /// \param result The data expression where the \@fset2fbag expression is put.
      /// \param s A sort expression.
      /// \param arg0 A data expression.
      inline
      void make_fset2fbag(data_expression& result, const sort_expression& s, const data_expression& arg0)
      {
        make_application(result, sort_bag::fset2fbag(s),arg0);
      }

      /// \brief Recogniser for application of \@fset2fbag.
      /// \param e A data expression.
      /// \return true iff e is an application of function symbol fset2fbag to a
      ///     number of arguments.
      inline
      bool is_fset2fbag_application(const atermpp::aterm& e)
      {
        return is_application(e) && is_fset2fbag_function_symbol(atermpp::down_cast<application>(e).head());
      }
      /// \brief Give all system defined mappings for bag
      /// \param s A sort expression
      /// \return All system defined mappings for bag
      inline
      function_symbol_vector bag_generate_functions_code(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_bag::bag_fbag(s));
        result.push_back(sort_bag::bag_comprehension(s));
        result.push_back(sort_bag::count(s, s, bag(s)));
        result.push_back(sort_bag::in(s, s, bag(s)));
        result.push_back(sort_bag::union_(s, bag(s), bag(s)));
        result.push_back(sort_bag::intersection(s, bag(s), bag(s)));
        result.push_back(sort_bag::intersection(s, sort_fbag::fbag(s), bag(s)));
        result.push_back(sort_bag::intersection(s, bag(s), sort_fbag::fbag(s)));
        result.push_back(sort_bag::difference(s, bag(s), bag(s)));
        result.push_back(sort_bag::difference(s, sort_fbag::fbag(s), bag(s)));
        result.push_back(sort_bag::bag2set(s));
        result.push_back(sort_bag::set2bag(s));
        result.push_back(sort_bag::zero_function(s));
        result.push_back(sort_bag::one_function(s));
        result.push_back(sort_bag::add_function(s));
        result.push_back(sort_bag::min_function(s));
        result.push_back(sort_bag::monus_function(s));
        result.push_back(sort_bag::nat2bool_function(s));
        result.push_back(sort_bag::bool2nat_function(s));
        result.push_back(sort_bag::join(s));
        result.push_back(sort_bag::fbag_intersect(s));
        result.push_back(sort_bag::fbag_difference(s));
        result.push_back(sort_bag::fbag2fset(s));
        result.push_back(sort_bag::fset2fbag(s));
        return result;
      }
      
      /// \brief Give all system defined mappings and constructors for bag
      /// \param s A sort expression
      /// \return All system defined mappings for bag
      inline
      function_symbol_vector bag_generate_constructors_and_functions_code(const sort_expression& s)
      {
        function_symbol_vector result=bag_generate_functions_code(s);
        for(const function_symbol& f: bag_generate_constructors_code(s))
        {
          result.push_back(f);
        }
        return result;
      }
      
      /// \brief Give all system defined mappings that can be used in mCRL2 specs for bag
      /// \param s A sort expression
      /// \return All system defined mappings for that can be used in mCRL2 specificationis bag
      inline
      function_symbol_vector bag_mCRL2_usable_mappings(const sort_expression& s)
      {
        function_symbol_vector result;
        result.push_back(sort_bag::bag_fbag(s));
        result.push_back(sort_bag::bag_comprehension(s));
        result.push_back(sort_bag::count(s, s, bag(s)));
        result.push_back(sort_bag::in(s, s, bag(s)));
        result.push_back(sort_bag::union_(s, bag(s), bag(s)));
        result.push_back(sort_bag::intersection(s, bag(s), bag(s)));
        result.push_back(sort_bag::intersection(s, sort_fbag::fbag(s), bag(s)));
        result.push_back(sort_bag::intersection(s, bag(s), sort_fbag::fbag(s)));
        result.push_back(sort_bag::difference(s, bag(s), bag(s)));
        result.push_back(sort_bag::difference(s, sort_fbag::fbag(s), bag(s)));
        result.push_back(sort_bag::bag2set(s));
        result.push_back(sort_bag::set2bag(s));
        result.push_back(sort_bag::zero_function(s));
        result.push_back(sort_bag::one_function(s));
        result.push_back(sort_bag::add_function(s));
        result.push_back(sort_bag::min_function(s));
        result.push_back(sort_bag::monus_function(s));
        result.push_back(sort_bag::nat2bool_function(s));
        result.push_back(sort_bag::bool2nat_function(s));
        result.push_back(sort_bag::join(s));
        result.push_back(sort_bag::fbag_intersect(s));
        result.push_back(sort_bag::fbag_difference(s));
        result.push_back(sort_bag::fbag2fset(s));
        result.push_back(sort_bag::fset2fbag(s));
        return result;
      }


      // The typedef is the sort that maps a function symbol to an function that rewrites it as well as a string of a function that can be used to implement it
      using implementation_map = std::map<function_symbol,std::pair<std::function<void(data_expression&, const data_expression&)>, std::string> >;
      /// \brief Give all system defined mappings that are to be implemented in C++ code for bag
      /// \param s A sort expression
      /// \return A mapping from C++ implementable function symbols to system defined mappings implemented in C++ code for bag
      inline
      implementation_map bag_cpp_implementable_mappings(const sort_expression& s)
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
        assert(is_constructor_application(e) || is_count_application(e) || is_in_application(e) || is_union_application(e) || is_intersection_application(e) || is_difference_application(e) || is_add_function_application(e) || is_min_function_application(e) || is_monus_function_application(e) || is_fbag2fset_application(e));
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
        assert(is_constructor_application(e) || is_count_application(e) || is_in_application(e) || is_union_application(e) || is_intersection_application(e) || is_difference_application(e) || is_add_function_application(e) || is_min_function_application(e) || is_monus_function_application(e) || is_fbag2fset_application(e));
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
        assert(is_bag_fbag_application(e) || is_bag_comprehension_application(e) || is_bag2set_application(e) || is_set2bag_application(e) || is_zero_function_application(e) || is_one_function_application(e) || is_nat2bool_function_application(e) || is_bool2nat_function_application(e) || is_fset2fbag_application(e));
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
        assert(is_join_application(e) || is_fbag_intersect_application(e) || is_fbag_difference_application(e));
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
        assert(is_join_application(e) || is_fbag_intersect_application(e) || is_fbag_difference_application(e));
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
        assert(is_join_application(e) || is_fbag_intersect_application(e) || is_fbag_difference_application(e));
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
        assert(is_join_application(e) || is_fbag_intersect_application(e) || is_fbag_difference_application(e));
        return atermpp::down_cast<application>(e)[3];
      }

      /// \brief Give all system defined equations for bag
      /// \param s A sort expression
      /// \return All system defined equations for sort bag
      inline
      data_equation_vector bag_generate_equations_code(const sort_expression& s)
      {
        variable vb("b",sort_fbag::fbag(s));
        variable vc("c",sort_fbag::fbag(s));
        variable vd("d",s);
        variable ve("e",s);
        variable vf("f",make_function_sort_(s, sort_nat::nat()));
        variable vg("g",make_function_sort_(s, sort_nat::nat()));
        variable vh("h",make_function_sort_(s, sort_bool::bool_()));
        variable vp("p",sort_pos::pos());
        variable vq("q",sort_pos::pos());
        variable vs("s",sort_fset::fset(s));
        variable vx("x",bag(s));
        variable vy("y",bag(s));
        variable vz("z",sort_fbag::fbag(s));
        variable vw("w",sort_fbag::fbag(s));

        data_equation_vector result;
        result.push_back(data_equation(variable_list({vb}), bag_fbag(s, vb), constructor(s, zero_function(s), vb)));
        result.push_back(data_equation(variable_list({vf}), sort_bag::bag_comprehension(s, vf), constructor(s, vf, sort_fbag::empty(s))));
        result.push_back(data_equation(variable_list({vb, ve, vf}), count(s, ve, constructor(s, vf, vb)), sort_nat::swap_zero(vf(ve), count(s, ve, vb))));
        result.push_back(data_equation(variable_list({ve, vx}), in(s, ve, vx), greater(count(s, ve, vx), sort_nat::c0())));
        result.push_back(data_equation(variable_list({vb, vc, vf, vg}), equal_to(constructor(s, vf, vb), constructor(s, vg, vc)), if_(equal_to(vf, vg), equal_to(vb, vc), forall(variable_list({vd}), equal_to(count(s, vd, constructor(s, vf, vb)), count(s, vd, constructor(s, vg, vc)))))));
        result.push_back(data_equation(variable_list({vx, vy}), less(vx, vy), sort_bool::and_(less_equal(vx, vy), not_equal_to(vx, vy))));
        result.push_back(data_equation(variable_list({vx, vy}), less_equal(vx, vy), equal_to(intersection(s, vx, vy), vx)));
        result.push_back(data_equation(variable_list({vb, vc, vf, vg}), union_(s, constructor(s, vf, vb), constructor(s, vg, vc)), constructor(s, add_function(s, vf, vg), join(s, vf, vg, vb, vc))));
        result.push_back(data_equation(variable_list({vx}), intersection(s, vx, vx), vx));
        result.push_back(data_equation(variable_list({vx, vy}), intersection(s, vx, intersection(s, vx, vy)), intersection(s, vx, vy)));
        result.push_back(data_equation(variable_list({vx, vy}), intersection(s, vx, intersection(s, vy, vx)), intersection(s, vy, vx)));
        result.push_back(data_equation(variable_list({vx, vy}), intersection(s, intersection(s, vx, vy), vx), intersection(s, vx, vy)));
        result.push_back(data_equation(variable_list({vx, vy}), intersection(s, intersection(s, vy, vx), vx), intersection(s, vy, vx)));
        result.push_back(data_equation(variable_list({vb, vc, vf, vg}), intersection(s, constructor(s, vf, vb), constructor(s, vg, vc)), constructor(s, min_function(s, vf, vg), fbag_intersect(s, vf, vg, vb, vc))));
        result.push_back(data_equation(variable_list({vx}), intersection(s, sort_fbag::empty(s), vx), sort_fbag::empty(s)));
        result.push_back(data_equation(variable_list({vb, vd, vp, vx}), intersection(s, sort_fbag::cons_(s, vd, vp, vb), vx), if_(in(s, vd, vx), sort_fbag::cons_(s, vd, sort_nat::minimum(vp, sort_nat::nat2pos(count(s, vd, vx))), intersection(s, vb, vx)), intersection(s, vb, vx))));
        result.push_back(data_equation(variable_list({vb, vx}), intersection(s, vx, vb), intersection(s, vb, vx)));
        result.push_back(data_equation(variable_list({vb, vc, vf, vg}), difference(s, constructor(s, vf, vb), constructor(s, vg, vc)), constructor(s, monus_function(s, vf, vg), fbag_difference(s, vf, vg, vb, vc))));
        result.push_back(data_equation(variable_list({vx}), difference(s, sort_fbag::empty(s), vx), sort_fbag::empty(s)));
        result.push_back(data_equation(variable_list({vb, vd, vp, vx}), difference(s, sort_fbag::cons_(s, vd, vp, vb), vx), if_(greater(sort_nat::pos2nat(vp), count(s, vd, vx)), sort_fbag::cons_(s, vd, sort_nat::nat2pos(sort_nat::monus(sort_nat::pos2nat(vp), count(s, vd, vx))), difference(s, vb, vx)), difference(s, vb, vx))));
        result.push_back(data_equation(variable_list({vb, vf}), bag2set(s, constructor(s, vf, vb)), sort_set::constructor(s, nat2bool_function(s, vf), fbag2fset(s, vf, vb))));
        result.push_back(data_equation(variable_list({vh, vs}), set2bag(s, sort_set::constructor(s, vh, vs)), constructor(s, bool2nat_function(s, vh), fset2fbag(s, vs))));
        result.push_back(data_equation(variable_list({ve}), zero_function(s, ve), sort_nat::c0()));
        result.push_back(data_equation(variable_list({ve}), one_function(s, ve), sort_nat::most_significant_digit_nat(sort_machine_word::one_word())));
        result.push_back(data_equation(variable_list(), equal_to(zero_function(s), one_function(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list(), equal_to(one_function(s), zero_function(s)), sort_bool::false_()));
        result.push_back(data_equation(variable_list({ve, vf, vg}), add_function(s, vf, vg)(ve), sort_nat::auxiliary_plus_nat(vf(ve), vg(ve))));
        result.push_back(data_equation(variable_list({vf}), add_function(s, vf, zero_function(s)), vf));
        result.push_back(data_equation(variable_list({vf}), add_function(s, zero_function(s), vf), vf));
        result.push_back(data_equation(variable_list({ve, vf, vg}), min_function(s, vf, vg)(ve), sort_nat::minimum(vf(ve), vg(ve))));
        result.push_back(data_equation(variable_list({vf}), min_function(s, vf, vf), vf));
        result.push_back(data_equation(variable_list({vf}), min_function(s, vf, zero_function(s)), zero_function(s)));
        result.push_back(data_equation(variable_list({vf}), min_function(s, zero_function(s), vf), zero_function(s)));
        result.push_back(data_equation(variable_list({ve, vf, vg}), monus_function(s, vf, vg)(ve), sort_nat::monus(vf(ve), vg(ve))));
        result.push_back(data_equation(variable_list({vf}), monus_function(s, vf, vf), zero_function(s)));
        result.push_back(data_equation(variable_list({vf}), monus_function(s, vf, zero_function(s)), vf));
        result.push_back(data_equation(variable_list({vf}), monus_function(s, zero_function(s), vf), zero_function(s)));
        result.push_back(data_equation(variable_list({ve, vf}), nat2bool_function(s, vf)(ve), greater(vf(ve), sort_nat::c0())));
        result.push_back(data_equation(variable_list(), nat2bool_function(s, zero_function(s)), sort_set::false_function(s)));
        result.push_back(data_equation(variable_list(), nat2bool_function(s, one_function(s)), sort_set::true_function(s)));
        result.push_back(data_equation(variable_list({ve, vh}), bool2nat_function(s, vh)(ve), if_(vh(ve), sort_nat::most_significant_digit_nat(sort_machine_word::one_word()), sort_nat::c0())));
        result.push_back(data_equation(variable_list(), bool2nat_function(s, sort_set::false_function(s)), zero_function(s)));
        result.push_back(data_equation(variable_list(), bool2nat_function(s, sort_set::true_function(s)), one_function(s)));
        result.push_back(data_equation(variable_list({vw, vz}), join(s, zero_function(s), zero_function(s), vw, vz), union_(s, vw, vz)));
        result.push_back(data_equation(variable_list({vf, vg}), join(s, vf, vg, sort_fbag::empty(s), sort_fbag::empty(s)), sort_fbag::empty(s)));
        result.push_back(data_equation(variable_list({vb, vd, vf, vg, vp}), join(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::empty(s)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_add(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::c0()), join(s, vf, vg, vb, sort_fbag::empty(s)))));
        result.push_back(data_equation(variable_list({vc, ve, vf, vg, vq}), join(s, vf, vg, sort_fbag::empty(s), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, ve, sort_nat::swap_zero_add(vf(ve), vg(ve), sort_nat::c0(), sort_nat::pos2nat(vq)), join(s, vf, vg, sort_fbag::empty(s), vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, vf, vg, vp, vq}), join(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, vd, vq, vc)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_add(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::pos2nat(vq)), join(s, vf, vg, vb, vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vf, vg, vp, vq}), less(vd, ve), join(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_add(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::c0()), join(s, vf, vg, vb, sort_fbag::cons_(s, ve, vq, vc)))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vf, vg, vp, vq}), less(ve, vd), join(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, ve, sort_nat::swap_zero_add(vf(ve), vg(ve), sort_nat::c0(), sort_nat::pos2nat(vq)), join(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), vc))));
        result.push_back(data_equation(variable_list({vf, vg}), fbag_intersect(s, vf, vg, sort_fbag::empty(s), sort_fbag::empty(s)), sort_fbag::empty(s)));
        result.push_back(data_equation(variable_list({vb, vd, vf, vg, vp}), fbag_intersect(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::empty(s)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_min(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::c0()), fbag_intersect(s, vf, vg, vb, sort_fbag::empty(s)))));
        result.push_back(data_equation(variable_list({vc, ve, vf, vg, vq}), fbag_intersect(s, vf, vg, sort_fbag::empty(s), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, ve, sort_nat::swap_zero_min(vf(ve), vg(ve), sort_nat::c0(), sort_nat::pos2nat(vq)), fbag_intersect(s, vf, vg, sort_fbag::empty(s), vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, vf, vg, vp, vq}), fbag_intersect(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, vd, vq, vc)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_min(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::pos2nat(vq)), fbag_intersect(s, vf, vg, vb, vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vf, vg, vp, vq}), less(vd, ve), fbag_intersect(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_min(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::c0()), fbag_intersect(s, vf, vg, vb, sort_fbag::cons_(s, ve, vq, vc)))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vf, vg, vp, vq}), less(ve, vd), fbag_intersect(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, ve, sort_nat::swap_zero_min(vf(ve), vg(ve), sort_nat::c0(), sort_nat::pos2nat(vq)), fbag_intersect(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), vc))));
        result.push_back(data_equation(variable_list({vf, vg}), fbag_difference(s, vf, vg, sort_fbag::empty(s), sort_fbag::empty(s)), sort_fbag::empty(s)));
        result.push_back(data_equation(variable_list({vb, vd, vf, vg, vp}), fbag_difference(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::empty(s)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_monus(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::c0()), fbag_difference(s, vf, vg, vb, sort_fbag::empty(s)))));
        result.push_back(data_equation(variable_list({vc, ve, vf, vg, vq}), fbag_difference(s, vf, vg, sort_fbag::empty(s), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, ve, sort_nat::swap_zero_monus(vf(ve), vg(ve), sort_nat::c0(), sort_nat::pos2nat(vq)), fbag_difference(s, vf, vg, sort_fbag::empty(s), vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, vf, vg, vp, vq}), fbag_difference(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, vd, vq, vc)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_monus(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::pos2nat(vq)), fbag_difference(s, vf, vg, vb, vc))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vf, vg, vp, vq}), less(vd, ve), fbag_difference(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, vd, sort_nat::swap_zero_monus(vf(vd), vg(vd), sort_nat::pos2nat(vp), sort_nat::c0()), fbag_difference(s, vf, vg, vb, sort_fbag::cons_(s, ve, vq, vc)))));
        result.push_back(data_equation(variable_list({vb, vc, vd, ve, vf, vg, vp, vq}), less(ve, vd), fbag_difference(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), sort_fbag::cons_(s, ve, vq, vc)), sort_fbag::cinsert(s, ve, sort_nat::swap_zero_monus(vf(ve), vg(ve), sort_nat::c0(), sort_nat::pos2nat(vq)), fbag_difference(s, vf, vg, sort_fbag::cons_(s, vd, vp, vb), vc))));
        result.push_back(data_equation(variable_list({vf}), fbag2fset(s, vf, sort_fbag::empty(s)), sort_fset::empty(s)));
        result.push_back(data_equation(variable_list({vb, vd, vf, vp}), fbag2fset(s, vf, sort_fbag::cons_(s, vd, vp, vb)), sort_fset::cinsert(s, vd, equal_to(equal_to(vf(vd), sort_nat::pos2nat(vp)), greater(vf(vd), sort_nat::c0())), fbag2fset(s, vf, vb))));
        result.push_back(data_equation(variable_list(), fset2fbag(s, sort_fset::empty(s)), sort_fbag::empty(s)));
        result.push_back(data_equation(variable_list({vd, vs}), fset2fbag(s, sort_fset::cons_(s, vd, vs)), sort_fbag::cinsert(s, vd, sort_nat::pos2nat(sort_pos::c1()), fset2fbag(s, vs))));
        return result;
      }

    } // namespace sort_bag

  } // namespace data

} // namespace mcrl2

#endif // MCRL2_DATA_BAG64_H
