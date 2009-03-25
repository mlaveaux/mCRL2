// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/new_data/expression_traits.h
/// \brief Contains term traits for data_expression.

#ifndef MCRL2_DATA_TERM_TRAITS_H
#define MCRL2_DATA_TERM_TRAITS_H

#include <functional>
#include "mcrl2/core/term_traits.h"
#include "mcrl2/core/print.h"
#include "mcrl2/new_data/data_expression.h"
#include "mcrl2/new_data/variable.h"
#include "mcrl2/new_data/find.h"
#include "mcrl2/new_data/detail/data_sequence_algorithm.h"

namespace mcrl2 {

namespace core {

  /// \brief Contains type information for data expressions.
  template <>
  struct term_traits<new_data::data_expression>
  {
    /// \brief The term type
    typedef new_data::data_expression term_type;

    /// \brief The variable type
    typedef new_data::variable variable_type;

    /// \brief The variable sequence type
    typedef new_data::variable_list variable_sequence_type;

    /// \brief The value true
    /// \return The value true
    static inline
    term_type true_()
    { return atermpp::aterm_appl(core::detail::gsMakeDataExprTrue()); }

    /// \brief The value false
    /// \return The value false
    static inline
    term_type false_()
    { return atermpp::aterm_appl(core::detail::gsMakeDataExprFalse()); }

    /// \brief Operator not
    /// \param p A term
    /// \return Operator not applied to p
    static inline
    term_type not_(term_type p)
    { return atermpp::aterm_appl(core::detail::gsMakeDataExprNot(p)); }

    /// \brief Operator and
    /// \param p A term
    /// \param q A term
    /// \return Operator and applied to p and q
    static inline
    term_type and_(term_type p, term_type q)
    { return atermpp::aterm_appl(core::detail::gsMakeDataExprAnd(p,q)); }

    /// \brief Operator or
    /// \param p A term
    /// \param q A term
    /// \return Operator or applied to p and q
    static inline
    term_type or_(term_type p, term_type q)
    { return atermpp::aterm_appl(core::detail::gsMakeDataExprOr(p,q)); }

    /// \brief Test for value true
    /// \param t A term
    /// \return True if the term has the value true
    static inline
    bool is_true(term_type t)
    { return core::detail::gsIsDataExprTrue(t); }

    /// \brief Test for value false
    /// \param t A term
    /// \return True if the term has the value false
    static inline
    bool is_false(term_type t)
    { return core::detail::gsIsDataExprFalse(t); }

    /// \brief Test for operator not
    /// \param t A term
    /// \return True if the term is of type not
    static inline
    bool is_not(term_type t)
    { return core::detail::gsIsDataExprNot(t); }

    /// \brief Test for operator and
    /// \param t A term
    /// \return True if the term is of type and
    static inline
    bool is_and(term_type t)
    { return core::detail::gsIsDataExprAnd(t); }

    /// \brief Test for operator or
    /// \param t A term
    /// \return True if the term is of type or
    static inline
    bool is_or(term_type t)
    { return core::detail::gsIsDataExprOr(t); }

    /// \brief Test for implication
    /// \param t A term
    /// \return True if the term is an implication
    static inline
    bool is_imp(term_type t)
    { return core::detail::gsIsDataExprImp(t);; }

    /// \brief Test for universal quantification
    /// \param t A term
    /// \return True if the term is an universal quantification
    static inline
    bool is_forall(term_type t)
    { return core::detail::gsIsDataExprForall(t); }

    /// \brief Test for existential quantification
    /// \param t A term
    /// \return True if the term is an existential quantification
    static inline
    bool is_exists(term_type t)
    { return core::detail::gsIsDataExprExists(t); }

    /// \brief Test for lambda abstraction
    /// \param t A term
    /// \return True if the term is a lambda expression
    static inline
    bool is_lambda(term_type t)
    { return core::detail::gsIsDataExprLambda(t); }

    /// \brief Conversion from variable to term
    /// \param v A variable
    /// \return The converted variable
    static inline
    term_type variable2term(variable_type v)
    {
      return v;
    }

    /// \brief Test if a term is a variable
    /// \param t A term
    /// \return True if the term is a variable
    static inline
    bool is_variable(term_type t)
    {
      return t.is_variable();
    }

    /// \brief Returns the free variables of a term
    /// \param t A term
    /// \return The free variables of a term
    static inline
    variable_sequence_type free_variables(term_type t)
    {
      std::set<variable_type> v = new_data::find_all_variables(t);
      return variable_sequence_type(v.begin(), v.end());
    }

    /// \brief Returns the difference of two unordered sets of variables
    /// \param v A sequence of data variables
    /// \param w A sequence of data variables
    /// \return The difference of two sets.
    static inline
    variable_sequence_type set_difference(const variable_sequence_type& v, const variable_sequence_type& w)
    {
      return new_data::detail::set_difference(v, w);
    }

    /// \brief Test if a term is constant
    /// \param t A term
    /// \return True if the term is constant. N.B. It is unknown if the current implementation
    /// works for quantifier expressions.
    static inline
    bool is_constant(term_type t)
    {
      struct local {
        static bool caster(atermpp::aterm p) {
          return new_data::data_expression(p).is_variable();
        }
      };
      return atermpp::find_if(t, std::ptr_fun(&local::caster)) == atermpp::aterm();
    }

    /// \brief Pretty print function
    /// \param t A term
    /// \return A pretty print representation of the term
    static inline
    std::string pp(term_type t)
    {
      return core::pp(t);
    }
  };

} // namespace core

namespace new_data {
  /// \brief TODO replace term_traits by expression_traits
  template < typename Expression >
  struct expression_traits : public core::term_traits< Expression >
  { };
}

} // namespace mcrl2

#endif // MCRL2_DATA_TERM_TRAITS_H
