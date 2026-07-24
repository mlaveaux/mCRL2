// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/detail/bqnf_quantifier_rewriter.h
/// \brief Replaces universal quantifier over conjuncts with conjuncts over universal quantifiers:
/// rewrite_bqnf_expression(forall x . /\_i phi_i) =
///   /\_i forall (x intersection free(phi_i)) . phi_i.
/// This rewriter is experimental.
#ifndef MCRL2_PBES_DETAIL_BQNF_QUANTIFIER_REWRITER_H
#define MCRL2_PBES_DETAIL_BQNF_QUANTIFIER_REWRITER_H

#include "mcrl2/pbes/detail/bqnf_visitor.h"

namespace mcrl2::pbes_system::detail
{

struct bqnf_quantifier_rewriter: public bqnf_visitor
  {
    /// \brief Constructor.
    bqnf_quantifier_rewriter() = default;

    /// \brief Filters the expression e such that subexpressions that are data expression that
    /// do not refer to variables in the set d are discarded:
    /// filter(e, d) = e if e is a data expression and the intersection of free(e) and d is not empty;
    /// filter(e1 'op' e2, d) = filter(e1, d) 'op' filter(e2, d).
    /// \param e a simple PBES expression.
    /// \param d a set of variables.
    /// \return the filtered expression.
    virtual pbes_expression filter(const pbes_expression& e, const std::set<data::variable>& d)
    {
      assert(is_simple_expression(e, false));
      pbes_expression empty;
      if (is_data(e))
      {
        std::vector<data::variable> intersection;
        for (const data::variable& var: pbes_system::find_free_variables(e))
        {
          if (d.find(var) != d.end())
          {
            intersection.push_back(var);
          }
        }
        if (intersection.empty())
        {
          return e;
        }
        else
        {
          return empty;
        }
      }
      else if (is_and(e) || is_or(e) || is_imp(e))
      {
        pbes_expression l = filter(pbes_system::accessors::left(e), d);
        pbes_expression r = filter(pbes_system::accessors::right(e), d);
        if (l==empty && r==empty)
        {
          return empty;
        }
        else if (l==empty)
        {
          return r;
        }
        else if (r==empty)
        {
          if (is_imp(e))
          {
            return not_(l);
          }
          else
          {
            return l;
          }
        }
        else
        {
          if (is_and(e))
          {
            return and_(l, r);
          }
          else if (is_or(e))
          {
            return or_(l, r);
          }
          else // is_imp(e)
          {
            return imp(l, r);
          }
        }
      }
      else
      {
        std::clog << "filter: Unexpected expression: " << pp(e) << std::endl;
        throw(std::runtime_error("filter: Unexpected expression."));
      }
    }



    /// \brief Filters a 'guard' g with respect to a set of variables d and
    /// an expression phi_i such that parts of g that not relevant for phi_i
    /// are discarded.
    /// filter_guard(g, phi_i, d) = ( exists_{d intersects (free(g) - free(phi_i))} . filter(g, {d - free(phi_i)))
    ///   && filter(g, d intersects free(phi_i)).
    /// \param g a simple expression that functions as guard.
    /// \param phi_i a PBES expression.
    /// \param d a set of relevant variables.
    /// \return the expression g filtered with respect to phi_i and d.
    virtual pbes_expression filter_guard(const pbes_expression& g, const pbes_expression& phi_i, const data::variable_list& d)
    {
      pbes_expression result = true_();
      data::variable_list free_g = free_variables(g);
      std::set<data::variable> free_phi_i;
      for (const data::variable& v: pbes_system::find_free_variables(phi_i))
      {
        free_phi_i.insert(v);
      }
      std::set<data::variable> free_g_minus_free_phi_i;
      for (const data::variable& var: free_g)
      {
        if (free_phi_i.find(var)==free_phi_i.end()) { // !free_phi_i.contains(v)
          free_g_minus_free_phi_i.insert(var);
        }
      }
      std::vector<data::variable> d_intersects_free_g_minus_free_phi_i;
      std::set<data::variable> d_minus_free_phi_i;
      std::set<data::variable> d_intersects_free_phi_i;
      for (const data::variable& var: d)
      {
        if (free_g_minus_free_phi_i.find(var) != free_g_minus_free_phi_i.end()) // free_g_minus_free_phi_i.contains(v)
        {
          d_intersects_free_g_minus_free_phi_i.push_back(var);
        }
        if (free_phi_i.find(var) != free_phi_i.end()) // free_phi_i.contains(v)
        {
          d_intersects_free_phi_i.insert(var);
        }
        else
        {
          d_minus_free_phi_i.insert(var);
        }
      }
      // Now I have the sets I need; let's generate terms.
      pbes_expression e_1 = filter(g, d_minus_free_phi_i);
      if (!d_intersects_free_g_minus_free_phi_i.empty())
      {
        e_1 =
// data::exists(    // N.B. Removing this, since it does not make sense to convert a pbes_system::exists to a data::exists (Wieger).
            optimized_exists(
            data::variable_list(d_intersects_free_g_minus_free_phi_i.begin(), d_intersects_free_g_minus_free_phi_i.end()),
            e_1)
// )
           ;
      }
      pbes_expression e_2 = filter(g, d_intersects_free_phi_i);
      pbes_expression empty;
      if (e_1 == empty)
      {
        if (e_2 != empty)
        {
          result = e_2;
        }
      }
      else
      {
        if (e_2 == empty)
        {
          result = e_1;
        }
        else
        {
          result = and_(e_1, e_2);
        }
      }
      return result;
    }


    /// \brief Rewrites a bounded universal quantifier expression.
    /// If the subexpression is a conjunction, the quantifier over the conjunction
    /// is transformed to a conjunction of quantifier expressions.
    /// rewrite_bounded_forall(forall x . phi => /\_i psi_i) =
    ///   /\_i forall (x intersection free(psi_i)) .
    ///   filter_guard(phi, psi_i, x) => rewrite_bqnf_expression(psi_i).
    /// \param e a PBES expression
    /// \return the expression resulting from the transformation.
    virtual pbes_expression rewrite_bounded_forall(const pbes_expression& e)
    {
      assert(is_forall(e));
      data::variable_list qvars = quantifier_variables(e);
      pbes_expression qexpr = pbes_system::accessors::arg(e);
      while (is_forall(qexpr)) {
        qvars = qvars + quantifier_variables(qexpr);
        qexpr = pbes_system::accessors::arg(qexpr);
      }
      // forall qvars . qexpr
      pbes_expression result;
      if (is_propositional_variable_instantiation(qexpr) || is_simple_expression(qexpr, false)) {
        // forall d . phi | forall d . X(e)
        result = e;
      } else {
        pbes_expression phi = is_or(qexpr) ? static_cast<const pbes_expression&>(true_()) : static_cast<const pbes_expression&>(false_());
        pbes_expression psi = qexpr;
        if (is_or(qexpr) || is_imp(qexpr)) {
          pbes_expression l = pbes_system::accessors::left(qexpr);
          pbes_expression r = pbes_system::accessors::right(qexpr);
          if (is_simple_expression(l, false)) {
            phi = l;
            psi = r;
          }
        }
        if (is_propositional_variable_instantiation(psi)) {
          // forall d . phi => X(e)
          result = e;
        } else if (!qvars.empty() || !(is_or(qexpr) ? is_true(phi) : is_false(phi))) {
          // forall d . phi => psi
          std::vector<pbes_expression> conjuncts;
          if (is_and(psi)) {
            conjuncts = split_conjuncts(psi);
          } else {
            conjuncts.push_back(psi);
          }
          pbes_expression conjunction = true_();
          for (auto phi_i: conjuncts)
          {
            pbes_expression r = rewrite_bqnf_expression(phi_i);
            if (is_or(qexpr)) {
              if (!is_true(phi)) {
                phi = filter_guard(phi, r, qvars);
                r = or_(phi, r);
              }
            } else {
              if (!is_false(phi)) {
                phi = filter_guard(phi, r, qvars);
                r = imp(phi, r);
              }
            }
            std::vector<data::variable> qvars_i;
            std::set<data::variable> free_phi_i;
            for (const data::variable& v : pbes_system::find_free_variables(phi_i))
            {
              free_phi_i.insert(v);
            }
            for (const data::variable& var: qvars)
            {
              if (free_phi_i.find(var) != free_phi_i.end()) // free_phi_i.contains(v)
              {
                qvars_i.push_back(var);
              }
            }
            // qvars_i = qvars intersects free(phi_i)
            if (!qvars_i.empty()) {
              r = optimized_forall(data::variable_list(qvars_i.begin(), qvars_i.end()), r);
            }
            if (is_true(conjunction)) {
              conjunction = r;
            } else {
              conjunction = and_(conjunction, r);
            }
          }
          result = conjunction;
        } else {
          // qexpr not of the form phi => psi
          //std::clog << "rewrite_bounded_forall: unexpected qexpr = " << pp(qexpr) << std::endl;
          throw(std::runtime_error("rewrite_bounded_forall: unexpected qexpr."));
        }
      }
      return result;
    }



    /// \brief Rewrites a bounded existential quantifier expression.
    /// rewrite_bounded_exists(exists v . phi) = exists v. rewrite_bqnf_expression(phi)
    /// \param e a PBES expression
    /// \return the expression resulting from the transformation.
    virtual pbes_expression rewrite_bounded_exists(const pbes_expression& e)
    {
      assert(is_exists(e));
      pbes_expression qexpr = pbes_system::accessors::arg(e);
      data::variable_list qvars = quantifier_variables(e);
      while (is_exists(qexpr)) {
        qvars = qvars + quantifier_variables(qexpr);
        qexpr = pbes_system::accessors::arg(qexpr);
      }
      pbes_expression r = rewrite_bqnf_expression(qexpr);
      pbes_expression result = optimized_exists(qvars, r);
      return result;
    }


    /// \brief Rewrites a conjunctive expression.
    /// rewrite_imp(/\_i phi_i) = /\_i rewrite_bqnf_expression(phi_i).
    /// \param e a PBES expression
    /// \return the expression resulting from the transformation.
    virtual pbes_expression rewrite_and(const pbes_expression& e)
    {
      pbes_expression conjunction = true_();
      std::vector<pbes_equation> new_eqns;
      std::vector<pbes_expression> conjuncts = split_conjuncts(e);
      for (auto expr: conjuncts)
      {
        pbes_expression r = rewrite_bqnf_expression(expr);
        if (is_true(conjunction)) {
          conjunction = r;
        } else {
          conjunction = and_(conjunction, r);
        }
      }
      return conjunction;
    }


    /// \brief Rewrites a disjunctive expression.
    /// rewrite_imp(\/_i phi_i) = \/_i rewrite_bqnf_expression(phi_i).
    /// \param e a PBES expression
    /// \return the expression resulting from the transformation.
    virtual pbes_expression rewrite_or(const pbes_expression& e)
    {
      pbes_expression disjunction = false_();
      std::vector<pbes_expression> new_exprs;
      std::vector<pbes_expression> disjuncts = split_disjuncts(e);
      for (auto expr: disjuncts)
      {
        pbes_expression r = rewrite_bqnf_expression(expr);
        if (is_false(disjunction)) {
          disjunction = r;
        } else {
          disjunction = or_(disjunction, r);
        }
      }
      return disjunction;
    }


    /// \brief Rewrites an implication:
    /// rewrite_imp(phi => psi) =
    ///   rewrite_bqnf_expression(phi) => rewrite_bqnf_expression(psi).
    /// \param e a PBES expression
    /// \return the expression resulting from the transformation.
    virtual pbes_expression rewrite_imp(const pbes_expression& e)
    {
      pbes_expression l = rewrite_bqnf_expression(pbes_system::accessors::left(e));
      pbes_expression r = rewrite_bqnf_expression(pbes_system::accessors::right(e));
      pbes_expression result = imp(l, r);
      return result;
    }


    /// \brief Rewrites a BQNF expression.
    /// Replaces universal quantifier over conjuncts with conjuncts over universal quantifiers:
    /// rewrite_bqnf_expression(forall x . /\_i phi_i) =
    ///   /\_i forall (x intersection free(phi_i)) . phi_i
    /// For subexpressions the transformation is done in the usual way:
    /// rewrite_bqnf_expression(phi 'op' psi) =
    ///   rewrite_bqnf_expression(phi) 'op' rewrite_bqnf_expression(psi).
    /// \param e a PBES expression
    /// \return the expression resulting from the transformation.
    virtual pbes_expression rewrite_bqnf_expression(const pbes_expression& e)
    {
      pbes_expression result;
      if (is_propositional_variable_instantiation(e) || is_simple_expression(e, false)) {
        // Eqn of the form sigma X(d: D) = phi && Y(h(d, l)), with phi a simple formula.
        // Add sigma X(d) = e.
        result = e;
      } else if (is_forall(e)) {
        result = rewrite_bounded_forall(e);
      } else if (is_exists(e)) {
        result = rewrite_bounded_exists(e);
      } else if (is_or(e)) {
        result = rewrite_or(e);
      } else if (is_imp(e)) {
        result = rewrite_imp(e);
      } else if (is_and(e)) {
        result = rewrite_and(e);
      } else {
        std::clog << "Unexpected expression: " << pp(e) << std::endl;
        throw(std::runtime_error("Unexpected expression."));
      }
      return result;
    }

  };

  } // namespace mcrl2::pbes_system::detail

#endif // MCRL2_PBES_DETAIL_BQNF_QUANTIFIER_REWRITER_H
