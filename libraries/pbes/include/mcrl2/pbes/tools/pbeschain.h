// Author(s): Jore Booy
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/tools/pbeschain.h
/// \brief This file provides a tool that can simplify PBESs by
///        substituting PBES equations for variables in the rhs,
///        simplifying the result, and keeping it when it can
///        eliminate PBES variables.

#ifndef MCRL2_PBES_TOOLS_PBESCHAIN_H
#define MCRL2_PBES_TOOLS_PBESCHAIN_H

#include "mcrl2/data/detail/prover/bdd_prover.h"
#include "mcrl2/data/rewriter.h"
#include "mcrl2/pbes/algorithms.h"
#include "mcrl2/pbes/detail/iteration_builders.h"
#include "mcrl2/pbes/detail/stategraph_pbes.h"
#include "mcrl2/pbes/io.h"
#include "mcrl2/pbes/pbes_equation.h"
#include "mcrl2/pbes/pbes_expression.h"
#include "mcrl2/pbes/rewrite.h"
#include "mcrl2/pbes/rewriters/data2pbes_rewriter.h"
#include "mcrl2/pbes/rewriters/pbes2data_rewriter.h"
#include <cstddef>
#include <chrono>

namespace mcrl2::pbes_system
{

struct pbeschain_options
{
  data::rewrite_strategy rewrite_strategy = data::rewrite_strategy::jitty;
  bool use_bdd_simplifier = false;
  double bdd_timeout = 0.25;
  bool back_substitution = true;
  int max_depth = 12;
  bool count_unique_pvi = false;
  bool fill_pvi = false;
  double timeout = 0.0; // timeout in seconds, 0 = no timeout
};

// Substitutor to target specific path, replace our specific pvi with true/false
template<template<class> class Builder>
struct substitute_propositional_variables_for_true_false_builder
    : public Builder<substitute_propositional_variables_for_true_false_builder<Builder>>
{
  using super = Builder<substitute_propositional_variables_for_true_false_builder<Builder>>;
  using super::apply;

  simplify_quantifiers_data_rewriter<data::rewriter> m_pbes_rewriter;
  propositional_variable_instantiation m_pvi;
  pbes_expression m_replacement;

  explicit substitute_propositional_variables_for_true_false_builder(
      simplify_quantifiers_data_rewriter<data::rewriter>& r)
      : m_pbes_rewriter(r)
  {}

  void set_pvi(const propositional_variable_instantiation x) { m_pvi = x; }
  void set_replacement(const pbes_expression x) { m_replacement = x; }

  template<class T>
  void apply(T& result, const propositional_variable_instantiation& x)
  {
    if (x == m_pvi)
    {
      result = m_replacement;
    }
    else
    {
      result = x;
    }
  }
};

template<template<class> class Builder>
struct rewrite_if_builder : public Builder<rewrite_if_builder<Builder>>
{
  using super = Builder<rewrite_if_builder<Builder>>;
  using super::apply;

  simplify_data_rewriter<data::rewriter> m_pbes_rewriter;

  explicit rewrite_if_builder(simplify_data_rewriter<data::rewriter>& r)
      : m_pbes_rewriter(r)
  {}

  data::data_expression simpl_data(const data::data_expression& x)
  {
    data::data_expression res = x;
    if (data::is_application(x))
    {
      data::application da = atermpp::down_cast<data::application>(x);

      if (data::is_function_symbol(da.head()))
      {
        data::function_symbol fname = atermpp::down_cast<data::function_symbol>(da.head());
        if (fname.name() == atermpp::aterm_string("if") && (da.head().sort().target_sort() == data::bool_()))
        {
          res = data::or_(data::and_(da[0], simpl_data(da[1])), data::and_(data::not_(da[0]), simpl_data(da[2])));
          res = atermpp::down_cast<data::data_expression>(m_pbes_rewriter(res));
        }
      }
    }
    return res;
  }

  template<class T>
  void apply(T& result, const data::data_expression& x)
  {
    result = simpl_data(x);
  }
};

template<template<class> class Builder>
struct substitute_propositional_variables_builder : public Builder<substitute_propositional_variables_builder<Builder>>
{
  using super = Builder<substitute_propositional_variables_builder<Builder>>;
  using super::apply;

  pbes_equation m_eq;
  core::identifier_string name;
  simplify_data_rewriter<data::rewriter> m_pbes_rewriter;
  bool m_stable = false;

  explicit substitute_propositional_variables_builder(simplify_data_rewriter<data::rewriter>& r)
      : m_pbes_rewriter(r)
  {}

  void set_stable(bool b) { m_stable = b; }

  bool stable() const { return m_stable; }

  void set_equation(const pbes_equation& eq) { m_eq = eq; }

  void set_name(const core::identifier_string& s) { name = s; }

  template<class T>
  void apply(T& result, const propositional_variable_instantiation& x)
  {
    if (x.name() == m_eq.variable().name())
    {
      data::mutable_indexed_substitution sigma;
      data::data_expression_list pars = x.parameters();
      for (const data::variable& v: m_eq.variable().parameters())
      {
        data::data_expression par = pars.front();
        pars.pop_front();
        if (par.sort() != v.sort())
        {
          // Parameters do not match with variables. Ignore this substitution.
          result = x;
          return;
        }
        sigma[v] = par;
      }
      pbes_expression p = pbes_rewrite(m_eq.formula(), m_pbes_rewriter, sigma);
      std::set<propositional_variable_instantiation> set = find_propositional_variable_instantiations(p);
      if (std::all_of(set.begin(),
              set.end(),
              [this](const propositional_variable_instantiation& v) { return v.name() != m_eq.variable().name(); }))
      {
        // The result does not contain the variable m_eq.variable().name() and is therefore considered simpler.
        mCRL2log(log::verbose) << "Replaced in PBES equation for " << name << ":\n" << x << " --> " << p << "\n";
        result = p;
        m_stable = false;
        return;
      }
      mCRL2log(log::debug) << "No Replacement in PBES equation for " << name << ":\n" << x << " --> " << p << "\n";
      result = x;
      return;
    }
    result = x;
  }
};

template<template<class> class Builder>
struct substitute_propositional_variables_i_builder
    : public Builder<substitute_propositional_variables_i_builder<Builder>>
{
  typedef Builder<substitute_propositional_variables_i_builder<Builder>> super;
  using super::apply;

  simplify_data_rewriter<data::rewriter> m_pbes_rewriter;
  std::vector<detail::predicate_variable> predvars;
  int i = 0;

  explicit substitute_propositional_variables_i_builder(simplify_data_rewriter<data::rewriter>& r,
      std::vector<detail::predicate_variable> predvars)
      : m_pbes_rewriter(r),
        predvars(predvars)
  {}

  template<class T>
  void apply(T& result, const propositional_variable_instantiation& x)
  {
    result = pbes_rewrite(pbes_expression(x), m_pbes_rewriter, predvars[i].sigma());
    i++;
  }
};

/// \brief Returns all data variables that occur in a range of expressions
/// \param[in] container a container with expressions
/// \return All data variables that occur in the term t
template<typename Container>
std::vector<propositional_variable_instantiation> get_propositional_variable_instantiations(Container const& container)
{
  std::vector<propositional_variable_instantiation> result;
  pbes_system::find_propositional_variable_instantiations(container, std::inserter(result, result.end()));
  return result;
}

inline std::set<propositional_variable_instantiation> filter_pvis(const propositional_variable_instantiation& needle,
    const std::vector<propositional_variable_instantiation>& haystack)
{
  std::set<propositional_variable_instantiation> result;

  std::copy_if(haystack.begin(),
      haystack.end(),
      std::inserter(result, result.end()),
      [&](const propositional_variable_instantiation& v)
      {
        if (v.name() != needle.name())
        {
          return false;
        }

        const auto& v_params = as_vector(v.parameters());
        const auto& needle_params = as_vector(needle.parameters());

        if (v_params.size() != needle_params.size())
        {
          return false;
        }

        for (std::size_t i = 0; i < v_params.size(); ++i)
        {
          if (v_params[i] != needle_params[i])
          // if (!(data::is_variable(needle_params[i]) || v_params[i] == needle_params[i]))
          {
            return false;
          }
        }

        return true;
      });

  return result;
}

inline bool pvi_in_set(const propositional_variable_instantiation needle,
    const std::set<propositional_variable_instantiation> haystack)
{
  return std::any_of(haystack.begin(),
      haystack.end(),
      [&](const propositional_variable_instantiation& v)
      {
        if (v.name() != needle.name())
        {
          return false;
        }

        const auto& v_params = as_vector(v.parameters());
        const auto& needle_params = as_vector(needle.parameters());

        if (v_params.size() != needle_params.size())
        {
          return false;
        }

        for (size_t i = 0; i < v_params.size(); i++)
        {
          if (v_params[i] != needle_params[i])
          // if (!(data::is_variable(needle_params[i]) || v_params[i] == needle_params[i]))
          {
            return false;
          }
        }
        return true;
      });
}

inline data::data_expression pbestodata(pbes_expression& expr,
    detail::replace_other_propositional_variables_with_functions_builder<pbes_system::pbes_expression_builder>&
        replace_substituter)
{
  pbes_expression asdf;
  replace_substituter.reset_variable_list();
  replace_substituter.set_forward(true);
  replace_substituter.apply(asdf, expr);

  data::data_expression data_expr = detail::pbes2data(asdf);

  return data_expr;
}

inline pbes_expression datatopbes(data::data_expression& data_expr,
    rewrite_if_builder<pbes_system::pbes_expression_builder>& if_substituter,
    detail::replace_other_propositional_variables_with_functions_builder<pbes_system::pbes_expression_builder>&
        replace_substituter)
{
  pbes_expression expr;
  pbes_expression res;
  replace_substituter.set_forward(false);
  if_substituter.apply(expr, data_expr);
  replace_substituter.apply(res, detail::data2pbes<pbes_expression>(expr));
  replace_substituter.reset_variable_list();
  return res;
}

inline pbes_expression simplify_expr(pbes_expression& phi,
    pbeschain_options options,
    rewrite_if_builder<pbes_system::pbes_expression_builder>& if_substituter,
    detail::replace_other_propositional_variables_with_functions_builder<pbes_system::pbes_expression_builder>&
        replace_substituter,
    simplify_quantifiers_data_rewriter<data::rewriter>& pbes_rewriter,
    mcrl2::data::detail::BDD_Prover& f_bdd_prover)
{
  std::vector<propositional_variable_instantiation> phi_vector = get_propositional_variable_instantiations(phi);
    if (options.use_bdd_simplifier)
    {
      data::data_expression expr = pbestodata(phi, replace_substituter);
      f_bdd_prover.set_formula(expr);
      expr = f_bdd_prover.get_bdd();
      return datatopbes(expr, if_substituter, replace_substituter);
    }
    else
    {
      pbes_expression res = pbes_rewrite(phi, pbes_rewriter);
      if_substituter.apply(res, res);
      return res;
    }
  return phi;
}

inline void self_substitute(pbes_equation& equation,
    substitute_propositional_variables_for_true_false_builder<pbes_system::pbes_expression_builder>& pvi_substituter,
    rewrite_if_builder<pbes_system::pbes_expression_builder>& if_substituter,
    detail::replace_other_propositional_variables_with_functions_builder<pbes_system::pbes_expression_builder>&
        replace_substituter,
    simplify_quantifiers_data_rewriter<data::rewriter>& pbes_rewriter,
    simplify_data_rewriter<data::rewriter>& pbes_default_rewriter,
    mcrl2::data::detail::BDD_Prover& f_bdd_prover,
    pbeschain_options options)
{
  bool stable = false;
  auto start_time = std::chrono::high_resolution_clock::now();
  
  if (options.timeout > 0.0)
  {
    mCRL2log(log::verbose) << "Starting substitution for equation " << equation.variable().name() 
                           << " with timeout " << options.timeout << "s" << std::endl;
  }

  while (!stable)
  {
    // Check timeout
    if (options.timeout > 0.0)
    {
      auto current_time = std::chrono::high_resolution_clock::now();
      auto elapsed = std::chrono::duration<double>(current_time - start_time).count();
      if (elapsed >= options.timeout)
      {
        mCRL2log(log::verbose) << "Timeout reached (" << options.timeout << "s) for equation " 
                               << equation.variable().name() << " after " << elapsed << "s, stopping substitution" << std::endl;
        break;
      }
    }
    stable = true;
    std::set<propositional_variable_instantiation> stable_set = {}; // To record pvi that have reach a max depth
    std::vector<propositional_variable_instantiation> set
        = get_propositional_variable_instantiations(equation.formula());
    std::size_t previous_size = set.size();
    for (const propositional_variable_instantiation& x: set)
    {
      if (equation.variable().name() != x.name())
      {
        continue;
      }
      if (auto it = stable_set.find(x); it != stable_set.end())
      {
        continue;
      }
      // Check if during the substitution of the other pvi this one got cancelled out
      std::set<propositional_variable_instantiation> path = {x};

      propositional_variable_instantiation cur_x = x;
      pbes_expression result = x;

      mCRL2log(log::debug) << " -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n"
                           << equation.formula() << "\n\nStart " << (cur_x) << "\n";
      bool pvi_done = false;
      int depth = 0;
      while (!pvi_done)
      {
        depth = depth + 1;

        // (1) simplify
        data::mutable_indexed_substitution sigma;
        data::data_expression_list pars = cur_x.parameters();
        for (const data::variable& v: equation.variable().parameters())
        {
          data::data_expression par = pars.front();
          pars.pop_front();
          if (par.sort() != v.sort())
          {
            // Parameters do not match with variables. Ignore this substitution.
            result = x;
            break;
          }
          sigma[v] = par;
        }
        pbes_expression phi = pbes_rewrite(equation.formula(), pbes_default_rewriter, sigma);

        std::vector<propositional_variable_instantiation> phi_vector = get_propositional_variable_instantiations(phi);

        // (2) replace all reoccuring with true (nu) and false (mu)
        auto gauss_set = filter_pvis(cur_x, phi_vector);
        for (const auto& gauss_pvi: gauss_set)
        {
          mCRL2log(log::debug) << "Need to replace this with true/false " << pp(gauss_pvi) << "\n";
          mCRL2log(log::debug) << phi << "\n";
          mCRL2log(log::debug) << equation.formula() << "\n";

          // pbes_expression p_;
          pvi_substituter.set_pvi(gauss_pvi);
          pvi_substituter.set_replacement(equation.symbol().is_nu() ? true_() : false_());
          pvi_substituter.apply(phi, phi);

          mCRL2log(log::debug) << phi << "\n";

          mCRL2log(log::debug) << "- - - - - - - - - - - - - - - - - - - - \n";
        }

        // Simplify

        
        phi = simplify_expr(phi, options, if_substituter, replace_substituter, pbes_rewriter, f_bdd_prover);
        phi_vector = get_propositional_variable_instantiations(phi);
        int size = phi_vector.size();
        if (options.count_unique_pvi)
        {
          size = std::set(phi_vector.begin(), phi_vector.end()).size();
        }

        // (3) check if simpler
        if (size == 1 && (*phi_vector.begin()).name() == equation.variable().name())
        {
          propositional_variable_instantiation new_x = *phi_vector.begin();

          mCRL2log(log::debug) << "Trying loop " << new_x << " in path with \n";
          for (const auto& itr: path)
          {
            mCRL2log(log::debug) << itr << "\n";
          }

          if (pvi_in_set(new_x, path))
          {
            // We have already seen this, so we are in a loop.
            mCRL2log(log::debug) << "Loop, seen " << new_x << " in path after " << cur_x << "    " << phi << "\n";
            for (const auto& itr: path)
            {
              mCRL2log(log::debug) << itr << "\n";
            }
            pvi_substituter.set_pvi(new_x);
            pvi_substituter.set_replacement(equation.symbol().is_nu() ? true_() : false_());
            pvi_substituter.apply(result, phi);
            pvi_substituter.set_pvi(cur_x);
            pvi_substituter.set_replacement(result);
            pvi_substituter.apply(equation.formula(), equation.formula());

            mCRL2log(log::debug) << "new_phi " << result << "\n";
            stable = false;
            pvi_done = true;
          }
          else
          {
            // The result does not contain the variable m_eq.variable().name() and is therefore considered simpler.
            mCRL2log(log::debug) << "Replaced in PBES equation for " << cur_x << "\n-->\n"
                                 << phi << "\n[" << new_x << "]\n";
            pvi_substituter.set_pvi(cur_x);
            pvi_substituter.set_replacement(phi);
            pvi_substituter.apply(equation.formula(), equation.formula());
            cur_x = new_x;
            path.insert(new_x);
          }
        }
        else if (size == 0)
        {
          pvi_substituter.set_pvi(cur_x);
          pvi_substituter.set_replacement(phi);
          pvi_substituter.apply(equation.formula(), equation.formula());
          stable = false;
          mCRL2log(log::debug) << "Replaced in PBES equation for " << cur_x << ":\n" << x << " \n-->\n " << phi << "\n";
          pvi_done = true;
        }
        else
        {
          mCRL2log(log::debug) << "Not simpler: " << cur_x << " \n-->\n " << phi << " and size " << phi_vector.size()
                               << "\n";
          pvi_done = true;
          if (depth > 1)
          {
            stable = false;
          }
        }
        if (pvi_done)
        {
          mCRL2log(log::debug) << "Nothing further to do\n";
          break;
        }
        else if (depth > options.max_depth)
        {
          stable_set.insert(x);
          stable_set.insert(cur_x);
          pvi_done = true;
          // pvi_substituter.set_pvi(cur_x);
          // pvi_substituter.set_replacement(x);
          // pvi_substituter.apply(equation.formula(), equation.formula());
          break;
        }
      }

      std::vector<propositional_variable_instantiation> set
          = get_propositional_variable_instantiations(equation.formula());

      std::size_t current_size = set.size();
      mCRL2log(log::verbose) << "\rNew number of pvi: " << current_size << "";

      // Simplify
      if (!options.use_bdd_simplifier && (current_size == 0 || (previous_size >= current_size + 10)))
      {
        previous_size = current_size;
        equation.formula() = simplify_expr(equation.formula(),
            options,
            if_substituter,
            replace_substituter,
            pbes_rewriter,
            f_bdd_prover);
      }
    }
  }
  mCRL2log(log::verbose) << "\n" << equation.variable() << "\n is stable! \n\n\n";
}

inline void substitute(pbes_equation& into,
    const pbes_equation& by,
    substitute_propositional_variables_builder<pbes_system::pbes_expression_builder>& substituter)
{
  substituter.set_equation(by);
  substituter.set_name(into.variable().name());
  pbes_expression p;
  substituter.apply(p, into.formula());
  into.formula() = p;
}

// Fill in the parameters of the pvi based on the guard of each pvi
inline pbes fill_pvi(pbes& p, data::rewriter data_rewriter)
{
  detail::stategraph_pbes stategraph(p, data_rewriter);
  simplify_data_rewriter<data::rewriter> pbes_default_rewriter(data_rewriter);

  // Preparation
  for (detail::stategraph_equation& equation: stategraph.equations())
  {
    for (detail::predicate_variable& predvar: equation.predicate_variables())
    {
      predvar.simplify_guard();
    }
  }
  stategraph.compute_source_target_copy();

  // Reduction
  std::vector<pbes_equation> eqn = {};
  for (detail::stategraph_equation& eq: stategraph.equations())
  {
    substitute_propositional_variables_i_builder<pbes_system::pbes_expression_builder> substituter(
        pbes_default_rewriter,
        eq.predicate_variables());

    pbes_expression new_formula;
    substituter.apply(new_formula, eq.formula());
    pbes_equation new_eq(eq);
    new_eq.formula() = new_formula;
    eqn.push_back(new_eq);
  }

  // Back to pbes
  pbes res(p.data(), p.global_variables(), eqn, p.initial_state());
  return res;
}

struct pbeschain_pbes_backward_substituter
{
  void run(pbes& p, pbeschain_options options)
  {
    data::rewriter data_rewriter(p.data(), options.rewrite_strategy);
    data::rewriter data_default_rewriter(p.data());
    simplify_quantifiers_data_rewriter<data::rewriter> pbes_rewriter(data_rewriter);
    simplify_data_rewriter<data::rewriter> pbes_rewriter2(data_rewriter);
    simplify_data_rewriter<data::rewriter> pbes_default_rewriter(data_default_rewriter);
    substitute_propositional_variables_builder<pbes_system::pbes_expression_builder> substituter(pbes_default_rewriter);
    rewrite_if_builder<pbes_system::pbes_expression_builder> if_rewriter(pbes_default_rewriter);
    substitute_propositional_variables_for_true_false_builder<pbes_system::pbes_expression_builder> pvi_substituter(
        pbes_rewriter);

    if (options.fill_pvi)
    {
      p = fill_pvi(p, data_rewriter);
    }

    mcrl2::data::detail::BDD_Prover f_bdd_prover(p.data(),
        data::used_data_equation_selector(p.data()),
        mcrl2::data::jitty,
        options.bdd_timeout);
    detail::replace_other_propositional_variables_with_functions_builder<pbes_system::pbes_expression_builder>
        replace_substituter(pbes_rewriter2);
    
    for (std::vector<pbes_equation>::reverse_iterator i = p.equations().rbegin(); i != p.equations().rend(); i++)
    {
      mCRL2log(log::verbose) << "Investigating the equation for " << i->variable().name() << "\n";
      self_substitute(*i,
          pvi_substituter,
          if_rewriter,
          replace_substituter,
          pbes_rewriter,
          pbes_default_rewriter,
          f_bdd_prover,
          options);

      std::set<propositional_variable_instantiation> pvi_set
          = find_propositional_variable_instantiations((*i).formula());
      mCRL2log(log::verbose) << "How many are left? " << pvi_set.size() << "\n";

      if (pvi_set.size() == 0 && options.back_substitution)
      {
        for (std::vector<pbes_equation>::reverse_iterator j = i + 1; j != p.equations().rend(); j++)
        {
          substitute(*j, *i, substituter);
        }
      }
    }
  }
};

inline void pbeschain(const std::string& input_filename,
    const std::string& output_filename,
    const utilities::file_format& input_format,
    const utilities::file_format& output_format,
    pbeschain_options options)
{
  pbes p;
  load_pbes(p, input_filename, input_format);
  algorithms::normalize(p);
  pbeschain_pbes_backward_substituter backward_substituter;
  backward_substituter.run(p, options);
  save_pbes(p, output_filename, output_format);
}

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_TOOLS_PBESCHAIN_H
