// Author(s): Muck van Weerdenburg
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/detail/rewrite.h

#ifndef MCRL2_DATA_DETAIL_REWRITE_H
#define MCRL2_DATA_DETAIL_REWRITE_H

#include "mcrl2/data/detail/enumerator_identifier_generator.h"
#include "mcrl2/data/rewrite_strategy.h"
#include "mcrl2/data/selection.h"
#include "mcrl2/data/substitutions/mutable_indexed_substitution.h"



namespace mcrl2::data::detail
{

/**
 * \brief Rewriter interface class.
 *
 * This is the interface class for the rewriters. To create a specific
 * rewriter, use createRewriter.
 *
 * Simple use of the rewriter would be as follows (with t a term in the mCRL2
 * internal format):
 *
 * \code
 *   std::shared_ptr<Rewriter> r = createRewriter(equations);
 *   t = r->rewrite(t);
 * \endcode
 **/
class Rewriter
{
  protected:
    enumerator_identifier_generator m_generator;  //name for variables.

    /** \brief The copy assignment operator is protected. Public copying is not allowed.
    **/
    Rewriter& operator=(const Rewriter& other) = default;

    /** \brief The copy constructor operator is protected. Public copying is not allowed.
    **/
    Rewriter(const Rewriter& other) = default;

  public:
    using substitution_type = mutable_indexed_substitution<>;

    used_data_equation_selector data_equation_selector;
    /**
     * \brief Constructor. Do not use directly; use createRewriter()
     *        function instead.
     * \sa createRewriter()
     **/
    Rewriter(const data_specification& data_spec, const used_data_equation_selector& eq_selector):
          data_equation_selector(eq_selector),
          m_data_specification_for_enumeration(data_spec)
    {
    }

    /** \brief Destructor. */
    virtual ~Rewriter() = default;

    /** \brief The fresh name generator of the rewriter */
    data::enumerator_identifier_generator& identifier_generator()
    {
      return m_generator;
    }

    /**
     * \brief Get rewriter strategy that is used.
     * \return Used rewriter strategy.
     **/
    virtual rewrite_strategy getStrategy() = 0;

    /**
     * \brief Rewrite an mCRL2 data term.
     * \param term The term to be rewritten. This term should be a data_term
     * \return The normal form of term.
     **/
    virtual data_expression rewrite(const data_expression& term, substitution_type& sigma) = 0;

    /**
     * \brief Rewrite an mCRL2 data term.
     * \param result The normal form of term.
     * \param term The term to be rewritten. This term should be a data_term
     * \post result contains the normal form of term.
     **/
    virtual void rewrite(data_expression& result, const data_expression& term, substitution_type& sigma) = 0;

    /**
     * \brief Provide the rewriter with a () operator, such that it can also
     *        rewrite terms using this operator.
     **/

    data_expression operator() (const data_expression& term, substitution_type& sigma)
    {
      return rewrite(term,sigma);
    }

    /**
     * \brief Clone a rewriter.
     * \return A (pointer to a) a clone of the rewriter.
     **/
    virtual std::shared_ptr<detail::Rewriter> clone() = 0;

  public:
  /* The functions below are public, because they are used in the compiling jitty rewriter */
    void existential_quantifier_enumeration(data_expression& result,
                                            const abstraction& t,
                                            substitution_type& sigma,
                                            bool t1_is_normal_form = false);

    void existential_quantifier_enumeration(data_expression& result,
                                            const variable_list& vl,
                                            const data_expression& t1,
                                            substitution_type& sigma,
                                            bool t1_is_normal_form = false);

    void universal_quantifier_enumeration(data_expression& result,
                                          const abstraction& t,
                                          substitution_type& sigma,
                                          bool t1_is_normal_form = false);

    void universal_quantifier_enumeration(data_expression& result,
                                          const variable_list& vl,
                                          const data_expression& t1,
                                          substitution_type& sigma,
                                          bool t1_is_normal_form = false);

    // Rewrite a where expression where the subdataexpressions are in internal format.
    // It yields a term without a where expression. The result is passed back in the variable result. 
    void rewrite_where(data_expression& result,
                       const where_clause& term,
                       substitution_type& sigma);

    // Rewrite an expression with a lambda as outermost symbol. The expression is in internal format.
    // Bound variables are replaced by new variables to avoid a clash with variables in the right hand sides
    // of sigma.

    void rewrite_single_lambda(data_expression& result,
      const variable_list& vl,
      const data_expression& body,
      substitution_type& sigma,
      bool body_in_normal_form);

    /// Rewrite t, assuming that the headsymbol of t, which can be nested, is a lambda term.
    void rewrite_lambda_application(data_expression& result,
      const data_expression& t,
      substitution_type& sigma,
      bool arguments_are_in_normal_form = false);

    void rewrite_lambda_application(data_expression& result,
      const abstraction& lambda_term,
      const application& t,
      substitution_type& sigma,
      bool arguments_are_in_normal_form = false);

    virtual void thread_initialise()
    {
    }

  protected:

    mcrl2::data::data_specification m_data_specification_for_enumeration;

    void quantifier_enumeration(data_expression& result,
      const variable_list& vl,
      const data_expression& t1,
      bool t1_is_normal_form,
      substitution_type& sigma,
      const binder_type& binder,
      data_expression (*lazy_op)(const data_expression&, const data_expression&),
      const data_expression& identity_element,
      const data_expression& absorbing_element);
};

/**
 * \brief Create a rewriter.
 * \param DataSpec A data specification.
 * \param Strategy The rewrite strategy to be used by the rewriter.
 * \return A (pointer to a) rewriter that uses the data specification DataSpec
 *         and strategy Strategy to rewrite.
 **/
std::shared_ptr<detail::Rewriter> createRewriter(const data_specification& DataSpec,
  const used_data_equation_selector& equations_selector,
  rewrite_strategy Strategy = jitty);

/**
 * \brief Check that an mCRL2 data equation is a valid rewrite rule. If not, an runtime_error is thrown indicating the problem.
 * \param DataEqn The mCRL2 data equation to be checked.
 * \throw std::runtime_error containing a reason why DataEqn is not a valid rewrite rule.
 **/
void CheckRewriteRule(const data_equation& data_eqn);

/**
 * \brief Check whether or not an mCRL2 data equation is a valid rewrite rule.
 * \param DataEqn The mCRL2 data equation to be checked.
 * \return Whether or not DataEqn is a valid rewrite rule.
 **/
bool isValidRewriteRule(const data_equation& data_eqn);


// This function calculates the cumulated length of all
// potential function arguments.
inline std::size_t getArity(const data::function_symbol& op)
{
  sort_expression sort = op.sort();
  std::size_t arity = 0;

  while (is_function_sort(sort))
  {
    const function_sort fsort(sort);
    const sort_expression_list& sort_dom = fsort.domain();
    arity += sort_dom.size();
    sort = fsort.codomain();
  }
  return arity;
}

// This function calculates the number of direct function arguments.
inline std::size_t get_direct_arity(const data::function_symbol& op)
{
  sort_expression sort = op.sort();
  if (is_function_sort(sort))
  {
    return atermpp::down_cast<function_sort>(sort).domain().size();
  }
  return 0;
}

} // namespace mcrl2::data::detail



#endif // MCRL2_DATA_DETAIL_REWRITE_H
