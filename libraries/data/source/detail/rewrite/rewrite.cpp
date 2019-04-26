// Author(s): Muck van Weerdenburg
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include "mcrl2/data/detail/rewrite/innermost.h"
#include "mcrl2/data/detail/rewrite/jitty.h"
#include "mcrl2/data/detail/rewrite/jitty_jittyc.h"
#ifdef MCRL2_JITTYC_AVAILABLE
#include "mcrl2/data/detail/rewrite/jittyc.h"
#endif

#include "mcrl2/data/detail/rewrite/with_prover.h"

#include "mcrl2/data/detail/rewriter_wrapper.h"
#include "mcrl2/data/enumerator_with_iterator.h"

using namespace mcrl2::core;
using namespace mcrl2::core::detail;
using namespace mcrl2::log;

namespace mcrl2
{
namespace data
{
namespace detail
{

static std::size_t npos()
{
  return std::size_t(-1);
}

#ifndef NDEBUG
// function object to test if it is an aterm_appl with function symbol "f"
struct is_a_variable
{
  bool operator()(const atermpp::aterm& t) const
  {
    return is_variable(atermpp::down_cast<atermpp::aterm_appl>(t));
  }
};

static
bool occur_check(const variable& v, const atermpp::aterm_appl& e)
{
  if (v==e)
  {
    // The variable is reset. This is ok.
    return true;
  }
  std::set<variable> s;
  find_all_if(e,is_a_variable(),std::inserter(s,s.begin()));
  if (s.count(v)>0)
  {
    return false; // Occur check failed.
  }
  return true;
}
#endif

data_expression Rewriter::rewrite_where(
                      const where_clause& term,
                      substitution_type& sigma)
{
  const assignment_list& assignments = term.assignments();
  const data_expression& body=term.body();

  mutable_map_substitution<std::map < variable,data_expression> > variable_renaming;
  for (const assignment& a: assignments)
  {
    const variable& v=a.lhs();
    const variable v_fresh(m_generator(), v.sort());
    variable_renaming[v]=v_fresh;
    sigma[v_fresh]=rewrite(a.rhs(),sigma);
  }
  const data_expression result=rewrite(replace_variables(body,variable_renaming),sigma);

  // Reset variables in sigma
  for (mutable_map_substitution<std::map < variable,data_expression> >::const_iterator it=variable_renaming.begin();
      it!=variable_renaming.end(); ++it)
  {
    sigma[atermpp::down_cast<variable>(it->second)]=it->second;
  }
  return result;
}

abstraction Rewriter::rewrite_single_lambda(
                      const variable_list& vl,
                      const data_expression& body,
                      const bool body_in_normal_form,
                      substitution_type& sigma)
{
  assert(vl.size()>0);
  // A lambda term without arguments; Take care that the bound variable is made unique with respect to
  // the variables occurring in sigma. But in case vl is empty, just rewrite...

  // First filter the variables in vl by those occuring as left/right hand sides in sigma.

  std::size_t number_of_renamed_variables=0;
  std::size_t count=0;
  std::vector<variable> new_variables(vl.size());
  {
    // Create new unique variables to replace the old and create storage for
    // storing old values for variables in vl.
    for(const variable& v: vl)
    {
      if (sigma(v)!=v || sigma.variable_occurs_in_a_rhs(v))
      {
        number_of_renamed_variables++;
        new_variables[count]=data::variable(m_generator(), v.sort());
        assert(occur_check(v, new_variables[count]));
      }
      else
      {
        new_variables[count]=v;
      }
      count++;
    }
  }

  if (number_of_renamed_variables==0)
  {
    return abstraction(lambda_binder(),vl,(body_in_normal_form?body:rewrite(body,sigma)));
  }

  data_expression result;
  variable_list::const_iterator v;
  if (body_in_normal_form)
  {
    // If the body is already in normal form, a simple replacement of the old variables
    // by the new ones will do
    mutable_map_substitution<std::map<variable,data_expression> > variable_renaming;
    for(v = vl.begin(), count = 0; v != vl.end(); ++v, ++count)
    {
      if (*v != new_variables[count])
      {
        variable_renaming[*v] = new_variables[count];
      }
    }
    result = replace_variables(body, variable_renaming);
  }
  else
  {
    // If the body is not in normal form, then we have to rewrite with an updated sigma.
    // We first change sigma and save the values in sigma we overwrote...
    std::vector<data_expression> saved_substitutions;
    for(v = vl.begin(), count = 0; v != vl.end(); ++v, ++count)
    {
      if (*v != new_variables[count])
      {
        saved_substitutions.push_back(sigma(*v));
        sigma[*v] = new_variables[count];
      }
    }
    // ... then we rewrite with the new sigma ...
    result = rewrite(body,sigma);
    // ... and then we restore sigma to its old state.
    std::size_t new_variable_count = 0;
    for(v = vl.begin(), count = 0; v != vl.end(); ++v, ++count)
    {
      if (*v != new_variables[count])
      {
        sigma[*v] = saved_substitutions[new_variable_count++];
      }
    }
  }
  variable_list new_variable_list(new_variables.begin(), new_variables.end());
  return abstraction(lambda_binder(),new_variable_list,result);
}


// The function rewrite_lambda_application assumes that t has the shape
// application(...application(lambda x:D...., arg1,...,argn),argn+1,...,argN).
// It applies the lambda term to its arguments, and rewrites the result to
// normal form.

data_expression Rewriter::rewrite_lambda_application(
                      const data_expression& t,
                      substitution_type& sigma)
{
  if (is_lambda(t))
  {
    const abstraction& ta(t);
    return rewrite_single_lambda(ta.variables(),ta.body(),false,sigma);
  }

  const application ta(t);
  if (is_lambda(ta.head()))
  {
    return rewrite_lambda_application(ta.head(),ta,sigma);
  }

  return rewrite(application(rewrite_lambda_application(ta.head(),sigma),ta.begin(),ta.end()),sigma);
}



// The function rewrite_lambda_application rewrites a lambda term to a set of
// arguments which are the arguments 1,...,n of t. If t has the shape
// #REWR#(t0,t1,....,tn) and the lambda term is L, we calculate the normal form
// in internal format for L(t1,...,tn). Note that the term t0 is ignored.
// Note that we assume that neither L, nor t is in normal form.

data_expression Rewriter::rewrite_lambda_application(
                      const abstraction& lambda_term,
                      const application& t,
                      substitution_type& sigma)
{
  assert(is_lambda(lambda_term));  // The function symbol in this position cannot be anything else than a lambda term.
  const variable_list& vl=lambda_term.variables();
  const data_expression& lambda_body=lambda_term.body();
  std::size_t arity=t.size();
  if (arity==0) // The term has shape application(lambda d..:D...t), i.e. without arguments.
  {
    data_expression r= rewrite_single_lambda(vl, lambda_body, true, sigma);
    return r;
  }
  assert(vl.size()<=arity);

  // The variable vl_backup will be used to first store the values to be substituted
  // for the variables in vl. Subsequently, it will be used to temporarily save the values in sigma
  // assigned to vl.
  data_expression* vl_backup = MCRL2_SPECIFIC_STACK_ALLOCATOR(data_expression,vl.size());

  // Calculate the values that must be substituted for the variables in vl and store these in vl_backup.
  for(std::size_t count=0; count<vl.size(); count++)
  {
    new (&vl_backup[count]) data_expression(rewrite(data_expression(t[count]),sigma));
  }

  // Swap the values assigned to variables in vl with those in vl_backup.
  std::size_t count=0;
  for(const variable& v: vl)
  {
    const data_expression temp=sigma(v);
    sigma[v]=vl_backup[count];
    vl_backup[count]=temp;
    count++;
  }

  const data_expression result=rewrite(lambda_body,sigma);

  // Reset variables in sigma and destroy the elements in vl_backup.
  count=0;
  for(const variable& v: vl)
  {
    sigma[v]=vl_backup[count];
    vl_backup[count].~data_expression();
    count++;
  }

  if (vl.size()==arity)
  {
    return result;
  }

  // There are more arguments than bound variables.
  // Rewrite the remaining arguments and apply the rewritten lambda term to them.
  return application(result,
                     t.begin()+vl.size(),
                     t.end(),
                     [this, &sigma](const data_expression& t) -> data_expression { return rewrite(t, sigma); });
}

data_expression Rewriter::existential_quantifier_enumeration(
     const abstraction& t,
     substitution_type& sigma)
{
  // This is a quantifier elimination that works on the existential quantifier as specified
  // in data types, i.e. without applying the implement function anymore.

  assert(is_exists(t));
  return existential_quantifier_enumeration(t.variables(), t.body(), false, sigma);
}

// Generate a term equivalent to exists vl.t1.
// The variable t1_is_normal_form indicates whether t1 is in normal
// form, but this information is not used as it stands.
data_expression Rewriter::existential_quantifier_enumeration(
      const variable_list& vl,
      const data_expression& t1,
      const bool t1_is_normal_form,
      substitution_type& sigma)
{
  return quantifier_enumeration(vl, t1, t1_is_normal_form, sigma, exists_binder(), &lazy::or_, sort_bool::false_(), sort_bool::true_());
}


data_expression Rewriter::universal_quantifier_enumeration(
     const abstraction& t,
     substitution_type& sigma)
{
  assert(is_forall(t));
  return universal_quantifier_enumeration(t.variables(),t.body(),false,sigma);
}

// Generate a term equivalent to forall vl.t1.
// The variable t1_is_normal_form indicates whether t1 is in normal
// form, but this information is not used as it stands.
data_expression Rewriter::universal_quantifier_enumeration(
      const variable_list& vl,
      const data_expression& t1,
      const bool t1_is_normal_form,
      substitution_type& sigma)
{
  return quantifier_enumeration(vl, t1, t1_is_normal_form, sigma, forall_binder(), &lazy::and_, sort_bool::true_(), sort_bool::false_());
}

data_expression Rewriter::quantifier_enumeration(
      const variable_list& vl,
      const data_expression& t1,
      const bool t1_is_normal_form,
      substitution_type& sigma,
      const binder_type& binder,
      data_expression (*lazy_op)(const data_expression&, const data_expression&),
      const data_expression& identity_element,
      const data_expression& absorbing_element
    )
{
  // Rename the bound variables to unique
  // variables, to avoid naming conflicts.

  mutable_map_substitution<std::map < variable,data_expression> > variable_renaming;
  variable_vector vl_new_v;
  for(const variable& v: vl)
  {
    if (sigma(v)!=v)
    {
      const variable v_fresh(m_generator(), v.sort());
      variable_renaming[v]=v_fresh;
      vl_new_v.push_back(v_fresh);
    }
    else
    {
      vl_new_v.push_back(v);
    }
  }

  const data_expression t2=replace_variables(t1,variable_renaming);
  const data_expression t3=(t1_is_normal_form?t2:rewrite(t2,sigma));

  // Check whether the bound variables occur free in the rewritten body
  std::set < variable > free_variables=find_free_variables(t3);
  variable_list vl_new_l_enum;
  variable_list vl_new_l_other;

  bool sorts_are_finite=true;
  for(variable_vector::const_reverse_iterator i=vl_new_v.rbegin(); i!=vl_new_v.rend(); ++i)
  {
    const variable v = *i;
    if (free_variables.count(v)>0)
    {
      if(is_enumerable(m_data_specification_for_enumeration, rewriter_wrapper(this), v.sort()))
      {
        vl_new_l_enum.push_front(v);
        sorts_are_finite = sorts_are_finite && m_data_specification_for_enumeration.is_certainly_finite(v.sort());
      }
      else
      {
        vl_new_l_other.push_front(v);
      }
    }
  }

  if (vl_new_l_enum.empty())
  {
    if(vl_new_l_other.empty())
    {
      return t3; // No quantified variables are bound.
    }
    else
    {
      return abstraction(binder, vl_new_l_other, t3);
    }
  }

  /* Find A solution*/
  rewriter_wrapper wrapped_rewriter(this);
  const std::size_t max_count = sorts_are_finite ? std::numeric_limits<std::size_t>::max() : data::detail::get_enumerator_iteration_limit();

  struct is_not
  {
    data_expression m_expr;
    is_not(const data_expression& expr)
    : m_expr(expr)
    {}

    bool operator()(const data_expression& expr)
    {
      return expr != m_expr;
    }
  };

  typedef enumerator_algorithm_with_iterator<rewriter_wrapper,
                                             enumerator_list_element<>,
                                             is_not,
                                             rewriter_wrapper,
                                             rewriter_wrapper::substitution_type> enumerator_type;
  try
  {
    enumerator_type enumerator(wrapped_rewriter, m_data_specification_for_enumeration,
      wrapped_rewriter, m_generator, max_count, true, is_not(identity_element));

    /* Create a list to store solutions */
    data_expression partial_result = identity_element;

    std::size_t loop_upperbound = (sorts_are_finite ? npos() : 10);
    data::enumerator_queue<enumerator_list_element<> > enumerator_solution_deque(enumerator_list_element<>(vl_new_l_enum, t3));

    enumerator_type::iterator sol = enumerator.begin(sigma, enumerator_solution_deque);
    for( ; loop_upperbound > 0 &&
           sol != enumerator.end();
           ++sol)
    {
      partial_result = lazy_op(partial_result, sol->expression());
      loop_upperbound--;
      if(partial_result == absorbing_element)
      {
        // We found a solution, so prevent the enumerator from doing any unnecessary work
        // Also prevents any further exceptions from the enumerator
        return absorbing_element;
      }
    }

    if (sol == enumerator.end())
    {
      if(vl_new_l_other.empty())
      {
        return partial_result;
      }
      else
      {
        return abstraction(binder, vl_new_l_other, partial_result);
      }
    }
    // One can consider to replace the variables by their original, in order to not show
    // internally generated variables in the output.
    assert(!sol->is_valid() || loop_upperbound == 0);
  }
  catch(const mcrl2::runtime_error&)
  {
    // It is not possible to enumerate one of the bound variables, so we just return
    // the simplified expression
  }

  return abstraction(binder,vl_new_l_enum+vl_new_l_other,rewrite(t3,sigma));
}


std::shared_ptr<Rewriter> createRewriter(
            const data_specification& data_spec,
            const used_data_equation_selector& equations_selector,
            const rewrite_strategy strategy)
{
  switch (strategy)
  {
    case innermost:
      return std::make_shared<InnermostRewriter>(data_spec,equations_selector);
    case jitty:
      return std::shared_ptr<Rewriter>(new RewriterJitty(data_spec,equations_selector));
#ifdef MCRL2_JITTYC_AVAILABLE
    case jitty_compiling:
      return std::shared_ptr<Rewriter>(new RewriterCompilingJitty(data_spec,equations_selector));
#endif
    case jitty_prover:
      return std::shared_ptr<Rewriter>(new RewriterProver(data_spec,jitty,equations_selector));
#ifdef MCRL2_JITTYC_AVAILABLE
    case jitty_compiling_prover:
      return std::shared_ptr<Rewriter>(new RewriterProver(data_spec,jitty_compiling,equations_selector));
#endif
    default: throw mcrl2::runtime_error("Cannot create a rewriter using strategy " + pp(strategy) + ".");
  }
}

//Prototype
static void check_vars(const data_expression& expr, const std::set <variable>& vars, std::set <variable>& used_vars);

static void check_vars(application::const_iterator begin,
                       const application::const_iterator& end,
                       const std::set <variable>& vars,
                       std::set <variable>& used_vars)
{
  while (begin != end)
  {
    check_vars(*begin++, vars, used_vars);
  }
}

static void check_vars(const data_expression& expr, const std::set <variable>& vars, std::set <variable>& used_vars)
{
  if (is_application(expr))
  {
    const application& a=atermpp::down_cast<application>(expr);
    check_vars(a.head(),vars,used_vars);
    check_vars(a.begin(),a.end(),vars,used_vars);
  }
  else if (is_variable(expr))
  {
    const variable& v=atermpp::down_cast<variable>(expr);
    used_vars.insert(v);

    if (vars.count(v)==0)
    {
      throw v;
    }
  }
}

//Prototype
static void checkPattern(const data_expression& p);

static void checkPattern(application::const_iterator begin,
                         const application::const_iterator& end)
{
  while (begin != end)
  {
    checkPattern(*begin++);
  }
}

static void checkPattern(const data_expression& p)
{
  if (is_application(p))
  {
    if (is_variable(atermpp::down_cast<application>(p).head()))
    {
      throw mcrl2::runtime_error(std::string("variable ") + data::pp(application(p).head()) +
               " is used as head symbol in an application, which is not supported");
    }
    const application& a=atermpp::down_cast<application>(p);
    checkPattern(a.head());
    checkPattern(a.begin(),a.end());
  }
}

void CheckRewriteRule(const data_equation& data_eqn)
{
  const variable_list& rule_var_list = data_eqn.variables();
  const std::set <variable> rule_vars(rule_var_list.begin(),rule_var_list.end());

  // collect variables from lhs and check that they are in rule_vars
  std::set <variable> lhs_vars;
  try
  {
    check_vars(data_eqn.lhs(),rule_vars,lhs_vars);
  }
  catch (variable& var)
  {
    // This should never occur if data_eqn is a valid data equation
    mCRL2log(log::error) << "Data equation: " << data_expression(data_eqn) << std::endl;
    assert(0);
    throw runtime_error("variable " + pp(var) + " occurs in left-hand side of equation but is not defined (in equation: " + pp(data_eqn) + ")");
  }

  // check that variables from the condition occur in the lhs
  try
  {
    std::set <variable> dummy;
    check_vars(data_eqn.condition(),lhs_vars,dummy);
  }
  catch (variable& var)
  {
    throw runtime_error("variable " + pp(var) + " occurs in condition of equation but not in left-hand side (in equation: " +
                    pp(data_eqn) + "); equation cannot be used as rewrite rule");
  }

  // check that variables from the rhs are occur in the lhs
  try
  {
    std::set <variable> dummy;
    check_vars(data_eqn.rhs(),lhs_vars,dummy);
  }
  catch (variable& var)
  {
    throw runtime_error("variable " + pp(var) + " occurs in right-hand side of equation but not in left-hand side (in equation: " +
                pp(data_eqn) + "); equation cannot be used as rewrite rule");
  }

  // check that the lhs is a supported pattern
  if (is_variable(data_eqn.lhs()))
  {
    throw runtime_error("left-hand side of equation is a variable; this is not allowed for rewriting");
  }
  try
  {
    checkPattern(data_eqn.lhs());
  }
  catch (mcrl2::runtime_error& s)
  {
    throw runtime_error(std::string(s.what()) + " (in equation: " + pp(data_eqn) + "); equation cannot be used as rewrite rule");
  }
}

bool isValidRewriteRule(const data_equation& data_eqn)
{
  try
  {
    CheckRewriteRule(data_eqn);
    return true;
  }
  catch (runtime_error&)
  {
    return false;
  }
  return false; // compiler warning
}

}
}
}
