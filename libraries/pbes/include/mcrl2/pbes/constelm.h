// Author(s): Wieger Wesselink
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/constelm.h
/// \brief The constelm algorithm.

#ifndef MCRL2_PBES_CONSTELM_H
#define MCRL2_PBES_CONSTELM_H

#include <ranges>

#include "mcrl2/pbes/algorithms.h"
#include "mcrl2/pbes/pbes_rewriter_type.h"
#include "mcrl2/pbes/print.h"
#include "mcrl2/pbes/replace.h"
#include "mcrl2/pbes/rewriters/enumerate_quantifiers_rewriter.h"
#include "mcrl2/pbes/detail/pbes_remove_counterexample_info.h"



namespace mcrl2::pbes_system
{

namespace detail
{

/// \cond INTERNAL_DOCS
inline
void make_constelm_substitution(const std::map<data::variable, data::data_expression>& m, data::rewriter::substitution_type& result)
{
  for (const auto& i : m)
  {
    result[i.first] = i.second;
  }
}

class quantified_variable
{
protected:
  bool m_is_forall;
  data::variable m_var;

public:
  quantified_variable(bool is_forall, const data::variable& var)
  : m_is_forall(is_forall)
  , m_var(var)
  {}

  bool is_forall() const
  {
    return m_is_forall;
  }

  const data::variable& variable() const
  {
    return m_var;
  }

  bool operator==(const quantified_variable& other) const
  {
    return m_is_forall == other.m_is_forall && m_var == other.m_var;
  }

  bool operator!=(const quantified_variable& other) const
  {
    return !(*this == other);
  }

  bool operator<(const quantified_variable& other) const
  {
    return m_is_forall < other.m_is_forall || (m_is_forall == other.m_is_forall && m_var < other.m_var);
  }

  pbes_expression make_expr(const pbes_expression& expr) const
  {
    return m_is_forall ? pbes_expression(forall(data::variable_list({m_var}), expr)) : pbes_expression(exists(data::variable_list({m_var}), expr));
  }

  std::string to_string() const
  {
    std::ostringstream out;
    out << (is_forall() ? "forall " : "exists ") << variable() << ": " << variable().sort() << ". ";
    return out.str();
  }
};

/// \brief A quantified predicate variable instantiation
struct QPVI
{
  std::list<quantified_variable> Q;
  propositional_variable_instantiation X_e;

  bool operator<(const QPVI& other) const
  {
    return std::tie(Q, X_e) < std::tie(other.Q, other.X_e);
  }
};

struct edge_details
{
  /// \brief Contains expressions that characterise when an edge is enabled.
  /// The conjunction of these expressions is a guard for some PVI.
  std::set<data::data_expression> conditions;
  /// \brief The set of free variables that occur on the other side of the conjunctions this
  /// PVI occurs in. Can be used to determine whether the quantifier inside rewriter
  /// can manage to push an existential quantifier all the way to this PVI.
  std::set<data::variable> conjunctive_context_FV;
  /// \brief The set of free variables that occur on the other side of the conjunctions this
  /// PVI occurs in. Can be used to determine whether the quantifier inside rewriter
  /// can manage to push a universal quantifier all the way to this PVI.
  std::set<data::variable> disjunctive_context_FV;
};

struct edge_traverser_stack_elem
{
  using edge_map = std::multimap<QPVI, edge_details>;

  data::data_expression Cpos;
  data::data_expression Cneg;
  std::set<data::variable> FV;
  edge_map edges;

  edge_traverser_stack_elem(const data::data_expression& cond_pos, const data::data_expression& cond_neg, std::set<data::variable>&& free_vars)
    : Cpos(cond_pos), Cneg(cond_neg), FV(std::move(free_vars))
  {}
};

struct edge_condition_traverser: public pbes_expression_traverser<edge_condition_traverser>
{
  using super = pbes_expression_traverser<edge_condition_traverser>;
  using super::enter;
  using super::leave;
  using super::apply;

  using stack_elem = edge_traverser_stack_elem;
  using edge_map = stack_elem::edge_map;
  using qvar_list = std::list<detail::quantified_variable>;

  std::vector<stack_elem> condition_fv_stack;
  std::list<pbes_expression> quantified_context;

  template <class... Args>
  void emplace(Args&&... args)
  {
    condition_fv_stack.emplace_back(std::forward<Args>(args)...);
  }

  void push(const stack_elem& x)
  {
    condition_fv_stack.push_back(x);
  }

  stack_elem& top()
  {
    return condition_fv_stack.back();
  }

  const stack_elem& top() const
  {
    return condition_fv_stack.back();
  }

  stack_elem pop()
  {
    // As the last element of the stack is removed anyway, we move the element
    // out of the stack. This avoids (rather expensive) allocation and deallocation
    // of a temporary stack_elem.
    // This cost was observed during profiling (Jeroen Keiren, 27/6/2025)
    stack_elem result = std::move(condition_fv_stack.back());
    condition_fv_stack.pop_back();
    return result;
  }

  // Merge conditions of stack elements ec1 and ec2 into ec.
  // As construction and destruction of edges (in particular edge details)
  // is expensive (observed by Jeroen Keiren, 27/6/2025), and ec1 and ec2 are
  // not used in the calling context, we here explicitly accept them as
  // rvalue reference; this allows us to move the edges out of ec1 and ec2.
  // NOLINTBEGIN(cppcoreguidelines-rvalue-reference-param-not-moved) ec1 and ec2 are consumed: their edges are moved into ec below.
  void merge_conditions(stack_elem&& ec1, bool negate1,
                        stack_elem&& ec2, bool negate2,
                        stack_elem& ec, bool is_conjunctive
                       )
  // NOLINTEND(cppcoreguidelines-rvalue-reference-param-not-moved)
  {
    for (auto& i: ec1.edges)
    {
      auto& [Q_X_e, details] = i;
      details.conditions.insert(negate2 ? ec2.Cneg : ec2.Cpos);
      (is_conjunctive ? details.conjunctive_context_FV : details.disjunctive_context_FV)
        .insert(ec2.FV.begin(), ec2.FV.end());
      ec.edges.insert(std::move(i));
    }
    for (auto& i: ec2.edges)
    {
      auto& [Q_X_e, details] = i;
      details.conditions.insert(negate1 ? ec1.Cneg : ec1.Cpos);
      (is_conjunctive ? details.conjunctive_context_FV : details.disjunctive_context_FV)
        .insert(ec1.FV.begin(), ec1.FV.end());
      ec.edges.insert(std::move(i));
    }
  }

  // enter functions related to maintaining the quantfier scope
  void enter(const not_&)
  {
    quantified_context.clear();
  }

  void enter(const and_&)
  {
    quantified_context.clear();
  }

  void enter(const or_&)
  {
    quantified_context.clear();
  }

  void enter(const imp&)
  {
    quantified_context.clear();
  }

  void enter(const forall& x)
  {
    quantified_context.push_back(x);
  }

  void enter(const exists& x)
  {
    quantified_context.push_back(x);
  }

  // leave functions, mostly used to build conditions and gather free variables
  void leave(const data::data_expression& x)
  {
    data::data_expression cond_not;
    data::make_optimized_not(cond_not, x);

    emplace(x, cond_not, data::find_free_variables(x));
  }

  void leave(const not_&)
  {
    std::swap(top().Cpos, top().Cneg);
  }

  void leave(const and_&)
  {
    stack_elem ec_right = pop();
    stack_elem ec_left = pop();
    data::data_expression cond_and;
    data::data_expression cond_or;
    data::make_optimized_and(cond_and, ec_left.Cpos, ec_right.Cpos);
    data::make_optimized_or(cond_or, ec_left.Cneg, ec_right.Cneg);

    emplace(cond_and, cond_or,
      utilities::detail::set_union(ec_left.FV, ec_right.FV));
    merge_conditions(std::move(ec_left), false, std::move(ec_right), false, top(), true);
  }

  void leave(const or_&)
  {
    stack_elem ec_right = pop();
    stack_elem ec_left = pop();
    data::data_expression cond_and;
    data::data_expression cond_or;

    data::make_optimized_and(cond_and, ec_left.Cneg, ec_right.Cneg);
    data::make_optimized_or(cond_or, ec_left.Cpos, ec_right.Cpos);

    emplace(cond_or, cond_and,
      utilities::detail::set_union(ec_left.FV, ec_right.FV));
    merge_conditions(std::move(ec_left), true, std::move(ec_right), true, top(), false);
  }

  void leave(const imp&)
  {
    stack_elem ec_right = pop();
    stack_elem ec_left = pop();
    data::data_expression cond_or;
    data::data_expression cond_and;

    data::make_optimized_or(cond_or, ec_left.Cneg, ec_right.Cpos);
    data::make_optimized_and(cond_and, ec_left.Cpos, ec_right.Cneg);

    emplace(cond_or, cond_and,
      utilities::detail::set_union(ec_left.FV, ec_right.FV));;
    merge_conditions(std::move(ec_left), false, std::move(ec_right), true, top(), false);
  }

  void leave(const forall& x)
  {
    // build conditions and free variable sets
    stack_elem ec = pop();
    for (auto& [X_e, details]: ec.edges)
    {
      auto& [cond_set, conj_FV, disj_FV] = details;

      // Update the conditions
      std::set<data::data_expression> new_conditions;
      for(const data::data_expression& e: cond_set)
      {
        data::data_expression t;
        data::make_optimized_exists(t, x.variables(), e, true);
        new_conditions.insert(t);
      }
      cond_set = std::move(new_conditions);
      data::data_expression forall;
      data::make_optimized_forall(forall, x.variables(), ec.Cpos, true);

      cond_set.insert(forall);

      // Update FV
      conj_FV = ec.FV;
    }
    data::make_optimized_forall(ec.Cpos, x.variables(), ec.Cpos, true);
    data::make_optimized_exists(ec.Cneg, x.variables(), ec.Cneg, true);
    std::set<data::variable> bound_vars{x.variables().begin(), x.variables().end()};
    ec.FV = utilities::detail::set_difference(ec.FV, bound_vars);
    push(ec);

    // maintain quantifier scope
    if(!quantified_context.empty() && quantified_context.back() == x)
    {
      quantified_context.pop_back();
    }
  }

  void leave(const exists& x)
  {
    // build conditions and free variable sets
    stack_elem ec = pop();
    for (auto& [X_e, details]: ec.edges)
    {
      auto& [cond_set, conj_FV, disj_FV] = details;

      // Update the conditions
      std::set<data::data_expression> new_conditions;
      for(const data::data_expression& e: cond_set)
      {
        data::data_expression t;
        data::make_optimized_exists(t, x.variables(), e, true);
        new_conditions.insert(t);
      }
      cond_set = std::move(new_conditions);

      data::data_expression forall;
      data::make_optimized_forall(forall, x.variables(), ec.Cneg, true);
      cond_set.insert(forall);

      // Update FV
      disj_FV = ec.FV;
    }
    data::make_optimized_exists(ec.Cpos, x.variables(), ec.Cpos, true);
    data::make_optimized_forall(ec.Cneg, x.variables(), ec.Cneg, true);
    std::set<data::variable> bound_vars{x.variables().begin(), x.variables().end()};
    ec.FV = utilities::detail::set_difference(ec.FV, bound_vars);
    push(ec);

    // maintain quantifier scope
    if(!quantified_context.empty() && quantified_context.back() == x)
    {
      quantified_context.pop_back();
    }
  }

  void leave(const propositional_variable_instantiation& x)
  {
    // Build list of qvars from quantifier scope
    qvar_list qvars;
    for(const pbes_expression& expr: quantified_context)
    {
      assert(is_forall(expr) || is_exists(expr));
      data::variable_list vars(is_forall(expr) ? atermpp::down_cast<forall>(expr).variables() : atermpp::down_cast<exists>(expr).variables());
      for(const data::variable& v: vars)
      {
        qvars.emplace_back(is_forall(expr), v);
      }
    }
    QPVI Q_X_e{.Q = qvars, .X_e = x};

    // Store the QPVI and the condition true
    emplace(data::sort_bool::true_(), data::sort_bool::true_(), data::find_free_variables(x.parameters()));
    top().edges.emplace(Q_X_e,
        edge_details{.conditions = std::set<data::data_expression>{data::sort_bool::true_()},
            .conjunctive_context_FV = std::set<data::variable>{},
            .disjunctive_context_FV = std::set<data::variable>{}});
  }

  const edge_map& result() const
  {
    assert(condition_fv_stack.size() == 1);
    return top().edges;
  }
};

} // namespace detail
/// \endcond


/// \brief Algorithm class for the constelm algorithm
template <typename DataRewriter, typename PbesRewriter>
class pbes_constelm_algorithm
{
  protected:
    /// \brief A map with constraints on the vertices of the graph
    using constraint_map = std::map<data::variable, data::data_expression>;
    using qvar_list = std::list<detail::quantified_variable>;

    /// \brief Compares data expressions for equality.
    const DataRewriter& m_data_rewriter;

    /// \brief Compares data expressions for equality.
    const PbesRewriter& m_pbes_rewriter;

    class vertex;

    /// \brief Represents an edge of the dependency graph. The assignments are stored
    /// implicitly using the 'right' parameter. The condition determines under
    /// what circumstances the influence of the edge is propagated to its target
    /// vertex.
    //
    // N.B. The attribute condition "pbes_expression condition;" needs to be protected.
    // This is achieved by deriving from pbes_expression. This is very ugly, but AFAIK
    // this is the least destructive solution to garbage collection problems.
    // Note that source and target are protected elsewhere.
    class edge: public data::data_expression
    {
      protected:
        /// \brief The propositional variable at the source of the edge
        const propositional_variable m_source;

        /// \brief The quantifiers in whose direct context the target PVI occurs
        const qvar_list m_qvars;

        /// \brief The propositional variable instantiation that determines the target of the edge
        const propositional_variable_instantiation m_target;

        const std::set<data::variable> m_conj_context;
        const std::set<data::variable> m_disj_context;

      public:
        /// \brief Constructor
        edge() = default;

        /// \brief Constructor
        /// \param src A propositional variable declaration
        /// \param tgt A propositional variable
        /// \param c A term
        edge(
          const propositional_variable& src,
          const qvar_list& qvars,
          const propositional_variable_instantiation& tgt,
          const std::set<data::variable>& conj_context,
          const std::set<data::variable>& disj_context,
          data::data_expression c = data::sort_bool::true_()
        )
        : data::data_expression(c)
        , m_source(src)
        , m_qvars(qvars)
        , m_target(tgt)
        , m_conj_context(conj_context)
        , m_disj_context(disj_context)
        {}

        /// \brief Returns a string representation of the edge.
        /// \return A string representation of the edge.
        std::string to_string() const
        {
          std::ostringstream out;
          out << "(" << m_source.name() << ", " << m_target.name() << ")  label = ";
          for(const detail::quantified_variable& qv: m_qvars)
          {
            out << qv.to_string();
          }
          out << m_target << "  condition = " << condition();
          out << "; conjunctive context = {";
          for (const data::variable& var: m_conj_context)
          {
            out << var << ": " << var.sort() << ", ";
          }
          out << "}; disjunctive context = {";
          for (const data::variable& var: m_disj_context)
          {
            out << var << ": " << var.sort() << ", ";
          }
          out << "}";
          return out.str();
        }

        /// \brief The propositional variable at the source of the edge
        const propositional_variable& source() const
        {
          return m_source;
        }

        const qvar_list& quantified_variables() const
        {
          return m_qvars;
        }

        /// \brief The propositional variable instantiation that determines the target of the edge
        const propositional_variable_instantiation& target() const
        {
          return m_target;
        }

        /// \brief The condition of the edge
        const data::data_expression& condition() const
        {
          return *this;
        }

        /// \brief Try to guess which quantifiers of Q can end up directly
        /// before target, when the quantifier inside rewriter is applied.
        qvar_list quantifier_inside_approximation(const vertex& source, const DataRewriter& rewr) const
        {
          const qvar_list& Q = source.quantified_variables();
          const data::variable_list& par_def = source.variable().parameters();
          const constraint_map& constraints = source.constraints();

          data::rewriter::substitution_type sigma;
          detail::make_constelm_substitution(constraints, sigma);

          qvar_list result;
          for (const auto& it: std::ranges::reverse_view(Q))
          {
            bool is_forall = it.is_forall();
            const data::variable& var = it.variable();
            // Variable of a universal quantifier cannot occur in the disjunctive context
            // Variable of an existential quantifier cannot occur in the conjunctive context
            const std::set<data::variable>& context = is_forall ? m_disj_context : m_conj_context;
            bool none_occurs_in_context = true;

            auto it_pvi = m_target.parameters().begin();
            auto it_def = par_def.begin();
            for (; it_pvi != m_target.parameters().end(); ++it_pvi, ++it_def)
            {
              if (context.find(*it_def) != context.end() && data::search_free_variable(rewr(*it_pvi, sigma), var))
              {
                none_occurs_in_context = false;
                break;
              }
            }

            if (none_occurs_in_context)
            {
              result.push_front(it);
            }
            else
            {
              break;
            }
          }
          return result;
        }
    };

    /// \brief Represents a vertex of the dependency graph.
    class vertex
    {
      protected:
        /// \brief The propositional variable that corresponds to the vertex
        propositional_variable m_variable;

        /// \brief The list of quantified variables that occur in the constraints
        qvar_list m_qvars;

        /// \brief Maps data variables to data expressions. If a parameter is not
        // present in this map, it means that it represents NaC ("not a constant").
        constraint_map m_constraints;

        /// \brief Indicates whether this vertex has been visited at least once.
        bool m_visited = false;

        /// \brief Returns true if the parameter v has been assigned a constant expression.
        /// \param v A parameter of this->variable()
        /// \return True if the data parameter v has been assigned a constant expression.
        bool is_constant(const data::variable& v) const
        {
          auto i = m_constraints.find(v);
          return i != m_constraints.end();
        }

        /// \brief Returns true iff all free variables in e are bound in qvars
        bool bound_in_quantifiers(const qvar_list& qvars, const data::data_expression& e)
        {
          std::set<data::variable> free_vars = data::find_free_variables(e);
          return std::all_of(free_vars.begin(), free_vars.end(), [&](const data::variable& v)
            {
              return std::find_if(qvars.begin(), qvars.end(), [&](const detail::quantified_variable& qvar){ return qvar.variable() == v; }) != qvars.end();
            });
        }

        /// \brief Weaken the constraints so they satisfy
        /// - vars(m_constraints[d]) subset vars(m_qvars); and
        /// - vars(deleted_constraints) intersection vars(m_qvars) = {}
        void fix_constraints(std::vector<data::data_expression> deleted_constraints)
        {
          while (!deleted_constraints.empty())
          {
            std::set<data::variable> vars_deleted;
            for (const data::data_expression& fi: deleted_constraints)
            {
              data::find_free_variables(fi, std::inserter(vars_deleted, vars_deleted.end()));
            }
            deleted_constraints.clear();

            auto del_i = std::find_if(m_qvars.rbegin(), m_qvars.rend(), [&](const detail::quantified_variable& qv)
            {
              return vars_deleted.find(qv.variable()) != vars_deleted.end();
            });
            // Remove quantified variables up to and including del_i
            m_qvars.erase(m_qvars.begin(), del_i.base());

            for (const data::variable& par: m_variable.parameters())
            {
              auto k = m_constraints.find(par);
              if(k == m_constraints.end())
              {
                continue;
              }
              if(!bound_in_quantifiers(m_qvars, k->second))
              {
                deleted_constraints.push_back(k->second);
                m_constraints.erase(k);
              }
            }
          }
        }

      public:
        /// \brief Constructor
        vertex() = default;

        /// \brief Constructor
        /// \param x A propositional variable declaration
        vertex(propositional_variable x)
          : m_variable(x)
        {}

        /// \brief The propositional variable that corresponds to the vertex
        const propositional_variable& variable() const
        {
          return m_variable;
        }

        const qvar_list& quantified_variables() const
        {
          return m_qvars;
        }

        /// \brief Maps data variables to data expressions. If the right hand side is a data
        /// variable, it means that it represents NaC ("not a constant").
        const constraint_map& constraints() const
        {
          return m_constraints;
        }

        /// \brief Returns the indices of the constant parameters of this vertex.
        /// \return The indices of the constant parameters of this vertex.
        std::vector<std::size_t> constant_parameter_indices() const
        {
          std::vector<std::size_t> result;
          std::size_t index = 0;
          for (const data::variable& parameter: m_variable.parameters())
          {
            if (is_constant(parameter))
            {
              result.push_back(index);
            }
            index++;
          }
          return result;
        }

        /// \brief Returns a string representation of the vertex.
        /// \return A string representation of the vertex.
        std::string to_string() const
        {
          std::ostringstream out;
          out << m_variable << "  assertions = ";
          for(const detail::quantified_variable& v: quantified_variables())
          {
            out << v.to_string();
          }
          for (const auto& constraint: m_constraints)
          {
            out << "{" << constraint.first << " := " << constraint.second << "} ";
          }
          return out.str();
        }

        /// \brief Assign new values to the parameters of this vertex, and update the constraints accordingly.
        /// The new values have a number of constraints.
        bool update(const qvar_list& qvars, const data::data_expression_list& e, const constraint_map& e_constraints, const DataRewriter& datar)
        {
          bool changed = false;

          data::variable_list params = m_variable.parameters();
          data::rewriter::substitution_type sigma;
          detail::make_constelm_substitution(e_constraints, sigma);

          if (!m_visited)
          {
            m_visited = true;
            changed = true;

            m_qvars = qvars;
            // Partition expressions in e based on whether their free variables
            // are bound in m_qvars.
            std::vector<data::data_expression> deleted_constraints;
            auto par = params.begin();
            for (auto i = e.begin(); i != e.end(); ++i, ++par)
            {
              data::data_expression e1 = datar(*i, sigma);
              if (bound_in_quantifiers(m_qvars, e1))
              {
                m_constraints[*par] = e1;
              }
              else
              {
                deleted_constraints.push_back(e1);
              }
            }
            fix_constraints(deleted_constraints);
          }
          else
          {
            // Find longest common suffix of qvars
            auto mismatch_it = std::mismatch(m_qvars.rbegin(), m_qvars.rend(), qvars.rbegin(), qvars.rend()).first;
            changed |= mismatch_it != m_qvars.rend();
            // Remove the outer quantifiers, up to and including mismatch_it
            m_qvars.erase(m_qvars.begin(), mismatch_it.base());

            // Find constraints for which f[i] = e[i] and for which all free
            // variables are bound in m_qvars (which may have changed).
            std::vector<data::data_expression> deleted_constraints;
            auto i = e.begin();
            for (auto par = params.begin(); i != e.end(); ++i, ++par)
            {
              auto k = m_constraints.find(*par);
              if(k == m_constraints.end())
              {
                continue;
              }
              const data::data_expression& fi = k->second;
              data::data_expression ei = datar(*i, sigma);
              if (fi != ei || !bound_in_quantifiers(m_qvars, fi))
              {
                changed = true;
                deleted_constraints.push_back(fi);
                deleted_constraints.push_back(ei);
                m_constraints.erase(k);
              }
            }
            fix_constraints(deleted_constraints);
          }
          return changed;
        }
    };

    /// \brief The storage type for vertices
    using vertex_map = std::map<core::identifier_string, vertex>;

    /// \brief The storage type for edges
    using edge_map = std::map<core::identifier_string, std::vector<edge>>;

    /// \brief The vertices of the dependency graph. They are stored in a map, to
    /// support searching for a vertex.
    vertex_map m_vertices;

    /// \brief The edges of the dependency graph. They are stored in a map, to
    /// easily access all out-edges corresponding to a particular vertex.
    edge_map m_edges;

    /// \brief The redundant parameters.
    std::map<core::identifier_string, std::vector<std::size_t> > m_redundant_parameters;

    /// \brief Logs the vertices of the dependency graph.
    std::string print_vertices() const
    {
      std::ostringstream out;
      for (const auto& v: m_vertices)
      {
        out << v.second.to_string() << std::endl;
      }
      return out.str();
    }

    /// \brief Logs the edges of the dependency graph.
    std::string print_edges()
    {
      std::ostringstream out;
      for (const auto& [source, targets]: m_edges)
      {
        for (const edge& e: targets)
        {
          out << e.to_string() << std::endl;
        }
      }
      return out.str();
    }

    std::string print_todo_list(const std::deque<propositional_variable>& todo)
    {
      std::ostringstream out;
      out << "\n<todo list> [";
      for (auto i = todo.begin(); i != todo.end(); ++i)
      {
        if (i != todo.begin())
        {
          out << ", ";
        }
        out << core::pp(i->name());
      }
      out << "]" << std::endl;
      return out.str();
    }

    std::string print_edge_update(const edge& e, const vertex& u, const vertex& v)
    {
      std::ostringstream out;
      out << "\n<updating edge> " << e.to_string() << std::endl;
      out << "  <source vertex       > " << u.to_string() << std::endl;
      out << "  <target vertex before> " << v.to_string() << std::endl;
      return out.str();
    }

    std::string print_condition(const edge& e, const vertex& u, const pbes_expression& value)
    {
      std::ostringstream out;
      data::rewriter::substitution_type sigma;
      detail::make_constelm_substitution(u.constraints(), sigma);
      out << "  <condition           > " << e.condition() << sigma << " to " << value << std::endl;
      return out.str();
    }

    std::string print_evaluation_failure(const edge& e, const vertex& u)
    {
      std::ostringstream out;
      data::rewriter::substitution_type sigma;
      detail::make_constelm_substitution(u.constraints(), sigma);
      out << "\nCould not evaluate condition " << e.condition() << sigma << " to true or false";
      return out.str();
    }

    template <typename E>
    std::list<E> concat(const std::list<E> a, const std::list<E> b)
    {
      std::list<E> result(a);
      result.insert(result.end(), b.begin(), b.end());
      return result;
    }

  public:

    /// \brief Constructor.
    /// \param datar A data rewriter
    /// \param pbesr A PBES rewriter
    pbes_constelm_algorithm(const DataRewriter& datar, const PbesRewriter& pbesr)
      : m_data_rewriter(datar), m_pbes_rewriter(pbesr)
    {}

    /// \brief Returns the parameters that have been removed by the constelm algorithm
    /// \return The removed parameters
    std::map<propositional_variable, std::vector<data::variable> > redundant_parameters() const
    {
      std::map<propositional_variable, std::vector<data::variable> > result;
      for (const std::pair<const core::identifier_string, std::vector<std::size_t>>& red_pair: m_redundant_parameters)
      {
        const vertex& v = m_vertices.find(red_pair.first)->second;
        std::vector<data::variable>& variables = result[v.variable()];
        for (const std::size_t par: red_pair.second)
        {
          typename data::variable_list::iterator k = v.variable().parameters().begin();
          std::advance(k, par);
          variables.push_back(*k);
        }
      }
      return result;
    }

    /// \brief Runs the constelm algorithm
    /// \param p A pbes
    /// \param compute_conditions If true, propagation conditions are computed. Note
    /// that the currently implementation has exponential behavior.
    void run(pbes& p, bool compute_conditions = false, bool check_quantifiers = true)
    {
      m_vertices.clear();
      m_edges.clear();
      m_redundant_parameters.clear();

      // compute the vertices and edges of the dependency graph
      for (pbes_equation& eqn: p.equations())
      {
        core::identifier_string name = eqn.variable().name();
        m_vertices[name] = vertex(eqn.variable());

        // use an edge_condition_traverser to compute the edges
        detail::edge_condition_traverser f;
        f.apply(eqn.formula());

        std::vector<edge>& edges = m_edges[name];
        for (const auto& [Q_X_e, details]: f.result())
        {
          const auto& [Q, X_e] = Q_X_e;
          const auto& [conditions, conj_FV, disj_FV] = details;

          // check options for quantifiers and conditions.
          qvar_list quantifier_list = check_quantifiers ? Q : qvar_list();
          data::data_expression condition = compute_conditions
            ? data::lazy::join_and(conditions.begin(), conditions.end())
            : data::data_expression(data::sort_bool::true_());

          edges.emplace_back(eqn.variable(), quantifier_list, X_e, conj_FV, disj_FV, condition);
        }
      }

      // initialize the todo list of vertices that need to be processed
      propositional_variable_instantiation init = p.initial_state();
      std::deque<propositional_variable> todo;
      const data::data_expression_list& e_init = init.parameters();
      vertex& u_init = m_vertices[init.name()];
      u_init.update(qvar_list(), e_init, constraint_map(), m_data_rewriter);
      todo.push_back(u_init.variable());

      mCRL2log(log::debug) << "\n--- initial vertices ---\n" << print_vertices();
      mCRL2log(log::debug) << "\n--- edges ---\n" << print_edges();

      // propagate constraints over the edges until the todo list is empty
      while (!todo.empty())
      {
        mCRL2log(log::debug) << print_todo_list(todo);
        propositional_variable var = todo.front();

        // remove all occurrences of var from todo
        todo.erase(std::remove(todo.begin(), todo.end(), var), todo.end());

        const vertex& u = m_vertices[var.name()];
        const std::vector<edge>& u_edges = m_edges[var.name()];

        for (const edge& e: u_edges)
        {
          vertex& v = m_vertices[e.target().name()];
          mCRL2log(log::debug) << print_edge_update(e, u, v);

          data::rewriter::substitution_type sigma;
          detail::make_constelm_substitution(u.constraints(), sigma);
          pbes_expression needs_update = m_pbes_rewriter(atermpp::down_cast<pbes_expression>(e.condition()), sigma);
          mCRL2log(log::debug) << print_condition(e, u, needs_update);

          if (!is_false(needs_update) && !is_true(needs_update))
          {
            mCRL2log(log::debug) << print_evaluation_failure(e, u);
          }
          if (!is_false(needs_update))
          {
            bool changed = v.update(
                              concat(e.quantifier_inside_approximation(u, m_data_rewriter), e.quantified_variables()),
                              e.target().parameters(),
                              u.constraints(),
                              m_data_rewriter);
            if (changed)
            {
              todo.push_back(v.variable());
            }
          }
          mCRL2log(log::debug) << "  <target vertex after > " << v.to_string() << "\n";
        }
      }

      mCRL2log(log::debug) << "\n--- final vertices ---\n" << print_vertices();

      // compute the redundant parameters and the redundant equations
      for (const pbes_equation& eqn: p.equations())
      {
        core::identifier_string name = eqn.variable().name();
        const vertex& v = m_vertices[name];
        if (!v.constraints().empty())
        {
          std::vector<std::size_t> r = v.constant_parameter_indices();
          if (!r.empty())
          {
            m_redundant_parameters[name] = r;
          }
        }
      }

      // Apply the constraints to the equations.
      for (pbes_equation& eqn: p.equations())
      {
        core::identifier_string name = eqn.variable().name();
        const vertex& v = m_vertices[name];

        if (!v.constraints().empty())
        {
          data::rewriter::substitution_type sigma;
          detail::make_constelm_substitution(v.constraints(), sigma);
          pbes_expression body = pbes_system::replace_free_variables(eqn.formula(), sigma);
          for (auto i = v.quantified_variables().crbegin(); i != v.quantified_variables().crend(); ++i)
          {
            body = i->make_expr(body);
          }
          eqn = pbes_equation(
                   eqn.symbol(),
                   eqn.variable(),
                   body
                 );
        }
      }

      // remove the redundant parameters
      pbes_system::algorithms::remove_parameters(p, m_redundant_parameters);

      // print the parameters and equation that are removed
      if (mCRL2logEnabled(log::verbose))
      {
        mCRL2log(log::verbose) << "\nremoved the following constant parameters:" << std::endl;
        for (const std::pair<const propositional_variable, std::vector<data::variable>>& i: redundant_parameters())
        {
          for (const data::variable& var: i.second)
          {
            mCRL2log(log::verbose) << "  (" << mcrl2::core::pp(i.first.name()) << ", " << data::pp(var) << ")" << std::endl;
          }
        }
      }
    }
};

/// \brief Apply the constelm algorithm
/// \param p A PBES to which the algorithm is applied
/// \param rewrite_strategy A data rewrite strategy
/// \param rewriter_type A PBES rewriter type
/// \param compute_conditions If true, conditions for the edges of the dependency graph are used N.B. Very inefficient!
/// \param remove_redundant_equations If true, unreachable equations will be removed.
inline
void constelm(pbes& p,
              data::rewrite_strategy rewrite_strategy,
              pbes_rewriter_type rewriter_type,
              bool compute_conditions = false,
              bool remove_redundant_equations = true,
              bool check_quantifiers = true
             )
{
  const bool has_counter_example = pbes_system::detail::has_counter_example_information(p);
  if (has_counter_example)
  {
    mCRL2log(log::warning) << "Warning: the PBES has counter example information, which may not be preserved by constant elimination." << std::endl;
  }
  // data rewriter
  data::rewriter datar(p.data(), rewrite_strategy);

  // pbes rewriter
  switch (rewriter_type)
  {
    case pbes_rewriter_type::simplify:
    {
      using pbes_rewriter = simplify_data_rewriter<data::rewriter>;
      pbes_rewriter pbesr(datar);
      pbes_constelm_algorithm<data::rewriter, pbes_rewriter> algorithm(datar, pbesr);
      algorithm.run(p, compute_conditions, check_quantifiers);
      if (remove_redundant_equations)
      {
        std::vector<propositional_variable> V = algorithms::remove_unreachable_variables(p);
        mCRL2log(log::verbose) << algorithms::print_removed_equations(V);
      }
      break;
    }
    case pbes_rewriter_type::quantifier_all:
    case pbes_rewriter_type::quantifier_finite:
    {
      const enumerate_quantifiers_mode enum_mode = (rewriter_type == pbes_rewriter_type::quantifier_all?
                                                         expand_infinite_sorts_and_use_data_rewriter:
                                                         expand_finite_sorts);
      enumerate_quantifiers_rewriter pbesr(datar, p.data(), enum_mode);
      pbes_constelm_algorithm<data::rewriter, enumerate_quantifiers_rewriter> algorithm(datar, pbesr);
      algorithm.run(p, compute_conditions, check_quantifiers);
      if (remove_redundant_equations)
      {
        std::vector<propositional_variable> V = algorithms::remove_unreachable_variables(p);
        mCRL2log(log::verbose) << algorithms::print_removed_equations(V);
      }
      break;
    }
    default:
    { }
  }
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_CONSTELM_H
