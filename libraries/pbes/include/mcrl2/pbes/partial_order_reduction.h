// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/partial_order_reduction.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PARTIAL_ORDER_REDUCTION_H
#define MCRL2_PBES_PARTIAL_ORDER_REDUCTION_H

#include <iomanip>
#include <boost/dynamic_bitset.hpp>
#include "mcrl2/data/rewriters/one_point_rule_rewriter.h"
#include "mcrl2/data/rewriters/quantifiers_inside_rewriter.h"
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/data/substitutions/maintain_variables_in_rhs.h"
#include "mcrl2/pbes/replace_capture_avoiding_with_an_identifier_generator.h"
#include "mcrl2/pbes/rewriters/enumerate_quantifiers_rewriter.h"
#include "mcrl2/pbes/pbes_equation_index.h"
#include "mcrl2/pbes/unify_parameters.h"
#include "mcrl2/smt/solver.h"
#include "mcrl2/utilities/skip.h"

namespace mcrl2::pbes_system
{

namespace detail {

inline
data::data_expression make_and(const data::data_expression& x1, const data::data_expression& x2, const data::data_expression& x3)
{
  return data::and_(x1, data::and_(x2, x3));
}

inline
data::data_expression equal_to(const data::data_expression_list& x, const data::data_expression_list& y)
{
  data::data_expression result = data::sort_bool::true_();
  auto xi = x.begin();
  auto yi = y.begin();
  for (; xi != x.end(); ++xi, ++yi)
  {
    result = data::lazy::and_(result, data::equal_to(*xi, *yi));
  }
  return result;
}

} // namespace detail

using summand_set = boost::dynamic_bitset<>;

inline
std::string print_summand_set(const summand_set& s)
{
  bool first = true;
  std::ostringstream buf;
  buf << "{ ";
  for (std::size_t k = s.find_first(); k != summand_set::npos; k = s.find_next(k))
  {
    buf << (first ? "" : ", ") << k;
    first = false;
  }
  buf << " }";
  return buf.str();
}

struct summand_class
{
  data::variable_list e;
  data::data_expression f;
  data::data_expression_list g;
  /// \brief Encodes the dependency relation belonging to this summand_class
  /// \details nxt[i] contains j iff X_i --this--> X_j
  std::vector<std::set<std::size_t>> nxt;
  summand_set NES;
  summand_set DNA;
  summand_set DNS;
  summand_set DNL;
  bool is_deterministic = false;

  summand_class() = default;

  // n is the number of PBES equations
  summand_class(data::variable_list  e_, data::data_expression f_, data::data_expression_list g_, std::size_t n)
   : e(std::move(e_)), f(std::move(f_)), g(std::move(g_))
  {
    nxt.resize(n);
  }

  void set_num_summands(const std::size_t N)
  {
    NES.resize(N);
    DNA.resize(N);
    DNS.resize(N);
    DNL.resize(N);
  }

  // returns X_i -k->
  [[nodiscard]]
  bool depends(std::size_t i) const
  {
    return !nxt[i].empty();
  }

  // returns X_i -k-> j
  [[nodiscard]]
  bool depends(std::size_t i, std::size_t j) const
  {
    using utilities::detail::contains;

    return contains(nxt[i], j);
  }

  void print(std::ostream& out, const std::set<std::size_t>& s, const std::size_t N) const
  {
    using utilities::detail::contains;

    for (std::size_t i = 0; i < N; i++)
    {
      out << (contains(s, i) ? '1' : '0');
    }
  }

  void print(std::ostream& out) const
  {
    out << "deterministic = " << std::boolalpha << is_deterministic << std::endl;
    out << "NES = " << print_summand_set(NES) << std::endl;
    out << "DNA = " << print_summand_set(DNA) << std::endl;
    out << "DNS = " << print_summand_set(DNS) << std::endl;
    out << "DNL = " << print_summand_set(DNL) << std::endl;
  }
};

// The part of a summand used for determining equivalence classes
struct summand_equivalence_key
{
  data::variable_list e;
  data::data_expression f;
  data::data_expression_list g;

  explicit summand_equivalence_key(data::variable_list e_, data::data_expression f_, data::data_expression_list g_)
  : e(std::move(e_)), f(std::move(f_)), g(std::move(g_))
  {}

  explicit summand_equivalence_key(const summand_class& summand)
   : e(summand.e), f(summand.f), g(summand.g)
  {}

  explicit summand_equivalence_key(const srf_summand& summand)
   : e(summand.parameters()), f(summand.condition()), g(summand.variable().parameters())
  {}

  bool operator<(const summand_equivalence_key& other) const
  {
    return std::tie(e, f, g) < std::tie(other.e, other.f, other.g);
  }

  bool operator==(const summand_equivalence_key& other) const
  {
    return std::tie(e, f, g) == std::tie(other.e, other.f, other.g);
  }
};

} // namespace mcrl2::pbes_system

namespace std {

/// \brief specialization of the standard std::hash function.
template<>
struct hash<mcrl2::pbes_system::summand_equivalence_key>
{
  std::size_t operator()(const mcrl2::pbes_system::summand_equivalence_key& x) const
  {
    std::size_t seed = std::hash<atermpp::aterm>()(x.f);
    if (!x.e.empty())
    {
      seed = std::hash<atermpp::aterm>()(x.e) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    if (!x.g.empty())
    {
      seed = std::hash<atermpp::aterm>()(x.g) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

} // namespace std

namespace mcrl2::pbes_system
{

enum tribool
{
  no, maybe, yes
};
static inline bool operator&&(tribool a, tribool b)
{
  return a == yes || b == yes || (a == maybe && b == maybe);
}
// Short-circuit version of the operator && for tribools
// The second function will be told whether a 'yes' answer is required to satisfy
// the expression
static inline bool operator&&(const std::function<tribool()>& a, const std::function<tribool(bool)>& b)
{
  tribool a_result = a();
  if(a_result == yes)
  {
    return true;
  }
  return a_result && b(a_result == no);
}

struct pbespor_options
{
  // if true use alternative A3 for maybe clauses in accordance conditions
  bool compute_weak_conditions = false;
  // if true, compute which summands are deterministic
  bool compute_determinism = true;
  // if true, compute triangle accordance
  bool compute_triangle_accordance = true;
  // if true, compute left accordance
  bool compute_left_accordance = true;
  // if true, compute NES relation
  bool compute_NES = true;
  // if true, perform static analysis and reduce state space
  bool reduction = true;
  // if true, apply condition L
  bool use_condition_L = true;

  data::rewrite_strategy rewrite_strategy = data::rewrite_strategy::jitty;

  bool use_smt_solver = false;
  std::chrono::milliseconds smt_timeout = std::chrono::milliseconds::zero();
};

class partial_order_reduction_algorithm
{
  protected:
    using enumerator_element = data::enumerator_list_element_with_substitution<data::data_expression>;

    struct invis_pair
    {
      summand_set Twork;
      summand_set Ts;
      invis_pair(summand_set Twork_, summand_set Ts_)
       : Twork(std::move(Twork_)), Ts(std::move(Ts_))
      {}

      bool operator<(const invis_pair& other) const
      {
        return std::tie(Twork, Ts) < std::tie(other.Twork, other.Ts);
      }
    };

    struct parameter_info
    {
      std::set<std::size_t> Ts; // test set
      std::set<std::size_t> Ws; // write set
      std::set<std::size_t> Rs; // read set
      std::set<std::size_t> Vs; // variable set
    };

    data::rewriter m_rewr;
    data::enumerator_identifier_generator m_id_generator;
    data::enumerator_algorithm<> m_enumerator;
    srf_pbes m_pbes;
    pbes_equation_index m_equation_index;
    data::mutable_indexed_substitution<> m_sigma;
    std::size_t m_largest_equation_size = 0;

    // the parameters of the PBES equations
    std::vector<data::variable> m_parameters;

    // maps parameters to their corresponding index
    std::map<data::variable, std::size_t> m_parameter_positions;

    // maps summands to the index of the corresponding summand class
    std::unordered_map<summand_equivalence_key, std::size_t> m_summand_index;

    // X_j \in nxt_k(X_i) <=> j \in m_summand_classes[k].nxt[i]
    // (X_i |- k -> X_j) <=> m_summand_classes[k].NES[i][j]
    std::vector<summand_class> m_summand_classes;

    summand_set m_invis; // invisible summand classes
    summand_set m_vis;   // visible summand classes

    // One NES for every predicate variable X_i that can be used for summand
    // class k when !depends(i,k).
    std::vector<summand_set> m_dependency_nes;

    std::chrono::high_resolution_clock::duration m_static_analysis_duration{};
    std::chrono::high_resolution_clock::duration m_exploration_duration{};

    std::unique_ptr<smt::smt_solver> m_solver;

    pbespor_options m_options;

    class summand_relations_data
    {
    private:
      partial_order_reduction_algorithm& parent;
      bool compute_weak_conditions;
      data::set_identifier_generator id_gen;

      data::variable_list qvars1_k;
      data::data_expression condition1_k;
      data::data_expression_list updates1_k;

      data::variable_list qvars1_k1;
      data::data_expression condition1_k1;
      data::data_expression_list updates1_k1;

      data::variable_list qvars2_k;
      data::data_expression condition2_k;
      data::data_expression_list updates2_k;

      data::variable_list qvars2_k1;
      data::data_expression condition2_k1;
      data::data_expression_list updates2_k1;

      data::mutable_indexed_substitution<> sigma_k;
      data::mutable_indexed_substitution<> sigma_k1;

      data::variable_list combined_quantified_vars;

      // Depending on whether the weak (A3) or strong condition (A4) is used, wrap the consequent in
      // an existential quantifier
      data::data_expression make_exists_if_strong(const data::variable_list& vars, const data::data_expression& body)
      {
        return compute_weak_conditions ? body : optimized_exists(vars, body);
      }

      data::data_expression left_accords_antecedent()
      {
        return data::sort_bool::and_(condition1_k1, data::replace_variables_capture_avoiding(condition1_k, sigma_k1, id_gen));
      }

      data::data_expression left_accords_consequent()
      {
        data::data_expression parameters_equal = detail::equal_to(data::replace_variables_capture_avoiding(updates2_k, sigma_k1, id_gen),
                                                                  data::replace_variables_capture_avoiding(updates2_k1, sigma_k, id_gen));
        data::data_expression body = detail::make_and(
          condition2_k,
          data::replace_variables_capture_avoiding(condition2_k1, sigma_k, id_gen),
          parameters_equal
        );
        return make_exists_if_strong(qvars2_k + qvars2_k1, body);
      }

      data::data_expression coenabled_antecedent()
      {
        return data::sort_bool::and_(condition1_k, condition1_k1);
      }

      data::data_expression square_accords_consequent()
      {
        data::data_expression parameters_equal = detail::equal_to(data::replace_variables_capture_avoiding(updates2_k, sigma_k1, id_gen),
                                                                  data::replace_variables_capture_avoiding(updates2_k1, sigma_k, id_gen));
        data::data_expression body = detail::make_and(
          data::replace_variables_capture_avoiding(condition2_k, sigma_k1, id_gen),
          data::replace_variables_capture_avoiding(condition2_k1, sigma_k, id_gen),
          parameters_equal
        );
        return make_exists_if_strong(qvars2_k + qvars2_k1, body);
      }

      data::data_expression triangle_accords_consequent_weak()
      {
        data::data_expression parameters_equal = detail::equal_to(updates1_k1, data::replace_variables_capture_avoiding(updates1_k1, sigma_k, id_gen));
        return data::sort_bool::and_(
          data::replace_variables_capture_avoiding(condition1_k1, sigma_k, id_gen),
          parameters_equal
        );
      }

      data::data_expression triangle_accords_consequent()
      {
        data::data_expression parameters_equal = detail::equal_to(updates2_k1, data::replace_variables_capture_avoiding(updates2_k1, sigma_k, id_gen));
        data::data_expression body = data::sort_bool::and_(
          data::replace_variables_capture_avoiding(condition2_k1, sigma_k, id_gen),
          parameters_equal
        );
        return make_exists_if_strong(qvars2_k1, body);
      }

      tribool accords_data(bool affect_set, bool needs_yes,
                           const std::function<data::data_expression()>& make_antecedent,
                           const std::function<data::data_expression()>& make_consequent)
      {
        // Check whether the maybe clause is satisfied by affect sets and it is sufficient to return maybe
        if (affect_set && !needs_yes)
        {
          return maybe;
        }

        data::data_expression antecedent = make_antecedent();
        data::data_expression yes_condition = data::optimized_forall(combined_quantified_vars, data::sort_bool::not_(antecedent));
        if (parent.is_true(yes_condition))
        {
          return yes;
        }
        if (needs_yes)
        {
          // we were not able to return yes, now it doesn't matter what we return
          return no;
        }

        data::data_expression consequent = make_consequent();
        data::data_expression condition = data::optimized_forall(combined_quantified_vars, data::sort_bool::implies(antecedent, consequent));

        return parent.is_true(condition) ? maybe : no;
      }

    public:
      summand_relations_data(partial_order_reduction_algorithm& p, const std::size_t k, const std::size_t k1)
      : parent(p)
      , compute_weak_conditions(p.m_options.compute_weak_conditions)
      {
        const summand_class& summand_k = parent.m_summand_classes[k];
        const summand_class& summand_k1 = parent.m_summand_classes[k1];

        const data::variable_list& parameters = parent.m_pbes.equations()[0].variable().parameters();
        for(const data::variable& v: parameters)
        {
          id_gen.add_identifier(v.name());
        }

        // For both summands, create a copy with fresh variables
        // These will be used when constructing accordance and NES conditions
        summand_equivalence_key new1_k = rename_duplicate_variables(id_gen, summand_equivalence_key(summand_k));
        summand_equivalence_key new1_k1 = rename_duplicate_variables(id_gen, summand_equivalence_key(summand_k1));
        qvars1_k  = new1_k.e;  condition1_k  = new1_k.f;  updates1_k  = new1_k.g;
        qvars1_k1 = new1_k1.e; condition1_k1 = new1_k1.f; updates1_k1 = new1_k1.g;
        data::add_assignments(sigma_k, parameters, updates1_k);
        data::add_assignments(sigma_k1, parameters, updates1_k1);

        if (!compute_weak_conditions)
        {
          // When using the stronger condition A4, create another fresh copy
          summand_equivalence_key new2_k = rename_duplicate_variables(id_gen, summand_equivalence_key(summand_k));
          summand_equivalence_key new2_k1 = rename_duplicate_variables(id_gen, summand_equivalence_key(summand_k1));
          qvars2_k  = new2_k.e;  condition2_k  = new2_k.f;  updates2_k  = new2_k.g;
          qvars2_k1 = new2_k1.e; condition2_k1 = new2_k1.f; updates2_k1 = new2_k1.g;
        }
        else
        {
          // In the weak case (A3), the two copies are identical
          qvars2_k  = qvars1_k;  condition2_k  = condition1_k;  updates2_k  = updates1_k;
          qvars2_k1 = qvars1_k1; condition2_k1 = condition1_k1; updates2_k1 = updates1_k1;
        }

        combined_quantified_vars = parameters + qvars1_k + qvars1_k1;
      }

      bool can_enable()
      {
        data::data_expression cannot_enable = data::optimized_forall(
          combined_quantified_vars,
          data::sort_bool::not_(
            detail::make_and(
              data::sort_bool::not_(condition1_k),
              condition1_k1,
              data::replace_variables_capture_avoiding(condition1_k, sigma_k1, id_gen)
            )
          )
        );

        // The condition is constructed in a negated way, so the approximation of the decision
        // procedure works the right way. Note that the result of this function is negated as well.
        return !parent.is_true(cannot_enable);
      }

      tribool left_accords_data(bool affect_set, bool needs_yes)
      {
        return accords_data(affect_set, needs_yes,
                            [this](){ return left_accords_antecedent(); },
                            [this](){ return left_accords_consequent(); });
      }

      tribool square_accords_data(bool affect_set, bool needs_yes)
      {
        return accords_data(affect_set, needs_yes,
                            [this](){ return coenabled_antecedent(); },
                            [this](){ return square_accords_consequent(); });
      }

      tribool triangle_accords_data(bool affect_set, bool needs_yes)
      {
        return accords_data(affect_set, needs_yes,
                            [this](){ return coenabled_antecedent(); },
                            [this](){ return triangle_accords_consequent(); });
      }
    };

    std::size_t summand_index(const srf_summand& summand) const
    {
      auto i = m_summand_index.find(summand_equivalence_key(summand));
      assert(i != m_summand_index.end());
      return i->second;
    }

    std::size_t parameter_position(const data::variable& v) const
    {
      auto i = m_parameter_positions.find(v);
      return i->second;
    }

    const summand_set& DNA(std::size_t k) const
    {
      return m_summand_classes[k].DNA;
    }

    summand_set& DNA(std::size_t k)
    {
      return m_summand_classes[k].DNA;
    }

    const summand_set& DNS(std::size_t k) const
    {
      return m_summand_classes[k].DNS;
    }

    summand_set& DNS(std::size_t k)
    {
      return m_summand_classes[k].DNS;
    }

    const summand_set& DNL(std::size_t k) const
    {
      return m_summand_classes[k].DNL;
    }

    summand_set& DNL(std::size_t k)
    {
      return m_summand_classes[k].DNL;
    }

    const summand_set& NES(std::size_t k) const
    {
      return m_summand_classes[k].NES;
    }

    summand_set& NES(std::size_t k)
    {
      return m_summand_classes[k].NES;
    }

    summand_set en(const propositional_variable_instantiation& X_e)
    {
      std::size_t N = m_summand_classes.size();

      summand_set result(N);
      std::size_t i = m_equation_index.index(X_e.name());
      const data::variable_list& d = m_pbes.equations()[i].variable().parameters();
      const data::data_expression_list& e = X_e.parameters();
      for (std::size_t k = 0; k < N; k++)
      {
        if (!depends(i, k))
        {
          continue;
        }
        const summand_class& summand_k = m_summand_classes[k];
        const data::variable_list& e_k = summand_k.e;
        const data::data_expression& f_k = summand_k.f;

        // Add assignments to d in every iteration, becuase they might have been
        // overwritten if a variable in e_k coincides with a parameter in d.
        add_assignments(m_sigma, d, e);
        m_enumerator.enumerate(enumerator_element(e_k, f_k),
                               m_sigma,
                               [&](const enumerator_element&) {
                                 result.set(k);
                                 return true;
                               },
                               data::is_false,
                               data::is_true
        );
        remove_assignments(m_sigma, e_k);
      }
      remove_assignments(m_sigma, d);
      return result;
    }

    summand_set invis(const summand_set& K)
    {
      return m_invis & K;
    }

    // Choose a NES according to whether k is present in X
    const summand_set& choose_minimal_NES(std::size_t k, const propositional_variable_instantiation& X_e) const
    {
      if(!depends(m_equation_index.index(X_e.name()), k))
      {
        return m_dependency_nes[m_equation_index.index(X_e.name())];
      }
      else
      {
        return m_summand_classes[k].NES;
      }
    }

    const summand_set& DNX(std::size_t k,
                                   const summand_set& Twork,
                                   const summand_set& Ts,
                                   const summand_set& en_X_e
                                  ) const
    {
      using utilities::detail::set_difference;
      using utilities::detail::set_intersection;
      using utilities::detail::set_union;

      const summand_class& summand_k = m_summand_classes[k];
      if (!summand_k.is_deterministic)
      {
        return DNL(k);
      }
      if (!m_options.compute_left_accordance)
      {
        return DNA(k);
      }

      summand_set Twork_Ts = Twork | Ts;

      summand_set T1 = Twork_Ts | en_X_e;
      summand_set T2 = Twork_Ts & en_X_e;

      auto h = [&](const summand_set& A)
      {
        return (A - T1).count() + m_largest_equation_size * (A - T2).count();
      };

      return h(DNS(k)) <= h(DNL(k)) ? DNS(k) : DNL(k);
    }

    summand_set stubborn_set(const propositional_variable_instantiation& X_e, const summand_set& en_X_e)
    {
      using utilities::detail::contains;
      using utilities::detail::has_empty_intersection;
      using utilities::detail::set_difference;
      using utilities::detail::set_includes;
      using utilities::detail::set_intersection;
      using utilities::detail::set_union;

      std::size_t N = m_summand_classes.size();

      struct compare_invis_pair
      {
        // Prevent allocations and destructions by using temp_set
        mutable summand_set temp_set;
        const summand_set& m_en_X_e;

        explicit compare_invis_pair(const summand_set& en)
        : m_en_X_e(en)
        {}

        std::size_t size(const invis_pair& p) const
        {
          temp_set = p.Twork;
          temp_set |= p.Ts;
          temp_set &= m_en_X_e;
          return temp_set.count();
        }

        bool operator()(const invis_pair& x, const invis_pair& y) const
        {
          std::size_t sizex = size(x);
          std::size_t sizey = size(y);
          return std::tie(sizex, x) < std::tie(sizey, y);
        };
      };

      // C will always be sorted according to the size of each element.
      std::set<invis_pair, compare_invis_pair> C{compare_invis_pair(en_X_e)};
      auto invis_en_X_e = invis(en_X_e);
      for (std::size_t k = invis_en_X_e.find_first(); k != summand_set::npos; k = invis_en_X_e.find_next(k))
      {
        // If DNL is not computed, only attempt to construct stubborn set for deterministic
        // summands.
        if(m_summand_classes[k].is_deterministic || m_options.compute_left_accordance)
        {
          invis_pair pair_k{summand_set(N), summand_set(N)};
          pair_k.Twork.set(k);
          C.insert(pair_k);
        }
      }
      if (C.empty())
      {
        return en_X_e;
      }

      while (!C.empty())
      {
        // The smallest element is the first element in the set
        auto p = C.extract(C.begin());
        auto& Twork = p.value().Twork;
        auto& Ts = p.value().Ts;

        if (Twork.none())
        {
          summand_set T = Ts & en_X_e;
          for (std::size_t k = T.find_first(); k != summand_set::npos; k = T.find_next(k))
          {
            if (DNA(k).is_subset_of(Ts))
            {
              return Ts;
            }
          }
          std::size_t k = T.find_first(); // TODO: choose k according to D2t
          Twork |= (DNA(k) - Ts);
        }
        else
        {
          std::size_t k = Twork.find_first();
          Twork.reset(k);
          Ts.set(k);
          if (en_X_e.test(k))
          {
            if (m_vis.test(k))
            {
              // To satisfy condition P, we fully expand in case a visible
              // transition occurs in the stubborn set.
              // Thus, do not add this candidate set back to C.
              // Effectively, we implement a stronger version of V:
              // If r(s) contains an enabled visible transition,
              // then it contains all transitions.
              continue;
            }
            auto& DNS_or_DNL = DNX(k, Twork, Ts, en_X_e);
            Twork |= (DNS_or_DNL - Ts);
          }
          else
          {
            auto& NES = choose_minimal_NES(k, X_e);
            Twork |= (NES - Ts);
          }
        }

        C.insert(std::move(p));
      }
      return en_X_e;
    }

    std::set<propositional_variable_instantiation> succ(const propositional_variable_instantiation& X_e, const summand_set& K)
    {
      const auto& d = m_parameters;
      const auto& e = X_e.parameters();

      std::set<propositional_variable_instantiation> result;
      std::size_t i = m_equation_index.index(X_e.name());
      for (std::size_t k = K.find_first(); k != summand_set::npos; k = K.find_next(k))
      {
        const summand_class& summand_k = m_summand_classes[k];
        const data::variable_list& e_k = summand_k.e;
        const data::data_expression& f_k = summand_k.f;
        const data::data_expression_list& g_k = summand_k.g;
        const auto& J = summand_k.nxt[i];

        // Add assignments for parameters during every iteration,
        // because they might have been removed on the previous one if
        // a parameter coincides with a quantified variable.
        data::add_assignments(m_sigma, d, e);
        data::remove_assignments(m_sigma, e_k);
        m_enumerator.enumerate(enumerator_element(e_k, f_k),
                               m_sigma,
                               [&](const enumerator_element& p) {
                                 p.add_assignments(e_k, m_sigma, m_rewr);
                                 data::data_expression_list g(g_k.begin(), g_k.end(), [&](const data::data_expression& x) { return m_rewr(x, m_sigma); });
                                 for (std::size_t j: J)
                                 {
                                   const core::identifier_string& X_j = m_pbes.equations()[j].variable().name();
                                   result.insert(propositional_variable_instantiation(X_j, g));
                                 }
                                 return false;
                               },
                               data::is_false
        );
        data::remove_assignments(m_sigma, e_k);
      }
      return result;
    }

    void compute_nxt()
    {
      std::size_t n = m_pbes.equations().size();
      for (std::size_t i = 0; i < n; i++)
      {
        const srf_equation& eqn = m_pbes.equations()[i];
        for (const srf_summand& summand: eqn.summands())
        {
          std::size_t j = m_equation_index.index(summand.variable().name());
          std::size_t k = summand_index(summand);
          m_summand_classes[k].nxt[i].insert(j);
        }
      }
    }

    // returns the indices of the parameters that occur freely in x
    std::set<std::size_t> FV(const pbes_expression& x) const
    {
      std::set<std::size_t> result;
      for (const data::variable& v: find_free_variables(x))
      {
        result.insert(parameter_position(v));
      }
      return result;
    }

    bool is_true(data::data_expression expr)
    {
      if(m_solver != nullptr)
      {
        bool negate = false;
        if(data::is_forall(expr))
        {
          negate = true;
          const data::forall& f = atermpp::down_cast<data::forall>(expr);
          expr = data::optimized_exists(f.variables(), data::sort_bool::not_(f.body()));
        }
        switch(m_solver->solve(data::variable_list(), expr, m_options.smt_timeout))
        {
          case smt::answer::SAT: return negate ^ true;
          case smt::answer::UNSAT: return negate ^ false;
          case smt::answer::UNKNOWN: return false;
        }
      }
      else
      {
        data::data_expression result = m_rewr(data::one_point_rule_rewrite((m_rewr(expr))));
        if (result != data::sort_bool::true_() && result != data::sort_bool::false_())
        {
          mCRL2log(log::verbose) << "Cannot rewrite " << result << " any further" << std::endl;
        }
        return result == data::sort_bool::true_();
      }
      // This code is unreachable, but necessary to silence a compiler warning
      return false;
    }

    /// \brief Return true iff k1 can never happen after k happens, as deduced from
    /// predicate dependencies.
    bool dependency_permanently_disables(const std::size_t k, const std::size_t k1) const
    {
      std::size_t N = m_summand_classes.size();
      const summand_class& summand_k = m_summand_classes[k];
      std::set<std::size_t> reachable_after_k;

      // Check to which equations k can lead
      for (std::size_t i = 0; i < m_pbes.equations().size(); i++)
      {
        reachable_after_k.insert(summand_k.nxt[i].begin(), summand_k.nxt[i].end());
      }

      // Explore the rest of the dependency relation
      std::list<std::size_t> todo(reachable_after_k.begin(), reachable_after_k.end());
      while (!todo.empty())
      {
        std::size_t i = todo.front();
        todo.pop_front();
        for (std::size_t k2 = 0; k2 < N; k2++)
        {
          for (const std::size_t j: m_summand_classes[k2].nxt[i])
          {
            if (reachable_after_k.count(j) == 0)
            {
              todo.push_back(j);
              reachable_after_k.insert(j);
            }
          }
        }
      }

      return !std::any_of(reachable_after_k.begin(), reachable_after_k.end(),
        [&](const std::size_t i) { return depends(i, k1); });
    }

    void compute_dependency_NES()
    {
      using utilities::detail::set_union;
      using utilities::detail::has_empty_intersection;

      std::size_t n = m_pbes.equations().size();
      std::size_t N = m_summand_classes.size();

      for (std::size_t i = 0; i < n; i++)
      {
        m_dependency_nes[i].resize(N);
        for (std::size_t k = 0; k < N; k++)
        {
          const std::set<std::size_t>& J = m_summand_classes[k].nxt[i];
          if (J.size() > 1 || (J.size() == 1 && *J.begin() != i))
          {
            m_dependency_nes[i].set(k);
          }
        }
      }
    }

    // returns X_i |--k--> X_j
    bool depends(std::size_t i, std::size_t k, std::size_t j) const
    {
      return m_summand_classes[k].depends(i, j);
    };

    // TODO: precompute this function
    // returns X_i |--k-->
    bool depends(std::size_t i, std::size_t k) const
    {
      std::size_t n = m_pbes.equations().size();
      for (std::size_t j = 0; j < n; j++)
      {
        if (depends(i, k, j))
        {
          return true;
        }
      }
      return false;
    };

    static summand_equivalence_key rename_duplicate_variables(data::set_identifier_generator& id_gen, const summand_equivalence_key& summ)
    {
      std::vector<data::variable> new_variables;
      data::maintain_variables_in_rhs< data::mutable_map_substitution<> > sigma;
      for (const data::variable& var: summ.e)
      {
        core::identifier_string new_name = id_gen(var.name());
        if (new_name != var.name())
        {
          sigma[var] = data::variable(new_name, var.sort());
        }
        new_variables.emplace_back(new_name, var.sort());
      }

      auto replace_vars = [&](const data::data_expression& e)
      {
        return data::replace_variables_capture_avoiding_with_an_identifier_generator(e, sigma, id_gen);
      };

      return summand_equivalence_key(
        data::variable_list(new_variables.begin(), new_variables.end()),
        replace_vars(summ.f),
        data::data_expression_list(summ.g.begin(), summ.g.end(), replace_vars)
      );
    }

    tribool left_accords_equations(std::size_t k, std::size_t k1) const
    {
      std::size_t n = m_pbes.equations().size();
      tribool result = yes;

      for (std::size_t i = 0; i < n; i++)
      {
        for (std::size_t i1 = 0; i1 < n; i1++)
        {
          bool X_k1_X1 = depends(i, k1, i1);
          for (std::size_t i_prime = 0; i_prime < n; i_prime++)
          {
            bool X1_k_Xprime = depends(i1, k, i_prime);
            if (X_k1_X1 && X1_k_Xprime)
            {
              result = maybe;
              bool found = false;
              for (std::size_t i2 = 0; i2 < n; i2++)
              {
                bool X_k_X2 = depends(i, k, i2);
                bool X2_k1_Xprime = depends(i2, k1, i_prime);
                if (X_k_X2 && X2_k1_Xprime)
                {
                  found = true;
                }
              }
              if (!found)
              {
                return no;
              }
            }
          }
        }
      }
      return result;
    }

    tribool square_accords_equations(std::size_t k, std::size_t k1) const
    {
      std::size_t n = m_pbes.equations().size();
      tribool result = yes;

      for (std::size_t i = 0; i < n; i++)
      {
        for (std::size_t i1 = 0; i1 < n; i1++)
        {
          bool X_k1_X1 = depends(i, k1, i1);
          for (std::size_t i2 = 0; i2 < n; i2++)
          {
            bool X_k_X2 = depends(i, k, i2);
            if (X_k1_X1 && X_k_X2)
            {
              result = maybe;
              bool found = false;
              for (std::size_t i_prime = 0; i_prime < n; i_prime++)
              {
                bool X1_k_Xprime = depends(i1, k, i_prime);
                bool X2_k1_Xprime = depends(i2, k1, i_prime);
                if (X1_k_Xprime && X2_k1_Xprime)
                {
                  found = true;
                  break;
                }
              }
              if (!found)
              {
                return no;
              }
            }
          }
        }
      }
      return result;
    }

    tribool triangle_accords_equations(std::size_t k, std::size_t k1) const
    {
      std::size_t n = m_pbes.equations().size();
      tribool result = yes;

      for (std::size_t i = 0; i < n; i++)
      {
        for (std::size_t i1 = 0; i1 < n; i1++)
        {
          bool X_k1_X1 = depends(i, k1, i1);
          for (std::size_t i2 = 0; i2 < n; i2++)
          {
            bool X_k_X2 = depends(i, k, i2);
            bool X2_k1_X1 = depends(i2, k1, i1);
            if (X_k1_X1 && X_k_X2)
            {
              result = maybe;
              if(!X2_k1_X1)
              {
                return no;
              }
            }
          }
        }
      }
      return result;
    }

    void compute_DNA_DNL_NES(const std::vector<parameter_info>& info)
    {
      using utilities::detail::has_empty_intersection;
      using utilities::detail::set_includes;
      using utilities::detail::set_intersection;
      using utilities::detail::set_union;

      std::size_t N = m_summand_classes.size();

      auto Rs = [&](const std::size_t k) { return info[k].Rs; };
      auto Ts = [&](const std::size_t k) { return info[k].Ts; };
      auto Vs = [&](const std::size_t k) { return info[k].Vs; };
      auto Ws = [&](const std::size_t k) { return info[k].Ws; };

      for (std::size_t k = 0; k < N; k++)
      {
        mCRL2log(log::verbose) << std::setw(3) << k << " = ";
        for (std::size_t k1 = 0; k1 < N; k1++)
        {
          if (k == k1)
          {
            mCRL2log(log::verbose) << ". ";
            continue;
          }
          bool DNL_DNS_affect_sets = has_empty_intersection(set_intersection(Vs(k), Vs(k1)), set_union(Ws(k), Ws(k1)));
          bool DNT_affect_sets = has_empty_intersection(Ws(k), Rs(k1)) && has_empty_intersection(Ws(k), Ts(k1)) && set_includes(Ws(k1), Ws(k));

          summand_relations_data summand_data(*this, k, k1);
          // Use lambda lifting for short-circuiting the && operator on tribools
          bool left_accords     = m_options.compute_left_accordance &&
                                  ([&]{ return left_accords_equations(k, k1); } &&
                                   [&](bool needs_yes) { return summand_data.left_accords_data(DNL_DNS_affect_sets, needs_yes); });
          // The DNS relation is symmetric
          bool square_accords   = (k1 < k && !DNS(k1).test(k)) ||
                                  (k1 > k &&
                                      ([&]{ return square_accords_equations(k, k1); } &&
                                       [&](bool needs_yes) { return summand_data.square_accords_data(DNL_DNS_affect_sets, needs_yes); }));
          bool accords          = square_accords ||
                                  (m_options.compute_triangle_accordance && ([&]{ return triangle_accords_equations(k, k1); } &&
                                                          [&](bool needs_yes) { return summand_data.triangle_accords_data(DNT_affect_sets, needs_yes); }));
          bool can_enable       = !m_options.compute_NES ||
                                  (!dependency_permanently_disables(k1, k) && !has_empty_intersection(Ts(k), Ws(k1)) && summand_data.can_enable());

          if (!left_accords)
          {
            DNL(k).set(k1);
          }
          if (!square_accords)
          {
            DNS(k).set(k1);
          }
          if (!accords)
          {
            DNA(k).set(k1);
            mCRL2log(log::verbose) << "- ";
          }
          else
          {
            mCRL2log(log::verbose) << (DNL_DNS_affect_sets ? ": " : "+ ");
          }
          mCRL2log(log::verbose) << std::flush;
          if (can_enable)
          {
            NES(k).set(k1);
          }
        }
        mCRL2log(log::verbose) << "\n";
      }
    }

    void compute_NES_DNA_DNL()
    {
      if (!m_options.reduction)
      {
        return;
      }

      using utilities::detail::set_union;

      std::size_t N = m_summand_classes.size();
      std::vector<parameter_info> info(N);
      const std::vector<data::variable>& d = m_parameters;

      auto compute_parameter_info = [&](summand_class& summand, parameter_info& info)
      {
        // compute Ts
        std::set<data::variable> FV = find_free_variables(summand.f);
        for (const data::variable& v: summand.e)
        {
          FV.erase(v);
        }
        for (const data::variable& v: FV)
        {
          info.Ts.insert(parameter_position(v));
        }

        // compute Ws and Rs
        auto gi = summand.g.begin();
        auto di = d.begin();
        for ( ; di != d.end(); ++di, ++gi)
        {
          if (*di != *gi)
          {
            std::size_t i = di - d.begin();
            info.Ws.insert(i);

            for (const data::variable& v: find_free_variables(*gi))
            {
              info.Rs.insert(parameter_position(v));
            }
          }
        }

        // compute Vs
        info.Vs = set_union(info.Ts, set_union(info.Ws, info.Rs));
      };

      for (std::size_t k = 0; k < N; k++)
      {
        compute_parameter_info(m_summand_classes[k], info[k]);
      }

      compute_dependency_NES();
      compute_DNA_DNL_NES(info);
    }

    bool compute_deterministic_equations(std::size_t k)
    {
      const summand_class& summand_k = m_summand_classes[k];

      std::size_t n = m_pbes.equations().size();
      for (std::size_t i = 0; i < n; i++)
      {
        if (summand_k.nxt[i].size() >= 2)
        {
          return false;
        }
      }
      return true;
    }

    bool compute_deterministic_data(std::size_t k)
    {
      const summand_class& summand_k = m_summand_classes[k];

      const data::variable_list& parameters = m_pbes.equations()[0].variable().parameters();
      data::set_identifier_generator id_gen;
      for(const data::variable& v: parameters)
      {
        id_gen.add_identifier(v.name());
      }

      summand_equivalence_key new1_k = rename_duplicate_variables(id_gen, summand_equivalence_key(summand_k));
      summand_equivalence_key new2_k = rename_duplicate_variables(id_gen, summand_equivalence_key(summand_k));
      data::variable_list qvars1_k = new1_k.e; data::data_expression condition1_k = new1_k.f; data::data_expression_list updates1_k = new1_k.g;
      data::variable_list qvars2_k = new2_k.e; data::data_expression condition2_k = new2_k.f; data::data_expression_list updates2_k = new2_k.g;

      data::data_expression antecedent = data::sort_bool::and_(condition1_k, condition2_k);
      data::data_expression consequent = data::sort_bool::true_();
      auto it1_k = updates1_k.begin();
      auto it2_k = updates2_k.begin();
      while (it1_k != updates1_k.end())
      {
        consequent = data::lazy::and_(consequent, data::equal_to(*it1_k, *it2_k));
        ++it1_k; ++it2_k;
      }
      data::data_expression condition = data::optimized_forall(parameters + qvars1_k + qvars2_k, data::sort_bool::implies(antecedent, consequent));

      return is_true(condition);
    }

    void compute_deterministic()
    {
      if (!m_options.compute_determinism || !m_options.reduction)
      {
        return;
      }
      std::size_t N = m_summand_classes.size();
      for (std::size_t k = 0; k < N; k++)
      {
        m_summand_classes[k].is_deterministic = compute_deterministic_equations(k) && compute_deterministic_data(k);
      }
    }

    void compute_summand_classes()
    {
      std::size_t n = m_pbes.equations().size();

      for (const srf_equation& eqn: m_pbes.equations())
      {
        for (const srf_summand& summand: eqn.summands())
        {
          summand_equivalence_key key(summand);
          auto i = m_summand_index.find(key);
          if (i == m_summand_index.end())
          {
            std::size_t k = m_summand_index.size();
            m_summand_index[key] = k;
            m_summand_classes.emplace_back(summand.parameters(), summand.condition(), summand.variable().parameters(), n);
          }
        }
      }
      for(summand_class& s: m_summand_classes)
      {
        s.set_num_summands(m_summand_classes.size());
      }
      compute_nxt();
      // optional steps
      compute_NES_DNA_DNL();
      compute_deterministic();
    }

    void compute_vis_invis()
    {
      using utilities::detail::contains;

      std::size_t n = m_pbes.equations().size();
      std::size_t N = m_summand_classes.size();

      m_vis.resize(N);
      for (std::size_t i = 0; i < n; i++)
      {
        const srf_equation& eqn = m_pbes.equations()[i];
        const core::identifier_string& X_i = eqn.variable().name();
        bool op_i = eqn.is_conjunctive();
        std::size_t rank_i = m_equation_index.rank(X_i);

        for (const srf_summand& summand: eqn.summands())
        {
          const core::identifier_string& X_j = summand.variable().name();
          std::size_t j = m_equation_index.index(X_j);
          std::size_t rank_j = m_equation_index.rank(X_j);
          bool op_j = m_pbes.equations()[j].is_conjunctive();
          bool is_invisible = op_i == op_j && rank_i == rank_j;
          if (!is_invisible)
          {
            std::size_t k = summand_index(summand);
            m_vis.set(k);
          }
        }
      }

      // Invis is the opposite of vis
      m_invis = m_vis;
      m_invis.flip();
    }

    std::string print_variables(const data::variable_list& v) const
    {
      std::ostringstream out;
      for (auto i = v.begin(); i != v.end(); ++i)
      {
        if (i != v.begin())
        {
          out << ", ";
        }
        out << *i << ": " << i->sort();
      }
      return out.str();
    }

    void print_summand(const srf_summand& summand, bool is_conjunctive) const
    {
      std::size_t k = summand_index(summand);
      mCRL2log(log::verbose) << "   (" << k << ") ";
      if (!summand.parameters().empty())
      {
        mCRL2log(log::verbose) << (is_conjunctive ? "forall " : "exists ") << print_variables(summand.parameters()) << ". ";
      }
      mCRL2log(log::verbose) << summand.condition()
                << (is_conjunctive ? " => " : " && ")
                << summand.variable()
                << std::endl;
    }

    void print_pbes() const
    {
      mCRL2log(log::verbose) << "srf_pbes" << std::endl;
      for (const srf_equation& eqn: m_pbes.equations())
      {
        mCRL2log(log::verbose) << eqn.symbol() << " " << eqn.variable() << " = " << (eqn.is_conjunctive() ? "conjunction" : "disjunction") << " of summands\n";
        for (const srf_summand& summand: eqn.summands())
        {
          print_summand(summand, eqn.is_conjunctive());
        }
        mCRL2log(log::verbose) << std::endl;
      }
    }

    void print_summand_classes() const
    {
      using utilities::detail::contains;

      if(mCRL2logEnabled(log::verbose))
      {
        std::size_t N = m_summand_classes.size();
        for (std::size_t k = 0; k < N; k++)
        {
          const summand_class& summand = m_summand_classes[k];
          mCRL2log(log::verbose) << "\n--- summand class " << k << " ---" << std::endl;
          mCRL2log(log::verbose) << "visible = " << std::boolalpha << m_vis.test(k) << "\n";
          summand.print(log::logger(log::verbose).get());
        }
        for (std::size_t i = 0; i < m_pbes.equations().size(); i++)
        {
          mCRL2log(log::verbose) << "dependency NES[" << std::setw(3) << i << "]  " << print_summand_set(m_dependency_nes[i]) << std::endl;
        }
      }
    }

  public:
    explicit partial_order_reduction_algorithm(const pbes& p,
          pbespor_options options
        )
     : m_rewr(p.data(),
              //TODO temporarily disabled used_data_equation_selector so the rewriter can rewrite accordance conditions
              // data::used_data_equation_selector(p.data(), pbes_system::find_function_symbols(p), p.global_variables()),
              options.rewrite_strategy),
       m_enumerator(m_rewr, p.data(), m_rewr, m_id_generator, true),
       m_pbes(pbes2srf(p)),
       m_equation_index(m_pbes),
       m_dependency_nes(m_pbes.equations().size()),
       m_solver(options.use_smt_solver ? std::make_unique<smt::smt_solver>(p.data()) : nullptr),
       m_options(options)
    {
      unify_parameters(m_pbes, false, true);

      // initialize m_parameters and m_parameter_positions
      const data::variable_list& parameters = m_pbes.equations().front().variable().parameters();
      m_parameters = std::vector<data::variable>{parameters.begin(), parameters.end()};
      for (std::size_t m = 0; m < m_parameters.size(); m++)
      {
        m_parameter_positions[m_parameters[m]] = m;
      }

      const std::chrono::time_point<std::chrono::high_resolution_clock> t_start =
        std::chrono::high_resolution_clock::now();
      compute_summand_classes();
      compute_vis_invis();
      m_static_analysis_duration = std::chrono::high_resolution_clock::now() - t_start;

      // initialize m_largest_equation_size;
      for (const srf_equation& eq: m_pbes.equations())
      {
        m_largest_equation_size = std::max(m_largest_equation_size, eq.summands().size());
      }

      print_pbes();
    }

    const propositional_variable_instantiation& initial_state() const
    {
      return m_pbes.initial_state();
    }

    const std::vector<data::variable>& parameters() const
    {
      return m_parameters;
    }

    const fixpoint_symbol& symbol(const core::identifier_string& X) const
    {
      std::size_t i = m_equation_index.index(X);
      return m_pbes.equations()[i].symbol();
    }

    void print() const
    {
      print_summand_classes();
    }

    template <
      typename EmitNode = utilities::skip,
      typename EmitEdge = utilities::skip
    >
    void explore(
      const propositional_variable_instantiation& X_init,
      EmitNode emit_node = EmitNode(),
      EmitEdge emit_edge = EmitEdge()
    )
    {
      using utilities::detail::contains;
      using utilities::detail::has_empty_intersection;
      using utilities::detail::set_difference;
      using utilities::detail::set_includes;
      using utilities::detail::set_intersection;
      using utilities::detail::set_union;

      const std::chrono::time_point<std::chrono::high_resolution_clock> t_start =
        std::chrono::high_resolution_clock::now();

      enum todo_state
      {
        NEW,            ///< Will be partially expanded
        DONE_PARTIALLY, ///< Has been partially expanded
        STARTS_CYCLE,   ///< Needs to be fully expanded, because it starts a cycle
        DONE            ///< Has been fully expanded
      };
      using todo_pair = std::pair<propositional_variable_instantiation, todo_state>;

      // The set seen also stores for each node an index and a boolean that expresses whether
      // the node is currently in the DFS stack and not fully explored, ie, its state in the
      // stack is DONE_PARTIALLY.
      std::unordered_map<propositional_variable_instantiation, std::pair<std::size_t, bool>> seen;
      std::deque<todo_pair> todo{todo_pair(X_init, NEW)};
      // Each state is given unique index, based on the order of discovery.
      // This means that the indices in the DFS stack are sorted from low to high.
      std::size_t index = 0;

      {
        std::size_t rank = m_equation_index.rank(X_init.name());
        std::size_t i = m_equation_index.index(X_init.name());
        bool is_conjunctive = m_pbes.equations()[i].is_conjunctive();
        emit_node(X_init, is_conjunctive, rank);
        seen.insert(std::make_pair(X_init, std::make_pair(index, true)));
        index++;
      }

      std::size_t iteration = 0;
      while (!todo.empty())
      {
        todo_pair& p = todo.back();
        const propositional_variable_instantiation X_e = p.first;
        todo_state& s = p.second;
        mCRL2log(log::debug) << "choose X_e = " << X_e << std::endl;

        if (s == DONE || s == DONE_PARTIALLY)
        {
          todo.pop_back();
          seen[X_e].second = false;
          continue;
        }

        std::set<propositional_variable_instantiation> next;
        summand_set en_X_e = en(X_e);

        if (s == NEW)
        {
          summand_set stubborn_set_X_e = stubborn_set(X_e, en_X_e);
          mCRL2log(log::debug) << "stubborn_set(X_e) = " << print_summand_set(stubborn_set_X_e) << std::endl;
          next = succ(X_e, stubborn_set_X_e & en_X_e);

          bool vis_expanded = m_vis.is_subset_of(stubborn_set_X_e);
          s = vis_expanded ? DONE : DONE_PARTIALLY;
          if (!vis_expanded)
          {
            seen[X_e].second = true;
          }

          if (m_options.use_condition_L)
          {
            // Check if a cycle is closed
            // At the same time, check whether some node on the stack is fully expanded
            // If both are true, some node will be fully expanded
            std::size_t num_cycles = 0;
            propositional_variable_instantiation min_node;
            for (const propositional_variable_instantiation& Y_f: next)
            {
              auto node = seen.find(Y_f);
              if (node == seen.end())
              {
                continue;
              }
              std::size_t node_instack = node->second.second;
              if (node_instack)
              {
                num_cycles++;
                min_node = Y_f;
                if (num_cycles > 1)
                {
                  break;
                }
              }
            }
            if (num_cycles == 1)
            {
              auto it = std::find_if(todo.rbegin(), todo.rend(), [&](const auto& pair){ return pair.first == min_node; });
              assert(it != todo.rend());
              auto& [Y_f, Y_f_state] = *it;
              if(Y_f_state == DONE_PARTIALLY)
              {
                Y_f_state = STARTS_CYCLE;
                seen[Y_f].second = false;
              }
            }
            else if (num_cycles > 1 && s == DONE_PARTIALLY)
            {
              // Generate the remainder of the successors
              for (auto& Y_f: succ(X_e, en_X_e - stubborn_set_X_e))
              {
                next.insert(Y_f);
              }
              s = DONE;
              seen[X_e].second = false;
            }
          }
        }
        else
        {
          assert(s == STARTS_CYCLE);
          next = succ(X_e, en_X_e);
          s = DONE;
        }

        mCRL2log(log::debug) << "next = " << core::detail::print_set(next) << std::endl;
        for (const propositional_variable_instantiation& Y_f: next)
        {
          if (seen.find(Y_f) == seen.end())
          {
            std::size_t rank = m_equation_index.rank(Y_f.name());
            std::size_t i = m_equation_index.index(Y_f.name());
            bool is_conjunctive = m_pbes.equations()[i].is_conjunctive();
            emit_node(Y_f, is_conjunctive, rank);
            seen.insert(std::make_pair(Y_f, std::make_pair(index, false)));
            index++;
            todo.emplace_back(Y_f, NEW);
          }
        }
        for (const propositional_variable_instantiation& Y_f: next)
        {
          emit_edge(X_e, Y_f);
        }

        iteration++;
        if(iteration == 100)
        {
          mCRL2log(log::status) << "Found " << seen.size() << " nodes. Todo set contains " << todo.size() << " nodes.\n";
          iteration = 0;
        }
      }
      mCRL2log(log::verbose) << "Finished exploration, found " << seen.size() << " nodes." << std::endl;

      m_exploration_duration = std::chrono::high_resolution_clock::now() - t_start;
      mCRL2log(log::info) << "timing pbespor (wall clock time in seconds):"
        "\n  static analysis: " << std::chrono::duration<double>(m_static_analysis_duration).count() <<
        "\n  exploration:     " << std::chrono::duration<double>(m_exploration_duration).count() << std::endl;
    }

    template <
      typename EmitNode = utilities::skip,
      typename EmitEdge = utilities::skip
    >
    void explore_full(
      const propositional_variable_instantiation& X_init,
      EmitNode emit_node = EmitNode(),
      EmitEdge emit_edge = EmitEdge()
    )
    {
      const std::chrono::time_point<std::chrono::high_resolution_clock> t_start =
        std::chrono::high_resolution_clock::now();

      std::unordered_set<propositional_variable_instantiation> seen;
      std::deque<propositional_variable_instantiation> todo{ X_init };

      {
        std::size_t rank = m_equation_index.rank(X_init.name());
        std::size_t i = m_equation_index.index(X_init.name());
        bool is_conjunctive = m_pbes.equations()[i].is_conjunctive();
        emit_node(X_init, is_conjunctive, rank);
        seen.insert(X_init);
      }

      std::size_t N = m_summand_classes.size();
      summand_set summands_X(N);

      std::size_t iteration = 0;
      while (!todo.empty())
      {
        const propositional_variable_instantiation X_e = todo.back();
        todo.pop_back();
        mCRL2log(log::debug) << "choose X_e = " << X_e << std::endl;

        std::size_t X_index = m_equation_index.index(X_e.name());
        for(std::size_t i = 0; i < N; i++)
        {
          if (depends(X_index, i))
          {
            summands_X.set(i);
          }
        }
        mCRL2log(log::debug) << "enabled according to dependencies = " << print_summand_set(summands_X) << std::endl;
        std::set<propositional_variable_instantiation> next = succ(X_e, summands_X);
        mCRL2log(log::debug) << "next = " << core::detail::print_set(next) << std::endl;
        summands_X.reset();

        for (const propositional_variable_instantiation& Y_f: next)
        {
          if (seen.find(Y_f) == seen.end())
          {
            std::size_t rank = m_equation_index.rank(Y_f.name());
            std::size_t i = m_equation_index.index(Y_f.name());
            bool is_conjunctive = m_pbes.equations()[i].is_conjunctive();
            emit_node(Y_f, is_conjunctive, rank);
            seen.insert(Y_f);
            todo.emplace_back(Y_f);
          }
        }
        for (const propositional_variable_instantiation& Y_f: next)
        {
          emit_edge(X_e, Y_f);
        }

        iteration++;
        if(iteration == 100)
        {
          mCRL2log(log::status) << "Found " << seen.size() << " nodes. Todo set contains " << todo.size() << " nodes.\n";
          iteration = 0;
        }
      }
      mCRL2log(log::verbose) << "Finished exploration, found " << seen.size() << " nodes." << std::endl;

      m_exploration_duration = std::chrono::high_resolution_clock::now() - t_start;
      mCRL2log(log::info) << "timing pbespor (wall clock time in seconds):"
        "\n  exploration:     " << std::chrono::duration<double>(m_exploration_duration).count() << std::endl;
    }
};

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_PARTIAL_ORDER_REDUCTION_H
