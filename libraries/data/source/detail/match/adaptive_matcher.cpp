// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/data/detail/match/adaptive_matcher.h"

#include "mcrl2/atermpp/aterm_io.h"
#include "mcrl2/utilities/stopwatch.h"
#include "mcrl2/core/index_traits.h"
#include "mcrl2/data/detail/enumerator_identifier_generator.h"
#include "mcrl2/data/detail/match/linearise.h"
#include "mcrl2/data/substitutions/sequence_substitution.h"

/// \brief Print the intermediate matches that succeeded.
constexpr bool PrintMatchSteps = false;

/// \brief Print the intermediate steps performed during construction.
constexpr bool PrintConstructionSteps = false;

// Various optimizations.

/// \brief Compute index positions in restrict and return a set of positions when available.
constexpr bool EnableIndexPositions = true;

/// \brief Remove positions where every pattern contains a variable.
constexpr bool EnableRemoveVariables = true;

/// \brief Stop whenever the prefix matches on of the left-hand sides already.
constexpr bool EnableGreedyMatching = true;

/// \brief Enable brute force minimization of the resulting automaton.
constexpr bool EnableBruteForce = false;

/// \brief Enable to reuse position indices based on the observation that each position is inspected at most once.
constexpr bool EnableReusableIndex = false;

using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

/// \brief Represents the unnamed variable omega.
class dontcare : public function_symbol
{
public:
  dontcare()
    : function_symbol(mcrl2::core::identifier_string("@@dont_care"), mcrl2::data::sort_expression())
  {}
};

class equal : public function_symbol
{
public:
  equal()
    : function_symbol(mcrl2::core::identifier_string("@@equal"), mcrl2::data::sort_expression())
  {}
};

class not_equal : public function_symbol
{
public:
  not_equal()
    : function_symbol(mcrl2::core::identifier_string("@@not_equal"), mcrl2::data::sort_expression())
  {}
};

/// \returns True iff the given data expression is of type omega.
inline bool is_dontcare(const data_expression& expression)
{
  return expression == dontcare();
}

/// \brief Decides whether the left and right terms unify, assuming that vars(left) and vars(right) are disjoint.
bool unify(const atermpp::aterm_appl& left, const atermpp::aterm_appl& right)
{
  if (is_variable(left) || is_variable(right))
  {
    return true;
  }
  else if (is_function_symbol(left) && is_function_symbol(right))
  {
    return left == right;
  }
  else
  {
    // The term and lhs are applications.
    const application& lhs_appl  = static_cast<const application&>(right);
    const application& term_appl = static_cast<const application&>(left);

    // Both must have the same arity, the head symbol must match and their arguments must match.
    if (lhs_appl.size() != term_appl.size())
    {
      return false;
    }

    if (!unify(term_appl.head(), lhs_appl.head()))
    {
      return false;
    }

    for (std::size_t i = 0; i < term_appl.size(); i++)
    {
      if (!unify(term_appl[i], lhs_appl[i]))
      {
        return false;
      }
    }

    return true;
  }
}

/// \brief Return true iff the given term matches the left-hand side.
bool matches(const data_expression& term,  const data_expression& lhs)
{
  if (is_function_symbol(lhs))
  {
    return term == lhs;
  }
  else if (is_variable(lhs))
  {
    return true;
  }
  else
  {
    // The term and lhs are applications.
    const application& lhs_appl  = static_cast<const application&>(lhs);
    const application& term_appl = static_cast<const application&>(term);

    // Both must have the same arity, the head symbol must match and their arguments must match.
    if (lhs_appl.size() != term_appl.size())
    {
      return false;
    }

    if (!matches(term_appl.head(), lhs_appl.head()))
    {
      return false;
    }

    for (std::size_t i = 0; i < term_appl.size(); i++)
    {
      if (!matches(term_appl[i], lhs_appl[i]))
      {
        return false;
      }
    }

    return true;
  }
}

/// \brief Lexicographical comparisons of left < right.
bool less_than(const position& left, const position& right)
{
  // Every element of left and right up to the length of smallest position should be less or equal.
  for (std::size_t index = 0; index <= std::min(left.size(), right.size()); ++index)
  {
    if (left[index] > right[index])
    {
      return false;
    }
  }

  return true;
}

/// \returns The number of arguments for the given data expressions.
inline std::size_t get_arity(const data_expression& t)
{
  if (!is_application(t) || is_dontcare(t))
  {
    return 0;
  }

  const application& ta = atermpp::down_cast<application>(t);
  return ta.size();
}

/// \returns Checks whether t is an actual application.
/// \details Does not fail whenever t is a newly introduced data_expression or whenever it is a sort_expression and whatever.
inline bool is_application_robust(const data_expression& t)
{
  return std::find_if(mcrl2::core::detail::function_symbols_DataAppl.begin(),
    mcrl2::core::detail::function_symbols_DataAppl.end(),
    [&](const std::unique_ptr<atermpp::function_symbol>& ptr)
    {
      return *ptr.get() == t.function();
    }
    ) != mcrl2::core::detail::function_symbols_DataAppl.end();
}

// Return the head symbol, nested within applications.
inline const data_expression& get_head(const data_expression& t)
{
  if (!is_application(t) || is_dontcare(t))
  {
    return t;
  }

  const application& ta = atermpp::down_cast<application>(t);
  return get_head(ta.head());
}

/// \returns A unique index for the head symbol that the given term starts with.
inline std::size_t get_head_index(const data_expression& term)
{
  return mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(get_head(term)));
}

// Public functions

template<typename Substitution>
AdaptiveMatcher<Substitution>::AdaptiveMatcher(const data_equation_vector& equations)
  : m_automaton()
{
  // Keep track of the construction time.
  mcrl2::utilities::stopwatch construction;

  // Preprocess the term rewrite system to be linear.
  std::size_t nof_nonlinear_equations = 0;
  std::vector<indexed_linear_data_equation> linear_equations;
  for (const data_equation& old_equation : equations)
  {
    // Rename the variables in the equation
    auto [equation, partition] = make_linear(old_equation);

    // If the partition is not empty then this equation is nonlinear.
    if (!partition.empty())
    {
      ++nof_nonlinear_equations;
    }

    // Add the index of the equation
    linear_equations.emplace_back(indexed_linear_data_equation(equation, partition, m_positions));
  }

  // Determine the index of dontcare.
  m_dontcare_index = mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(dontcare()));
  m_not_equal_index = mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(not_equal()));
  m_equal_index = mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(equal()));

  // Construct the automaton.
  try
  {
    const std::set<std::pair<std::size_t, std::size_t> > E;
    const std::set<std::pair<std::size_t, std::size_t> > N;
    
    m_automaton = construct_apma(position_variable(position()), linear_equations, E, N);

    // Metrics for the automaton.
    std::size_t nof_ambiguous_matches = 0; ///< The number of final states with multiple matches.
    std::size_t nof_final_states = 0; ///< The number of final states.
    for (std::size_t state = 1; state < m_automaton.states(); ++state)
    {
      //if (m_automaton.transition())

    }

    mCRL2log(debug) << "AdaptiveMatcher (EnableIndexPositions = " << EnableIndexPositions
      << ", EnableRemoveVariables = " << EnableRemoveVariables
      << ", EnableGreedyMatching = " << EnableGreedyMatching
      << ", EnableBruteForce = " << EnableBruteForce
      << "):\n";
    mCRL2log(debug) << "  There are " << nof_nonlinear_equations << " nonlinear left-hand sides out of " << linear_equations.size() << " rewrite rules.\n";

    mCRL2log(debug) << "  Automaton (states: " << m_automaton.states() << ", transitions: " << m_automaton.transitions() << ") construction took " << construction.time() << " milliseconds.\n";

    mCRL2log(debug) << "  There are " << m_positions.size() << " positions indexed.\n";

    mCRL2log(debug) << "  There are " << nof_ambiguous_matches << " ambiguous match sets out of " << nof_final_states << " final states.\n";

    // We might consider to make the position indices dense again.

    // Ensure that the subterm indexing can store all possible terms in the required places.
    m_subterms.resize(m_positions.size());
  }
  catch (const std::exception& ex)
  {
    mCRL2log(debug) << "AdaptiveMatcher: construction failed with " << ex.what() << ".\n";
  }
}

template<typename Substitution>
typename AdaptiveMatcher<Substitution>::const_iterator AdaptiveMatcher<Substitution>::match(const data_expression& term, Substitution& matching_sigma)
{
  // Start with the root state.
  std::size_t s = m_automaton.root();
  std::size_t s_old;

  // Store the root position where it is expected.
  m_subterms[m_automaton.label(m_automaton.root()).position] = term;

  while (true)
  {
    // Store the labelling of the current state for convenience.
    const apma_state& state = m_automaton.label(s);

    // Keep track of the current state.
    s_old = s;

    if (state.is_matching_state())
    {
      // t[p] is given by the subterm table at the current position.
      const data_expression& subterm = static_cast<const data_expression&>(m_subterms[state.position]);
      assert (subterm.defined());

      if (PrintMatchSteps) { mCRL2log(info) << "Matching m_subterms[" << state.position << "] = " << subterm << "\n"; }

      // The number of arguments must match the current state.
      if (is_application(subterm))
      {
        const auto& appl = static_cast<const application&>(subterm);

        // If delta(s, f) = (s', update)
        std::size_t s_prime = m_automaton.transition(s, get_head_index(appl.head()));
        if (s_prime != 0)
        {
          if (PrintMatchSteps) { mCRL2log(info) << "Took transition from " << s << " to " << s_prime << " with label " << appl.head() << "\n"; }
          s = s_prime;

          // Insert the subterms into the position mapping.
          auto it = state.argument_positions.begin();
          for (const atermpp::aterm& argument : appl)
          {
            assert(it != state.argument_positions.end());
            assert(*it < m_subterms.size());
            m_subterms[*it] = argument;
            ++it;
          }
        }
        else
        {
          // If delta(s, \neq)
          std::size_t s_prime = m_automaton.transition(s, m_dontcare_index);
          if (s_prime != 0)
          {
            if (PrintMatchSteps) { mCRL2log(info) << "Took transition from " << s << " to " << s_prime << " with label don't care.\n"; }

            s = s_prime;
          }
        }
      }
      // This transition can only be taken for constant function symbols.
      else if (is_function_symbol(subterm))
      {
        // If delta(s0, a) is defined for some term a followed by suffix t'.
        std::size_t s_prime = m_automaton.transition(s, get_head_index(subterm));
        if (s_prime != 0)
        {
          if (PrintMatchSteps) { mCRL2log(info) << "Took transition " << s << " to " << s_prime << " with label " << static_cast<function_symbol>(subterm) << "\n"; }
          s = s_prime;
        }
        else
        {
          // If delta(s, \dc)
          std::size_t s_prime = m_automaton.transition(s, m_dontcare_index);
          if (s_prime != 0)
          {
            if (PrintMatchSteps)
            {
              mCRL2log(info) << "Took transition from " << s << " to " << s_prime << " with label don't care.\n";
            }

            s = s_prime;
          }
        }
      }

      // If we are stuck in the current state then matching failed.
      if (s == s_old)
      {
        if (PrintMatchSteps) { mCRL2log(info) << "Matching failed.\n"; }
        return const_iterator(matching_sigma);
      }
    }
    else if (state.is_consistency_state())
    {
      assert(m_subterms[state.compare.first].defined());
      assert(m_subterms[state.compare.second].defined());
      if (m_subterms[state.compare.first] == m_subterms[state.compare.second])
      {
        std::size_t s_prime = m_automaton.transition(s, m_equal_index);
        assert(s_prime != 0);
        if (PrintMatchSteps) { mCRL2log(info) << "Took transition from " << s << " to " << s_prime << " with label equal.\n"; }
        s = s_prime;
      }
      else
      {
        std::size_t s_prime = m_automaton.transition(s, m_not_equal_index);
        assert(s_prime != 0);
        if (PrintMatchSteps) { mCRL2log(info) << "Took transition from " << s << " to " << s_prime << " with label not_equal.\n"; }
        s = s_prime;
      }
    }
    // 1. If s in Sfin then return Lfin(s).
    else
    {
      // 3.1 (R, P) := Lfin(s)
      // 3.2 sigma := id

      // 3.3 for p in P do sigma := sigma[x_p -> t[p]]
      for (const auto& [var, pos] : m_automaton.label(s).assignments)
      {
        assert(m_subterms[pos].defined());
        auto& expression = static_cast<const data_expression&>(m_subterms[pos]);
        if (PrintMatchSteps) { mCRL2log(info) << "sigma(" << var << ") := " << expression << ".\n"; }
        matching_sigma[var] = expression;
      }

      // 3.5 Return (R, sigma)
      if (PrintMatchSteps) { mCRL2log(info) << "Matching succeeded.\n"; }

      return const_iterator(&m_automaton.label(s).match_set, m_subterms, matching_sigma);
    }
  }
}

// Private functions

template<typename Substitution>
typename AdaptiveMatcher<Substitution>::Automaton AdaptiveMatcher<Substitution>::construct_apma(
  data_expression pref,
  const std::vector<indexed_linear_data_equation>& old_L,
  const std::set<std::pair<std::size_t, std::size_t> >& E,
  const std::set<std::pair<std::size_t, std::size_t> >& N)
{
  //if (PrintConstructionSteps) { mCRL2log(info) << "state = " << state << "("; }

  // 2. L := { l in L | l unifies with pref }
  std::vector<indexed_linear_data_equation> L;
  for (const indexed_linear_data_equation& equation : old_L)
  {
    bool can_be_consistent = true;
    // not exists C in P : exists x, y in C : {x, y} in N
    for (const std::vector<std::size_t>& consistency : equation.partition())
    {
      for (std::size_t x : consistency)
      {
        for (std::size_t y : consistency)
        {
          if (x < y)
          {
            if (N.find(std::make_pair(x, y)) != N.end())
            {
              can_be_consistent = false;
              break;
            }
          }
        }

        if (can_be_consistent == false) { break; }
      }
      if (can_be_consistent == false) { break; }
    }

    if (unify(equation.equation().lhs(), pref) && can_be_consistent)
    {
      L.emplace_back(equation);
    }
  }

  std::set<std::pair<std::size_t, std::size_t>> workC = {};

  for (const indexed_linear_data_equation& equation : L)
  {
    for (const std::vector<std::size_t>& consistency : equation.partition())
    {
      if (!consistency.empty())
      {
        std::size_t x = *consistency.begin();

        for (std::size_t y : consistency)
        {
          if (x != y)
          {
            std::size_t small = std::min(x, y);
            std::size_t large = std::max(x, y);

            if (E.find(std::make_pair(small, large)) == E.end()
              && at_position(pref, m_positions[x])
              && at_position(pref, m_positions[y]))
            {
              workC.insert(std::make_pair(small, large));
            }
          }
        }
      }
    }
  }

  // 3. F := restrict(fringe(pref), pref).
  std::set<position> workF = restrict(fringe(pref), L);

  // Greedy: If any of the elements in L match then we can chose that one.
  if constexpr (EnableGreedyMatching)
  {
    // If the prefix matches some left-hand side that was already linear.
    auto it = std::find_if(L.begin(),
      L.end(),
      [&pref](const indexed_linear_data_equation& equation)
      {
        return matches(pref, equation.equation().lhs()) && equation.partition().empty() && equation.equation().condition() == sort_bool::true_();
      });

    if (it != L.end())
    {
      // Change L to only be this equation.
      workF.clear();
      workC.clear();
      L = {*it};
    }
  }

  // The resulting automaton.
  Automaton M;

  // if (workC = emptyset and workF = emptyset) or L = emptyset
  if ((workF.empty() && workC.empty()) || L.empty())
  {
    // The labelling for the current state.
    apma_state& root = M.label(M.root());

    // Postprocessing: R := { r_i | l_i in L(s) }
    root.match_set.insert(root.match_set.begin(), L.begin(), L.end());

    // Postprocessing: P := union r_i in R : fringe(r_i) \ { L(s') in path(s) }
    if (PrintConstructionSteps) { mCRL2log(info) << "P = "; }

    // positions_rhs := { fringe(r_i) | l_i in L }.
    std::set<position> P;

    for (const indexed_linear_data_equation& equation : L)
    {
      // Find the variables of the rhs and condition.
      std::set<variable> vars = data::find_free_variables(equation.equation().rhs());
      std::set<variable> condition_vars = data::find_free_variables(equation.equation().condition());
      vars.insert(condition_vars.begin(), condition_vars.end());

      // Find where these variables occur in the left-hand side (equivalently in the prefix)
      std::set<position> lhs_fringe = fringe(equation.equation().lhs());

      // Find the corresponding position in the lhs for each variable.
      for (const variable& var : vars)
      {
        for(const position& pos : lhs_fringe)
        {
          std::optional<data_expression> expr = at_position(equation.equation().lhs(), pos);
          if (expr)
          {
            assert(is_variable(expr.value()));
            if (expr.value() == var)
            {
              P.insert(pos);
            }
          }
        }
      }

      assert(vars.size() <= P.size());
    }

    // Convert P to a vector for faster access.
    for (const position& pos : P)
    {
      if (PrintConstructionSteps) { mCRL2log(info) << pos << ", "; }
      root.assignments.emplace_back(std::make_pair(position_variable(pos), position_index(pos)));
    }

    // L' := L'[s -> (R, P)]
    if (PrintConstructionSteps) { mCRL2log(info) << ") \n"; }
  }
  else if (!workF.empty())
  {
    // Find the smallest automaton for all positions, only perform consistency checks when matching is finished (for now).
    while (!workF.empty())
    {
      // pos := select(F).
      position pos = select(workF, L);

      // A temporary automaton for the current choice.
      Automaton M_prime;

      // The labelling for the current state.
      apma_state& root = M_prime.label(M_prime.root());

      // M := M[L := L[s -> pos]]
      root.position = position_index(pos);
      if (PrintConstructionSteps) { mCRL2log(info) << "L = " << pos << ")\n"; }

      // for f in F s.t. exist l in L : head(l[pos]) = f
      std::set<std::pair<mcrl2::data::function_symbol, std::size_t>> symbols;

      for (const indexed_linear_data_equation& equation : L)
      {
        // t := l[pos]
        std::optional<data_expression> t = at_position(equation.equation().lhs(), pos);

        if (t && (is_function_symbol(t.value()) || is_application_robust(t.value())))
        {
          // head(l[pos]) in F.
          symbols.insert(std::make_pair(static_cast<const function_symbol&>(get_head(t.value())), get_arity(t.value())));
        }
      }

      std::size_t max_arity = 0;
      for (const auto& [symbol, arity] : symbols)
      {
        // M := M[S := (S cup {s'})], where s' is a fresh state.
        std::size_t s_prime = M_prime.add_state();

        // M := M[delta := (delta := delta(s, f) -> s')
        M_prime.add_transition(M_prime.root(), mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(symbol), s_prime);

        //if (PrintConstructionSteps) { mCRL2log(info) << "added transition from "
        //  << "0 with label " << static_cast<atermpp::aterm>(symbol) << " to state " << s_prime << ".\n"; }

        // pref := f(omega^ar(f)), but we store in the correct position variables directly in the prefix.
        std::vector<data_expression> arguments;

        if (arity > 0)
        {
          // Name the arguments of the new prefix.
          position var_position = pos;
          var_position.push_back(0);

          // In this case there are multiple variables are introduced, their positions are derived from current position.i where i is their
          // position in the application f(x_0, ..., x_n).
          for (std::size_t index = 0; index < arity; ++index)
          {
            var_position.back() = index + 1;
            arguments.push_back(position_variable(var_position));
          }

          M_prime.merge(s_prime, construct_apma(assign_at_position(pref, pos, application(symbol, arguments.begin(), arguments.end())), L, E, N));
          max_arity = std::max(max_arity, arity);
        }
        else
        {
          // In this case the symbol is just a function symbol.
          M_prime.merge(s_prime, construct_apma(assign_at_position(pref, pos, symbol), L, E, N));
        }
      }

      // Here, we also ensure that the arguments of this symbol can be stored in the subterm table during evaluation.
      position var_position = pos;
      var_position.push_back(0);
      for (std::size_t index = 0; index < max_arity; ++index)
      {
        var_position.back() = index + 1;
        root.argument_positions.emplace_back(position_index(var_position));
      }

      // if exists l in L : (exists pos' <= pos : head(l[pos']) in V then
      if (std::find_if(L.begin(),
        L.end(),
        [&](const indexed_linear_data_equation& equation)
        {
          return has_variable_higher(equation.equation().lhs(), pos);
        }) != L.end())
      {
        // M := M[S := (S cup {s'})], where s' is a fresh state.
        std::size_t s_prime = M_prime.add_state();

        // M := M[delta := (delta := delta(s, \neq) -> s')
        M_prime.add_transition(M_prime.root(), m_dontcare_index, s_prime);
        //if (PrintConstructionSteps) { mCRL2log(info) << "added transition from " << state << " with label " << static_cast<atermpp::aterm>(dontcare()) << " to state " << s_prime << ".\n"; }

        M_prime.merge(s_prime, construct_apma(assign_at_position(pref, pos, dontcare()), L, E, N));
      }

      if (EnableBruteForce)
      {
        // Try all possible positions until the smallest automaton has been found.
        if (M_prime.size() < (M.size() == 1 ? std::numeric_limits<std::size_t>::max() : M.size()))
        {
          if (PrintConstructionSteps) { mCRL2log(debug) << "Created smallest automaton of size " << M_prime.size() << " for position " << pos << "\n"; }
          M = M_prime;
        }

        // Remove the selection position from the possibilies.
        workF.erase(pos);
      }
      else
      {
        // Accept the first selected position and exit the loop.
        M = M_prime;
        break;
      }
    }
  }
  else
  {
    // Add the consistency checks
    assert(!workC.empty());
    std::pair<std::size_t, std::size_t> compare = *workC.begin();
    workC.clear();

    // The labelling for the current state.
    apma_state& root = M.label(M.root());
    root.compare = compare;

    if (PrintConstructionSteps) { mCRL2log(info) << "compare = (" << compare.first << ", " << compare.second << ")\n"; }

    // M := M[S := (S cup {s'})], where s' is a fresh state.
    std::size_t s_cmark = M.add_state();

    // M := M[delta := (delta := delta(s, cmark) -> s')
    M.add_transition(M.root(), m_equal_index, s_cmark);
    {
      std::set<std::pair<std::size_t, std::size_t>> new_E = E;
      new_E.insert(compare);
      M.merge(s_cmark, construct_apma(pref, L, new_E, N));
    }

    // M := M[S := (S cup {s'})], where s' is a fresh state.
    std::size_t s_xmark = M.add_state();

    // M := M[delta := (delta := delta(s, cmark) -> s')
    M.add_transition(M.root(), m_not_equal_index, s_xmark);
    {
      std::set<std::pair<std::size_t, std::size_t>> new_N = N;
      new_N.insert(compare);
      M.merge(s_xmark, construct_apma(pref, L, E, new_N));
    }
  }

  return M;
}

template<typename Substitution>
std::set<position> AdaptiveMatcher<Substitution>::restrict(const std::set<position>& F, const std::vector<indexed_linear_data_equation>& L)
{
  // A restriction on F.
  std::set<position> result;

  if constexpr (EnableRemoveVariables)
  {
    // intersection := intersection_{l in L} fringe(l)
    std::set<position> intersection;
    bool first = true;
    for (const indexed_linear_data_equation& equation : L)
    {
      std::set<position> local;
      std::set<position> fringe_l = fringe(equation.equation().lhs());

      if (first)
      {
        intersection = fringe_l;
      }

      std::set_intersection(intersection.begin(), intersection.end(), fringe_l.begin(), fringe_l.end(), std::inserter(local, local.begin()));
      intersection = local;
      first = false;
    }

    // F := F \ intersection_{l in L} fringe(l)
    std::set_difference(F.begin(), F.end(), intersection.begin(), intersection.end(), std::inserter(result, result.begin()));
  }
  else
  {
    // Set F for further computations.
    result = F;
  }

  return result;
}

template<typename Substitution>
position AdaptiveMatcher<Substitution>::select(const std::set<position>& F, const std::vector<indexed_linear_data_equation>& L)
{
  // Compute a set of indices.
  if constexpr(EnableIndexPositions)
  {
    std::set<position> indices;
    for (const position& position : F)
    {
      if (std::all_of(L.begin(), L.end(),
        [&](const indexed_linear_data_equation& equation)
        {
          std::optional<data_expression> result = at_position(equation.equation().lhs(), position);
          if (result)
          {
            return !is_variable(result.value());
          }
          else
          {
            return true;
          }
        }))
      {
        indices.insert(position);
      }
    }

    // If there is an index then return it.
    if (!indices.empty())
    {
      return *indices.begin();
    }
  }

  // This corresponds to the left-to-right depth-first traversal order.
  return *std::min_element(F.begin(), F.end(), less_than);
}

template<typename Substitution>
std::size_t AdaptiveMatcher<Substitution>::position_index(const position& pos)
{
  return m_positions.insert(pos).first;
}

template class mcrl2::data::detail::AdaptiveMatcher<sequence_substitution>;
template class mcrl2::data::detail::AdaptiveMatcher<mutable_map_substitution<>>;
